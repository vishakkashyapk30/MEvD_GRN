"""RegDiffusion baseline (diffusion-based GRN inference, RNA-only; plan 9.2).

Install:  pip install regdiffusion
Fastest DL GRN method — establishes the runtime lower bound. Produces a gene×gene
influence matrix; we score TF->target as |A[tf, tg]|.
"""
from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np

from src.baselines.common import DenseMatrixScorer, load_expression_df, universe_gene_names


def run_regdiffusion(rna_path: str, gene_index: Dict[str, int], pcfg: dict,
                     device: str = "cuda:0",
                     tf_names: Optional[List[str]] = None,
                     max_genes: int = 3500) -> DenseMatrixScorer:
    try:
        import regdiffusion as rd
    except ImportError as e:  # pragma: no cover
        raise ImportError("regdiffusion not installed — `pip install regdiffusion`") from e

    genes = universe_gene_names(gene_index)                    # order == gene_index
    expr_full = load_expression_df(rna_path, gene_index, pcfg)[genes].to_numpy(np.float32)

    # RegDiffusion learns a dense, fully-parameterized (n_gene x n_gene)
    # adjacency matrix plus Adam state over it -- O(n_gene^2) memory that
    # cannot fit ~23k genes on an 8GB GPU regardless of memory_efficient
    # (that flag only avoids materializing extra n^2 buffers during the L1
    # loss, not the parameter/optimizer state itself). Restrict to the
    # union of highly-variable genes (same n_hvg threshold MEvD-GRN uses for
    # its own message-passing graph) and all TFs, so every TF->target edge
    # touching a TF is still scorable; the remaining (typically low-signal)
    # genes score 0, which is reported alongside the results.
    n_hvg = min(int(pcfg.get("n_hvg", 3000)), max_genes)
    variance = expr_full.var(axis=0)
    hvg_idx = set(np.argsort(variance)[::-1][:n_hvg].tolist())
    if tf_names:
        name_to_idx = {g: i for i, g in enumerate(genes)}
        hvg_idx |= {name_to_idx[t] for t in tf_names if t in name_to_idx}
    if len(hvg_idx) > max_genes:
        hvg_idx = set(sorted(hvg_idx, key=lambda i: -variance[i])[:max_genes])
    universe_mask = np.zeros(len(genes), dtype=bool)
    universe_mask[sorted(hvg_idx)] = True
    n_restricted = int((~universe_mask).sum())
    if n_restricted:
        print(f"[regdiffusion] restricting to {int(universe_mask.sum())}/{len(genes)} genes "
              f"(top-{n_hvg} HVG ∪ TFs) -- RegDiffusion's O(n_gene^2) adjacency parameter "
              "does not fit the full gene universe on this GPU; excluded genes score 0",
              flush=True)
    genes_full = genes
    expr_full = expr_full[:, universe_mask]
    genes = [g for g, k in zip(genes_full, universe_mask) if k]

    # RegDiffusion internally min-max normalizes PER CELL, (x - cell_min) /
    # cell_range, then rejects genes whose column is constant afterward. A
    # gene that is always the single highest- (or lowest-) expressed gene in
    # every cell (e.g. a dominant housekeeping gene like EEF1A1) becomes an
    # exact constant (1.0 or 0.0) under that transform for every cell -- a
    # genuine rank property, not a floating-point artifact, so no amount of
    # jitter fixes it. Drop such genes before training (as RegDiffusion's own
    # docs recommend) and reinsert them as all-zero rows/cols afterward so the
    # full gene_index alignment is preserved for scoring.
    cell_min = expr_full.min(axis=1, keepdims=True)
    cell_max = expr_full.max(axis=1, keepdims=True)
    cell_range = np.where(cell_max > cell_min, cell_max - cell_min, 1.0)
    normalized = (expr_full - cell_min) / cell_range
    keep_mask = normalized.std(axis=0) > 0
    n_drop = int((~keep_mask).sum())
    if n_drop:
        dropped = [genes[i] for i in np.where(~keep_mask)[0][:5]]
        print(f"[regdiffusion] excluding {n_drop} gene(s) RegDiffusion's own per-cell "
              f"normalization renders zero-variance (e.g. {dropped}); scored as 0 for "
              "any edge touching them", flush=True)
    kept_genes = [g for g, k in zip(genes, keep_mask) if k]
    expr = expr_full[:, keep_mask]
    print(f"[regdiffusion] {expr.shape[0]} cells × {expr.shape[1]} genes on {device}", flush=True)

    dev = device if str(device).startswith("cuda") else "cpu"
    trainer = _make_trainer(rd, expr, kept_genes, dev)
    trainer.train()
    A_kept = _extract_adjacency(trainer, kept_genes)
    A_kept = np.abs(np.asarray(A_kept, dtype=np.float64))
    np.fill_diagonal(A_kept, 0.0)

    # Reinsert into the FULL gene_index-aligned matrix: first the restricted
    # (HVG ∪ TF) universe, then the RegDiffusion-internal zero-variance drops
    # within it. Genes outside the restricted universe score 0.
    n_genes_full = len(genes_full)
    full_idx = {g: i for i, g in enumerate(genes_full)}
    restricted_pos = np.array([full_idx[g] for g in genes], dtype=np.int64)
    kept_pos = restricted_pos[np.where(keep_mask)[0]]
    A = np.zeros((n_genes_full, n_genes_full), dtype=np.float64)
    A[np.ix_(kept_pos, kept_pos)] = A_kept
    return DenseMatrixScorer(A)


def _make_trainer(rd, expr, genes, device):
    # ~23k genes makes the default RegDiffusion model's dense n_gene x n_gene
    # adjacency (plus optimizer state / L1-reg buffers) exceed an 8GB GPU;
    # memory_efficient=True (RegDiffusionME) avoids materializing that full
    # matrix and is the officially supported large-gene-count path.
    for kwargs in ({"gene_names": genes, "device": device, "memory_efficient": True},
                   {"device": device, "memory_efficient": True},
                   {"gene_names": genes, "device": device}, {"device": device}, {}):
        try:
            return rd.RegDiffusionTrainer(expr, **kwargs)
        except TypeError:
            continue
    return rd.RegDiffusionTrainer(expr)


def _extract_adjacency(trainer, genes) -> np.ndarray:
    # Preferred: GRN object with an adjacency matrix.
    try:
        grn = trainer.get_grn(genes)
        for attr in ("adj_matrix", "adj", "A"):
            if hasattr(grn, attr):
                return getattr(grn, attr)
    except Exception:
        pass
    for attr in ("adj_matrix", "adj", "get_adj"):
        obj = getattr(trainer, attr, None)
        if obj is not None:
            return obj() if callable(obj) else obj
    raise RuntimeError("Could not extract adjacency from RegDiffusion trainer; "
                       "check the installed regdiffusion API.")
