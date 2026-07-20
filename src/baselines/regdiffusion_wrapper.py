"""RegDiffusion baseline (diffusion-based GRN inference, RNA-only; plan 9.2).

Install:  pip install regdiffusion
Fastest DL GRN method — establishes the runtime lower bound. Produces a gene×gene
influence matrix; we score TF->target as |A[tf, tg]|.
"""
from __future__ import annotations

from typing import Dict

import numpy as np

from src.baselines.common import DenseMatrixScorer, load_expression_df, universe_gene_names


def run_regdiffusion(rna_path: str, gene_index: Dict[str, int], pcfg: dict,
                     device: str = "cuda:0") -> DenseMatrixScorer:
    try:
        import regdiffusion as rd
    except ImportError as e:  # pragma: no cover
        raise ImportError("regdiffusion not installed — `pip install regdiffusion`") from e

    genes = universe_gene_names(gene_index)                    # order == gene_index
    expr = load_expression_df(rna_path, gene_index, pcfg)[genes].to_numpy(np.float32)
    # RegDiffusion rejects zero-variance genes; add tiny noise to constant columns
    # (keeps gene_index alignment so we can score any edge).
    zero_var = expr.std(axis=0) == 0
    if zero_var.any():
        rng = np.random.default_rng(0)
        expr[:, zero_var] += rng.normal(0, 1e-6, size=(expr.shape[0], int(zero_var.sum()))).astype(np.float32)
        print(f"[regdiffusion] denoised {int(zero_var.sum())} zero-variance genes", flush=True)
    print(f"[regdiffusion] {expr.shape[0]} cells × {expr.shape[1]} genes on {device}", flush=True)

    dev = device if str(device).startswith("cuda") else "cpu"
    trainer = _make_trainer(rd, expr, genes, dev)
    trainer.train()
    A = _extract_adjacency(trainer, genes)
    A = np.abs(np.asarray(A, dtype=np.float64))
    np.fill_diagonal(A, 0.0)
    return DenseMatrixScorer(A)


def _make_trainer(rd, expr, genes, device):
    # API has shifted across versions; try the documented signatures in turn.
    for kwargs in ({"gene_names": genes, "device": device}, {"device": device}, {}):
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
