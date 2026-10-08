"""PBMC10k Multiome -> MeVD-GRN processed dirs, one per evaluated cell type.

Everything here is LABEL-FREE: node features, the co-expression kNN and the
TF-candidate graph are built from the cell type's own cells only, never from
any label set (the 02_preprocess.py label-dependent candidate graph is not
reproduced; see docs/experiments/scmultiomegrn_generalization.md s10).

Layout of a processed dir (read by src.data.dataset.load_celltype_data):
  rna_features_aligned.npy  atac_features_aligned.npy  rna_signature.npy
  openness.npy  gene_index.json  tf_indices.json  coexpr_edges.pt
  tf_candidate_edges.pt  motif_edges.pt (empty)  negative_pool.pt (empty;
  replaced per run by the degree-matched train pool)  evidence_label.pt (empty)
  + tf_candidate_edges_rnaonly.pt (candidate graph without the accessibility
  gate, for the RNA-only ablation)  + summary.json
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch

from src.data import graph_builder as gb
from src.data.preprocessing import (_col_detection_rate, _col_mean_var, _load_tss, _normalize_log,
                                    _parse_peaks, _peak_gene_weighted_incidence,
                                    _regulatory_potential_descriptor, canon,
                                    coexpression_signatures)
from src.utils.io import save_json


def load_multiome(h5: str):
    """(rna cells x genes CSR, gene symbols UPPER, atac cells x peaks CSR, peak names, barcodes)."""
    import h5py
    with h5py.File(h5, "r") as f:        # single read of the CSC matrix (127M nnz)
        m = f["matrix"]
        nf, nb = int(m["shape"][0]), int(m["shape"][1])
        M = sp.csc_matrix((m["data"][:].astype(np.float32), m["indices"][:].astype(np.int32),
                           m["indptr"][:]), shape=(nf, nb))
        ft = m["features"]["feature_type"][:]
        names = [x.decode() for x in m["features"]["name"][:]]
        inter = [x.decode() for x in m["features"]["interval"][:]]
        bc = [b.decode() for b in m["barcodes"][:]]
    M = M.T.tocsr()                       # barcodes x features
    g_rows = np.where(ft == b"Gene Expression")[0]
    p_rows = np.where(ft == b"Peaks")[0]
    rna = M[:, g_rows].tocsr()
    atac = M[:, p_rows].tocsr()
    del M
    return rna, [canon(names[i]) for i in g_rows], atac, [inter[i] for i in p_rows], bc


def dedup_columns(X: sp.csr_matrix, names: List[str]) -> Tuple[sp.csr_matrix, List[str]]:
    """Collapse duplicate symbols by summing counts (10x has a handful)."""
    s = pd.Series(range(len(names)), index=names)
    if s.index.is_unique:
        return X, names
    uniq = list(dict.fromkeys(names))
    pos = {g: i for i, g in enumerate(uniq)}
    M = sp.csr_matrix((np.ones(len(names)), (np.arange(len(names)), [pos[g] for g in names])),
                      shape=(len(names), len(uniq)))
    return (X @ M).tocsr(), uniq


def gene_universe(rna: sp.csr_matrix, genes: List[str], min_cells: int,
                  restrict_to: Optional[Sequence[str]] = None,
                  always: Sequence[str] = ()) -> List[str]:
    """Genes detected in >= min_cells cells (all cells), optionally intersected with
    a protocol gene list, plus `always` (evaluation TFs) if present in the matrix."""
    det = np.asarray((rna > 0).sum(axis=0)).ravel()
    keep = [g for g, d in zip(genes, det) if d >= min_cells]
    if restrict_to is not None:
        r = {canon(g) for g in restrict_to}
        keep = [g for g in keep if g in r]
    have = set(genes)
    extra = [canon(t) for t in always if canon(t) in have and canon(t) not in set(keep)]
    return sorted(set(keep) | set(extra))


def build_processed(out_dir: Path, rna: sp.csr_matrix, genes: List[str], atac: sp.csr_matrix,
                    peaks: List[str], cell_mask: np.ndarray, universe: List[str],
                    tf_names: Sequence[str], gtf_path: str, dcfg: dict, acfg: dict,
                    cell_type: str) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    gene_index = {g: i for i, g in enumerate(universe)}
    n = len(universe)
    col = {g: i for i, g in enumerate(genes)}
    sel = np.array([col[g] for g in universe])

    # --- RNA (this cell type's cells only)
    Xr = _normalize_log(rna[cell_mask][:, sel], dcfg)
    mean, var = _col_mean_var(Xr)
    det = _col_detection_rate(Xr)
    rna_feats = np.stack([mean, var, det], axis=1).astype(np.float32)
    sig = coexpression_signatures(Xr, int(dcfg.get("signature_dim", 50))).astype(np.float32)
    if sig.shape[0] != n:
        raise RuntimeError("signature / universe size mismatch")

    # --- ATAC: RP-weighted gene activity + locus descriptor
    Xa = _normalize_log(atac[cell_mask], dcfg)
    pk = _parse_peaks(peaks)
    tss = _load_tss("hg38", universe, gtf_path)
    inc, per_gene = _peak_gene_weighted_incidence(pk, tss, gene_index, len(peaks),
                                                  int(dcfg.get("atac_window_bp", 100_000)),
                                                  float(acfg.get("rp_decay_bp", 10_000)))
    act = Xa @ inc
    am, av = _col_mean_var(act)
    ad = _col_detection_rate(act)
    atac_feats = np.stack([am, av, ad], axis=1).astype(np.float32)
    openness = _regulatory_potential_descriptor(per_gene, n, int(acfg.get("rp_top_k", 10)))

    # --- label-free graphs
    tf_idx = sorted(gene_index[canon(t)] for t in set(tf_names) if canon(t) in gene_index)
    coexpr = gb.build_coexpression_graph(sig, k=int(dcfg["coexpr_knn_k"]))
    cand = gb.build_prior_graph(gene_index, tf_idx, atac_feats, rna_feats,
                                top_k=int(dcfg["tf_candidate_topk"]), signatures=sig,
                                openness=openness, evidence=None, exclude_positives=False)
    cand_rna = gb.build_prior_graph(gene_index, tf_idx, np.zeros_like(atac_feats), rna_feats,
                                    top_k=int(dcfg["tf_candidate_topk"]), signatures=sig,
                                    openness=None, evidence=None, exclude_positives=False)

    np.save(out_dir / "rna_features_aligned.npy", rna_feats)
    np.save(out_dir / "atac_features_aligned.npy", atac_feats)
    np.save(out_dir / "rna_signature.npy", sig)
    np.save(out_dir / "openness.npy", openness)
    save_json(gene_index, out_dir / "gene_index.json")
    save_json(tf_idx, out_dir / "tf_indices.json")
    torch.save(coexpr, out_dir / "coexpr_edges.pt")
    torch.save(cand, out_dir / "tf_candidate_edges.pt")
    torch.save(cand_rna, out_dir / "tf_candidate_edges_rnaonly.pt")
    empty = torch.zeros((2, 0), dtype=torch.long)
    torch.save(empty, out_dir / "motif_edges.pt")
    torch.save(empty, out_dir / "negative_pool.pt")
    torch.save(empty, out_dir / "evidence_label.pt")
    # dense log-normalised expression of the regulators is NOT stored; baselines
    # recompute from the h5 (keeps processed dirs small)
    summ = {"cell_type": cell_type, "n_cells": int(cell_mask.sum()), "n_genes": n,
            "n_tfs_in_candidate_graph": len(tf_idx), "coexpr_edges": int(coexpr.shape[1]),
            "tf_candidate_edges": int(cand.shape[1]), "tf_candidate_edges_rnaonly": int(cand_rna.shape[1]),
            "genes_with_tss": int(len(tss)), "atac_nonzero_frac": float((atac_feats[:, 0] > 0).mean()),
            "proximal_peak_frac": float((openness[:, 1] > 0.1).mean()),
            "label_free": True, "prior_exclude_positives": False}
    save_json(summ, out_dir / "summary.json")
    return summ
