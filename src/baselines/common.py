"""Shared baseline infrastructure.

Every baseline must (a) consume the SAME RNA input and gene universe as MEvD-GRN
and (b) expose a uniform per-edge scorer, so all methods are evaluated on the
identical test splits with the identical metric code (fair comparison).
"""
from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from src.data.preprocessing import _normalize_log, canon, extract_symbol, load_features_by_cells


def load_expression_df(rna_path: str, gene_index: Dict[str, int], pcfg: dict) -> pd.DataFrame:
    """Return a DENSE (cells × universe-genes) log-normalized DataFrame.

    Columns are gene names in gene_index order (missing genes -> all-zero column).
    Uses the same CP10K+log1p transform MEvD-GRN derives its features from.
    """
    X, raw_names = load_features_by_cells(rna_path)            # cells × genes (CSR)
    names = [extract_symbol(g) for g in raw_names]
    # light cell QC (match preprocessing)
    min_genes = int(pcfg.get("min_genes_per_cell", 200))
    keep = np.asarray((X > 0).sum(axis=1)).ravel() >= min_genes
    if keep.sum() == 0:
        keep[:] = True
    X = _normalize_log(X[keep], pcfg).tocsc()

    col_of: Dict[str, int] = {}
    for j, g in enumerate(names):                              # first occurrence wins
        col_of.setdefault(g, j)
    universe = sorted(gene_index, key=lambda g: gene_index[g])
    n_cells = X.shape[0]
    mat = np.zeros((n_cells, len(universe)), dtype=np.float32)
    for k, g in enumerate(universe):
        j = col_of.get(g)
        if j is not None:
            mat[:, k] = X[:, j].toarray().ravel()
    return pd.DataFrame(mat, columns=universe)


class DictScorer:
    """Sparse per-edge scorer (default 0 for unseen edges); e.g. GRNBoost2."""

    def __init__(self, lookup: Dict[Tuple[int, int], float]):
        self.lookup = lookup

    def score(self, tf_idx: np.ndarray, tg_idx: np.ndarray) -> np.ndarray:
        return np.array([self.lookup.get((int(a), int(b)), 0.0)
                         for a, b in zip(tf_idx, tg_idx)], dtype=np.float64)


class DenseMatrixScorer:
    """Dense gene×gene score matrix scorer; e.g. RegDiffusion / GMF."""

    def __init__(self, matrix: np.ndarray):
        self.M = matrix                                        # (n_genes, n_genes)

    def score(self, tf_idx: np.ndarray, tg_idx: np.ndarray) -> np.ndarray:
        return np.asarray(self.M[np.asarray(tf_idx), np.asarray(tg_idx)], dtype=np.float64)


class EmbeddingScorer:
    """Score(tf, tg) = <Z[tf], Z[tg]>; for embedding/MF methods (e.g. GMF-GAE)."""

    def __init__(self, Z: np.ndarray):
        self.Z = np.asarray(Z, dtype=np.float64)

    def score(self, tf_idx: np.ndarray, tg_idx: np.ndarray) -> np.ndarray:
        a = self.Z[np.asarray(tf_idx)]
        b = self.Z[np.asarray(tg_idx)]
        return (a * b).sum(axis=1)


def edgelist_to_dict(df: pd.DataFrame, gene_index: Dict[str, int],
                     tf_col: str = "TF", tg_col: str = "target",
                     w_col: str = "importance") -> Dict[Tuple[int, int], float]:
    """Map a (TF, target, importance) edgelist to {(tf_idx, tg_idx): weight}."""
    out: Dict[Tuple[int, int], float] = {}
    for tf, tg, w in zip(df[tf_col], df[tg_col], df[w_col]):
        i, j = gene_index.get(canon(tf)), gene_index.get(canon(tg))
        if i is not None and j is not None:
            out[(i, j)] = float(w)
    return out


def universe_gene_names(gene_index: Dict[str, int]) -> List[str]:
    return sorted(gene_index, key=lambda g: gene_index[g])
