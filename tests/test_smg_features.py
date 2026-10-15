"""Joint-histogram port == official compute_edge_attr (on the official demo, if cached)."""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.benchmarks.scmgrn_features import joint_histogram_pair_feats, svd_reduce


def _official(rna_genes_by_cells, names, i, j, bins=16):
    a = np.log10(rna_genes_by_cells.loc[names[i]].values + 1e-2)
    b = np.log10(rna_genes_by_cells.loc[names[j]].values + 1e-2)
    H = np.histogram2d(a, b, bins=bins)[0].T
    return (np.log10(H / len(a) + 1e-4) + 4) / 4


def test_hist_matches_numpy_official_formula():
    rng = np.random.RandomState(0)
    X = rng.poisson(1.5, size=(7, 300)).astype(float)
    X[3] = 2.0                                           # constant TF (degenerate range)
    names = [f"T{i}" for i in range(7)]
    df = pd.DataFrame(X, index=names)
    F = joint_histogram_pair_feats(X).astype(np.float64)
    iu = np.triu_indices(16)
    for i, j in [(0, 1), (2, 5), (1, 1), (3, 4), (6, 0)]:
        E = _official(df, names, i, j)
        S = (E + E.T) / 2
        np.testing.assert_allclose(F[i, j], S[iu], atol=2e-3)       # float16 storage
    np.testing.assert_allclose(F, F.transpose(1, 0, 2), atol=1e-6)  # symmetric


def test_svd_reduce_shape_and_determinism():
    M = np.random.RandomState(1).randn(40, 500)
    a, b = svd_reduce(M, 64), svd_reduce(M, 64)
    assert a.shape == (40, 39) and np.allclose(a, b)
