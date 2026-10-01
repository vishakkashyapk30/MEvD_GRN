"""Baselines for the PBMC10k / LINGER-Cistrome benchmark.

Every baseline returns a dense score matrix S (n_regulators x n_genes) for a
list of regulator indices, over the SAME gene universe as MeVD-GRN, so the
same evaluation code (src/benchmarks/pbmc_eval.py) scores everything.

Trivial baselines (mandatory, docs/experiments/pbmc10k_linger_benchmark.md s4.7):
  degree   score(t, g) = in-degree of g among the TRAIN positives (TF-agnostic)
  geneid   InfoSEM-style logistic regression on one-hot(TF) (+) one-hot(target),
           fit on the same train positives + degree-matched negatives
  pearson  |r| (and signed r) between TF and target log-expression in the cell type
  grnboost2 arboreto SGBM per target (seed 0), regulators = label TFs + eval TFs
"""
from __future__ import annotations

from typing import Dict, List, Sequence

import numpy as np
import pandas as pd
import scipy.sparse as sp


def degree_scores(train_pos: np.ndarray, n_genes: int, regs: Sequence[int]) -> np.ndarray:
    indeg = np.bincount(train_pos[:, 1], minlength=n_genes).astype(np.float32)
    return np.tile(indeg, (len(regs), 1))


def geneid_scores(train_pos: np.ndarray, train_neg: np.ndarray, n_genes: int,
                  regs: Sequence[int], C: float = 1.0, seed: int = 0) -> np.ndarray:
    """One-hot(TF) || one-hot(target) -> logistic regression. Unseen TFs get only
    the intercept + target weight (their TF column never fires in training)."""
    from sklearn.linear_model import LogisticRegression
    pairs = np.concatenate([train_pos, train_neg])
    y = np.concatenate([np.ones(len(train_pos)), np.zeros(len(train_neg))])
    m = len(pairs)
    rows = np.repeat(np.arange(m), 2)
    cols = np.stack([pairs[:, 0], n_genes + pairs[:, 1]], axis=1).ravel()
    X = sp.csr_matrix((np.ones(2 * m, dtype=np.float32), (rows, cols)), shape=(m, 2 * n_genes))
    lr = LogisticRegression(C=C, max_iter=2000, solver="liblinear", random_state=seed)
    lr.fit(X, y)
    w = lr.coef_.ravel()
    w_tf, w_tg = w[:n_genes], w[n_genes:]
    return (w_tf[np.asarray(regs)][:, None] + w_tg[None, :] + lr.intercept_[0]).astype(np.float32)


def pearson_scores(X_log: sp.csr_matrix, regs: Sequence[int], signed: bool = False) -> np.ndarray:
    """X_log: cells x genes (log-normalised, the universe order). Returns r or |r|."""
    X = X_log.tocsc().astype(np.float64)
    C = X.shape[0]
    mean = np.asarray(X.mean(axis=0)).ravel()
    sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    std = np.sqrt(np.maximum(sq - mean ** 2, 0.0))
    std[std == 0] = np.inf                                   # constant genes -> r = 0
    R = X[:, list(regs)].toarray()                            # cells x regs
    num = (X.T @ R).T / C - np.outer(mean[list(regs)], mean)  # cov (regs x genes)
    r = num / np.outer(std[list(regs)], std)
    r = np.nan_to_num(r, nan=0.0, posinf=0.0, neginf=0.0)
    return (r if signed else np.abs(r)).astype(np.float32)


def grnboost2_scores(X_log: sp.csr_matrix, genes: List[str], regulators: Sequence[str],
                     eval_regs: Sequence[str], n_jobs: int = 8, seed: int = 0,
                     targets: Sequence[str] | None = None) -> np.ndarray:
    """GRNBoost2 with the full regulator list; returns rows for `eval_regs` only."""
    from src.benchmarks.scmgrn_features import grnboost2_features
    gi = {g: i for i, g in enumerate(genes)}
    use = sorted(set(targets) | set(regulators)) if targets is not None else genes
    cols = [gi[g] for g in use]
    df = pd.DataFrame(X_log[:, cols].toarray().T, index=use)
    net = grnboost2_features(df, [r for r in regulators if r in gi], seed=seed, n_jobs=n_jobs)
    S = np.zeros((len(eval_regs), len(genes)), dtype=np.float32)
    for k, r in enumerate(eval_regs):
        if r in net.index:
            row = net.loc[r]
            S[k, [gi[g] for g in row.index]] = row.to_numpy(dtype=np.float32)
    return S
