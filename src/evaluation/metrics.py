"""Evaluation metrics (plan Part 8): AUROC, AUPR, Early Precision, EPR."""
from __future__ import annotations

from typing import Dict

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score


def compute_auroc(labels: np.ndarray, scores: np.ndarray) -> float:
    if len(np.unique(labels)) < 2:
        return float("nan")
    return float(roc_auc_score(labels, scores))


def compute_aupr(labels: np.ndarray, scores: np.ndarray) -> float:
    if len(np.unique(labels)) < 2:
        return float("nan")
    return float(average_precision_score(labels, scores))


def compute_early_precision(labels: np.ndarray, scores: np.ndarray, k: int | None = None) -> float:
    """Precision among the top-k scored items (k defaults to #positives)."""
    labels = np.asarray(labels)
    scores = np.asarray(scores)
    n_pos = int(labels.sum())
    if k is None:
        k = n_pos
    k = min(k, len(labels))
    if k <= 0:
        return float("nan")
    top = np.argsort(-scores)[:k]
    return float(labels[top].sum() / k)


def compute_epr(early_precision: float, n_positives: int, n_total: int) -> float:
    """Enrichment over random = EP / (n_pos / n_total)."""
    base = n_positives / n_total if n_total > 0 else 0.0
    if base == 0:
        return float("nan")
    return float(early_precision / base)


def compute_all_metrics(labels: np.ndarray, scores: np.ndarray,
                        n_total_candidates: int | None = None) -> Dict[str, float]:
    labels = np.asarray(labels)
    scores = np.asarray(scores)
    n_pos = int(labels.sum())
    n_total = n_total_candidates if n_total_candidates is not None else len(labels)
    ep = compute_early_precision(labels, scores)
    return {
        "auroc": compute_auroc(labels, scores),
        "aupr": compute_aupr(labels, scores),
        "early_precision": ep,
        "epr": compute_epr(ep, n_pos, n_total),
        "n_pos": n_pos,
        "n_total": int(len(labels)),
    }
