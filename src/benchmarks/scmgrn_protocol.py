"""scMultiomeGRN evaluation protocol: exact split replica, fold loading, metrics.

Split = re-implementation of `GraphDataModule.split_dataset` (official
`src/dataset.py`; defaults dataset_seed=666, n_splits=10, train_split_rate=0.8,
rest_all_train=True). Lightning's `seed_everything` is replaced by the three
seeds it sets (random / numpy / torch); only numpy's global RNG and
sklearn's KFold(random_state) are consumed, so the files are byte-identical
(verified against the official function by scripts/16 --verify_official).

Metrics = re-implementation of official `src/utils.py: metric_fn`:
  auroc = auc(roc_curve), aupr = trapezoidal auc(recall, precision),
  ap = average_precision_score, acc/precision/recall at threshold = median of
  the scores (`pred >= thr`) unless a threshold is passed in.
Two evaluation sets, as logged by the official Trainer:
  "test": the held-out fold (pos + equal number of random non-edges), each
          unordered pair in both directions
  "all" : every N x N entry (incl. diagonal, label 0 there), label = full
          ground-truth adjacency, threshold = the test-set median
The paper's reported AUROC/AUPR correspond to "all" (see the experiment doc).
"""
from __future__ import annotations

import copy
import os
import random
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import torch
from sklearn import metrics as skm
from sklearn.model_selection import KFold

SPLIT_TAG = "seed-{seed}_fold-{n}_rate-{rate}_all-train-{rest}"


def load_adj(adj_file: str) -> Tuple[np.ndarray, list]:
    data = pd.read_csv(adj_file, sep="\t", index_col=0)
    adj = data.values.astype(int)
    adj = adj - np.diag(np.diag(adj))
    return adj, [str(x) for x in data.index]


def split_dir(root: str, seed: int = 666, n_splits: int = 10, rate: float = 0.8,
              rest_all_train: bool = True) -> str:
    return os.path.join(root, SPLIT_TAG.format(seed=seed, n=n_splits, rate=rate, rest=rest_all_train))


def make_splits(adj_file: str, save_root: str, seed: int = 666, rate: float = 0.8,
                n_splits: int = 10, rest_all_train: bool = True) -> str:
    """Write split-XXX_{train,val,test}-edge-<npos>_<N>xN.csv exactly as the
    official split_dataset does. Returns the split directory."""
    save_dir = split_dir(save_root, seed, n_splits, rate, rest_all_train)
    os.makedirs(save_dir, exist_ok=True)
    adj, _ = load_adj(adj_file)
    edge_pos = np.vstack(np.nonzero(np.triu(adj, k=1))).T
    neg_mask = np.triu(np.ones_like(adj) ^ adj, k=1)
    edge_neg = np.vstack(np.nonzero(neg_mask)).T
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)   # = pl.seed_everything
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for i, (trval_idx, test_pos_idx) in enumerate(kf.split(edge_pos)):
        train_num = int(len(trval_idx) * rate)
        np.random.shuffle(trval_idx)
        train_pos_idx, val_pos_idx = trval_idx[:train_num], trval_idx[train_num:]
        neg_perm = np.random.permutation(len(edge_neg))
        test_neg_idx = neg_perm[:len(test_pos_idx)]
        m = copy.deepcopy(neg_mask)
        m[edge_neg[test_neg_idx][:, 0], edge_neg[test_neg_idx][:, 1]] = 0
        np.fill_diagonal(m, 1)
        trval_neg = np.vstack(np.nonzero(m)).T
        perm = np.random.permutation(len(trval_neg))
        n_val = len(trval_idx) - train_num
        val_neg_idx = perm[:n_val]
        train_neg_idx = perm[n_val:] if rest_all_train else perm[n_val:len(trval_idx)]
        parts = {
            "train": (edge_pos[train_pos_idx], trval_neg[train_neg_idx]),
            "val": (edge_pos[val_pos_idx], trval_neg[val_neg_idx]),
            "test": (edge_pos[test_pos_idx], edge_neg[test_neg_idx]),
        }
        for name, (p, n) in parts.items():
            e = np.concatenate([p, n])
            lab = np.concatenate([np.ones(len(p)), np.zeros(len(n))]).astype(int)
            pd.DataFrame({"row": e[:, 0], "col": e[:, 1], "label": lab}).to_csv(
                os.path.join(save_dir, f"split-{i+1:03d}_{name}-edge-{lab.sum()}_"
                                       f"{adj.shape[0]}x{adj.shape[1]}.csv"), index=False)
    return save_dir


def fold_files(sdir: str, fold: int) -> Dict[str, str]:
    pat = f"split-{fold:03d}_"
    files = [f for f in os.listdir(sdir) if f.startswith(pat)]
    return {k: os.path.join(sdir, next(f for f in files if f"_{k}-edge" in f))
            for k in ("train", "val", "test")}


def _both_dirs(e: np.ndarray) -> np.ndarray:
    return np.concatenate([e, e[:, ::-1]], axis=0)


def load_fold(sdir: str, fold: int, n_nodes: int) -> Dict[str, Dict[str, np.ndarray]]:
    """Directed (both-direction) pos/neg pairs per part; diagonal pairs (which
    the official loader relabels as positives, GAE's A+I convention) are dropped
    from train/val -- they are not TF-TF interactions and never appear in test.
    Also returns the official test evaluation set: unique directed pairs of the
    test file with labels from the test adjacency (as GraphDataModule.load_data)."""
    out: Dict[str, Dict[str, np.ndarray]] = {}
    for part, path in fold_files(sdir, fold).items():
        df = pd.read_csv(path)
        e = df[["row", "col"]].to_numpy()
        lab = df["label"].to_numpy()
        off = e[:, 0] != e[:, 1]
        out[part] = {"pos": _both_dirs(e[(lab == 1) & off]), "neg": _both_dirs(e[(lab == 0) & off])}
        if part == "test":
            A = np.zeros((n_nodes, n_nodes), dtype=int)
            A[e[:, 0], e[:, 1]] = lab
            A = A + A.T
            np.fill_diagonal(A, 1)
            ei = np.unique(np.concatenate([e, e[:, ::-1]]).T, axis=1)
            out["test_eval"] = {"edge_index": ei, "label": A[ei[0], ei[1]]}
    return out


def metric_fn(predict: np.ndarray, label: np.ndarray, threshold: float | None = None) -> Dict[str, float]:
    predict = np.asarray(predict, dtype=np.float64)
    label = np.asarray(label).astype(int)
    fpr, tpr, _ = skm.roc_curve(y_true=label, y_score=predict)
    prec, rec, _ = skm.precision_recall_curve(label, predict)
    if threshold is None:
        threshold = float(np.nanmedian(predict))
    pb = predict >= threshold
    return {
        "auroc": float(skm.auc(fpr, tpr)),
        "aupr": float(skm.auc(rec, prec)),
        "ap": float(skm.average_precision_score(label, predict)),
        "acc": float(((pb.astype(int) + label) % 2 == 0).mean()),
        "precision": float(skm.precision_score(label, pb, zero_division=0)),
        "recall": float(skm.recall_score(label, pb, zero_division=0)),
        "threshold": float(threshold),
        "n": int(len(label)), "n_pos": int(label.sum()),
    }


def evaluate_score_matrix(S: np.ndarray, adj: np.ndarray, test_eval: dict) -> Dict[str, dict]:
    """S: (N, N) symmetric score matrix. Returns {"test": ..., "all": ...,
    "all_offdiag": ...} computed exactly like the official best_test/best_all."""
    ei, lab = test_eval["edge_index"], test_eval["label"]
    test = metric_fn(S[ei[0], ei[1]], lab)
    all_ = metric_fn(S.reshape(-1), adj.reshape(-1), threshold=test["threshold"])
    off = ~np.eye(adj.shape[0], dtype=bool)
    all_off = metric_fn(S[off], adj[off], threshold=test["threshold"])
    return {"test": test, "all": all_, "all_offdiag": all_off}
