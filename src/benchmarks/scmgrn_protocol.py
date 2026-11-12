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


def degree_matched_test(adj: np.ndarray, fold: dict, seed: int = 0, n_near: int = 3) -> Dict[str, np.ndarray]:
    """Second, degree-controlled test set ("Hidden sampling biases inflate performance in GRN inference", 2025).

    Positives = the fold's test positives (unordered, u < v). Negatives = true non-links (adj == 0) that are not
    in the fold's train/val files, one per positive: drawn at random among the n_near candidate pairs closest to
    the positive in (lower, higher) log1p TRAINING degree of the two TFs. A pure degree rule is then near chance.
    Deterministic for a given (fold, seed), so every method is scored on the same pairs. Returns {"u", "v", "label"}.
    """
    n = adj.shape[0]
    tr = fold["train"]["pos"]
    ld = np.log1p(np.bincount(tr[:, 0], minlength=n).astype(float))
    seen = np.zeros((n, n), dtype=bool)
    for part in ("train", "val"):
        for k in ("pos", "neg"):
            e = fold[part][k]
            seen[e[:, 0], e[:, 1]] = True
    tp = fold["test"]["pos"]
    tp = tp[tp[:, 0] < tp[:, 1]]
    iu, ju = np.triu_indices(n, 1)
    ok = (adj[iu, ju] == 0) & ~seen[iu, ju] & ~seen[ju, iu]
    cu, cv = iu[ok], ju[ok]
    c_lo, c_hi = np.minimum(ld[cu], ld[cv]), np.maximum(ld[cu], ld[cv])
    rng = np.random.default_rng(seed)
    neg = []
    for u, v in tp:
        d = np.abs(c_lo - min(ld[u], ld[v])) + np.abs(c_hi - max(ld[u], ld[v]))
        near = np.argpartition(d, min(n_near, len(d) - 1))[:n_near]
        j = near[rng.integers(len(near))]
        neg.append((cu[j], cv[j]))
    neg = np.array(neg, dtype=int).reshape(-1, 2)
    return {"u": np.concatenate([tp[:, 0], neg[:, 0]]), "v": np.concatenate([tp[:, 1], neg[:, 1]]),
            "label": np.concatenate([np.ones(len(tp), int), np.zeros(len(neg), int)])}


def degree_caliper_test(adj: np.ndarray, fold: dict, seed: int = 0, caliper: float = 0.1) -> Dict[str, np.ndarray]:
    """Official split, strict version of degree_matched_test: a test link (u, v) is kept only if some test-pool
    non-link (adj == 0, in no train/val file) has BOTH sorted endpoint log1p training degrees within `caliper` of
    the link's; one such non-link is drawn at random. Unmatched links are dropped, so the degree rule is at chance."""
    n = adj.shape[0]
    tr = fold["train"]["pos"]
    ld = np.log1p(np.bincount(tr[:, 0], minlength=n).astype(float))
    seen = np.zeros((n, n), dtype=bool)
    for part in ("train", "val"):
        for k in ("pos", "neg"):
            e = fold[part][k]
            seen[e[:, 0], e[:, 1]] = True
    tp = fold["test"]["pos"]
    tp = tp[tp[:, 0] < tp[:, 1]]
    iu, ju = np.triu_indices(n, 1)
    ok = (adj[iu, ju] == 0) & ~seen[iu, ju] & ~seen[ju, iu]
    cu, cv = iu[ok], ju[ok]
    c_lo, c_hi = np.minimum(ld[cu], ld[cv]), np.maximum(ld[cu], ld[cv])
    rng = np.random.default_rng(seed)
    us, vs, labs = [], [], []
    for u, v in tp:
        lo, hi = min(ld[u], ld[v]), max(ld[u], ld[v])
        cand = np.where((np.abs(c_lo - lo) <= caliper) & (np.abs(c_hi - hi) <= caliper))[0]
        if len(cand) == 0:
            continue
        j = cand[rng.integers(len(cand))]
        us += [u, cu[j]]; vs += [v, cv[j]]; labs += [1, 0]
    return {"u": np.array(us, int), "v": np.array(vs, int), "label": np.array(labs, int)}


def evaluate_degree_matched(S: np.ndarray, dm: Dict[str, np.ndarray]) -> Dict[str, float]:
    S = np.asarray(S, dtype=np.float64)
    return metric_fn((S[dm["u"], dm["v"]] + S[dm["v"], dm["u"]]) / 2.0, dm["label"])


# ---------------------------------------------------------------------------
# Held-out-TF split (ours, not scMultiomeGRN's): no model can count a hidden TF's links.
# Pair classes follow Park & Marcotte (Nat Methods 2012): C2 = one TF of the pair hidden, C3 = both hidden.
# ---------------------------------------------------------------------------
SPLIT_KINDS = ("official", "tfholdout")
TFHOLDOUT_TAG = "tfholdout_seed-{seed}_fold-{n}"


def resolve_split_dir(ct_split_root: str, split_cfg: dict, kind: str = "official") -> str:
    """Directory of one cell type's fold files for the given split kind."""
    if kind == "official":
        return split_dir(ct_split_root, int(split_cfg["seed"]), int(split_cfg["n_splits"]),
                         float(split_cfg["train_split_rate"]), bool(split_cfg["rest_all_train"]))
    if kind == "tfholdout":
        return os.path.join(ct_split_root, TFHOLDOUT_TAG.format(seed=int(split_cfg["seed"]),
                                                                n=int(split_cfg["n_splits"])))
    raise ValueError(f"unknown split kind {kind!r}")


def results_subdir(results_dir: str, kind: str = "official") -> str:
    """Results of non-official splits live next to the official ones: results/<ds>_<kind>/..."""
    return results_dir if kind == "official" else f"{results_dir}_{kind}"


def make_tf_holdout_splits(adj_file: str, save_root: str, seed: int = 666, n_splits: int = 10) -> str:
    """TFs are shuffled (seed) into n_splits groups. Fold k: test TFs = group k, validation TFs = group k+1,
    training TFs = the rest.
      train: every pair of two training TFs (links label 1, non-links label 0) plus the diagonal, as the official
             split also lists diagonal pairs among training negatives (the official loader relabels them).
      val:   all links with >= 1 validation TF and no test TF, plus as many random non-links of the same pairs.
      test:  all links with >= 1 test TF and no validation TF (C2: test-train, C3: test-test), plus, per class,
             as many random non-links. Pairs of a test TF with a validation TF are used nowhere.
    Written in the official file format (row, col, label; split-XXX_<part>-edge-<npos>_<N>xN.csv), so the official
    scMultiomeGRN code runs on it unchanged. A sidecar split-XXX_tfs.json lists the test / validation TFs.
    """
    import json
    save_dir = os.path.join(save_root, TFHOLDOUT_TAG.format(seed=seed, n=n_splits))
    os.makedirs(save_dir, exist_ok=True)
    adj, _ = load_adj(adj_file)
    n = adj.shape[0]
    rng = np.random.default_rng(seed)
    groups = np.array_split(rng.permutation(n), n_splits)
    iu, ju = np.triu_indices(n, 1)
    for k in range(n_splits):
        test_tfs, val_tfs = groups[k], groups[(k + 1) % n_splits]
        role = np.zeros(n, dtype=int)                    # 0 train, 1 val, 2 test
        role[val_tfs], role[test_tfs] = 1, 2
        ru, rv = role[iu], role[ju]
        lab = adj[iu, ju]
        parts = {}
        tr = (ru == 0) & (rv == 0)
        diag = np.arange(n)
        parts["train"] = (np.concatenate([np.stack([iu[tr], ju[tr]], 1), np.stack([diag, diag], 1)]),
                          np.concatenate([lab[tr], np.zeros(n, int)]))

        def balanced(mask_list):
            rows, labs = [], []
            for m in mask_list:
                pos = np.where(m & (lab == 1))[0]
                neg = np.where(m & (lab == 0))[0]
                neg = rng.choice(neg, size=min(len(pos), len(neg)), replace=False) if len(neg) else neg
                idx = np.concatenate([pos, neg])
                rows.append(np.stack([iu[idx], ju[idx]], 1))
                labs.append(lab[idx])
            return np.concatenate(rows), np.concatenate(labs)

        has_val = (ru == 1) | (rv == 1)
        has_test = (ru == 2) | (rv == 2)
        parts["val"] = balanced([has_val & ~has_test])
        c2 = has_test & ~has_val & ((ru == 0) | (rv == 0))
        c3 = (ru == 2) & (rv == 2)
        parts["test"] = balanced([c2, c3])
        for name, (e, l) in parts.items():
            pd.DataFrame({"row": e[:, 0], "col": e[:, 1], "label": l.astype(int)}).to_csv(
                os.path.join(save_dir, f"split-{k+1:03d}_{name}-edge-{int(l.sum())}_{n}x{n}.csv"), index=False)
        with open(os.path.join(save_dir, f"split-{k+1:03d}_tfs.json"), "w") as fh:
            json.dump({"test_tfs": sorted(int(x) for x in test_tfs), "val_tfs": sorted(int(x) for x in val_tfs)}, fh)
    return save_dir


def _tfholdout_meta(sdir: str, fold: int):
    import json
    p = os.path.join(sdir, f"split-{fold:03d}_tfs.json")
    return json.load(open(p)) if os.path.exists(p) else None


def tfholdout_degree_matched_test(adj: np.ndarray, fold: dict, test_tfs, val_tfs, seed: int = 0,
                                  caliper: float = 0.1):
    """Hub-controlled C2 test with caliper matching. For each C2 link (hidden t, training TF r), the negative is a
    non-link (t, w) of the SAME hidden TF whose training-TF partner w has |log1p deg(w) - log1p deg(r)| <= caliper
    (drawn at random among those). Links with no such partner are dropped, so the degree of both endpoints is
    matched on every kept pair. The hidden TF has no training links and C2 pairs are in no training file, so the
    pool is every C2 non-link of that TF."""
    n = adj.shape[0]
    tr = fold["train"]["pos"]
    ld = np.log1p(np.bincount(tr[:, 0], minlength=n).astype(float))
    hidden = np.zeros(n, bool)
    hidden[list(test_tfs)] = True
    val = np.zeros(n, bool)
    val[list(val_tfs)] = True
    partners = np.where(~hidden & ~val)[0]
    rng = np.random.default_rng(seed)
    us, vs, labs = [], [], []
    for t in np.where(hidden)[0]:
        links = partners[adj[t, partners] == 1]
        nonlinks = partners[adj[t, partners] == 0]
        for r in links:
            ok = nonlinks[np.abs(ld[nonlinks] - ld[r]) <= caliper]
            if len(ok) == 0:
                continue
            w = ok[rng.integers(len(ok))]
            us += [t, t]; vs += [r, w]; labs += [1, 0]
    return {"u": np.array(us, int), "v": np.array(vs, int), "label": np.array(labs, int)}


def evaluate_all(S: np.ndarray, adj: np.ndarray, f: dict, sdir: str, fold: int) -> Dict[str, dict]:
    """Every metric set for one fold, for either split kind.
      official : test, all, all_offdiag, test_degmatched (degree_matched_test, nearest-3),
                 test_degcaliper (degree_caliper_test: unmatched links dropped, degree rule at chance)
      tfholdout: test (C2 + C3), all, all_offdiag, test_c2, test_c3,
                 test_degmatched = test_degcaliper = caliper-matched C2 (same hidden TF, matched partner degree)
    """
    S = np.asarray(S, dtype=np.float64)
    S = np.nan_to_num((S + S.T) / 2.0, nan=0.0, posinf=0.0, neginf=0.0)
    m = evaluate_score_matrix(S, adj, f["test_eval"])
    meta = _tfholdout_meta(sdir, fold)
    if meta is None:
        m["test_degmatched"] = evaluate_degree_matched(S, degree_matched_test(adj, f, seed=fold))
        dc = degree_caliper_test(adj, f, seed=fold)
        if 0 < dc["label"].sum() < len(dc["label"]):
            m["test_degcaliper"] = evaluate_degree_matched(S, dc)
        return m
    hidden = np.zeros(adj.shape[0], bool)
    hidden[meta["test_tfs"]] = True
    ei, lab = f["test_eval"]["edge_index"], f["test_eval"]["label"]
    keep = ei[0] < ei[1]
    u, v, y = ei[0][keep], ei[1][keep], lab[keep]
    both = hidden[u] & hidden[v]
    for name, sel in (("test_c2", ~both), ("test_c3", both)):
        if sel.sum() and 0 < y[sel].sum() < sel.sum():
            m[name] = metric_fn(S[u[sel], v[sel]], y[sel])
    dm = tfholdout_degree_matched_test(adj, f, meta["test_tfs"], meta["val_tfs"], seed=fold)
    if len(dm["label"]):
        m["test_degmatched"] = evaluate_degree_matched(S, dm)
        m["test_degcaliper"] = m["test_degmatched"]
    return m
