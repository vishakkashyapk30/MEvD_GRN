"""LINGER-protocol Cistrome evaluation (+ internal label-set test metrics).

`bm_trans_scores` is LINGER's `LingerGRN/Benchmk.py::bm_trans` (LingerGRN==1.110)
with plotting removed and the inferred GRN passed in memory. KEGNI's
benchmark_cistrome.ipynb uses the identical code. Verbatim semantics:

    data0 = pd.read_csv(Groundtruth, sep='\\t', skiprows=5, header=0)
    data1 = data0.groupby(['symbol'])['score'].max()
    label = data1.sort_values(axis=0, ascending=False).reset_index()
    label = label['symbol'].iloc[:1000]                 # N = 1000
    d1 = np.zeros(len(TGset)); d1[np.isin(TGset, label)] = 1
    auc   = roc_auc_score(d1, Score)
    auprr = average_precision_score(d1, Score) * len(d1) / sum(d1)

Candidate space = the TG rows handed in. For a like-for-like comparison every
method is scored over the SAME gene set (runbook amendment A3):
  linger_tg  LINGER re-run's cell-type-specific TG set (<root>/linger/tg_<ct>.txt)
  expressed  the benchmark gene universe (LINGER's filter_genes(min_cells=3))
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

N_TOP = 1000


# ----------------------------------------------------------------------------- ground truth
def chip_top_targets(path: str, n_top: int = N_TOP) -> List[str]:
    data0 = pd.read_csv(path, sep="\t", skiprows=5, header=0)
    data1 = data0.groupby(["symbol"])["score"].max()
    label = data1.sort_values(axis=0, ascending=False).reset_index()
    return [str(s) for s in label["symbol"].iloc[:n_top]]


def load_groundtruth(dcfg: dict, root: Path) -> dict:
    """{'datasets': [{id, tf, cell_type, file, top (list), in_linger_table}]}"""
    gt = dcfg["groundtruth"]
    gdir = Path(root) / gt["dir"]
    out = []
    for d in gt["datasets"]:
        f = gdir / f"{d['id']}_{d['tf']}_gene_score_5fold.txt"
        out.append({"id": str(d["id"]), "tf": d["tf"].upper(), "cell_type": d["cell_type"],
                    "file": str(f), "top": chip_top_targets(str(f)),
                    "in_linger_table": bool(d.get("in_linger_table", True))})
    return {"datasets": out}


# ----------------------------------------------------------------------------- LINGER metric
def bm_trans_scores(TGset: Sequence[str], Score: np.ndarray, label: Sequence[str]) -> dict:
    TGset = np.asarray(TGset)
    d1 = np.zeros(len(TGset))
    loc = np.where(np.isin(TGset, np.asarray(label)))[0]
    d1[loc] = 1
    if d1.sum() == 0 or d1.sum() == len(d1):
        return {"auc": float("nan"), "aupr_ratio": float("nan"), "n_pos": int(d1.sum()), "n": len(d1)}
    auc = roc_auc_score(d1, Score)
    aupr = average_precision_score(d1, Score)
    auprr = aupr * len(d1) / sum(d1)
    return {"auc": float(auc), "aupr_ratio": float(auprr), "aupr": float(aupr),
            "n_pos": int(d1.sum()), "n": int(len(d1))}


def candidate_genes(dcfg: dict, root: Path, cell_type: str, genes: Sequence[str],
                    space: str) -> Optional[List[str]]:
    if space == "expressed":
        return list(genes)
    if space == "linger_tg":
        f = Path(root) / "linger" / f"tg_{cell_type}.txt"
        if not f.exists():
            return None
        return [l.strip() for l in open(f) if l.strip()]
    raise ValueError(space)


def evaluate_chip(score_rows: Dict[str, np.ndarray], genes: Sequence[str], cell_type: str,
                  gt: dict, dcfg: dict, root: Optional[Path] = None,
                  spaces: Sequence[str] = ("expressed", "linger_tg")) -> dict:
    """score_rows: TF -> score vector over `genes`. Evaluates every Cistrome dataset
    whose cell type == cell_type. Genes of the candidate space missing from `genes`
    are given the minimum score (a method that cannot score a gene ranks it last)."""
    root = Path(root) if root is not None else Path(dcfg.get("_root", "."))
    gi = {g: i for i, g in enumerate(genes)}
    res = {}
    for space in spaces:
        cand = candidate_genes(dcfg, root, cell_type, genes, space)
        if cand is None:
            continue
        idx = np.array([gi.get(g, -1) for g in cand])
        rows = []
        for d in gt["datasets"]:
            if d["cell_type"] != cell_type or d["tf"] not in score_rows:
                continue
            v = np.asarray(score_rows[d["tf"]], dtype=np.float64)
            s = np.where(idx >= 0, v[np.clip(idx, 0, None)], np.nanmin(v) - 1.0)
            m = bm_trans_scores(cand, s, d["top"])
            m.update({"id": d["id"], "tf": d["tf"], "in_linger_table": d["in_linger_table"],
                      "missing_genes": int((idx < 0).sum())})
            rows.append(m)
        res[space] = rows
    return {"per_dataset": res, "summary": summarize(res.get("expressed", []))}


def summarize(rows: List[dict]) -> dict:
    r19 = [r for r in rows if r["in_linger_table"] and np.isfinite(r["auc"])]
    rall = [r for r in rows if np.isfinite(r["auc"])]
    f = lambda rs, k: float(np.mean([r[k] for r in rs])) if rs else float("nan")
    return {"n": len(r19), "auroc_mean": f(r19, "auc"), "aupr_ratio_mean": f(r19, "aupr_ratio"),
            "n_all": len(rall), "auroc_mean_all": f(rall, "auc"),
            "aupr_ratio_mean_all": f(rall, "aupr_ratio")}


# ----------------------------------------------------------------------------- internal test
def internal_test_metrics(S: np.ndarray, regs: np.ndarray, arrays: dict, n: int, regime: str) -> dict:
    """Held-out part of the TRAINING label set (not the headline).
    per_tf_auroc: for each test TF, rank all candidate targets (test positives = 1;
    that TF's train/val positives excluded from the ranking; target regime: only
    test targets); pooled uniform / degree-matched AUROC+AUPR on the stored sets."""
    row = {int(r): k for k, r in enumerate(regs)}
    te = arrays["test"]
    known = {}
    for p in ("train", "val"):
        for t, g in arrays[p]["pos"]:
            known.setdefault(int(t), set()).add(int(g))
    allowed = np.unique(te["pos"][:, 1]) if regime in ("target", "target_all") else np.arange(n)
    per_tf = []
    pos_by_tf = {}
    for t, g in te["pos"]:
        pos_by_tf.setdefault(int(t), set()).add(int(g))
    for t, pos in pos_by_tf.items():
        cand = np.array([g for g in allowed if g != t and g not in known.get(t, ())])
        y = np.isin(cand, list(pos)).astype(int)
        if 0 < y.sum() < len(y):
            per_tf.append(roc_auc_score(y, S[row[t], cand]))

    def pooled(neg_key):
        neg = te.get(neg_key)
        if neg is None or len(neg) == 0:
            return float("nan"), float("nan")
        pairs = np.concatenate([te["pos"], neg])
        y = np.r_[np.ones(len(te["pos"])), np.zeros(len(neg))]
        s = S[[row[int(t)] for t in pairs[:, 0]], pairs[:, 1]]
        return float(roc_auc_score(y, s)), float(average_precision_score(y, s))

    ua, up = pooled("neg")
    da, dp = pooled("neg_dm")
    return {"per_tf_auroc_mean": float(np.mean(per_tf)) if per_tf else float("nan"),
            "n_test_tfs": len(per_tf), "uniform_auroc": ua, "uniform_aupr": up,
            "dm_auroc": da, "dm_aupr": dp}
