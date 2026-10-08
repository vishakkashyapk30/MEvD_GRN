"""Python port of BEAR-GRN's scoring code (UzunLab/BEAR-GRN v1.0.0, R package
`BEARGRN`, Zenodo 10.5281/zenodo.22015920), so any GRN edge list can be scored
exactly as the paper scores its methods.

Ported functions (R file -> here):
  R/benchmark_new_method_filtered.R      benchmark_new_method()        -> score_roc_pr()
  R/benchmark_new_method_early_metrics.R benchmark_new_method_early_metrics() -> score_early()
  R/benchmark_new_method_stability.R     read_and_sort_network_internal() +
                                         compute_pairwise_jaccard_internal() -> stability_ji()

Semantics reproduced line by line (see the R source for the originals):
  * GT: columns Source, Target, upper-cased, distinct. Tested TFs / targets =
    unique GT sources / targets.
  * Edge list: |score|, upper-cased, duplicates collapsed by max(score),
    filtered to (tf in tested TFs) & (target in tested targets).
  * AUROC: on the method's INFERRED (filtered) edges only; balanced 1:1 by
    down-sampling the larger class; pROC::roc(direction "<") AUC, which equals
    the Mann-Whitney AUC with ties counted 1/2 (sklearn roc_auc_score).
  * AUPRC: on the COMPLETE universe tested_TFs x tested_targets, unrecovered
    pairs scored 0; all positives + a random 10x-positives sample of negatives;
    PRROC::pr.curve(...)$auc.integral (continuous interpolation, ported below).
  * Random AUPRC: the same, with runif scores over the complete universe.
  * Early metrics: |score| sorted descending, filtered to GT TFs/targets,
    head(max_edges = 10,000), precision = TP/(TP+FP), recall = TP/|GT|, F1.
  * Stability: per network, distinct (Source, Target, Score), drop NA/inf,
    top round(10%) by Score (and a random 10% control); pairwise Jaccard over
    all network pairs; median.

The only intended difference is the random-number generator: R's
`set.seed(42 + i); sample_n()` cannot be reproduced bit-for-bit in numpy, so
the 1:1 (AUROC) and 1:10 (AUPRC) subsamples differ. `n_rep` repeats the
subsampling and reports the mean (n_rep=1 = the paper's single draw). For
byte-exact numbers run the official R code: scripts/24_bear_official_score.R.
"""
from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

# method -> (tf_col, target_col, score_col) as used by BEAR-GRN's ROC/PR code
# (benchmark_new_method_filtered.R, existing_method_info default).
ROC_PR_COLUMNS = {
    "CellOracle": ("source", "target", "coef_mean"),
    "SCENIC+": ("Source", "Target", "Score"),
    "Pando": ("tf", "target", "estimate"),
    "LINGER": ("Source", "Target", "Score"),
    "FigR": ("Motif", "DORC", "Score"),
    "TRIPOD": ("TF", "gene", "abs_coef"),
    "GRaNIE": ("TF.name", "gene.name", "TF_gene.r"),
    "Pando_xgb": ("tf", "target", "corr"),       # reproduce_ROC_PR_full.R adds this row
}
# early-metrics mapping (benchmark_new_method_early_metrics.R): note the
# different CellOracle score column and the extra Pando_xgb / DIRECTNET rows.
EARLY_COLUMNS = {
    "CellOracle": ("source", "target", "coef_abs"),
    "FigR": ("Motif", "DORC", "Score"),
    "GRaNIE": ("TF.name", "gene.name", "TF_gene.r"),
    "LINGER": ("Source", "Target", "Score"),
    "Pando": ("tf", "target", "estimate"),
    "Pando_xgb": ("tf", "target", "corr"),
    "SCENIC+": ("Source", "Target", "Score"),
    "TRIPOD": ("TF", "gene", "abs_coef"),
    "DIRECTNET": ("TF", "Gene", "Score"),
}
# dataset -> primary GT file (identical in both R functions)
DATASET_GT = {
    "K562": "filtered_RN117_K562.tsv",
    "Macrophage_S1": "filtered_RN204_Buffer1.tsv",
    # reproduce_ROC_PR_full.R (the paper's numbers) uses Buffer2 for S2;
    # benchmark_new_method*.R map S2 to Buffer1 (a bug). We follow the paper.
    "Macrophage_S2": "filtered_RN204_Buffer2.tsv",
    "mESC_E7.5_rep1": "filtered_RN111_E7.5_rep1.tsv",
    "mESC_E7.5_rep2": "filtered_RN111_E7.5_rep2.tsv",
    "mESC_E8.5_rep1": "filtered_RN111_E8.5_rep1.tsv",
    "mESC_E8.5_rep2": "filtered_RN111_E8.5_rep2.tsv",
    "iPS": "filtered_RN000_iPS.tsv",
    "Naive_mESC": "filtered_RN111_Naive_mESC.tsv",
}


# ----------------------------------------------------------------------------- IO
def read_table(path: str) -> pd.DataFrame:
    """readr-like read. Gene symbols such as "NA" / "NAN" must stay strings:
    in R, read_tsv turns "NA" into NA, and paste() then turns it back into
    the string "NA" on both the GT and the method side, so it behaves as a
    gene name. Only empty fields are treated as missing."""
    sep = "," if str(path).endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, low_memory=False, keep_default_na=False, na_values=[""])


def load_gt(path: str) -> pd.DataFrame:
    gt = read_table(path)
    gt = pd.DataFrame({"Source": gt["Source"].astype(str).str.upper(),
                       "Target": gt["Target"].astype(str).str.upper()}).drop_duplicates()
    return gt.reset_index(drop=True)


def standardise_edges(df: pd.DataFrame, tf_col: str, target_col: str,
                      score_col: Optional[str]) -> pd.DataFrame:
    """select -> |score| -> upper -> group_by(tf, target) max(score)."""
    out = pd.DataFrame({"tf": df[tf_col].astype(str).str.upper(),
                        "target": df[target_col].astype(str).str.upper(),
                        "score": np.abs(pd.to_numeric(df[score_col], errors="coerce"))
                        if score_col is not None else 1.0})
    return out.groupby(["tf", "target"], as_index=False, sort=False)["score"].max()


# ----------------------------------------------------------------------------- PRROC port
def prroc_auc_integral(scores_pos: np.ndarray, scores_neg: np.ndarray) -> float:
    """PRROC::pr.curve(scores.class0=pos, scores.class1=neg)$auc.integral
    (PRROC 1.4, compute.pr, unweighted two-vector call; continuous
    interpolation between the operating points of distinct scores)."""
    sp = np.asarray(scores_pos, dtype=np.float64)
    sn = np.asarray(scores_neg, dtype=np.float64)
    scores = np.concatenate([sp, sn])
    wpos = np.concatenate([np.ones_like(sp), np.zeros_like(sn)])
    wneg = 1.0 - wpos
    o = np.argsort(-scores, kind="stable")
    scores, wpos, wneg = scores[o], wpos[o], wneg[o]
    cpos, cneg = np.cumsum(wpos), np.cumsum(wneg)
    use = np.append(scores[:-1] != scores[1:], True)
    tp, fp = cpos[use], cneg[use]
    r = wpos.sum()
    tp_prev = np.concatenate([[0.0], tp[:-1]])
    fp_prev = np.concatenate([[0.0], fp[:-1]])
    with np.errstate(divide="ignore", invalid="ignore"):
        h = (fp - fp_prev) / (tp - tp_prev)
        a = 1.0 + h
        b = (fp_prev - h * tp_prev) / r
        same = tp == tp_prev
        h[same] = 1.0
        a[same] = 1.0
        b[same] = 0.0
        v = (tp / r - tp_prev / r - b / a * (np.log(a * tp / r + b) - np.log(a * tp_prev / r + b))) / a
        v2 = (tp / r - tp_prev / r) / a
    v = np.where(b == 0, v2, v)
    return float(np.sum(v))


# ----------------------------------------------------------------------------- ROC / PR
def score_roc_pr(edges: pd.DataFrame, gt: pd.DataFrame, rng: np.random.Generator,
                 n_rep: int = 1, filtered: bool = True) -> Dict[str, float]:
    """edges: standardised (tf, target, score). Returns AUROC, AUPRC,
    AUPRC_random (+ counts), each the mean over `n_rep` subsampling draws."""
    tested_tfs = pd.unique(gt["Source"])
    tested_targets = pd.unique(gt["Target"])
    gt_keys = set(zip(gt["Source"], gt["Target"]))
    n_total_edges = len(edges)
    inf = edges
    if filtered:
        inf = inf[inf["tf"].isin(tested_tfs) & inf["target"].isin(tested_targets)]
    lab = np.fromiter(((t, g) in gt_keys for t, g in zip(inf["tf"], inf["target"])),
                      dtype=bool, count=len(inf))
    sc = inf["score"].to_numpy(dtype=np.float64)
    pos_s, neg_s = sc[lab], sc[~lab]
    res = {"n_total_edges": int(n_total_edges), "n_evaluable_edges": int(len(inf)),
           "n_positives": int(lab.sum()), "n_negatives": int((~lab).sum())}
    if lab.sum() == 0 or (~lab).sum() == 0:
        res.update({"AUROC": float("nan"), "AUPRC": float("nan"), "AUPRC_random": float("nan")})
        return res

    # complete universe: tested_tfs x tested_targets, unrecovered = 0
    tf_ix = {t: i for i, t in enumerate(tested_tfs)}
    tg_ix = {t: i for i, t in enumerate(tested_targets)}
    n_tf, n_tg = len(tested_tfs), len(tested_targets)
    full = np.zeros(n_tf * n_tg, dtype=np.float64)
    ki = np.fromiter((tf_ix[t] * n_tg + tg_ix[g] for t, g in zip(inf["tf"], inf["target"])),
                     dtype=np.int64, count=len(inf))
    full[ki] = np.nan_to_num(sc, nan=0.0)          # R: is.na(score) -> 0 only for unmatched;
    gt_ki = np.fromiter((tf_ix[s] * n_tg + tg_ix[t] for s, t in gt_keys), dtype=np.int64,
                        count=len(gt_keys))
    is_pos = np.zeros(n_tf * n_tg, dtype=bool)
    is_pos[gt_ki] = True
    pos_idx = np.where(is_pos)[0]
    neg_idx = np.where(~is_pos)[0]
    n_neg10 = min(len(pos_idx) * 10, len(neg_idx))

    aurocs, auprcs, rands = [], [], []
    # When the 10x cap is not binding, the AUPRC sample is the whole negative set,
    # so AUPRC is deterministic: compute it once (identical value), saving time.
    pr_fixed = n_neg10 == len(neg_idx)
    for rep in range(max(1, n_rep)):
        if len(neg_s) < len(pos_s):
            ps, ns = rng.choice(pos_s, size=len(neg_s), replace=False), neg_s
        else:
            ps, ns = pos_s, rng.choice(neg_s, size=len(pos_s), replace=False)
        y = np.concatenate([np.ones(len(ps)), np.zeros(len(ns))])
        x = np.concatenate([ps, ns])
        ok = np.isfinite(x)                       # pROC::roc(na.rm = TRUE)
        aurocs.append(float(roc_auc_score(y[ok], x[ok])))
        if pr_fixed and rep > 0:
            auprcs.append(auprcs[0])
        else:
            nsel = neg_idx if pr_fixed else rng.choice(neg_idx, size=n_neg10, replace=False)
            auprcs.append(prroc_auc_integral(full[pos_idx], full[nsel]))
        if rep < 3:                     # random baseline: 3 draws are plenty (sd ~1e-4)
            ru = rng.random(n_tf * n_tg)
            nsel_r = neg_idx if pr_fixed else rng.choice(neg_idx, size=n_neg10, replace=False)
            rands.append(prroc_auc_integral(ru[pos_idx], ru[nsel_r]))
    res.update({"AUROC": float(np.mean(aurocs)), "AUPRC": float(np.mean(auprcs)),
                "AUPRC_random": float(np.mean(rands)), "n_rep": int(max(1, n_rep)),
                "universe_size": int(n_tf * n_tg), "universe_positives": int(len(pos_idx)),
                "AUROC_sd": float(np.std(aurocs)), "AUPRC_sd": float(np.std(auprcs))})
    return res


# ----------------------------------------------------------------------------- early metrics
def score_early(df_sorted: pd.DataFrame, gt: pd.DataFrame, max_edges: int = 10_000) -> Dict[str, float]:
    """df_sorted: columns source, target, already sorted by |score| desc (R:
    arrange(desc(abs(score)))). Duplicated rows are NOT collapsed (as in R)."""
    gt_pairs = set(gt["Source"] + "_" + gt["Target"])
    gt_tfs, gt_tgs = set(gt["Source"]), set(gt["Target"])
    s = df_sorted["source"].astype(str).str.upper()
    t = df_sorted["target"].astype(str).str.upper()
    keep = s.isin(gt_tfs) & t.isin(gt_tgs)
    s, t = s[keep], t[keep]
    n = min(max_edges, len(s))
    if n == 0:
        return {"Precision": float("nan"), "Recall": float("nan"), "F1": float("nan"),
                "Final_network_size": 0}
    pairs = set((s.iloc[:n] + "_" + t.iloc[:n]).tolist())
    tp = len(pairs & gt_pairs)
    fp = len(pairs - gt_pairs)
    fn = len(gt_pairs - pairs)
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    return {"Precision": prec, "Recall": rec, "F1": f1, "TP": tp, "FP": fp, "FN": fn,
            "Filtered_network_size": int(keep.sum()), "Final_network_size": int(n)}


def early_from_edges(edges: pd.DataFrame, gt: pd.DataFrame, max_edges: int = 10_000) -> Dict[str, float]:
    """Convenience: standardised (tf, target, score) -> sort desc -> score_early."""
    d = edges.sort_values("score", ascending=False, kind="stable")
    return score_early(d.rename(columns={"tf": "source"}), gt, max_edges)


# ----------------------------------------------------------------------------- stability
def top_fraction(edges: pd.DataFrame, frac: float = 0.1) -> set:
    """read_and_sort_network_internal: distinct, drop NA/inf, top round(frac*n)."""
    d = edges.drop_duplicates()
    d = d[np.isfinite(d["score"].to_numpy(dtype=np.float64))]
    n_top = max(1, min(int(round(frac * len(d))), len(d)))
    d = d.sort_values("score", ascending=False, kind="stable").head(n_top)
    return set(zip(d["tf"], d["target"]))


def jaccard(a: set, b: set) -> float:
    u = len(a | b)
    return float("nan") if u == 0 else len(a & b) / u


def stability_ji(networks: Sequence[pd.DataFrame], frac: float = 0.1) -> Dict[str, float]:
    tops = [top_fraction(n, frac) for n in networks]
    ji = [jaccard(tops[i], tops[j]) for i in range(len(tops)) for j in range(i + 1, len(tops))]
    return {"median_JI": float(np.nanmedian(ji)) if ji else float("nan"), "n_pairs": len(ji),
            "JI": ji}
