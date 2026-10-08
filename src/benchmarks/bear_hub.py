"""Hub control for BEAR-GRN (pre-registered s12.7): does a method carry signal
beyond the TF-disjoint target in-degree prior? Ported from the K562 analysis in
s10.3 (~/.cache/local_runs/scripts/bear_hub_control.py, bear_per_tf.py) and
generalised to any dataset / method / seed.

In-degree = the score of results/<ds>/baseline_indegree/seed<s>/grn.tsv.gz, i.e. the
target's in-degree among the training TFs of the SAME TF-disjoint fold (same
fold function and seed as scripts/22, checked edge by edge).
  raw            the method's score
  resid_lin      per fold: logit(score) minus its OLS fit on a cubic in log1p(in-degree)
  within_stratum per fold: percentile of the score within its in-degree stratum
                 (0 alone, then 20 quantile bins of the non-zero in-degrees)
Pooled with BEAR's metric (bear_metrics.score_roc_pr), plus per-TF AUROC of each
GT TF's own row (>= 10 pos and >= 10 neg): raw and within 20 in-degree bins.

  python -m src.benchmarks.bear_hub --root ~/.cache/bear/data --dataset K562 \
      --method mevd_fm_h384_hub --seeds 42 43 --n_rep 20
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

from src.benchmarks import bear_metrics as bm
from src.utils.io import load_config

REPO = Path(__file__).resolve().parents[2]


def _read(p: Path) -> pd.DataFrame:
    return pd.read_csv(p, sep="\t", keep_default_na=False, na_values=[""])


def strata(indeg: np.ndarray, nbins: int = 20) -> np.ndarray:
    s = np.zeros(len(indeg), dtype=np.int64)
    nz = indeg > 0
    if nz.any():
        q = np.unique(np.quantile(indeg[nz], np.linspace(0, 1, nbins + 1)[1:-1]))
        s[nz] = 1 + np.searchsorted(q, indeg[nz], side="right")
    return s


def resid_lin(x: np.ndarray, indeg: np.ndarray) -> np.ndarray:
    eps = 1e-6
    y = np.log(np.clip(x, eps, 1 - eps) / (1 - np.clip(x, eps, 1 - eps))) if x.max() <= 1 and x.min() >= 0 else x
    z = np.log1p(indeg)
    X = np.c_[np.ones_like(z), z, z ** 2, z ** 3]
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def hub_control(root: Path, ds: str, method: str, seed: int, n_rep: int = 20,
                gts_filter=None) -> list:
    R = root / "results" / ds
    m = _read(R / method / f"seed{seed}" / "grn.tsv.gz")
    d = _read(R / "baseline_indegree" / f"seed{seed}" / "grn.tsv.gz")
    ok = (len(m) == len(d) and (m["Source"].values == d["Source"].values).all()
          and (m["Target"].values == d["Target"].values).all()
          and ("fold" not in m or (m["fold"].values == d["fold"].values).all()))
    if not ok:     # align by pair (e.g. L2 single-fold outputs)
        d = d.set_index(["Source", "Target"]).reindex(pd.MultiIndex.from_arrays([m.Source, m.Target])).reset_index()
        d["fold"] = m["fold"].values if "fold" in m else 0
    df = pd.DataFrame({"Source": m["Source"], "Target": m["Target"],
                       "fold": m["fold"] if "fold" in m else 0,
                       "x": m["Score"].astype(float), "indeg": d["Score"].astype(float).fillna(0.0)})
    df["stratum"] = -1
    rl, ws = np.zeros(len(df)), np.zeros(len(df))
    for _f, g in df.groupby("fold"):
        ix = g.index.to_numpy()
        st = strata(g["indeg"].to_numpy())
        df.loc[ix, "stratum"] = st
        rl[ix] = resid_lin(g["x"].to_numpy(), g["indeg"].to_numpy())
        ws[ix] = pd.Series(g["x"].to_numpy()).groupby(st).rank(pct=True).to_numpy()
    variants = {"raw": df["x"].to_numpy(), "resid_lin": rl - rl.min() + 1e-3, "within_stratum": ws + 1e-3}
    cfg = load_config(REPO / "configs/bear/datasets.yaml")["datasets"][ds]
    rows = []
    for gname, pth in cfg["gts"].items():
        if gts_filter and gname not in gts_filter:
            continue
        if not (root / pth).exists():
            continue
        gt = bm.load_gt(str(root / pth))
        for vname, sc in variants.items():
            e = bm.standardise_edges(pd.DataFrame({"S": df.Source, "T": df.Target, "x": sc}), "S", "T", "x")
            r = bm.score_roc_pr(e, gt, np.random.default_rng(43), n_rep=n_rep)
            rows.append({"dataset": ds, "method": method, "seed": seed, "gt": gname, "variant": vname,
                         "AUROC": r["AUROC"], "AUPRC": r["AUPRC"], "AUPRC_random": r["AUPRC_random"]})
        # per-TF AUROC (raw and within in-degree bins) over in-space pairs
        keys = set(zip(gt.Source, gt.Target))
        S_u, T_u = df.Source.str.upper(), df.Target.str.upper()
        ins = df[S_u.isin(set(gt.Source)) & T_u.isin(set(gt.Target))]
        y = np.fromiter(((a.upper(), b.upper()) in keys for a, b in zip(ins.Source, ins.Target)), bool, len(ins))
        raw_a, win_a, ind_a = [], [], []
        for _tf, idx in pd.Series(np.arange(len(ins))).groupby(ins["Source"].to_numpy()):
            ii = idx.to_numpy()
            yy = y[ii]
            if yy.sum() < 10 or (~yy).sum() < 10:
                continue
            x = ins["x"].to_numpy()[ii]
            dg = ins["indeg"].to_numpy()[ii]
            raw_a.append(roc_auc_score(yy, x))
            ind_a.append(roc_auc_score(yy, dg))
            win = pd.Series(x).groupby(strata(dg)).rank(pct=True).to_numpy()
            win_a.append(roc_auc_score(yy, win))
        rows.append({"dataset": ds, "method": method, "seed": seed, "gt": gname, "variant": "per_TF",
                     "n_tfs": len(raw_a), "perTF_raw": float(np.mean(raw_a)) if raw_a else np.nan,
                     "perTF_indegree": float(np.mean(ind_a)) if ind_a else np.nan,
                     "perTF_within_bins": float(np.mean(win_a)) if win_a else np.nan,
                     "perTF_share_within_gt_0.5": float(np.mean(np.array(win_a) > 0.5)) if win_a else np.nan})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--method", required=True, help="results dir name, e.g. mevd_fm_h384_hub")
    ap.add_argument("--seeds", nargs="+", type=int, default=[42])
    ap.add_argument("--n_rep", type=int, default=20)
    ap.add_argument("--gts", nargs="*", default=None)
    args = ap.parse_args()
    root = Path(args.root)
    out = root / "results" / args.dataset / "hub_control"
    out.mkdir(parents=True, exist_ok=True)
    for s in args.seeds:
        rows = hub_control(root, args.dataset, args.method, s, args.n_rep, args.gts)
        pd.DataFrame(rows).to_csv(out / f"{args.method}__seed{s}.csv", index=False)
        for r in rows:
            print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}), flush=True)


if __name__ == "__main__":
    main()
