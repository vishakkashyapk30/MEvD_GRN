#!/usr/bin/env python
"""A2 (P2), A3 (P7) and the per-method statistics for 04_report.py.

For every (dataset, GT, method[, seed]):
  dr        doubly-ranked AUROC, same-fold TF pairs only (README D1)
  dr_null   mean DR over 4 marginal-null GTs (fixed per dataset x GT, shared by all methods)
  dr_excess dr - dr_null
  rowm      row-matched AUROC (pre-registered; fails its controls on sparse GTs, see D1)
  full      full-universe AUROC (unscored pairs = 0)
Released (unsupervised) methods have no folds: they are evaluated under each seed's fold partition (42-46).
Per item, C/N matrices (real and null, same-fold masked) go to results/stats/<ds>/<gt>/<method>__<seed>.npz.

  python 02_audit.py [--datasets K562 ...] [--only mid_] [--no_released]
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

import mdh_lib as L

OURS = ["baseline_pearson", "baseline_grnboost2", "baseline_indegree", "mevd_fm_h384_hubmotif"]
SEEDS = [42, 43, 44, 45, 46]
N_NULL = 4


def our_runs(ds: str, only):
    out = []
    base = L.ROOT / "results" / ds
    for p in sorted(base.iterdir()):
        if not p.is_dir() or not (p.name in OURS or p.name.startswith(("mid_", "mid2_"))):
            continue
        if only and not p.name.startswith(only):
            continue
        for f in sorted(p.glob("seed*/grn.tsv.gz")) + sorted(p.glob("seed*/score_matrix.npz")):
            out.append((p.name, int(f.parent.name[4:]), f))
    return out


def score_matrix(D, rows, path):
    if path.name == "score_matrix.npz":
        z = np.load(path)
        assert np.array_equal(z["tf"], D["tf"])
        return np.nan_to_num(z["S"][rows].astype(np.float64), nan=0.0)
    return L.edges_to_matrix(D, rows, L.read_ours(path))


def evaluate(D, ds, gt, blk, name, S, seeds, recs):
    rows, Y, Z, M, nulls = blk
    C, N = L.dr_pairs(S, Y, M)
    CN = [L.dr_pairs(S, y, M) for y in nulls]
    U, Dn = L.matched_stats(S, Y, Z, M)
    full = L.full_auroc(S, Y, M)
    d = L.OUT / "stats" / ds / gt
    d.mkdir(parents=True, exist_ok=True)
    for s in seeds:
        F = L.fold_pairs(D, rows, s)
        np.savez_compressed(d / f"{name}__{s}.npz", C=(C * F).astype(np.float32), N=(N * F).astype(np.float32),
                            Cn=np.stack([c * F for c, _ in CN]).astype(np.float32),
                            Nn=np.stack([n * F for _, n in CN]).astype(np.float32), rows=rows)
        dr = L.dr_auroc(C, N, F=F)
        drn = float(np.mean([L.dr_auroc(c, n, F=F) for c, n in CN]))
        rec = {"dataset": ds, "gt": gt, "method": name, "seed": s, "dr": dr, "dr_null": drn,
               "dr_excess": dr - drn, "rowm": float(U.sum() / Dn.sum()), "full": full, "n_tfs": int(len(rows))}
        recs.append(rec)
    print(json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in recs[-1].items()}), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="*", default=L.DATASETS)
    ap.add_argument("--only", default=None)
    ap.add_argument("--no_released", action="store_true")
    args = ap.parse_args()
    csv = L.OUT / "a2_metrics.csv"
    recs = []
    for ds in args.datasets:
        D = L.load_ds(ds, features=False)
        blocks = {}
        for gt in L.ds_gts(D):
            rows, Y, Z, M = L.eval_block(D, gt)
            if Y.sum() == 0:
                continue
            blocks[gt] = (rows, Y, Z, M, L.marginal_null_draws(Y, M, N_NULL, seed=0))
        if not (args.no_released or args.only):
            for m in L.RELEASED:
                e = L.read_released(ds, m)
                if e is None:
                    continue
                for gt, blk in blocks.items():
                    evaluate(D, ds, gt, blk, m, L.edges_to_matrix(D, blk[0], e), SEEDS, recs)
                del e
        for name, s, path in our_runs(ds, args.only):
            for gt, blk in blocks.items():
                seeds = SEEDS if name in ("baseline_pearson", "baseline_grnboost2") else [s]
                evaluate(D, ds, gt, blk, name, score_matrix(D, blk[0], path), seeds, recs)
        new = pd.DataFrame(recs)
        if csv.exists():
            old = pd.read_csv(csv)
            key = ["dataset", "gt", "method", "seed"]
            old = old.merge(new[key].drop_duplicates(), on=key, how="left", indicator=True)
            old = old[old["_merge"] == "left_only"].drop(columns="_merge")
            new = pd.concat([old, new], ignore_index=True)
        new.sort_values(["dataset", "gt", "method", "seed"]).to_csv(csv, index=False)
        recs = []

    # A3 (P7): across released methods, pooled BEAR ChIP AUROC vs Intersection AUROC (cached port scores)
    from scipy.stats import spearmanr
    rows = []
    for ds in L.DATASETS:
        for m in L.RELEASED:
            r = {}
            for gt in ("ChIP", "Intersection"):
                f = L.ROOT / "results" / ds / "scores" / f"{m}__{gt}.json"
                if f.exists():
                    r[gt] = json.load(open(f)).get("AUROC")
            if len(r) == 2 and all(v is not None and np.isfinite(v) for v in r.values()):
                rows.append({"dataset": ds, "method": m, **r})
    a3 = pd.DataFrame(rows)
    if len(a3):
        rho, p = spearmanr(a3["ChIP"], a3["Intersection"])
        per = {ds: float(spearmanr(g["ChIP"], g["Intersection"])[0]) for ds, g in a3.groupby("dataset")}
        a3.to_csv(L.OUT / "a3_chip_vs_intersection.csv", index=False)
        json.dump({"spearman_rho_pooled": float(rho), "p": float(p), "n": len(a3), "per_dataset_rho": per},
                  open(L.OUT / "a3_summary.json", "w"), indent=1)
        print("A3", rho, p, per)


if __name__ == "__main__":
    main()
