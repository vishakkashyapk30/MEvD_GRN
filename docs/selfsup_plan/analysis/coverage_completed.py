"""Coverage-completed released networks: append every measured (TF node x gene) pair the method did
not emit, with score = 0.5 x the method's minimum positive score, then re-score with BEAR's rules.
Puts every method at the same coverage as a dense label-free method. Usage: ds [inf_root]."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, "/home/vishak/research/mevd_grn/docs/selfsup_plan/analysis")
import labelfree_baselines as lb
from src.benchmarks import bear_metrics as bm
ds = sys.argv[1]
inf_root = Path(sys.argv[2]) if len(sys.argv) > 2 else lb.ROOT
names, tf, *_ = lb.load_proc(ds)
n, T = len(names), len(tf)
src = np.repeat(np.arange(T), n); dst = np.tile(np.arange(n), T)
keep = tf[src] != dst
cand = pd.DataFrame({"tf": names[tf[src[keep]]], "target": names[dst[keep]]})
rel_root = inf_root / "INFERRED.GRNS" / lb.CFG["inferred_dir"][ds]
rows = []
for mdir, (tc, gc, scol) in bm.ROC_PR_COLUMNS.items():
    d = rel_root / mdir
    files = sorted(f for f in d.iterdir() if f.is_file()) if d.exists() else []
    if not files:
        continue
    df = bm.read_table(str(files[0]))
    if scol not in df.columns:
        continue
    e = bm.standardise_edges(df, tc, gc, scol); del df
    e = e[np.isfinite(e["score"].to_numpy(np.float64))]
    mn = e.loc[e["score"] > 0, "score"].min()
    m = cand.merge(e, on=["tf", "target"], how="left", indicator=True)
    add = m[m["_merge"] == "left_only"][["tf", "target"]].assign(score=0.5 * mn)
    full = pd.concat([e, add], ignore_index=True)
    for g in ["ChIP", "Union"]:
        gt = lb.gts_of(ds).get(g)
        if gt is None:
            continue
        sc = lb.Scorer(gt)
        for tag, ee in (("original", e), ("coverage_completed", full)):
            ok, key = sc.keys(ee["tf"].to_numpy(), ee["target"].to_numpy())
            r = sc.score(key, ee["score"].to_numpy(np.float64)[ok], np.random.default_rng(43), n_rep=3)
            rows.append({"dataset": ds, "gt": g, "method": mdir, "version": tag, **r})
            print(ds, g, mdir, tag, round(r["AUROC"], 3), round(r["AUPRC"], 3), r["n_eval"], flush=True)
pd.DataFrame(rows).to_csv(lb.OUT / f"coverage_completed__{ds}.csv", index=False)
