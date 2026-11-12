"""Can a LABEL-FREE TF-level quantity reproduce the TF-marginal part of BEAR AUROC?
TF-only scores (same score for all targets of a TF), then x target expression.
Uses labelfree_baselines helpers; writes tf_marginal_probe.csv."""
import sys
import numpy as np, pandas as pd
sys.path.insert(0, "/home/vishak/research/mevd_grn/docs/selfsup_plan/analysis")
import labelfree_baselines as lb
rows = []
for ds in sys.argv[1:]:
    names, tf, rna, atac, op, F = lb.load_proc(ds)
    n, T = len(names), len(tf)
    src = np.repeat(np.arange(T), n); dst = np.tile(np.arange(n), T)
    keep = tf[src] != dst; src, dst = src[keep], dst[keep]
    tfn, tgn = names[tf[src]], names[dst]
    hit = (F[:, :, 1] > 0).mean(1)                 # fraction of genes whose top-5 peaks carry a motif hit
    has = F[:, :, 2].max(1) > 0
    hit_f = np.where(has, hit, np.median(hit[has]) if has.any() else 0)
    t_expr = lb.pct(rna[tf, 0]); t_det = lb.pct(rna[tf, 2]); t_hit = lb.pct(hit_f)
    g_expr = lb.pct(rna[:, 0])
    V = {"tf_expr_only": 1 + t_expr[src], "tf_detection_only": 1 + t_det[src],
         "tf_motif_hitrate_only": 1 + t_hit[src],
         "tf_motif_hitrate_x_tgt_expr": 1 + t_hit[src] * g_expr[dst],
         "tf_expr_x_tgt_expr": 1 + t_expr[src] * g_expr[dst]}
    for g, gt in lb.gts_of(ds).items():
        if g not in ("ChIP", "Union"):
            continue
        sc = lb.Scorer(gt); ok, key = sc.keys(tfn, tgn)
        rng = np.random.default_rng(43)
        for v, x in V.items():
            r = sc.score(key, x[ok], rng, n_rep=3)
            rows.append({"dataset": ds, "gt": g, "variant": v, "AUROC": r["AUROC"], "AUPRC": r["AUPRC"]})
            print(ds, g, v, round(r["AUROC"], 3), round(r["AUPRC"], 3), flush=True)
pd.DataFrame(rows).to_csv(lb.OUT / "tf_marginal_probe.csv", index=False)
