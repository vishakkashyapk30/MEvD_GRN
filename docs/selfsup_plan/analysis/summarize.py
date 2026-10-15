"""Compile labelfree__*.csv + BEAR published numbers into labelfree_summary.md (ChIP / Union / KO)."""
import glob
import pandas as pd
P = pd.read_csv("/home/vishak/.cache/bear/paper/bear_grn_auroc_auprc_tidy.csv")
P["dataset"] = P["dataset"].replace({"E7.5_rep1": "mESC_E7.5_rep1", "E7.5_rep2": "mESC_E7.5_rep2",
                                     "E8.5_rep1": "mESC_E8.5_rep1", "E8.5_rep2": "mESC_E8.5_rep2"})
D = pd.concat([pd.read_csv(f) for f in sorted(glob.glob("labelfree__*.csv"))])
LF = ["coverage", "tgt_expr_mean", "tgt_detection", "tgt_atac_activity", "tgt_rp_proximal", "tgt_n_peaks",
      "tfexpr_x_tgtatac", "motif_x_atac"]
order = ["K562", "Macrophage_S1", "Macrophage_S2", "iPS", "mESC_E7.5_rep1", "mESC_E7.5_rep2",
         "mESC_E8.5_rep1", "mESC_E8.5_rep2", "Naive_mESC"]
out = ["# Label-free baselines vs published BEAR methods (computed 2026-10-10)", "",
       "Published = BEAR Supp Data 2-5 (tidy CSV). Label-free rows: scored with BEAR's rules (fast port,",
       "AUROC mean of 3 draws). `oracle_*` rows are GT-derived diagnostics, NOT methods.", ""]
for gt in ["ChIP", "Union", "KO", "Intersection", "Core", "CellTypeExclusive"]:
    rows = []
    for ds in order:
        d = D[(D["dataset"] == ds) & (D["gt"] == gt)].set_index("variant")
        if d.empty:
            continue
        p = P[(P["dataset"] == ds) & (P["gt"] == gt)].set_index("method")
        lf = d.loc[[v for v in LF if v in d.index]]
        bA, bP = lf["AUROC"].idxmax(), lf["AUPRC"].idxmax()
        r = {"dataset": ds, "random": d["AUPRC_random"].iloc[0],
             "cov AUPRC": d.loc["coverage", "AUPRC"],
             "tgt_expr AUROC/AUPRC": f"{d.loc['tgt_expr_mean','AUROC']:.3f}/{d.loc['tgt_expr_mean','AUPRC']:.3f}",
             "motif_x_atac AUROC/AUPRC": f"{d.loc['motif_x_atac','AUROC']:.3f}/{d.loc['motif_x_atac','AUPRC']:.3f}",
             "best label-free AUROC": f"{lf.loc[bA,'AUROC']:.3f} ({bA})",
             "best label-free AUPRC": f"{lf.loc[bP,'AUPRC']:.3f} ({bP})",
             "motif_x_atac swap gap": d.loc["motif_x_atac", "swap_gap"],
             "motif TFs / TF nodes": f"{int(d['motif_tfs'].iloc[0])}/{int(d['tf_nodes'].iloc[0])}",
             "oracle TF-outdeg AUROC": d.loc["oracle_tf_outdegree", "AUROC"],
             "oracle LOO in-deg AUROC/AUPRC": f"{d.loc['oracle_loo_indegree','AUROC']:.3f}/{d.loc['oracle_loo_indegree','AUPRC']:.3f}"}
        if len(p):
            r["LINGER AUROC/AUPRC"] = f"{p.loc['LINGER','AUROC']:.3f}/{p.loc['LINGER','AUPRC']:.3f}" if "LINGER" in p.index else ""
            r["best published AUROC"] = f"{p['AUROC'].max():.3f} ({p['AUROC'].idxmax()})"
            r["best published AUPRC"] = f"{p['AUPRC'].max():.3f} ({p['AUPRC'].idxmax()})"
            r["label-free beats all published (AUROC, AUPRC)"] = \
                f"{'Y' if lf['AUROC'].max() > p['AUROC'].max() else 'n'}, {'Y' if lf['AUPRC'].max() > p['AUPRC'].max() else 'n'}"
        rows.append(r)
    if not rows:
        continue
    t = pd.DataFrame(rows)
    out += [f"## {gt}", "", t.to_markdown(index=False, floatfmt=".3f"), ""]
open("labelfree_summary.md", "w").write("\n".join(out))
print("\n".join(out))
