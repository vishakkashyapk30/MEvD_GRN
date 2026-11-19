#!/usr/bin/env python
"""Tables for the BEAR-GRN degree-controlled evaluation (src/benchmarks/bear_matched.py output).

Per ground truth: one row per method, mean over seeds within a dataset, then over datasets (with the number of
datasets); columns = in-space AUROC, degree-matched, caliper-matched (AUROC and pair-level win rate),
covariate-adjusted (TF x in-degree group) and degree-corrected (Aiyappa et al.) AUROC. Also a per-dataset table of
the caliper-matched AUROC.

  python scripts/bear/27_bear_matched_table.py --root ~/.cache/bear/data --out docs/results/bear/matched
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

COLS = [("AUROC_inspace", "AUROC"), ("AUROC_perTF", "AUROC per TF"), ("AUROC_degmatched", "degree-matched"), ("AUROC_caliper", "caliper-matched"),
        ("pairwin_caliper", "caliper, pair win rate"), ("AUROC_adjusted", "covariate-adjusted"),
        ("AUROC_aiyappa", "degree-corrected (Aiyappa)")]
LABEL = {
    "mevd_fm_h384_hubmotif": "MeVD-GRN", "mevd_fm_h384_hubmotif_L2": "MeVD-GRN, L2 (no cell-type labels)",
    "mevd_fm_h384_hubmotif_noatac": "MeVD-GRN, ATAC features zeroed (old ablation)",
    "mevd_fm_h384_hubmotif_atacnone": "MeVD-GRN, ATAC removed (all routes)",
    "mevd_fm_h384_hubmotif_atacshuffle": "MeVD-GRN, ATAC shuffled across genes",
    "mevd_fm_h384_hubmotif_nofm": "MeVD-GRN without Geneformer", "mevd_fm_h384_hubmotif_seq": "MeVD-GRN, sequential curriculum",
    "baseline_indegree": "In-degree", "baseline_pearson": "Pearson", "baseline_grnboost2": "GRNBoost2",
    "baseline_coverage": "All measured pairs", "lf_ggm_ridge": "Gaussian graphical model (ridge, label-free)",
    "lf_regdiffusion": "RegDiffusion (label-free)",
}
MODEL = {"lr": "logistic regression", "lda": "LDA", "qda": "QDA", "gnb": "Gaussian naive Bayes",
         "gmm": "Gaussian mixture classifier", "rbf": "Gaussian-kernel logistic regression", "mlp": "MLP (32 units)",
         "hgb": "gradient boosting"}


def label(m: str) -> str:
    if m in LABEL:
        return LABEL[m]
    if m.startswith("simple_"):
        _, mod, fs = m.split("_")
        return f"{MODEL.get(mod, mod)}, {'cell data only' if fs == 'data' else 'all features'}"
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", required=True, help="output prefix (writes <out>.md and <out>.csv)")
    a = ap.parse_args()
    files = sorted((Path(a.root).expanduser() / "results" / "matched").glob("*.csv"))
    d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    d = d.drop_duplicates(["dataset", "gt", "method", "seed"], keep="last")
    num = [c for c, _ in COLS]
    per_ds = d.groupby(["gt", "method", "kind", "dataset"], as_index=False)[num].mean()
    agg = per_ds.groupby(["gt", "method", "kind"], as_index=False).agg(
        **{c: (c, "mean") for c in num}, n_datasets=("dataset", "nunique"))
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    agg.to_csv(out.with_suffix(".csv"), index=False)
    per_ds.to_csv(out.with_name(out.name + "_per_dataset.csv"), index=False)
    lines = ["# BEAR-GRN: degree-controlled evaluation of every method", "",
             "Mean over seeds within a dataset, then over datasets. Pairs a method does not output score 0. Covariate = "
             "target in-degree among the training TFs of the fold holding the scored TF out (seed 42). In-space pairs = "
             "(ground-truth TF, measured ground-truth target). Matched sets are identical for every method.", ""]
    order = {"released": 0, "ours": 1}
    for gt in ["ChIP", "KO", "Union", "Intersection"]:
        t = agg[agg["gt"] == gt].copy()
        if t.empty:
            continue
        t["o"] = t["kind"].map(order)
        t = t.sort_values(["o", "AUROC_caliper"], ascending=[True, False])
        lines += [f"## {gt}", "", "| Method | Datasets | " + " | ".join(h for _, h in COLS) + " |",
                  "|---|---|" + "---|" * len(COLS)]
        for _, r in t.iterrows():
            lines.append(f"| {label(r['method'])} | {int(r['n_datasets'])} | "
                         + " | ".join(f"{r[c]:.3f}" for c, _ in COLS) + " |")
        lines.append("")
    out.with_suffix(".md").write_text("\n".join(lines) + "\n")
    print(f"wrote {out.with_suffix('.md')} ({len(agg)} rows)")


if __name__ == "__main__":
    main()
