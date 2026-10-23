#!/usr/bin/env python
"""Tables for the statistical baselines and the scMultiomeGRN ablations (scripts/22, scripts/18 --ablation).

For every method: held-out-test AUPR / AUROC (official test fold) and the same on the degree-matched test
(scmgrn_protocol.degree_matched_test), averaged over folds (and seeds) within a cell type, then over the cell
types where both the method and the reference are complete; mean fold-paired difference to the reference and a
Wilcoxon test over all cell-type x fold pairs.

Usage: python scripts/25_scmgrn_baseline_table.py --root ~/scmgrn --ref official_scmultiomegrn_full --out base_vs_smg
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

SECTIONS = [
    ("Reference models", [
        ("mevd_fm_h384_obs_dm", "MeVD-GRN"),
        ("official_scmultiomegrn_full", "scMultiomeGRN"),
    ]),
    ("Counting rules on the training links", [
        ("heur_degree_product", "Degree product (preferential attachment)"),
        ("heur_degree_sum", "Degree sum"),
        ("heur_common_neighbours", "Common neighbours"),
        ("heur_jaccard", "Jaccard index"),
        ("heur_adamic_adar", "Adamic-Adar index"),
        ("heur_resource_allocation", "Resource allocation index"),
        ("heur_katz", "Katz index"),
        ("heur_rwr", "Random walk with restart"),
        ("heur_svd_lowrank", "Low-rank SVD of the training links"),
        ("heur_spectral_lowrank", "Low-rank spectral reconstruction"),
        ("cf_lookalike", "Look-alike TFs (matrix completion of the training links)"),
        ("cf_lookalike_norm", "Look-alike TFs, normalised"),
        ("gf_family_top5", "Geneformer look-alikes of u linked to v (training links)"),
    ]),
    ("Scores from the cell data, no labels", [
        ("unsup_pearson_abs", "Pearson correlation of expression"),
        ("unsup_spearman_abs", "Spearman correlation of expression"),
        ("unsup_partial_corr_abs", "Partial correlation of expression"),
        ("unsup_mutual_information", "Mutual information of expression"),
        ("unsup_clr", "CLR on mutual information"),
        ("unsup_grnboost_tf_tf", "GRNBoost2 between the two TFs"),
        ("unsup_genie3", "GENIE3 between the two TFs"),
        ("unsup_ridge_regression", "Ridge regression coefficients"),
        ("unsup_lasso_regression", "Lasso regression coefficients"),
        ("unsup_maestro_pearson_abs", "Pearson correlation of MAESTRO accessibility"),
        ("unsup_grnboost_cosine", "Cosine similarity of GRNBoost2 profiles"),
        ("unsup_geneformer_cosine", "Cosine similarity of Geneformer embeddings"),
    ]),
    ("Logistic regression on features", [
        ("lr_maestro", "MAESTRO accessibility only"),
        ("lr_grnboost", "GRNBoost2 profiles only"),
        ("lr_geneformer", "Geneformer embeddings only"),
        ("lr_histogram", "Joint expression histograms only"),
        ("lr_promoter_openness", "Promoter openness of the two TFs only"),
        ("lr_maestro_grnboost_concat", "MAESTRO and GRNBoost2 side by side"),
        ("lr_maestro_x_grnboost", "MAESTRO multiplied by GRNBoost2 (per TF)"),
        ("lr_within_pair_product", "Product of the two TFs' features only"),
        ("linreg_maestro_grnboost_concat", "Linear regression, MAESTRO and GRNBoost2"),
        ("rf_maestro_grnboost_concat", "Random forest, MAESTRO and GRNBoost2"),
    ]),
    ("Logistic regression on training links plus features", [
        ("lr_degree_only", "Degree only"),
        ("lr_topology_only", "Degree and neighbourhood scores"),
        ("lr_degree_maestro_grnboost", "Degree + MAESTRO + GRNBoost2"),
        ("lr_degree_geneformer", "Degree + Geneformer"),
        ("lr_topology_maestro_grnboost", "Degree and neighbourhood + MAESTRO + GRNBoost2"),
        ("lr_topology_geneformer", "Degree and neighbourhood + Geneformer"),
    ]),
    ("scMultiomeGRN ablations", [
        ("official_scmultiomegrn_no_prior", "scMultiomeGRN without the TF-TF graph"),
        ("official_scmultiomegrn_string_prior", "scMultiomeGRN with STRING instead of the training links"),
        ("official_scmultiomegrn_no_edge", "scMultiomeGRN without the pair histograms"),
        ("official_scmultiomegrn_atac_only", "scMultiomeGRN with ATAC features only"),
        ("official_scmultiomegrn_scrna_only", "scMultiomeGRN with scRNA features only"),
        ("official_scmultiomegrn_plus_geneformer", "scMultiomeGRN + Geneformer embeddings"),
    ]),
]
EVALS = (("test", "official test"), ("test_degmatched", "degree-matched test"))
METRICS = ("aupr", "auroc")


def collect(root: Path, keep):
    """runs[(ds, ct)][method][fold] -> {(eval, metric): [per-seed values]}"""
    runs = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
    for p in root.glob("results/*/*/*/seed*/fold*.json"):
        ds, ct, method = p.parts[-5], p.parts[-4], p.parts[-3]
        if ds not in keep:
            continue
        d = json.load(open(p))
        if method.startswith("mevd_") and not (d.get("leak_status") or {}).get("zero_leak", False):
            raise SystemExit(f"{p}: not zero-leak")
        fold = int(p.stem.replace("fold", ""))
        for e, _ in EVALS:
            if e in d["metrics"]:
                for m in METRICS:
                    runs[(ds, ct)][method][fold][(e, m)].append(d["metrics"][e][m])
    return runs


def summarise(runs, method, ref, n_folds):
    per_ct, diffs, n_ct = defaultdict(list), defaultdict(list), 0
    for key, bym in runs.items():
        if method not in bym or ref not in bym or len(bym[method]) < n_folds or len(bym[ref]) < n_folds:
            continue
        n_ct += 1
        folds = sorted(set(bym[method]) & set(bym[ref]))[:n_folds]
        for e, _ in EVALS:
            for m in METRICS:
                k = (e, m)
                if not all(bym[method][f][k] and bym[ref][f][k] for f in folds):
                    continue
                v = [np.mean(bym[method][f][k]) for f in folds]
                h = [np.mean(bym[ref][f][k]) for f in folds]
                per_ct[k].append(np.mean(v))
                diffs[k] += [x - y for x, y in zip(v, h)]
    row = {"method": method, "n_cell_types": n_ct}
    for e, _ in EVALS:
        for m in METRICS:
            k = (e, m)
            if not per_ct[k]:
                continue
            dd = np.array(diffs[k])
            row[f"{e}_{m}"] = float(np.mean(per_ct[k]))
            row[f"{e}_{m}_diff"] = float(dd.mean())
            row[f"{e}_{m}_p"] = float(wilcoxon(dd).pvalue) if method != ref and len(dd) >= 5 and np.any(dd != 0) else None
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--ref", default="official_scmultiomegrn_full")
    ap.add_argument("--out", default="baselines")
    ap.add_argument("--n_folds", type=int, default=10)
    ap.add_argument("--datasets", default="lung,hema")
    a = ap.parse_args()
    runs = collect(Path(a.root).expanduser(), set(a.datasets.split(",")))
    out = {"ref": a.ref, "datasets": a.datasets, "sections": []}
    lines = [f"Reference for differences and p-values: {a.ref}", "",
             "| Method | Cell types | AUPR | diff | p | AUROC | diff | p | AUPR (degree-matched) | diff | AUROC (degree-matched) | diff |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|"]

    def f(r, k, fmt):
        return "–" if r.get(k) is None else format(r[k], fmt)
    for title, items in SECTIONS:
        sec = {"title": title, "rows": []}
        lines.append(f"| **{title}** | | | | | | | | | | | |")
        for method, label in items:
            r = summarise(runs, method, a.ref, a.n_folds)
            r["label"] = label
            sec["rows"].append(r)
            lines.append(f"| {label} | {r['n_cell_types']} | {f(r, 'test_aupr', '.4f')} | {f(r, 'test_aupr_diff', '+.4f')} | "
                         f"{f(r, 'test_aupr_p', '.2g')} | {f(r, 'test_auroc', '.4f')} | {f(r, 'test_auroc_diff', '+.4f')} | "
                         f"{f(r, 'test_auroc_p', '.2g')} | {f(r, 'test_degmatched_aupr', '.4f')} | "
                         f"{f(r, 'test_degmatched_aupr_diff', '+.4f')} | {f(r, 'test_degmatched_auroc', '.4f')} | "
                         f"{f(r, 'test_degmatched_auroc_diff', '+.4f')} |")
        out["sections"].append(sec)
    json.dump(out, open(a.out + ".json", "w"), indent=1)
    Path(a.out + ".md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
