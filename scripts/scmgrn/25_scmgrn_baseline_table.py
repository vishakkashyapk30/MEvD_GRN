#!/usr/bin/env python
"""Result tables for the scMultiomeGRN benchmark: MeVD-GRN, scMultiomeGRN, baselines and ablations, for either split.

For every method and evaluation set (official test; degree-matched; caliper-matched; held-out-TF C2 / C3):
AUPR / AUROC averaged over folds (and seeds) within a cell type, then over the cell types where both the method and
the reference are complete; a 95% CI from a cluster bootstrap over cell types (folds of one cell type are not
independent); the mean fold-paired difference to the reference with a Wilcoxon signed-rank test over all
cell-type x fold pairs; and the share of positives (= random AUPR) of the evaluation set.
Rows of methods with no complete cell type are left out of the Markdown table.

Usage: python scripts/scmgrn/25_scmgrn_baseline_table.py --root ~/scmgrn --ref official_scmultiomegrn_full --out tables/vs_smg
       python scripts/scmgrn/25_scmgrn_baseline_table.py --root ~/scmgrn --datasets lung_tfholdout,hema_tfholdout \
           --ref official_scmultiomegrn --evals test_c2,test_c3,test_degcaliper --out tables/tfholdout_vs_smg
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
        ("mevd_fm_h384_obs", "MeVD-GRN"),
        ("official_scmultiomegrn_full", "scMultiomeGRN"),
        ("official_scmultiomegrn", "scMultiomeGRN"),
    ]),
    ("Graph neural network baselines (scMultiomeGRN paper)", [
        ("deeptfni", "DeepTFni"),
        ("genelink", "GENELink"),
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
    ("MeVD-GRN ablations", [
        ("mevd_abl_nofm_obs_dm", "MeVD-GRN without Geneformer"),
        ("mevd_abl_noatac_obs_dm", "MeVD-GRN without the ATAC tower"),
        ("mevd_abl_nognn_obs_dm", "MeVD-GRN without the graph network"),
        ("mevd_abl_coexpr_obs_dm", "MeVD-GRN, co-expression graph only"),
        ("mevd_abl_tfcand_obs_dm", "MeVD-GRN, TF-candidate graph only"),
        ("mevd_abl_base_obs_dm", "MeVD-GRN, base size"),
        ("mevd_fm_h384_noobs_dm", "MeVD-GRN without the training-link relation"),
        ("mevd_abl_min_noobs_dm", "Floor model (no Geneformer, ATAC tower, graph network or training links)"),
        ("mevd_abl_fl_atac_noobs_dm", "Floor model + ATAC tower"),
        ("mevd_abl_fl_fm_noobs_dm", "Floor model + Geneformer"),
        ("mevd_abl_fl_gnn_noobs_dm", "Floor model + graph network"),
        ("mevd_abl_nofm_obs", "MeVD-GRN without Geneformer"),
        ("mevd_abl_noatac_obs", "MeVD-GRN without the ATAC tower"),
        ("mevd_abl_nognn_obs", "MeVD-GRN without the graph network"),
        ("mevd_abl_coexpr_obs", "MeVD-GRN, co-expression graph only"),
        ("mevd_abl_tfcand_obs", "MeVD-GRN, TF-candidate graph only"),
        ("mevd_abl_base_obs", "MeVD-GRN, base size"),
        ("mevd_fm_h384_noobs", "MeVD-GRN without the training-link relation"),
        ("mevd_abl_min_noobs", "Floor model (no Geneformer, ATAC tower, graph network or training links)"),
        ("mevd_abl_fl_atac_noobs", "Floor model + ATAC tower"),
        ("mevd_abl_fl_fm_noobs", "Floor model + Geneformer"),
        ("mevd_abl_fl_gnn_noobs", "Floor model + graph network"),
    ]),
]
EVALS = (("test", "official test"), ("test_degmatched", "degree-matched"), ("test_degcaliper", "caliper-matched"),
         ("test_c2", "held-out TF, C2"), ("test_c3", "held-out TF, C3"))
METRICS = ("aupr", "ap", "auroc")   # ap = average precision; unlike the official trapezoidal AUPR it
                                    # is not inflated for tied scores (all-tied -> share of positives)


def collect(root: Path, keep):
    """runs[(ds, ct)][method][fold] -> {(eval, metric): [per-seed values]}; also (eval, "base") = share of positives."""
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
                mm = d["metrics"][e]
                runs[(ds, ct)][method][fold][(e, "base")].append(mm["n_pos"] / max(mm["n"], 1))
    return runs


def cluster_ci(per_ct, n_boot=2000, seed=0):
    """95% CI of the mean over cell types, resampling cell types with replacement."""
    x = np.asarray(per_ct, float)
    if len(x) < 2:
        return None, None
    rng = np.random.default_rng(seed)
    bs = x[rng.integers(0, len(x), size=(n_boot, len(x)))].mean(1)
    return float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def summarise(runs, method, ref, n_folds):
    per_ct, diffs, base, n_ct = defaultdict(list), defaultdict(list), defaultdict(list), 0
    for key, bym in runs.items():
        if method not in bym or ref not in bym or len(bym[method]) < n_folds or len(bym[ref]) < n_folds:
            continue
        n_ct += 1
        folds = sorted(set(bym[method]) & set(bym[ref]))[:n_folds]
        for e, _ in EVALS:
            ok = [f for f in folds if bym[method][f][(e, "aupr")] and bym[ref][f][(e, "aupr")]]
            if not ok:
                continue
            base[e].append(np.mean([np.mean(bym[method][f][(e, "base")]) for f in ok]))
            for m in METRICS:
                v = [np.mean(bym[method][f][(e, m)]) for f in ok]
                h = [np.mean(bym[ref][f][(e, m)]) for f in ok]
                per_ct[(e, m)].append(np.mean(v))
                diffs[(e, m)] += [x - y for x, y in zip(v, h)]
    row = {"method": method, "n_cell_types": n_ct}
    for e, _ in EVALS:
        if base[e]:
            row[f"{e}_base"] = float(np.mean(base[e]))
        for m in METRICS:
            k = (e, m)
            if not per_ct[k]:
                continue
            dd = np.array(diffs[k])
            row[f"{e}_{m}"] = float(np.mean(per_ct[k]))
            row[f"{e}_{m}_ci"] = cluster_ci(per_ct[k])
            row[f"{e}_{m}_diff"] = float(dd.mean())
            row[f"{e}_{m}_p"] = float(wilcoxon(dd).pvalue) if method != ref and len(dd) >= 5 and np.any(dd != 0) else None
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--ref", default="official_scmultiomegrn_full")
    ap.add_argument("--out", default="baselines")
    ap.add_argument("--n_folds", type=int, default=10)
    ap.add_argument("--datasets", default="lung,hema", help="results/<name> folders, e.g. lung_tfholdout,hema_tfholdout")
    ap.add_argument("--evals", default="test,test_degmatched,test_degcaliper", help="evaluation sets shown in Markdown")
    a = ap.parse_args()
    runs = collect(Path(a.root).expanduser(), set(a.datasets.split(",")))
    evals = [e for e in a.evals.split(",") if e]
    names = dict(EVALS)
    out = {"ref": a.ref, "datasets": a.datasets, "sections": []}
    head = "| Method | Cell types | " + " | ".join(f"AUPR, {names[e]} [95% CI] | diff (p) | AP | AUROC, {names[e]}" for e in evals) + " |"
    lines = [f"Reference for differences and p-values: {a.ref}. Datasets: {a.datasets}.", "", head,
             "|---|---|" + "---|---|---|---|" * len(evals)]

    def cell(r, e):
        if r.get(f"{e}_aupr") is None:
            return ["–", "–", "–", "–"]
        lo, hi = r[f"{e}_aupr_ci"] or (None, None)
        ci = f" [{lo:.3f}, {hi:.3f}]" if lo is not None else ""
        p = r.get(f"{e}_aupr_p")
        d = "–" if r["method"] == a.ref else f"{r[f'{e}_aupr_diff']:+.3f}" + ("" if p is None else f" ({p:.1g})")
        return [f"{r[f'{e}_aupr']:.3f}{ci}", d, f"{r[f'{e}_ap']:.3f}", f"{r[f'{e}_auroc']:.3f}"]

    base_row = None
    for title, items in SECTIONS:
        sec = {"title": title, "rows": []}
        body, shown = [], set()
        for method, label in items:      # aliases share a label: the first one with results is shown
            r = summarise(runs, method, a.ref, a.n_folds)
            r["label"] = label
            sec["rows"].append(r)
            if r["n_cell_types"] and label not in shown:
                shown.add(label)
                body.append(f"| {label} | {r['n_cell_types']} | " + " | ".join(" | ".join(cell(r, e)) for e in evals) + " |")
                if base_row is None and all(r.get(f"{e}_base") is not None for e in evals):
                    base_row = r
        if body:
            lines.append(f"| **{title}** |" + " |" * (1 + 4 * len(evals)))
            lines += body
        out["sections"].append(sec)
    if base_row is not None:
        lines.append("| *Random AUPR (share of positives)* | | " + " | ".join(f"{base_row[f'{e}_base']:.3f} | | {base_row[f'{e}_base']:.3f} | " for e in evals) + "|")
    json.dump(out, open(a.out + ".json", "w"), indent=1)
    Path(a.out + ".md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
