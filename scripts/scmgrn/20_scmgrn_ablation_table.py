#!/usr/bin/env python
"""Ablation table for MeVD-GRN on the scMultiomeGRN benchmark (lung + hematopoiesis).

Reads the fold-level result JSONs written by scripts/17 (and 18 for the official re-run) under
<root>/results/<dataset>/<cell type>/<variant>/seed*/fold*.json. For every variant: the held-out-test AUROC / AUPR,
averaged over folds and seeds within a cell type, then over the cell types where BOTH the variant and the
headline variant H are complete; the mean difference to H; and a fold-paired Wilcoxon test against H
(per-fold values averaged over seeds, all cell-type x fold pairs). Rows are one change each (see configs/scmgrn/mevd_abl_*.yaml
and the variants A1-A3).

Usage: python scripts/scmgrn/20_scmgrn_ablation_table.py --root ~/scmgrn --out ablation
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

HEAD = "mevd_fm_h384_obs_smg"
ORDER = [
    (HEAD, "H: headline (FM + 384-wide GNN, undirected, training graph, scMultiomeGRN inputs)"),
    ("mevd_fm_h384_obs_smgnode", "only the scMultiomeGRN node features (no pair histograms)"),
    ("mevd_fm_h384_obs_smgpair", "only the pair histograms (no scMultiomeGRN node features)"),
    ("mevd_fm_h384_obs", "A1: no scMultiomeGRN inputs"),
    ("mevd_fm_h384_noobs", "A2: no scMultiomeGRN inputs, no training graph"),
    ("mevd_fm_h384_directed_obs", "A3: old directed model, scores averaged"),
    ("mevd_abl_nofm_obs_smg", "no Geneformer gene embedding"),
    ("mevd_abl_noatac_obs_smg", "no ATAC channel"),
    ("mevd_abl_nognn_obs_smg", "no graph network (features straight to the decoder)"),
    ("mevd_abl_coexpr_obs_smg", "co-expression graph only (also drops the training links)"),
    ("mevd_abl_tfcand_obs_smg", "TF-candidate graph only (also drops the training links)"),
    ("mevd_abl_base_obs_smg", "base size (128-wide, no Geneformer)"),
    ("mevd_abl_min_noobs", "floor: 128-wide, no Geneformer, no ATAC channel, no graph network, no training graph, no scMultiomeGRN inputs"),
    ("mevd_abl_fl_atac_noobs", "floor + only the ATAC channel"),
    ("mevd_abl_fl_fm_noobs", "floor + only the Geneformer embedding"),
    ("mevd_abl_fl_gnn_noobs", "floor + only the graph network (label-free graphs)"),
    ("mevd_abl_min_noobs_smgnode", "floor + only the scMultiomeGRN node features"),
    ("official_scmultiomegrn", "scMultiomeGRN, official code (reference)"),
    ("degree_product", "trivial baseline: product of the two TFs' training-link counts (no model, no features)"),
    ("degree_logreg", "trivial baseline: logistic regression on the two TFs' training-link counts"),
]
METRICS = ("aupr", "auroc")


def collect(root: Path):
    """runs[(dataset, ct)][variant][fold] -> {metric: [per-seed values]}"""
    runs = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
    for p in root.glob("results/*/*/*/seed*/fold*.json"):
        ds, ct, variant = p.parts[-5], p.parts[-4], p.parts[-3]      # results/<ds>/<ct>/<variant>/seed*/fold*.json
        fold = int(p.stem.replace("fold", ""))
        d = json.load(open(p))
        if variant.startswith("mevd_") and not (d.get("leak_status") or {}).get("zero_leak", False):
            raise SystemExit(f"{p}: not zero-leak")
        for m in METRICS:
            runs[(ds, ct)][variant][fold][m].append(d["metrics"]["test"][m])
    return runs


def main():
    global HEAD
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", default="ablation")
    ap.add_argument("--n_folds", type=int, default=10)
    ap.add_argument("--datasets", default="lung,hema")
    ap.add_argument("--head", default=HEAD, help="reference variant that differences and tests are computed against")
    a = ap.parse_args()
    HEAD = a.head
    keep = set(a.datasets.split(","))
    runs = {k: v for k, v in collect(Path(a.root)).items() if k[0] in keep}

    def complete(by_fold):
        return len(by_fold) >= a.n_folds

    rows = []
    for variant, label in ORDER:
        per_ct = {}
        diffs = {m: [] for m in METRICS}
        n_ct = 0
        for key, by_variant in runs.items():
            if variant not in by_variant or HEAD not in by_variant:
                continue
            if not (complete(by_variant[variant]) and complete(by_variant[HEAD])):
                continue
            n_ct += 1
            for m in METRICS:
                v = [np.mean(by_variant[variant][f][m]) for f in sorted(by_variant[variant])][: a.n_folds]
                h = [np.mean(by_variant[HEAD][f][m]) for f in sorted(by_variant[HEAD])][: a.n_folds]
                per_ct.setdefault(m, []).append((np.mean(v), np.mean(h)))
                diffs[m] += [x - y for x, y in zip(v, h)]
        if not n_ct:
            rows.append({"variant": variant, "label": label, "n_cell_types": 0})
            continue
        row = {"variant": variant, "label": label, "n_cell_types": n_ct}
        for m in METRICS:
            row[m] = float(np.mean([x for x, _ in per_ct[m]]))
            row[f"{m}_diff_vs_H"] = float(np.mean(diffs[m]))
            dd = np.array(diffs[m])
            row[f"{m}_p"] = float(wilcoxon(dd).pvalue) if variant != HEAD and len(dd) >= 5 and np.any(dd != 0) else None
            row[f"{m}_frac_worse_than_H"] = float(np.mean(dd < 0)) if variant != HEAD else None
        rows.append(row)

    out = Path(a.out)
    json.dump(rows, open(str(out) + ".json", "w"), indent=1)
    lines = ["| Variant | Cell types | test AUPR | diff vs H | p | test AUROC | diff vs H | p |", "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if not r["n_cell_types"]:
            lines.append(f"| {r['label']} | 0 | – | – | – | – | – | – |")
            continue
        def pv(x):
            return "–" if x is None else f"{x:.2g}"
        def df(x, base):
            return "–" if base else f"{x:+.4f}"
        is_h = r["variant"] == HEAD
        lines.append(f"| {r['label']} | {r['n_cell_types']} | {r['aupr']:.4f} | {df(r['aupr_diff_vs_H'], is_h)} | {pv(r['aupr_p'])} | "
                     f"{r['auroc']:.4f} | {df(r['auroc_diff_vs_H'], is_h)} | {pv(r['auroc_p'])} |")
    Path(str(out) + ".md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
