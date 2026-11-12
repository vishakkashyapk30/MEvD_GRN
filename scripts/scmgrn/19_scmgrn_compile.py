#!/usr/bin/env python
"""Aggregate scMultiomeGRN-benchmark results and compare with the paper.

Per method/variant and cell type: mean over the 10 folds for each seed, then mean +/- std across
seeds. Lung average = mean over the 9 cell types. Two evaluation sets:
  "test" : the held-out fold (1:1 pos/neg), no training edge involved -> the HEADLINE (zero leakage)
  "all"  : the paper's protocol, the full N x N matrix, which CONTAINS the training edges (the
           message-passing graph is the train adjacency); reported only to compare with the
           published 0.924 / 0.790 and labelled as such.
Also: fold-paired Wilcoxon tests of the headline MeVD-GRN variant against the official
scMultiomeGRN re-run (fold values averaged over seeds, 10 folds x 9 cell types), and a hard
failure if any MeVD-GRN fold record is not zero-leak.

Usage:
  python scripts/scmgrn/19_scmgrn_compile.py --config configs/scmgrn/lung.yaml --root $SCMGRN_ROOT
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.utils.io import load_config, save_json  # noqa: E402

METRICS = ("auroc", "aupr", "acc")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--root", default=os.environ.get("SCMGRN_ROOT", "data/scmgrn"))
    ap.add_argument("--out", default=None, help="default <root>/<results_dir>/summary.{json,md}")
    ap.add_argument("--n_folds", type=int, default=10)
    ap.add_argument("--head", default="mevd_fm_h384_obs_smg", help="headline MeVD-GRN variant dir name")
    ap.add_argument("--ref", default="official_scmultiomegrn", help="reference (official code re-run) dir name")
    args = ap.parse_args()
    dcfg = load_config(args.config)
    rdir = Path(args.root) / dcfg["results_dir"]
    cts = list(dcfg["cell_types"])

    # runs[method][ct][seed] = list of fold records
    runs = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for p in rdir.glob("*/*/seed*/fold*.json"):
        ct, method, seed = p.parts[-4], p.parts[-3], p.parts[-2]
        rec = json.load(open(p))
        if method.startswith("mevd_") and not (rec.get("leak_status") or {}).get("zero_leak", False):
            raise SystemExit(f"{p}: MeVD-GRN fold is not zero-leak ({rec.get('leak_status')})")
        runs[method][ct][seed].append(rec)

    summary = {}
    for method, by_ct in sorted(runs.items()):
        summary[method] = {}
        for ct in cts:
            seeds = by_ct.get(ct, {})
            per_seed = {es: {m: [] for m in METRICS} for es in ("all", "test")}
            fold_sd = {es: {m: [] for m in METRICS} for es in ("all", "test")}
            complete = 0
            for sd, recs in seeds.items():
                if len(recs) < args.n_folds:
                    continue
                complete += 1
                for es in ("all", "test"):
                    for m in METRICS:
                        v = [r["metrics"][es][m] for r in recs]
                        per_seed[es][m].append(float(np.mean(v)))
                        fold_sd[es][m].append(float(np.std(v)))
            if not complete:
                continue
            summary[method][ct] = {
                "n_seeds_complete": complete,
                **{f"{es}_{m}": float(np.mean(per_seed[es][m])) for es in ("all", "test") for m in METRICS},
                **{f"{es}_{m}_seed_sd": float(np.std(per_seed[es][m])) for es in ("all", "test") for m in METRICS},
                **{f"{es}_{m}_fold_sd": float(np.mean(fold_sd[es][m])) for es in ("all", "test") for m in METRICS},
            }
        done = [ct for ct in cts if ct in summary[method]]
        if len(done) == len(cts):
            summary[method]["_average"] = {f"{es}_{m}": float(np.mean([summary[method][c][f"{es}_{m}"] for c in cts]))
                                           for es in ("all", "test") for m in METRICS}

    NAMES = {"mevd_fm_h384_obs_smg": "MeVD-GRN H (undirected, +train adj, +smg inputs)",
             "mevd_fm_h384_obs": "A1 (no smg inputs)", "mevd_fm_h384_noobs": "A2 (label-free graphs only)",
             "mevd_fm_h384_directed_obs": "A3 (old directed model, symmetrised)",
             "official_scmultiomegrn": "scMultiomeGRN official code (re-run)"}
    order = [m for m in NAMES if m in summary] + [m for m in sorted(summary) if m not in NAMES]

    def fmt(d, es):
        if not d:
            return "–"
        return f"{d[f'{es}_auroc']:.3f} / {d[f'{es}_aupr']:.3f}"

    def table(es, title, note):
        out = [f"## {title}", "", note, "", "AUROC / AUPR per cell type, mean over 10 folds then over seeds.", ""]
        out += ["| Cell type | " + " | ".join(NAMES.get(m, m) for m in order) + " |", "|---|" + "---|" * len(order)]
        for ct in cts + ["_average"]:
            out.append("| " + ct.replace("_", " ").strip() + " | " +
                       " | ".join(fmt(summary[m].get(ct), es) for m in order) + " |")
        return out + [""]

    lines = [f"# scMultiomeGRN benchmark: {dcfg['dataset']}", ""]
    lines += table("test", "Headline: held-out test fold (zero leakage, same for both methods)",
                   "1:1 positives/negatives of the official fold; training edges are never scored.")
    pa = dcfg.get("paper_average") or {}
    lines += table("all", "Paper protocol: full N x N matrix (CONTAINS the training edges; not leakage-free)",
                   f"Only for comparison with the published numbers (paper average AUROC {pa.get('auroc')} / "
                   f"AUPR {pa.get('aupr')}, scMultiomeGRN; baselines in lung.yaml).")
    # fold-paired Wilcoxon: head vs ref
    stats = {}
    if args.head in runs and args.ref in runs:
        for es in ("test", "all"):
            for m in ("aupr", "auroc"):
                diffs, per_ct = [], {}
                for ct in cts:
                    def foldvals(method):
                        by_fold = defaultdict(list)
                        for recs in runs[method].get(ct, {}).values():
                            for r in recs:
                                by_fold[r["fold"]].append(r["metrics"][es][m])
                        return {f: float(np.mean(v)) for f, v in by_fold.items()}
                    h, r_ = foldvals(args.head), foldvals(args.ref)
                    common = sorted(set(h) & set(r_))
                    d = [h[f] - r_[f] for f in common]
                    if len(d) >= 5 and any(x != 0 for x in d):
                        per_ct[ct] = {"mean_diff": float(np.mean(d)), "n_folds": len(d),
                                      "wilcoxon_p": float(wilcoxon(d).pvalue)}
                    diffs += d
                if len(diffs) >= 5 and any(x != 0 for x in diffs):
                    stats[f"{es}_{m}"] = {"mean_diff": float(np.mean(diffs)), "n_pairs": len(diffs),
                                          "wilcoxon_p": float(wilcoxon(diffs).pvalue),
                                          "frac_head_better": float(np.mean(np.array(diffs) > 0)),
                                          "per_cell_type": per_ct}
    if stats:
        lines += [f"## Fold-paired comparison: {args.head} minus {args.ref}", "",
                  "| set / metric | mean diff | frac. folds better | Wilcoxon p | n pairs |", "|---|---|---|---|---|"]
        for k, v in stats.items():
            lines.append(f"| {k} | {v['mean_diff']:+.4f} | {v['frac_head_better']:.2f} | {v['wilcoxon_p']:.2g} | {v['n_pairs']} |")
        lines.append("")
    summary["_paired_vs_ref"] = stats
    md = "\n".join(lines) + "\n"
    out = Path(args.out) if args.out else rdir / "summary"
    save_json(summary, str(out) + ".json")
    Path(str(out) + ".md").write_text(md)
    print(md)
    print(f"-> {out}.json / {out}.md")


if __name__ == "__main__":
    main()
