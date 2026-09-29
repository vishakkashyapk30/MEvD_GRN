#!/usr/bin/env python
"""Aggregate scMultiomeGRN-benchmark results and compare with the paper.

Per method/variant and cell type: mean over the 10 folds for each seed, then
mean +/- std across seeds (fold-std also reported, which is what the paper's
"+/-" is). Lung average = mean over the 9 cell types. Both evaluation sets are
tabulated: "all" (full N x N matrix -- the set the paper's AUROC/AUPR come
from) and "test" (held-out fold, 1:1, leakage-free).

Usage:
  python scripts/19_scmgrn_compile.py --config configs/scmgrn/lung.yaml --root $SCMGRN_ROOT
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from src.utils.io import load_config, save_json  # noqa: E402

METRICS = ("auroc", "aupr", "acc")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--root", default=os.environ.get("SCMGRN_ROOT", "data/scmgrn"))
    ap.add_argument("--out", default=None, help="default <root>/<results_dir>/summary.{json,md}")
    ap.add_argument("--n_folds", type=int, default=10)
    args = ap.parse_args()
    dcfg = load_config(args.config)
    rdir = Path(args.root) / dcfg["results_dir"]
    cts = list(dcfg["cell_types"])

    # runs[method][ct][seed] = list of fold records
    runs = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    for p in rdir.glob("*/*/seed*/fold*.json"):
        ct, method, seed = p.parts[-4], p.parts[-3], p.parts[-2]
        runs[method][ct][seed].append(json.load(open(p)))

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

    lines = [f"# scMultiomeGRN benchmark: {dcfg['dataset']}", "",
             f"AUROC / AUPR / Acc, mean over {args.n_folds} folds (then over seeds). 'all' = full-matrix set "
             "(what the paper reports); 'test' = held-out fold 1:1 (leakage-free).", ""]
    hdr = "| Cell type | Paper scMultiomeGRN (all) | " + " | ".join(f"{m} all | {m} test" for m in sorted(summary)) + " |"
    lines += [hdr, "|" + "---|" * (2 + 2 * len(summary))]

    def fmt(d, es):
        if not d:
            return "–"
        return f"{d[f'{es}_auroc']:.3f} / {d[f'{es}_aupr']:.3f} / {d[f'{es}_acc']:.3f}"
    for ct in cts + ["_average"]:
        t = dcfg["cell_types"].get(ct, {}) if ct != "_average" else (dcfg.get("paper_average") or {})
        paper = (f"{t.get('auroc', float('nan')):.3f} / {t.get('aupr', float('nan')):.3f} / "
                 f"{t['acc']:.3f}" if t.get("acc") else
                 (f"{t['auroc']:.3f} / {t['aupr']:.3f}" if t.get("auroc") else "–"))
        row = [ct.replace("_", " ").strip(), paper]
        for m in sorted(summary):
            d = summary[m].get(ct)
            row += [fmt(d, "all"), fmt(d, "test")]
        lines.append("| " + " | ".join(row) + " |")
    md = "\n".join(lines) + "\n"
    out = Path(args.out) if args.out else rdir / "summary"
    save_json(summary, str(out) + ".json")
    Path(str(out) + ".md").write_text(md)
    print(md)
    print(f"-> {out}.json / {out}.md")


if __name__ == "__main__":
    main()
