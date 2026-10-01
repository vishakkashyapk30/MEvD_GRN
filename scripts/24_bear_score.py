#!/usr/bin/env python
"""BEAR-GRN step 4: score every GRN (released BEAR methods, MeVD-GRN variants,
our baselines) against every ground truth of a dataset with the BEAR-GRN
metric port (src/benchmarks/bear_metrics.py), and compile tables.

  # port validation: re-score the released outputs and compare with Supp Data 2-5
  python scripts/24_bear_score.py --root $BEAR_ROOT --datasets K562 Macrophage_S1 \
      --released LINGER CellOracle --gts ChIP --paper_csv ~/.cache/bear/paper/bear_grn_auroc_auprc_tidy.csv
  # everything for a dataset (released + results/<ds>/*/seed*/grn.tsv.gz)
  python scripts/24_bear_score.py --root $BEAR_ROOT --datasets K562 --released all --ours all

Outputs: $ROOT/results/<ds>/scores/<method>__<gt>.json (cached), and
$ROOT/results/summary_scores.csv / summary_tables.md.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from src.benchmarks import bear_metrics as bm   # noqa: E402
from src.utils.io import load_config            # noqa: E402

# released method dir -> (ROC/PR columns, early-metric columns); paper method names
PAPER_NAME = {"Pando": "Pando-GLM", "Pando_xgb": "Pando-XGB", "DIRECTNET": "DIRECT-NET"}


def released_files(root: Path, inf_dir: str):
    d = root / "INFERRED.GRNS" / inf_dir
    out = {}
    if not d.exists():
        return out
    for m in sorted(p.name for p in d.iterdir() if p.is_dir()):
        files = sorted(f for f in (d / m).iterdir() if f.is_file())   # R: list.files()[1]
        if files:
            out[m] = files[0]
    return out


def our_files(root: Path, ds: str):
    out = {}
    for f in sorted((root / "results" / ds).glob("*/seed*/grn.tsv.gz")):
        out[f"{f.parent.parent.name}/{f.parent.name}"] = f
    return out


def standardise_released(df: pd.DataFrame, method: str):
    """-> (edges for ROC/PR or None, sorted frame for early metrics)."""
    roc = bm.ROC_PR_COLUMNS.get(method)
    early = bm.EARLY_COLUMNS.get(method)
    e_roc = None
    if roc and roc[2] in df.columns and roc[0] in df.columns:
        e_roc = bm.standardise_edges(df, roc[0], roc[1], roc[2])
    e_early = None
    if early and early[0] in df.columns and early[1] in df.columns:
        d = df.rename(columns={early[0]: "source", early[1]: "target"})
        if early[2] in d.columns:
            d = d.assign(_s=np.abs(pd.to_numeric(d[early[2]], errors="coerce")))
            d = d.sort_values("_s", ascending=False, kind="stable", na_position="last")
        e_early = d[["source", "target"]]
    return e_roc, e_early


def standardise_ours(df: pd.DataFrame, rank_within_fold: bool = False):
    if rank_within_fold and "fold" in df.columns:
        df = df.assign(Score=df.groupby("fold")["Score"].rank(pct=True))
    e = bm.standardise_edges(df, "Source", "Target", "Score")
    d = df.assign(_s=np.abs(df["Score"])).sort_values("_s", ascending=False, kind="stable")
    return e, d.rename(columns={"Source": "source", "Target": "target"})[["source", "target"]]


def score_one(e_roc, e_early, gt, rng, n_rep):
    rec = {}
    if e_roc is not None:
        rec.update(bm.score_roc_pr(e_roc, gt, rng, n_rep=n_rep))
    if e_early is not None:
        rec.update({f"top10k_{k}": v for k, v in bm.score_early(e_early, gt).items()})
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.environ.get("BEAR_ROOT", "data/bear"))
    ap.add_argument("--config", default=str(REPO / "configs/bear/datasets.yaml"))
    ap.add_argument("--datasets", nargs="+", required=True)
    ap.add_argument("--gts", nargs="*", default=None, help="default: all GTs of the dataset")
    ap.add_argument("--released", nargs="*", default=[], help="'all' or method dir names")
    ap.add_argument("--ours", nargs="*", default=[], help="'all' or result names (e.g. mevd_fm_h384/seed42)")
    ap.add_argument("--n_rep", type=int, default=1, help="subsample draws averaged (1 = paper)")
    ap.add_argument("--rng_seed", type=int, default=43, help="R uses set.seed(42 + i)")
    ap.add_argument("--rank_within_fold", action="store_true")
    ap.add_argument("--sparse_match", action="store_true",
                    help="also score our GRNs truncated to the median released in-space edge count")
    ap.add_argument("--paper_csv", default=None, help="tidy Supp Data CSV for a reproduction check")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--compile_only", action="store_true",
                    help="only (re)build summary_tables.md from results/summary_scores.csv")
    args = ap.parse_args()
    if args.compile_only:
        compile_tables(Path(args.root), args.paper_csv)
        return
    root = Path(args.root)
    dcfg = load_config(args.config)
    rows = []
    for ds in args.datasets:
        c = dcfg["datasets"][ds]
        gts = {g: root / p for g, p in c["gts"].items() if (root / p).exists()
               and (args.gts is None or g in args.gts)}
        gt_df = {g: bm.load_gt(str(p)) for g, p in gts.items()}
        sdir = root / "results" / ds / "scores"
        sdir.mkdir(parents=True, exist_ok=True)
        rel = released_files(root, dcfg["inferred_dir"][ds])
        if args.released and args.released != ["all"]:
            rel = {m: f for m, f in rel.items() if m in args.released}
        elif not args.released:
            rel = {}
        ours = our_files(root, ds)
        if args.ours and args.ours != ["all"]:
            ours = {m: f for m, f in ours.items() if m in args.ours}
        elif not args.ours:
            ours = {}
        jobs = [("released", m, f) for m, f in rel.items()] + [("ours", m, f) for m, f in ours.items()]
        median_n = {}
        for kind, m, f in jobs:
            tag = m.replace("/", "__") + ("__rankfold" if (kind == "ours" and args.rank_within_fold) else "")
            def stale(g):
                fp = sdir / f"{tag}__{g}.json"
                if args.force or not fp.exists():
                    return True
                # our GRNs are re-scored if cached with a different number of draws
                return kind == "ours" and json.load(open(fp)).get("n_rep", 1) != args.n_rep
            todo = [g for g in gt_df if stale(g)]
            if todo:
                df = bm.read_table(str(f)) if not str(f).endswith(".gz") else \
                    pd.read_csv(f, sep="\t", keep_default_na=False, na_values=[""])
                if kind == "released":
                    e_roc, e_early = standardise_released(df, m)
                else:
                    e_roc, e_early = standardise_ours(df, args.rank_within_fold)
                del df
                for g in todo:
                    rng = np.random.default_rng(args.rng_seed)
                    rec = {"dataset": ds, "gt": g, "method": PAPER_NAME.get(m, m), "kind": kind,
                           "file": str(f)}
                    rec.update(score_one(e_roc, e_early, gt_df[g], rng, args.n_rep))
                    with open(sdir / f"{tag}__{g}.json", "w") as fh:
                        json.dump(rec, fh, indent=1)
                    print(f"[{ds} {g}] {rec['method']}: AUROC={rec.get('AUROC', float('nan')):.4f} "
                          f"AUPRC={rec.get('AUPRC', float('nan')):.4f} rnd={rec.get('AUPRC_random', float('nan')):.4f} "
                          f"P@10k={rec.get('top10k_Precision', float('nan')):.4f}", flush=True)
                    if kind == "ours" and args.sparse_match:
                        rel_n = [json.load(open(p)).get("n_evaluable_edges") for p in
                                 sdir.glob(f"*__{g}.json") if json.load(open(p)).get("kind") == "released"]
                        rel_n = [x for x in rel_n if x]
                        if rel_n:
                            n = int(np.median(rel_n))
                            tested = (e_roc["tf"].isin(gt_df[g]["Source"]) & e_roc["target"].isin(gt_df[g]["Target"]))
                            top = e_roc[tested].sort_values("score", ascending=False, kind="stable").head(n)
                            rs = {"dataset": ds, "gt": g, "method": PAPER_NAME.get(m, m) + "@sparse",
                                  "kind": "ours_sparse", "file": str(f), "sparse_n": n}
                            rs.update(bm.score_roc_pr(top, gt_df[g], np.random.default_rng(args.rng_seed),
                                                      n_rep=args.n_rep))
                            with open(sdir / f"{tag}__sparse__{g}.json", "w") as fh:
                                json.dump(rs, fh, indent=1)
        for p in sdir.glob("*.json"):
            rows.append(json.load(open(p)))
    if not rows:
        return
    allr = pd.DataFrame(rows)
    out = root / "results" / "summary_scores.csv"
    if out.exists():
        prev = pd.read_csv(out)
        prev = prev[~prev["dataset"].isin(args.datasets)]
        allr = pd.concat([prev, allr], ignore_index=True)
    allr.to_csv(out, index=False)
    print(f"[write] {out} ({len(allr)} rows)")
    if args.paper_csv:
        paper = pd.read_csv(args.paper_csv)
        paper["dataset"] = paper["dataset"].replace({"E7.5_rep1": "mESC_E7.5_rep1", "E7.5_rep2": "mESC_E7.5_rep2",
                                                     "E8.5_rep1": "mESC_E8.5_rep1", "E8.5_rep2": "mESC_E8.5_rep2"})
        mine = allr[allr["kind"] == "released"]
        j = mine.merge(paper, on=["dataset", "gt", "method"], suffixes=("_ours", "_paper"))
        cols = ["dataset", "gt", "method", "AUROC_ours", "AUROC_paper", "AUPRC_ours", "AUPRC_paper",
                "AUPRC_random_ours", "AUPRC_random_paper", "n_evaluable_edges", "evaluable_edges",
                "n_positives", "positives"]
        j = j[[c for c in cols if c in j.columns]]
        j["dAUROC"] = j["AUROC_ours"] - j["AUROC_paper"]
        j["dAUPRC"] = j["AUPRC_ours"] - j["AUPRC_paper"]
        pd.set_option("display.width", 250)
        print(j.round(4).to_string(index=False))
        j.to_csv(root / "results" / "port_validation.csv", index=False)


def _split_method(m: str):
    if "/seed" in m:
        v, sd = m.rsplit("/seed", 1)
        return v, sd.split("@")[0], ("@sparse" if m.endswith("@sparse") else "")
    return m, "", ""


def compile_tables(root: Path, paper_csv: str | None = None) -> None:
    """mean +- sd over seeds per (dataset, gt, method variant); markdown tables."""
    df = pd.read_csv(root / "results" / "summary_scores.csv")
    parts = df["method"].map(_split_method)
    df["variant"] = [a + c for a, b, c in parts]
    df["seed"] = [b for a, b, c in parts]
    metrics = [m for m in ["AUROC", "AUPRC", "AUPRC_random", "top10k_Precision", "top10k_F1"] if m in df]
    g = df.groupby(["dataset", "gt", "kind", "variant"])[metrics]
    agg = g.mean().join(g.std().add_suffix("_sd")).join(g.size().rename("n_seeds")).reset_index()
    agg.to_csv(root / "results" / "summary_agg.csv", index=False)
    lines = ["# MeVD-GRN under BEAR-GRN: compiled scores", "",
             "Mean over seeds (sd in parentheses; n = seeds). Released BEAR methods are",
             "re-scored with the verified port (single draw, as the paper).", ""]
    for gt in sorted(agg["gt"].unique()):
        sub = agg[agg["gt"] == gt]
        for met in ["AUROC", "AUPRC"]:
            if met not in sub:
                continue
            lines += [f"## {gt} - {met}", ""]
            piv = {}
            for _, r in sub.iterrows():
                v = r[met]
                if pd.isna(v):
                    continue
                cell = f"{v:.3f}" if r["n_seeds"] <= 1 or pd.isna(r[f"{met}_sd"]) else \
                    f"{v:.3f} ({r[f'{met}_sd']:.3f}, n={int(r['n_seeds'])})"
                piv.setdefault(r["variant"], {})[r["dataset"]] = cell
            dsets = sorted(sub["dataset"].unique())
            rnd = {d: sub[sub["dataset"] == d]["AUPRC_random"].mean() for d in dsets}
            lines.append("| method | " + " | ".join(dsets) + " |")
            lines.append("|---|" + "---|" * len(dsets))
            if met == "AUPRC":
                lines.append("| *random* | " + " | ".join(f"{rnd[d]:.3f}" for d in dsets) + " |")
            for v in sorted(piv, key=lambda x: (not x.startswith("mevd"), x)):
                lines.append(f"| {v} | " + " | ".join(piv[v].get(d, "") for d in dsets) + " |")
            lines.append("")
    (root / "results" / "summary_tables.md").write_text("\n".join(lines))
    print(f"[write] {root / 'results' / 'summary_tables.md'}")


if __name__ == "__main__":
    main()
