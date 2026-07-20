#!/usr/bin/env python
"""Compile all result JSONs into a CSV table + paper figures (plan Step 14-15).

Usage:  python scripts/07_compile_results.py
Outputs: results/summary_table.csv and results/figures/*.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.utils.io import ensure_dir, load_json

FIG = ensure_dir("results/figures")
TIERS = ["localization", "perturbation", "dual_evidence"]
BLUE, GREY = "#1E407C", "#9bb4d6"


def _rows_from_test(name, test_dict):
    return [{"method": name, "tier": t, **{k: test_dict[t].get(k) for k in
            ("aupr", "auroc", "early_precision", "epr")}}
            for t in test_dict if isinstance(test_dict[t], dict)]


def main():
    rows = []
    # MEvD-GRN main results
    for ct in ("K562", "ESC"):
        p = Path(f"results/{ct}_results.json")
        if p.exists():
            rows += [{"cell_type": ct, **r} for r in
                     _rows_from_test(f"MEvD-GRN", load_json(p)["test"])]
    # transfer
    for p in Path("results").glob("*_transfer.json"):
        obj = load_json(p)
        rows += [{"cell_type": p.stem, **r} for r in _rows_from_test("MEvD-GRN(transfer)", obj)]
    # baselines
    for p in Path("results/baselines").glob("*.json"):
        obj = load_json(p)
        if "results" in obj:
            rows += [{"cell_type": obj.get("cell_type", "?"), **r}
                     for r in _rows_from_test(obj["baseline"], obj["results"])]
    # ablations (K562)
    abl = {}
    for p in Path("results/ablations").glob("*.json"):
        obj = load_json(p)
        res = obj.get("results", {})
        key = obj.get("key_tier", "dual_evidence")
        if key in res:
            abl[obj["ablation"]] = res[key]["aupr"]
        rows += [{"cell_type": obj.get("cell_type", "?"),
                  "method": f"abl:{obj['ablation']}", **r}
                 for r in _rows_from_test(f"abl:{obj['ablation']}", res)]

    df = pd.DataFrame(rows)
    if df.empty:
        print("No results found yet.")
        return
    df.to_csv("results/summary_table.csv", index=False)
    print("Wrote results/summary_table.csv\n")
    with pd.option_context("display.max_rows", None, "display.width", 160):
        print(df.round(4).to_string(index=False))

    # Figure 1: MEvD-GRN per-tier AUPR (K562 vs ESC)
    main_df = df[df["method"] == "MEvD-GRN"]
    if not main_df.empty:
        piv = main_df.pivot_table(index="tier", columns="cell_type", values="aupr").reindex(TIERS)
        piv.plot(kind="bar", figsize=(7, 4), color=[BLUE, GREY])
        plt.ylabel("AUPR"); plt.title("MEvD-GRN AUPR by evidence tier")
        plt.xticks(rotation=15); plt.tight_layout()
        plt.savefig(FIG / "fig1_mevdgrn_per_tier.png", dpi=200); plt.close()

    # Figure 2: ablation bars (dual-evidence AUPR)
    if abl:
        order = sorted(abl, key=lambda k: -abl[k])
        colors = [BLUE if k == "full_curriculum" else GREY for k in order]
        plt.figure(figsize=(7, 4))
        plt.barh(order, [abl[k] for k in order], color=colors)
        plt.gca().invert_yaxis(); plt.xlabel("dual-evidence AUPR")
        plt.title("Ablations (K562)"); plt.tight_layout()
        plt.savefig(FIG / "fig2_ablation_dual_aupr.png", dpi=200); plt.close()

    # Figure 3: MEvD-GRN vs baselines on dual-evidence (per cell type)
    dual = df[(df["tier"] == "dual_evidence") & (~df["method"].str.startswith("abl:"))]
    if not dual.empty:
        piv = dual.pivot_table(index="method", columns="cell_type", values="aupr")
        piv.plot(kind="bar", figsize=(8, 4))
        plt.ylabel("dual-evidence AUPR"); plt.title("MEvD-GRN vs baselines (dual-evidence)")
        plt.xticks(rotation=25, ha="right"); plt.tight_layout()
        plt.savefig(FIG / "fig3_vs_baselines_dual.png", dpi=200); plt.close()

    print(f"\nFigures written to {FIG}/")


if __name__ == "__main__":
    main()
