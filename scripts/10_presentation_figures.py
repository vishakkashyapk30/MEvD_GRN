#!/usr/bin/env python
"""Generate supplementary figures for presentation.md (beyond 07_compile_results.py's
fig1-3): full per-tier baseline comparison, catastrophic-forgetting AUROC view,
full ablation heatmap, K562 training curves, data/evidence overview, and a
head-to-head MEvD-GRN-vs-scMultiomeGRN comparison across candidate curricula.

Usage:  python scripts/10_presentation_figures.py
Outputs: results/figures/fig4..fig9_*.png
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.utils.io import ensure_dir, load_json

FIG = ensure_dir("results/figures")
TIERS = ["localization", "perturbation", "dual_evidence"]
TIER_LABEL = {"localization": "Localization\n(test)", "perturbation": "Perturbation\n(test)",
              "dual_evidence": "Dual-Evidence\n(zero-shot)"}
BLUE, GREY, RED, GREEN = "#1E407C", "#9bb4d6", "#C0392B", "#2E8B57"


def main():
    df = pd.read_csv("results/summary_table.csv")

    # ---- Fig 4: full per-tier comparison, MEvD-GRN vs all baselines (AUPR) ----
    methods = ["MEvD-GRN", "scmultiomegrn", "grnboost2", "regdiffusion", "gmf_gae"]
    labels = ["MEvD-GRN\n(ours)", "scMultiomeGRN", "GRNBoost2", "RegDiffusion", "GMF-GAE"]
    sub = df[df["method"].isin(methods) & (df["cell_type"] == "K562")]
    piv = sub.pivot_table(index="method", columns="tier", values="aupr").reindex(methods)[TIERS]
    piv.index = labels
    fig, ax = plt.subplots(figsize=(9, 5))
    piv.plot(kind="bar", ax=ax, color=[BLUE, "#5A8FC7", GREY],
              edgecolor="white", width=0.75)
    ax.set_ylabel("AUPR")
    ax.set_title("MEvD-GRN vs. all baselines, per evidence tier (K562)")
    ax.legend([TIER_LABEL[t].replace("\n", " ") for t in TIERS], title=None, fontsize=9)
    ax.set_xticklabels(labels, rotation=15, ha="right")
    fig.tight_layout()
    fig.savefig(FIG / "fig4_full_baseline_comparison.png", dpi=200)
    plt.close(fig)

    # ---- Fig 5: catastrophic forgetting -- AUROC by tier for the main model ----
    main_df = df[(df["method"] == "MEvD-GRN") & (df["cell_type"] == "K562")]
    piv2 = main_df.set_index("tier").reindex(TIERS)["auroc"]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    colors = [RED if v < 0.5 else BLUE for v in piv2.values]
    bars = ax.bar([TIER_LABEL[t] for t in TIERS], piv2.values, color=colors, edgecolor="white")
    ax.axhline(0.5, color="black", linestyle="--", linewidth=1, label="random (AUROC=0.5)")
    for b, v in zip(bars, piv2.values):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.02 if v >= 0 else v - 0.05,
                f"{v:.3f}", ha="center", fontsize=10, fontweight="bold")
    ax.set_ylabel("AUROC")
    ax.set_ylim(0, 1.05)
    ax.set_title("MEvD-GRN default (2-stage, no dual training):\ncatastrophic forgetting of localization")
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "fig5_catastrophic_forgetting_auroc.png", dpi=200)
    plt.close(fig)

    # ---- Fig 6: full ablation heatmap (AUPR, all tiers) ----
    abl = df[df["method"].str.startswith("abl:")].copy()
    abl["ablation"] = abl["method"].str.replace("abl:", "", regex=False)
    order = ["full_curriculum", "loc_only", "pert_only", "dual_only", "all_at_once",
             "rna_only", "gated_fusion", "concat_fusion", "no_gnn", "with_replay"]
    piv3 = abl.pivot_table(index="ablation", columns="tier", values="aupr").reindex(order)[TIERS]
    fig, ax = plt.subplots(figsize=(7.5, 6))
    im = ax.imshow(piv3.values, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(TIERS)))
    ax.set_xticklabels([TIER_LABEL[t].replace("\n", " ") for t in TIERS], fontsize=9)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=10)
    for i in range(len(order)):
        for j in range(len(TIERS)):
            v = piv3.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.3f}", ha="center", va="center", fontsize=9,
                         color="black" if 0.3 < v < 0.75 else "white")
    ax.set_title("MEvD-GRN ablations \u2014 AUPR per tier (K562)")
    fig.colorbar(im, ax=ax, label="AUPR", fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(FIG / "fig6_ablation_heatmap.png", dpi=200)
    plt.close(fig)

    # ---- Fig 7: K562 training curves (Stage1 localization, Stage2 perturbation) ----
    stage1 = {"epoch": [5, 10, 15, 20, 25, 30],
              "val_aupr": [0.9293, 0.9345, 0.9371, 0.9381, 0.9391, 0.9396],
              "val_auroc": [0.9280, 0.9337, 0.9361, 0.9368, 0.9378, 0.9383],
              "loss": [0.3919, 0.3656, 0.3575, 0.3524, 0.3493, 0.3481]}
    stage2 = {"epoch": [5, 10, 15],
              "val_aupr": [0.5541, 0.5743, 0.5779],
              "val_auroc": [0.8662, 0.8755, 0.8764],
              "loss": [0.5893, 0.5538, 0.5469]}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=False)
    for ax, stage, name, lr in zip(axes, [stage1, stage2],
                                    ["Stage 1: Localization", "Stage 2: Perturbation"],
                                    ["lr=1e-3, 30 epochs", "lr=3e-4, 15 epochs"]):
        ax.plot(stage["epoch"], stage["val_aupr"], "o-", color=BLUE, label="val AUPR")
        ax.plot(stage["epoch"], stage["val_auroc"], "s-", color=GREEN, label="val AUROC")
        ax2 = ax.twinx()
        ax2.plot(stage["epoch"], stage["loss"], "^--", color=GREY, label="train loss", alpha=0.8)
        ax2.set_ylabel("train loss", color=GREY)
        ax.set_xlabel("epoch")
        ax.set_ylabel("validation score")
        ax.set_title(f"{name}\n({lr})")
        ax.set_ylim(0, 1.02)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc="center right", fontsize=8)
    fig.suptitle("MEvD-GRN main model training curves (K562, 2-stage curriculum + replay)")
    fig.tight_layout()
    fig.savefig(FIG / "fig7_training_curves.png", dpi=200)
    plt.close(fig)

    # ---- Fig 8: evidence tier sizes + nesting (data overview) ----
    summ = load_json("data/processed/K562/summary.json")
    sizes = summ["evidence_sizes"]
    nesting = summ["nesting"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    order2 = ["localization", "perturbation", "dual_evidence"]
    vals = [sizes[t] for t in order2]
    bars = ax.bar(["Localization", "Perturbation", "Dual-Evidence"], vals,
                  color=[BLUE, "#5A8FC7", GREY], edgecolor="white")
    ax.set_yscale("log")
    ax.set_ylim(top=max(vals) * 3)
    ax.set_ylabel("# positive TF\u2192gene edges (log scale)")
    ax.set_title("K562 evidence tier sizes\n(after gene-universe alignment)")
    ax.tick_params(axis="x", labelsize=9)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v * 1.3, f"{v:,}", ha="center", fontsize=9)

    ax = axes[1]
    keys = list(nesting.keys())
    nvals = [nesting[k] * 100 for k in keys]
    nice = [k.replace("\u2286", " \u2286 ") for k in keys]
    bars = ax.barh(nice, nvals, color=RED, edgecolor="white")
    ax.set_xlabel("% of edges also present in the other tier")
    ax.set_xlim(0, 105)
    ax.set_title("Tier nesting (the leakage hazard)")
    for b, v in zip(bars, nvals):
        ax.text(v + 2, b.get_y() + b.get_height() / 2, f"{v:.1f}%", va="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "fig8_data_overview_nesting.png", dpi=200)
    plt.close(fig)

    # ---- Fig 9: MEvD-GRN (3 configs) vs scMultiomeGRN, head-to-head, AUPR ----
    configs = {
        "scMultiomeGRN\n(4.08M params/tier)": "scmultiomegrn",
        "MEvD-GRN default\n(2-stage, shipped)": "MEvD-GRN",
        "MEvD-GRN all_at_once\n(candidate)": "abl:all_at_once",
        "MEvD-GRN with_replay\n(candidate)": "abl:with_replay",
    }
    rows = []
    for label, key in configs.items():
        sub = df[(df["method"] == key) & (df["cell_type"] == "K562")]
        piv = sub.set_index("tier").reindex(TIERS)["aupr"]
        rows.append(piv.values)
    mat = np.array(rows)  # (4 configs, 3 tiers)
    x = np.arange(len(TIERS))
    width = 0.2
    fig, ax = plt.subplots(figsize=(9.5, 5))
    colors = ["#888888", "#1E407C", "#5A8FC7", "#2E8B57"]
    for i, (label, _) in enumerate(configs.items()):
        bars = ax.bar(x + (i - 1.5) * width, mat[i], width, label=label, color=colors[i],
                       edgecolor="white")
        for b, v in zip(bars, mat[i]):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.015, f"{v:.2f}", ha="center", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([TIER_LABEL[t].replace("\n", " ") for t in TIERS])
    ax.set_ylabel("AUPR")
    ax.set_ylim(0, 1.08)
    ax.set_title("MEvD-GRN vs. scMultiomeGRN: head-to-head across candidate curricula (K562)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "fig9_mevdgrn_vs_scmultiomegrn.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote fig4-fig9 to {FIG}/")


if __name__ == "__main__":
    main()
