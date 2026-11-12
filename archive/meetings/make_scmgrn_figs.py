#!/usr/bin/env python
"""Figures for the MeVD-GRN vs scMultiomeGRN apples-to-apples deck, from scripts/19's summary.json.

  fig_test.png       headline: held-out test fold, AUPR and AUROC per cell type, MeVD-GRN H vs the
                     official scMultiomeGRN code re-run (zero leakage, same folds)
  fig_ablation.png   variants H / A1 / A2 / A3 on the same metric (test AUPR, mean over finished cell types)
  fig_paper.png      paper protocol ('all' = full matrix, contains training edges): our official re-run
                     vs the numbers published in the paper (fidelity check) and MeVD-GRN H

Usage: python docs/presentation/make_scmgrn_figs.py --summary <summary.json> --config configs/scmgrn/lung.yaml
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.utils.io import load_config  # noqa: E402

OFF, HEAD = "official_scmultiomegrn", "mevd_fm_h384_obs"
LAB = {HEAD: "MeVD-GRN", OFF: "scMultiomeGRN (official code, re-run)"}
COL = {HEAD: "#1b6ca8", OFF: "#c8602b"}


def short(ct):
    return ct.replace("_cells", "").replace("_and_alveolar_epithelial", "/alveolar").replace("_", " ")


def grouped(ax, cts, series, key, title):
    w = 0.8 / max(len(series), 1)
    for i, m in enumerate(series):
        v = [S[m][c][key] if c in S[m] else np.nan for c in cts]
        sd = [S[m][c].get(key + "_fold_sd", 0) if c in S[m] else 0 for c in cts]
        ax.bar(np.arange(len(cts)) + i * w - 0.4 + w / 2, v, w, yerr=sd, color=COL.get(m, "#555"), label=LAB.get(m, m),
               error_kw=dict(lw=0.8, capsize=2))
    ax.set_xticks(range(len(cts))); ax.set_xticklabels([short(c) for c in cts], rotation=30, ha="right", fontsize=8)
    ax.set_title(title, fontsize=10, loc="left"); ax.spines[["top", "right"]].set_visible(False)


ap = argparse.ArgumentParser()
ap.add_argument("--summary", required=True); ap.add_argument("--config", required=True)
ap.add_argument("--out", default=str(Path(__file__).parent / "figs_scmgrn"))
a = ap.parse_args()
S = json.load(open(a.summary)); S.pop("_paired_vs_ref", None)
cfg = load_config(a.config); out = Path(a.out); out.mkdir(exist_ok=True, parents=True)
cts = [c for c in cfg["cell_types"] if all(c in S.get(m, {}) for m in (HEAD, OFF) if m in S)]
cts = [c for c in cts if c in S.get(HEAD, {}) or c in S.get(OFF, {})]

fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
for ax, k, t in zip(axs, ("test_aupr", "test_auroc"), ("Held-out test fold: AUPR", "Held-out test fold: AUROC")):
    grouped(ax, cts, [m for m in (HEAD, OFF) if m in S], k, t); ax.set_ylim(0.4 if "auroc" in k else 0.3, 1.0)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, fontsize=9, frameon=False, loc="upper center", ncol=2)
fig.text(0.5, 0.01, "Same datasets, same 10 official folds, same metric code, zero leakage (bars: mean over folds, whiskers: fold sd)",
         fontsize=8.5, ha="center", va="bottom")
fig.tight_layout(rect=(0, 0.07, 1, 0.92)); fig.savefig(out / "fig_test.png", dpi=170); plt.close(fig)

fig, ax = plt.subplots(figsize=(6, 3.6)); ms = [m for m in LAB if m in S and m != OFF]
for i, m in enumerate(ms + ([OFF] if OFF in S else [])):
    vals = [S[m][c]["test_aupr"] for c in cts if c in S[m]]
    if vals: ax.barh(i, np.mean(vals), color=COL.get(m, "#555")); ax.text(np.mean(vals) + .005, i, f"{np.mean(vals):.3f}", va="center", fontsize=8)
ax.set_yticks(range(len(ms) + (OFF in S))); ax.set_yticklabels([LAB[m] for m in ms + ([OFF] if OFF in S else [])], fontsize=8)
ax.invert_yaxis(); ax.set_xlabel(f"test AUPR (mean over {len(cts)} cell types)"); ax.spines[["top", "right"]].set_visible(False)
ax.set_xlim(0.3, 1.0); fig.tight_layout(); fig.savefig(out / "fig_ablation.png", dpi=170); plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
for ax, k, pk, t in zip(axs, ("all_auroc", "all_aupr"), ("auroc", "aupr"),
                        ("Paper protocol (contains training edges): AUROC", "Paper protocol (contains training edges): AUPR")):
    paper = [cfg["cell_types"][c].get(pk, np.nan) for c in cts]
    ax.bar(np.arange(len(cts)) - 0.27, paper, 0.27, color="#d9d9d9", label="published (paper)")
    for j, m in enumerate([OFF, HEAD]):
        if m in S:
            ax.bar(np.arange(len(cts)) - 0.27 + 0.27 * (j + 1), [S[m][c][k] if c in S[m] else np.nan for c in cts], 0.27,
                   color=COL[m], label=LAB[m])
    ax.set_xticks(range(len(cts))); ax.set_xticklabels([short(c) for c in cts], rotation=30, ha="right", fontsize=8)
    ax.set_title(t, fontsize=10, loc="left"); ax.set_ylim(0.3 if "aupr" in k else 0.6, 1.0); ax.spines[["top", "right"]].set_visible(False)
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, fontsize=9, frameon=False, loc="upper center", ncol=3)
fig.tight_layout(rect=(0, 0, 1, 0.92)); fig.savefig(out / "fig_paper.png", dpi=170); plt.close(fig)
print("figures ->", out, "| cell types:", cts)
