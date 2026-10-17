#!/usr/bin/env python
"""One chart for the whole story: held-out test AUPR (mean over the 12 cell types) of every method and ablation, sorted.
Reads docs/experiments/scmgrn_results/a1_ablation.json (written by scripts/20_scmgrn_ablation_table.py)."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
rows = {r["variant"]: r for r in json.load(open(REPO / "docs/experiments/scmgrn_results/a1_ablation.json")) if r["n_cell_types"]}
spec = [  # (variant key, label, colour group)
    ("mevd_abl_min_noobs", "Stripped-down MeVD-GRN", "floor"),
    ("mevd_abl_fl_atac_noobs", "Stripped-down + ATAC channel", "route"),
    ("mevd_abl_fl_gnn_noobs", "Stripped-down + graph network", "route"),
    ("official_scmultiomegrn", "scMultiomeGRN (official code)", "official"),
    ("degree_product", "No model: product of the two TFs' training-link counts", "trivial"),
    ("degree_logreg", "No model: logistic regression on the two counts", "trivial"),
    ("mevd_abl_fl_fm_noobs", "Stripped-down + Geneformer embedding", "route"),
    ("mevd_fm_h384_obs", "MeVD-GRN", "ours"),
]
col = {"floor": "#8A94A0", "route": "#7fb3d5", "official": "#c8602b", "trivial": "#E69F00", "ours": "#1b6ca8"}
spec = [s for s in spec if s[0] in rows]
spec.sort(key=lambda s: rows[s[0]]["aupr"])
fig, ax = plt.subplots(figsize=(10, 5.0))
for i, (k, lab, g) in enumerate(spec):
    v = rows[k]["aupr"]
    ax.barh(i, v, color=col[g], height=0.72)
    ax.text(v + 0.006, i, f"{v:.3f}", va="center", fontsize=9)
ax.set_yticks(range(len(spec)))
ax.set_yticklabels([s[1] for s in spec], fontsize=9)
ax.set_xlim(0.5, 1.0)
ax.set_xlabel("held-out test AUPR, mean over 12 cell types (9 lung + 3 hematopoiesis)")
ax.spines[["top", "right"]].set_visible(False)
handles = [plt.Rectangle((0, 0), 1, 1, color=col[g]) for g in ("ours", "official", "trivial", "route", "floor")]
fig.legend(handles, ["MeVD-GRN", "scMultiomeGRN", "no-model baselines", "stripped-down + one part", "stripped-down"],
           fontsize=8.5, frameon=False, loc="lower center", ncol=5)
fig.tight_layout(rect=(0, 0.07, 1, 1))
out = REPO / "docs/presentation/figs_scmgrn/fig_ladder.png"
fig.savefig(out, dpi=170)
print("written", out)
