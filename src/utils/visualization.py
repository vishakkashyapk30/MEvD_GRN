"""Plotting helpers for paper figures (plan Step 14). Matplotlib only."""
from __future__ import annotations

from typing import Dict, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import average_precision_score, precision_recall_curve


def plot_pr_curves(method_results: Dict[str, Tuple[np.ndarray, np.ndarray]],
                   tier: str, output_path: str) -> None:
    """method_results: {method_name: (labels, scores)}."""
    plt.figure(figsize=(6, 5))
    for name, (labels, scores) in method_results.items():
        prec, rec, _ = precision_recall_curve(labels, scores)
        ap = average_precision_score(labels, scores)
        plt.plot(rec, prec, label=f"{name} (AUPR={ap:.3f})")
    plt.xlabel("Recall"); plt.ylabel("Precision")
    plt.title(f"PR curves — {tier}")
    plt.legend(loc="upper right", fontsize=8); plt.tight_layout()
    plt.savefig(output_path, dpi=200); plt.close()


def plot_ablation_bars(ablation_results: Dict[str, float], metric: str = "aupr",
                       output_path: str = "results/figures/ablation.pdf",
                       highlight: str = "full_curriculum") -> None:
    names = list(ablation_results.keys())
    vals = [ablation_results[n] for n in names]
    colors = ["#1E407C" if n == highlight else "#9bb4d6" for n in names]
    plt.figure(figsize=(7, 0.5 * len(names) + 1))
    plt.barh(names, vals, color=colors)
    plt.xlabel(metric.upper()); plt.title(f"Ablation — {metric.upper()}")
    plt.tight_layout(); plt.savefig(output_path, dpi=200); plt.close()


def plot_training_curves(history: Dict[str, list], output_path: str) -> None:
    """history: {stage_name: [ {epoch, train_loss, val_aupr?}, ... ]}."""
    fig, ax1 = plt.subplots(figsize=(8, 4))
    ax2 = ax1.twinx()
    off = 0
    for stage, recs in history.items():
        ep = [off + r["epoch"] for r in recs]
        ax1.plot(ep, [r["train_loss"] for r in recs], "-o", ms=3, label=f"{stage} loss")
        va = [(off + r["epoch"], r["val_aupr"]) for r in recs if "val_aupr" in r]
        if va:
            ax2.plot([x for x, _ in va], [y for _, y in va], "--s", ms=3, color="green")
        if recs:
            off += recs[-1]["epoch"]
            ax1.axvline(off + 0.5, color="grey", ls=":", lw=0.8)
    ax1.set_xlabel("epoch"); ax1.set_ylabel("train loss")
    ax2.set_ylabel("val AUPR", color="green")
    ax1.legend(fontsize=7, loc="upper right"); plt.tight_layout()
    plt.savefig(output_path, dpi=200); plt.close()
