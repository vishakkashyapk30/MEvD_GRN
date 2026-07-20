"""Baseline evaluation + result compilation (plan Step 9)."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import torch

from src.evaluation.metrics import compute_all_metrics


def evaluate_scorer_on_splits(scorer, splits_per_tier: Dict[str, dict],
                              n_total_candidates: int) -> Dict[str, dict]:
    """Score each tier's TEST split with a baseline `scorer` (has .score(tf,tg))
    and compute the same metrics used for MEvD-GRN — identical splits, fair."""
    out: Dict[str, dict] = {}
    for tier, sp in splits_per_tier.items():
        pos, neg = sp["test"]["pos"], sp["test"]["neg"]
        edges = torch.cat([pos, neg], dim=1).numpy()
        labels = np.concatenate([np.ones(pos.shape[1]), np.zeros(neg.shape[1])])
        scores = scorer.score(edges[0], edges[1])
        out[tier] = compute_all_metrics(labels, scores, n_total_candidates)
    return out


def compile_results_table(results_dir: str, pattern: str = "*.json") -> pd.DataFrame:
    rows: List[dict] = []
    for p in sorted(Path(results_dir).rglob(pattern)):
        try:
            with open(p) as f:
                obj = json.load(f)
        except Exception:
            continue
        _flatten(obj, prefix=p.stem, rows=rows, source=str(p))
    return pd.DataFrame(rows)


def _flatten(obj, prefix, rows, source):
    """Flatten {tier: {metric: val}} or {metric: val} into rows."""
    if isinstance(obj, dict) and obj and all(isinstance(v, dict) for v in obj.values()):
        for tier, metrics in obj.items():
            row = {"name": prefix, "tier": tier, "source": source}
            row.update({k: v for k, v in metrics.items()})
            rows.append(row)
    elif isinstance(obj, dict):
        row = {"name": prefix, "tier": "-", "source": source}
        row.update(obj)
        rows.append(row)
