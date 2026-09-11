"""Baseline evaluation + result compilation (plan Step 9)."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd
import torch

from src.evaluation.metrics import compute_all_metrics


def eval_edges_for_tier(sp: dict, tier: str, held_out_tiers: Optional[Iterable[str]] = None):
    """Pick the evaluation edge set for one tier.

    Tiers that are actually TRAINED on (the default: localization,
    perturbation) use the standard TEST split only. Tiers held out of
    training entirely (default: dual_evidence, since it is nested inside the
    other tiers and is treated as zero-shot inference -- see
    `configs/default.yaml: curriculum.main_curriculum_tiers`) use VAL+TEST
    combined: since nothing in this tier was ever trained on, both partitions
    are equally "unseen" and combining them gives a larger, less noisy
    evaluation set.
    """
    held_out_tiers = set(held_out_tiers) if held_out_tiers is not None else {"dual_evidence"}
    if tier not in held_out_tiers:
        return sp["test"]["pos"], sp["test"]["neg"]
    pos = torch.cat([sp["val"]["pos"], sp["test"]["pos"]], dim=1)
    neg = torch.cat([sp["val"]["neg"], sp["test"]["neg"]], dim=1)
    return pos, neg


def evaluate_scorer_on_splits(scorer, splits_per_tier: Dict[str, dict],
                              held_out_tiers: Optional[Iterable[str]] = None) -> Dict[str, dict]:
    """Score each tier's evaluation split with a baseline `scorer` (has
    .score(tf,tg)) and compute the same metrics used for MEvD-GRN — identical
    splits, fair. See `eval_edges_for_tier` for which partition is used per tier.

    EPR's base rate is taken from the eval split itself (pos+sampled neg),
    matching the pool early precision is ranked over. Using a genome-wide
    TF x gene count here instead would inflate EPR by orders of magnitude
    since the eval pool is deliberately positive-enriched, not genome-wide.
    """
    out: Dict[str, dict] = {}
    for tier, sp in splits_per_tier.items():
        pos, neg = eval_edges_for_tier(sp, tier, held_out_tiers)
        edges = torch.cat([pos, neg], dim=1).numpy()
        labels = np.concatenate([np.ones(pos.shape[1]), np.zeros(neg.shape[1])])
        if hasattr(scorer, "score_tier"):
            scores = scorer.score_tier(tier, edges[0], edges[1])
        else:
            scores = scorer.score(edges[0], edges[1])
        out[tier] = compute_all_metrics(labels, scores)
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
