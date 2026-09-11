"""Negative sampling (plan 6.3, 9.3): random pool draws + tier-aware hard negatives."""
from __future__ import annotations

from typing import Dict, List, Optional

import torch


def sample_negatives(neg_pool: torch.Tensor, n: int,
                     generator: Optional[torch.Generator] = None) -> torch.Tensor:
    """Randomly draw `n` edges (without replacement) from the negative pool."""
    total = neg_pool.shape[1]
    if total == 0 or n <= 0:
        return torch.zeros((2, 0), dtype=torch.long)
    n = min(n, total)
    idx = torch.randperm(total, generator=generator)[:n]
    return neg_pool[:, idx]


def get_hard_negatives(train_evidence: Dict[str, torch.Tensor],
                       full_evidence: Dict[str, torch.Tensor], current_tier: str,
                       hierarchy: List[str], n_hard: int,
                       generator: Optional[torch.Generator] = None,
                       exclude: Optional[torch.Tensor] = None) -> torch.Tensor:
    """Hard negatives = edges in the immediately-lower tier that are NOT in the
    current tier (e.g. binds but doesn't functionally regulate). Empty for tier 0.

    The lower-tier candidate pool (`lower_set`) MUST come from that tier's TRAIN
    split only (`train_evidence`), never its full evidence set. Otherwise the
    lower tier's own VAL/TEST positives get trained on as label=0 examples here,
    which directly contradicts the ground truth those edges are later scored
    against when evaluating that lower tier -- this is what was driving
    localization AUROC below random (not just "forgetting" toward it) after
    Stage 2, since it's active, systematic mislabeling of the eval set rather
    than mere signal decay.

    The current tier's exclusion set (`cur_set`) uses `full_evidence` (all
    splits) instead, since we want to guarantee NO current-tier positive --
    train, val, or test -- is ever mistakenly sampled as a negative for its
    own tier; that direction of over-exclusion is safe.

    `exclude` additionally removes any edge simultaneously being fed in as a
    label=1 replay example this epoch, so the same edge can't be label=0
    (hard negative) and label=1 (replay) in the same training pool.
    """
    if current_tier not in hierarchy:
        return torch.zeros((2, 0), dtype=torch.long)
    ci = hierarchy.index(current_tier)
    if ci == 0 or n_hard <= 0:
        return torch.zeros((2, 0), dtype=torch.long)
    lower = hierarchy[ci - 1]
    if lower not in train_evidence or current_tier not in full_evidence:
        return torch.zeros((2, 0), dtype=torch.long)
    lower_set = set(zip(train_evidence[lower][0].tolist(), train_evidence[lower][1].tolist()))
    cur_set = set(zip(full_evidence[current_tier][0].tolist(), full_evidence[current_tier][1].tolist()))
    exclude_set = set()
    if exclude is not None and exclude.shape[1]:
        exclude_set = set(zip(exclude[0].tolist(), exclude[1].tolist()))
    cands = sorted(lower_set - cur_set - exclude_set)
    if not cands:
        return torch.zeros((2, 0), dtype=torch.long)
    n = min(n_hard, len(cands))
    idx = torch.randperm(len(cands), generator=generator)[:n].tolist()
    chosen = [cands[i] for i in idx]
    src, dst = zip(*chosen)
    return torch.tensor([list(src), list(dst)], dtype=torch.long)
