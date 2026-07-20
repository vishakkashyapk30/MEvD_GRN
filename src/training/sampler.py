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


def get_hard_negatives(evidence: Dict[str, torch.Tensor], current_tier: str,
                       hierarchy: List[str], n_hard: int,
                       generator: Optional[torch.Generator] = None) -> torch.Tensor:
    """Hard negatives = edges in the immediately-lower tier that are NOT in the
    current tier (e.g. binds but doesn't functionally regulate). Empty for tier 0.
    """
    if current_tier not in hierarchy:
        return torch.zeros((2, 0), dtype=torch.long)
    ci = hierarchy.index(current_tier)
    if ci == 0 or n_hard <= 0:
        return torch.zeros((2, 0), dtype=torch.long)
    lower = hierarchy[ci - 1]
    if lower not in evidence or current_tier not in evidence:
        return torch.zeros((2, 0), dtype=torch.long)
    lower_set = set(zip(evidence[lower][0].tolist(), evidence[lower][1].tolist()))
    cur_set = set(zip(evidence[current_tier][0].tolist(), evidence[current_tier][1].tolist()))
    cands = sorted(lower_set - cur_set)
    if not cands:
        return torch.zeros((2, 0), dtype=torch.long)
    n = min(n_hard, len(cands))
    idx = torch.randperm(len(cands), generator=generator)[:n].tolist()
    chosen = [cands[i] for i in idx]
    src, dst = zip(*chosen)
    return torch.tensor([list(src), list(dst)], dtype=torch.long)
