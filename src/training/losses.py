"""Loss functions (plan 7.1). Operate on LOGITS for numerical stability."""
from __future__ import annotations

import torch
import torch.nn.functional as F


def bce_weighted_loss(logits: torch.Tensor, labels: torch.Tensor,
                      pos_weight: float = 5.0,
                      sample_weight: torch.Tensor | None = None) -> torch.Tensor:
    """Weighted BCE-with-logits. Positives upweighted by `pos_weight` to counter
    class imbalance; optional per-sample weight (used for memory replay)."""
    pw = torch.tensor(float(pos_weight), device=logits.device)
    loss = F.binary_cross_entropy_with_logits(logits, labels, pos_weight=pw, reduction="none")
    if sample_weight is not None:
        loss = loss * sample_weight
    return loss.mean()


def focal_loss(logits: torch.Tensor, labels: torch.Tensor,
               alpha: float = 0.75, gamma: float = 2.0,
               sample_weight: torch.Tensor | None = None) -> torch.Tensor:
    """Focal loss (Lin et al. 2017) on logits — down-weights easy examples."""
    p = torch.sigmoid(logits)
    ce = F.binary_cross_entropy_with_logits(logits, labels, reduction="none")
    p_t = p * labels + (1 - p) * (1 - labels)
    alpha_t = alpha * labels + (1 - alpha) * (1 - labels)
    loss = alpha_t * (1 - p_t).pow(gamma) * ce
    if sample_weight is not None:
        loss = loss * sample_weight
    return loss.mean()
