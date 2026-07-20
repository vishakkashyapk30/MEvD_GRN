"""Modality fusion (plan 5.4). Gated fusion is the default; concat is the ablation."""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class GatedFusion(nn.Module):
    """h = g ⊙ h_rna + (1−g) ⊙ h_atac, with g = σ(W[h_rna; h_atac]) per gene/dim."""

    def __init__(self, dim: int = 128):
        super().__init__()
        self.gate = nn.Sequential(nn.Linear(dim * 2, dim), nn.Sigmoid())

    def forward(self, h_rna: torch.Tensor, h_atac: torch.Tensor) -> torch.Tensor:
        g = self.gate(torch.cat([h_rna, h_atac], dim=-1))
        return g * h_rna + (1.0 - g) * h_atac


class ConcatFusion(nn.Module):
    """Ablation: concatenate then project (no gating)."""

    def __init__(self, dim: int = 128):
        super().__init__()
        self.proj = nn.Linear(dim * 2, dim)

    def forward(self, h_rna: torch.Tensor, h_atac: torch.Tensor) -> torch.Tensor:
        return F.gelu(self.proj(torch.cat([h_rna, h_atac], dim=-1)))
