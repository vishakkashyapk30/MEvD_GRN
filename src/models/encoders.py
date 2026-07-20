"""Per-gene modality encoders (plan 5.2, 5.3)."""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class RNAEncoder(nn.Module):
    """2-layer MLP: (N_G, in_dim) -> (N_G, out_dim). LayerNorm + GELU + dropout."""

    def __init__(self, in_dim: int = 2, hidden_dim: int = 64,
                 out_dim: int = 128, dropout: float = 0.2):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, out_dim),
            nn.LayerNorm(out_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class ATACEncoder(nn.Module):
    """Single projection: (N_G, in_dim) -> (N_G, out_dim). Linear -> LayerNorm -> GELU."""

    def __init__(self, in_dim: int = 2, out_dim: int = 128):
        super().__init__()
        self.proj = nn.Linear(in_dim, out_dim)
        self.norm = nn.LayerNorm(out_dim)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.norm(F.gelu(self.proj(x)))
