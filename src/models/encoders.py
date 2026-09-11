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


class FMEncoder(nn.Module):
    """2-layer MLP: (N_G, in_dim) -> (N_G, out_dim). LayerNorm + GELU + dropout.

    Projects a frozen, pretrained foundation-model gene embedding (e.g.
    Geneformer's input token embedding, in_dim=768 -- see
    scripts/11_extract_fm_embeddings.py) down to the model's hidden size.
    Unlike RNAEncoder's input, this embedding is CELL-TYPE-INVARIANT: the
    same gene gets the same vector regardless of which cell type is being
    modeled (plan.md Section 5). Same shape as RNAEncoder/ATACEncoder.
    """

    def __init__(self, in_dim: int = 768, hidden_dim: int = 128,
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
    """2-layer MLP: (N_G, in_dim) -> (N_G, out_dim). LayerNorm + GELU + dropout.

    Previously a single Linear+LayerNorm+GELU with no hidden layer -- fine
    when the input was a near-rank-1, 2-dim signal (see plan.md Section 2),
    but a real bottleneck now that the input is the richer RP-weighted
    [mean, var, detection] gene activity. Matches RNAEncoder's shape.
    """

    def __init__(self, in_dim: int = 3, hidden_dim: int = 64,
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
