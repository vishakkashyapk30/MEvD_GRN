"""MeVD-GRN with BEAR-specific decoder terms (pre-registered M2 / M3,
docs/experiments/bear_grn_benchmark.md s12.6). A subclass, so the shared
model code (src/models) is untouched.

  logit(tf, g) = MEvDGRN.decode(...)                                  (M0 / M1)
               + w_hub * hub[g]                                        (M2, M3)
               + Linear(pair_feats[tf, g])                             (M3)

hub[g] = log1p(in-degree of g among the CURRENT fold's training TFs) / log1p(max);
it is set per fold with `set_hub()` and never contains a held-out TF's label.
pair_feats = bear_motif.pair_features (label-free sequence features).
"""
from __future__ import annotations

from typing import Optional

import numpy as np
import torch
import torch.nn as nn

from src.models.mevd_grn import MEvDGRN


def hub_vector(train_pos: np.ndarray, n_genes: int) -> np.ndarray:
    indeg = np.bincount(np.asarray(train_pos)[1], minlength=n_genes).astype(np.float32)
    h = np.log1p(indeg)
    return h / max(float(h.max()), 1e-6)


class BearMEvDGRN(MEvDGRN):
    def __init__(self, *args, use_hub: bool = False, pair_feats: Optional[np.ndarray] = None,
                 tf_nodes: Optional[np.ndarray] = None, n_genes: int = 0, **kwargs):
        super().__init__(*args, **kwargs)
        self.use_hub = bool(use_hub)
        self.register_buffer("hub", torch.zeros(max(n_genes, 1)))
        if self.use_hub:
            self.w_hub = nn.Parameter(torch.zeros(1))
        self.use_pair = pair_feats is not None
        if self.use_pair:
            row = torch.full((max(n_genes, 1),), -1, dtype=torch.long)
            row[torch.as_tensor(np.asarray(tf_nodes), dtype=torch.long)] = torch.arange(len(tf_nodes))
            self.register_buffer("tf_row", row)
            self.register_buffer("pair", torch.as_tensor(np.asarray(pair_feats, dtype=np.float32)))
            self.pair_lin = nn.Linear(self.pair.shape[-1], 1)
            nn.init.zeros_(self.pair_lin.weight)
            nn.init.zeros_(self.pair_lin.bias)

    def set_hub(self, hub: np.ndarray) -> None:
        self.hub.copy_(torch.as_tensor(hub, dtype=torch.float32, device=self.hub.device))

    def decode(self, emb, tf_idx, target_idx, signatures=None, openness=None):
        out = super().decode(emb, tf_idx, target_idx, signatures, openness)
        if self.use_hub:
            out = out + self.w_hub * self.hub[target_idx]
        if self.use_pair:
            r = self.tf_row[tf_idx]
            f = self.pair[r.clamp(min=0), target_idx] * (r >= 0).unsqueeze(-1).float()
            out = out + self.pair_lin(f).squeeze(-1)
        return out
