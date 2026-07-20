"""MEvD-GRN main model: encoders -> (GNN over co-expression + TF-candidate graphs)
-> role-aware decoder (plan Part 5, 6.6).

Biological integration. RNA (expression / co-expression) and ATAC (chromatin
accessibility) are NOT interchangeable, so the default `role_aware` integration
keeps them in separate channels and combines them only at the edge: the target's
RNA state is gated by its ATAC accessibility (a regulatability precondition), TF
activity comes from the TF's RNA embedding, and explicit co-expression / locus
openness terms are added. See src/models/decoder.py.

Ablation hooks (set at construction):
  use_atac=False        -> rna_only (ATAC branch zeroed/ignored)
  use_gnn=False         -> no_gnn (features go straight to the decoder)
  integration='gated'   -> legacy convex gated fusion + bilinear decoder
  integration='concat'  -> concat fusion + bilinear decoder
decode/forward return LOGITS (use predict() for probabilities).
"""
from __future__ import annotations

from typing import Optional, Tuple

import torch
import torch.nn as nn

from src.models.decoder import BilinearDecoder, RoleAwareDecoder
from src.models.encoders import ATACEncoder, RNAEncoder
from src.models.fusion import ConcatFusion, GatedFusion
from src.models.gnn import GNNBackbone

Emb = Tuple[torch.Tensor, torch.Tensor]


class MEvDGRN(nn.Module):
    def __init__(self, rna_in_dim: int = 3, atac_in_dim: int = 2, hidden_dim: int = 128,
                 n_gnn_layers: int = 2, dropout: float = 0.2,
                 use_atac: bool = True, use_gnn: bool = True,
                 integration: str = "role_aware"):
        super().__init__()
        self.use_atac = use_atac
        self.use_gnn = use_gnn
        self.integration = integration
        self.rna_encoder = RNAEncoder(rna_in_dim, hidden_dim // 2, hidden_dim, dropout)
        self.atac_encoder = ATACEncoder(atac_in_dim, hidden_dim)

        if integration == "role_aware":
            self.gnn_rna = GNNBackbone(hidden_dim, hidden_dim, n_gnn_layers, dropout) if use_gnn else None
            self.gnn_atac = GNNBackbone(hidden_dim, hidden_dim, n_gnn_layers, dropout) if use_gnn else None
            self.fusion = None
            self.gnn = None
            self.decoder = RoleAwareDecoder(hidden_dim)
        else:                                                # legacy fused integrations (ablations)
            self.fusion = GatedFusion(hidden_dim) if integration == "gated" else ConcatFusion(hidden_dim)
            self.gnn = GNNBackbone(hidden_dim, hidden_dim, n_gnn_layers, dropout) if use_gnn else None
            self.gnn_rna = self.gnn_atac = None
            self.decoder = BilinearDecoder(hidden_dim)

    # -------------------------------------------------------------------------
    def encode(self, rna_features: torch.Tensor, atac_features: torch.Tensor,
               coexpr_edges: torch.Tensor, tf_candidate_edges: torch.Tensor) -> Emb:
        """Returns (h_rna, h_target-context). For role_aware the two channels stay
        separate; for fused integrations both entries are the shared fused embedding."""
        h_rna = self.rna_encoder(rna_features)
        h_atac = self.atac_encoder(atac_features) if self.use_atac else torch.zeros_like(h_rna)
        graphs = [coexpr_edges, tf_candidate_edges]
        if self.integration == "role_aware":
            if self.use_gnn:
                h_rna = self.gnn_rna(h_rna, graphs)
                h_atac = self.gnn_atac(h_atac, graphs)
            return h_rna, h_atac
        h = self.fusion(h_rna, h_atac)
        if self.use_gnn:
            h = self.gnn(h, graphs)
        return h, h

    def decode(self, emb: Emb, tf_idx: torch.Tensor, target_idx: torch.Tensor,
               signatures: Optional[torch.Tensor] = None,
               openness: Optional[torch.Tensor] = None) -> torch.Tensor:
        h_rna, h_atac = emb
        if self.integration != "role_aware":
            return self.decoder(h_rna[tf_idx], h_rna[target_idx])          # bilinear on fused
        if signatures is not None and signatures.numel():
            coexpr = (signatures[tf_idx] * signatures[target_idx]).sum(dim=-1)
        else:
            coexpr = torch.zeros(tf_idx.shape[0], device=h_rna.device)
        open_j = openness[target_idx] if openness is not None else \
            torch.zeros(tf_idx.shape[0], device=h_rna.device)
        return self.decoder(h_rna[tf_idx], h_rna[target_idx], h_atac[target_idx],
                            coexpr, open_j)

    def forward(self, rna_features, atac_features, coexpr_edges, tf_candidate_edges,
                tf_idx, target_idx, signatures=None, openness=None) -> torch.Tensor:
        emb = self.encode(rna_features, atac_features, coexpr_edges, tf_candidate_edges)
        return self.decode(emb, tf_idx, target_idx, signatures, openness)

    @torch.no_grad()
    def predict(self, *args, **kwargs) -> torch.Tensor:
        return torch.sigmoid(self.forward(*args, **kwargs))

    def count_parameters(self) -> int:
        return sum(p.numel() for p in self.parameters() if p.requires_grad)

    def set_encoder_frozen(self, frozen: bool) -> None:
        for p in self.rna_encoder.parameters():
            p.requires_grad = not frozen
        for p in self.atac_encoder.parameters():
            p.requires_grad = not frozen
