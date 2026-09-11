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
  graph_mode='both'|'coexpr'|'tf_candidate'  -> which prior graph(s) the GNN sees
  use_edge_mlp=True     -> nonlinear edge embedding in the role-aware decoder
decode/forward return LOGITS (use predict() for probabilities).
"""
from __future__ import annotations

from typing import List, Optional, Tuple

import torch
import torch.nn as nn

from src.models.decoder import BilinearDecoder, RoleAwareDecoder
from src.models.encoders import ATACEncoder, FMEncoder, RNAEncoder
from src.models.fusion import ConcatFusion, GatedFusion
from src.models.gnn import GNNBackbone

Emb = Tuple[torch.Tensor, torch.Tensor]
_GRAPH_MODES = ("both", "coexpr", "tf_candidate")


class MEvDGRN(nn.Module):
    def __init__(self, rna_in_dim: int = 3, atac_in_dim: int = 3, hidden_dim: int = 128,
                 n_gnn_layers: int = 2, dropout: float = 0.2,
                 use_atac: bool = True, use_gnn: bool = True,
                 integration: str = "role_aware",
                 graph_mode: str = "both",
                 use_edge_mlp: bool = False,
                 combine_mode: str = "sum",
                 use_fm: bool = False, fm_in_dim: int = 768):
        super().__init__()
        if graph_mode not in _GRAPH_MODES:
            raise ValueError(f"graph_mode must be one of {_GRAPH_MODES}, got {graph_mode!r}")
        self.use_atac = use_atac
        self.use_gnn = use_gnn
        self.integration = integration
        self.graph_mode = graph_mode
        self.use_edge_mlp = use_edge_mlp
        self.combine_mode = combine_mode
        self.use_fm = use_fm
        self.rna_encoder = RNAEncoder(rna_in_dim, hidden_dim // 2, hidden_dim, dropout)
        self.atac_encoder = ATACEncoder(atac_in_dim, hidden_dim // 2, hidden_dim, dropout)
        # Pretrained foundation-model gene embedding (plan.md Section 5), fused
        # additively into the RNA channel -- cell-type-invariant, unlike the
        # hand-built co-expression signature used at decode time.
        self.fm_encoder = FMEncoder(fm_in_dim, hidden_dim, hidden_dim, dropout) if use_fm else None

        if integration == "role_aware":
            self.gnn_rna = GNNBackbone(hidden_dim, hidden_dim, n_gnn_layers, dropout,
                                       combine_mode=combine_mode) if use_gnn else None
            self.gnn_atac = GNNBackbone(hidden_dim, hidden_dim, n_gnn_layers, dropout,
                                        combine_mode=combine_mode) if use_gnn else None
            self.fusion = None
            self.gnn = None
            self.decoder = RoleAwareDecoder(hidden_dim, use_edge_mlp=use_edge_mlp)
        else:                                                # legacy fused integrations (ablations)
            self.fusion = GatedFusion(hidden_dim) if integration == "gated" else ConcatFusion(hidden_dim)
            self.gnn = GNNBackbone(hidden_dim, hidden_dim, n_gnn_layers, dropout,
                                   combine_mode=combine_mode) if use_gnn else None
            self.gnn_rna = self.gnn_atac = None
            self.decoder = BilinearDecoder(hidden_dim)

    def get_relation_weights(self) -> dict:
        """Learned per-layer relation-importance weights (combine_mode="gated"
        only), keyed by which GNN backbone they came from. Empty dict if
        combine_mode="sum" or use_gnn=False."""
        out = {}
        for name, gnn in (("rna", self.gnn_rna), ("atac", self.gnn_atac), ("fused", self.gnn)):
            if gnn is not None:
                w = gnn.get_relation_weights()
                if w is not None:
                    out[name] = w
        return out

    def _select_graphs(self, coexpr_edges: torch.Tensor,
                       tf_candidate_edges: torch.Tensor) -> List[torch.Tensor]:
        """Pick which structural prior(s) feed the relational GraphSAGE."""
        empty = torch.zeros((2, 0), dtype=torch.long, device=coexpr_edges.device)
        if self.graph_mode == "coexpr":
            return [coexpr_edges, empty]
        if self.graph_mode == "tf_candidate":
            return [empty, tf_candidate_edges]
        return [coexpr_edges, tf_candidate_edges]

    # -------------------------------------------------------------------------
    def encode(self, rna_features: torch.Tensor, atac_features: torch.Tensor,
               coexpr_edges: torch.Tensor, tf_candidate_edges: torch.Tensor,
               fm_features: Optional[torch.Tensor] = None) -> Emb:
        """Returns (h_rna, h_target-context). For role_aware the two channels stay
        separate; for fused integrations both entries are the shared fused embedding."""
        h_rna = self.rna_encoder(rna_features)
        if self.use_fm and fm_features is not None and fm_features.numel():
            h_rna = h_rna + self.fm_encoder(fm_features)
        h_atac = self.atac_encoder(atac_features) if self.use_atac else torch.zeros_like(h_rna)
        graphs = self._select_graphs(coexpr_edges, tf_candidate_edges)
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
            torch.zeros(tf_idx.shape[0], 4, device=h_rna.device)
        return self.decoder(h_rna[tf_idx], h_rna[target_idx], h_atac[target_idx],
                            coexpr, open_j)

    def forward(self, rna_features, atac_features, coexpr_edges, tf_candidate_edges,
                tf_idx, target_idx, signatures=None, openness=None,
                fm_features=None) -> torch.Tensor:
        emb = self.encode(rna_features, atac_features, coexpr_edges, tf_candidate_edges,
                          fm_features)
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
        if self.fm_encoder is not None:
            for p in self.fm_encoder.parameters():
                p.requires_grad = not frozen
