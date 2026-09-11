"""TF -> Target scoring heads (plan 5.6).

Two heads:
  * BilinearDecoder: legacy asymmetric bilinear on a single fused embedding
    (used by the `gated`/`concat` integration ablations).
  * RoleAwareDecoder: the biologically-motivated head. RNA and ATAC play distinct,
    non-interchangeable roles. Regulation of target j by TF i requires that (a) the
    target locus is ACCESSIBLE (an ATAC property) and (b) TF and target are
    co-expressed / compatible (an RNA property). We therefore gate the target's
    RNA representation by an accessibility gate derived from its ATAC embedding
    (accessibility is a PRECONDITION, not an interchangeable signal), score the
    asymmetric TF->target compatibility bilinearly, and add explicit co-expression
    and locus-openness terms.

Both heads return LOGITS; apply sigmoid in the loss/eval for numerical stability.
"""
from __future__ import annotations

import torch
import torch.nn as nn


class BilinearDecoder(nn.Module):
    """score(TF_i, TG_j) = h_i . W . h_j + b  — asymmetric via learnable W."""

    def __init__(self, dim: int = 128):
        super().__init__()
        self.W = nn.Parameter(torch.empty(dim, dim))
        self.b = nn.Parameter(torch.zeros(1))
        nn.init.xavier_uniform_(self.W)

    def forward(self, h_tf: torch.Tensor, h_target: torch.Tensor) -> torch.Tensor:
        logits = (h_tf @ self.W * h_target).sum(dim=-1) + self.b
        return logits


class RoleAwareDecoder(nn.Module):
    """Role-aware TF->target score with ATAC-as-accessibility-gate.

        a_j    = sigmoid(gate(h_atac_j))            # per-dim regulatability gate in [0,1]
        z_tg   = h_rna_j * a_j                        # target regulatable-expression state
        raw    = (h_rna_i . W) . z_tg                 # asymmetric TF->target compatibility

    Default edge terms (linear):
        logit  = raw + w_coexp * coexpr(i,j) + open_proj(openness_j) + b

    With use_edge_mlp=True (MPNN-style edge embedding; Gilmer et al. 2017):
        e_ij   = EdgeMLP([coexpr(i,j), openness_j])   # learned edge embedding
        logit  = raw + EdgeToLogit(e_ij) + b

    coexpr(i,j) is the signature dot product (~ Pearson correlation) supplied by the
    caller; openness_j is the target locus's regulatory-potential descriptor (a
    small vector -- mean/max RP weight, mean peak distance, peak count -- see
    plan.md Section 2), not a single collapsed scalar.
    """

    def __init__(self, dim: int = 128, use_edge_mlp: bool = False,
                 edge_hidden: int = 32, openness_dim: int = 4):
        super().__init__()
        self.use_edge_mlp = use_edge_mlp
        self.W = nn.Parameter(torch.empty(dim, dim))
        self.acc_gate = nn.Linear(dim, dim)
        self.b = nn.Parameter(torch.zeros(1))
        nn.init.xavier_uniform_(self.W)
        if use_edge_mlp:
            # (1 + openness_dim) edge features -> nonlinear embedding -> scalar contribution
            self.edge_mlp = nn.Sequential(
                nn.Linear(1 + openness_dim, edge_hidden),
                nn.GELU(),
                nn.Linear(edge_hidden, edge_hidden),
                nn.GELU(),
            )
            self.edge_to_logit = nn.Linear(edge_hidden, 1)
            self.w_coexp = self.open_proj = None
        else:
            self.edge_mlp = self.edge_to_logit = None
            self.w_coexp = nn.Parameter(torch.ones(1))
            self.open_proj = nn.Linear(openness_dim, 1)

    def forward(self, h_tf_rna: torch.Tensor, h_tg_rna: torch.Tensor,
                h_tg_atac: torch.Tensor, coexpr: torch.Tensor,
                openness_j: torch.Tensor) -> torch.Tensor:
        a_j = torch.sigmoid(self.acc_gate(h_tg_atac))          # (B, d) regulatability gate
        z_tg = h_tg_rna * a_j
        raw = (h_tf_rna @ self.W * z_tg).sum(dim=-1)           # (B,)
        if self.use_edge_mlp:
            e = self.edge_mlp(torch.cat([coexpr.unsqueeze(-1), openness_j], dim=-1))
            return raw + self.edge_to_logit(e).squeeze(-1) + self.b
        return raw + self.w_coexp * coexpr + self.open_proj(openness_j).squeeze(-1) + self.b
