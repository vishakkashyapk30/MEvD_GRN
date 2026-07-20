"""GraphSAGE backbone over TWO biological relations (plan Part 5.5).

The GNN aggregates context from two graphs with distinct meaning:
  * co-expression kNN (gene-gene, undirected): co-regulation modules;
  * TF -> accessible/co-expressed target candidates (directed): regulatory context.

Each layer runs one SAGEConv per relation and sums the results, so a node's update
mixes co-expression neighbours and candidate-regulatory neighbours. The TF->target
relation is made bidirectional for representation learning; directionality is
reintroduced only at the decoder.
"""
from __future__ import annotations

from typing import List, Optional

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import SAGEConv


def to_bidirectional(edge_index: torch.Tensor) -> torch.Tensor:
    """Concatenate reverse edges (dedup not required for SAGE mean-aggregation)."""
    if edge_index.numel() == 0:
        return edge_index
    rev = edge_index.flip(0)
    return torch.cat([edge_index, rev], dim=1)


class GNNBackbone(nn.Module):
    """Relational GraphSAGE: per layer, one SAGEConv per relation, summed.

    `n_relations` defaults to 2 (co-expression + TF-candidate). forward() accepts a
    list of edge_index tensors of that length; empty relations contribute nothing.
    """

    def __init__(self, in_dim: int = 128, hidden_dim: int = 128,
                 num_layers: int = 2, dropout: float = 0.2, n_relations: int = 2):
        super().__init__()
        self.n_relations = n_relations
        self.layers = nn.ModuleList()
        self.norms = nn.ModuleList()
        for i in range(num_layers):
            d_in = in_dim if i == 0 else hidden_dim
            self.layers.append(nn.ModuleList(
                [SAGEConv(d_in, hidden_dim, aggr="mean") for _ in range(n_relations)]))
            self.norms.append(nn.LayerNorm(hidden_dim))
        self.dropout = nn.Dropout(dropout)
        self.residual = in_dim == hidden_dim

    def forward(self, x: torch.Tensor, edge_indices: List[torch.Tensor]) -> torch.Tensor:
        rels = self._prepare(edge_indices, x.device)
        h = x
        for i, (convs, norm) in enumerate(zip(self.layers, self.norms)):
            agg = None
            for conv, ei in zip(convs, rels):
                out = conv(h, ei)
                agg = out if agg is None else agg + out
            h_new = self.dropout(F.gelu(norm(agg)))
            if self.residual and i > 0:
                h_new = h_new + h
            h = h_new
        return h

    def _prepare(self, edge_indices: List[torch.Tensor], device) -> List[torch.Tensor]:
        """Normalize to exactly `n_relations` bidirectional edge_index tensors."""
        eis = list(edge_indices)[: self.n_relations]
        while len(eis) < self.n_relations:
            eis.append(torch.zeros((2, 0), dtype=torch.long, device=device))
        out = []
        for ei in eis:
            ei = ei.to(device)
            if ei.numel() == 0:                       # self-loops keep SAGE well-defined
                ei = torch.zeros((2, 0), dtype=torch.long, device=device)
            out.append(to_bidirectional(ei))
        return out
