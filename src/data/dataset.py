"""Processed-data container: loads everything a training run needs for one cell type."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

import numpy as np
import torch

from src.utils.io import load_json


@dataclass
class CellTypeData:
    cell_type: str
    rna_features: torch.Tensor          # (N_G, 3) [mean, var, detection]
    atac_features: torch.Tensor         # (N_G, 2) [mean, var]
    tf_candidate_edges: torch.Tensor    # (2, E) TF -> accessible/co-expressed target (message passing)
    coexpr_edges: torch.Tensor          # (2, E) gene-gene co-expression kNN (message passing)
    rna_signature: torch.Tensor         # (N_G, d) co-expression signatures (s_i.s_j ~ corr)
    openness: torch.Tensor              # (N_G,) locus accessibility proxy
    gene_index: Dict[str, int]
    tf_indices: torch.Tensor            # (N_TF,)
    evidence: Dict[str, torch.Tensor]   # tier -> (2, E)
    negative_pool: torch.Tensor         # (2, N_neg)

    @property
    def n_genes(self) -> int:
        return self.rna_features.shape[0]

    def to(self, device) -> "CellTypeData":
        self.rna_features = self.rna_features.to(device)
        self.atac_features = self.atac_features.to(device)
        self.tf_candidate_edges = self.tf_candidate_edges.to(device)
        self.coexpr_edges = self.coexpr_edges.to(device)
        self.rna_signature = self.rna_signature.to(device)
        self.openness = self.openness.to(device)
        self.tf_indices = self.tf_indices.to(device)
        return self


def load_celltype_data(processed_dir: str, evidence_tiers: List[str]) -> CellTypeData:
    d = Path(processed_dir)
    cell_type = d.name
    rna = torch.from_numpy(np.load(d / "rna_features_aligned.npy")).float()
    atac = torch.from_numpy(np.load(d / "atac_features_aligned.npy")).float()
    # TF-candidate message-passing graph (formerly "prior_edges").
    tf_cand_path = d / "tf_candidate_edges.pt"
    tf_cand = torch.load(tf_cand_path if tf_cand_path.exists() else d / "prior_edges.pt")
    coexpr_path = d / "coexpr_edges.pt"
    coexpr = torch.load(coexpr_path) if coexpr_path.exists() else torch.zeros((2, 0), dtype=torch.long)
    sig_path = d / "rna_signature.npy"
    signature = (torch.from_numpy(np.load(sig_path)).float() if sig_path.exists()
                 else torch.zeros((rna.shape[0], 1)))
    open_path = d / "openness.npy"
    openness = (torch.from_numpy(np.load(open_path)).float() if open_path.exists()
                else torch.zeros(rna.shape[0]))
    gene_index = load_json(d / "gene_index.json")
    tf_indices = torch.tensor(load_json(d / "tf_indices.json"), dtype=torch.long)
    neg_pool = torch.load(d / "negative_pool.pt")
    evidence = {}
    for tier in evidence_tiers:
        p = d / f"evidence_{tier}.pt"
        if p.exists():
            evidence[tier] = torch.load(p)
    return CellTypeData(cell_type, rna, atac, tf_cand, coexpr, signature, openness,
                        gene_index, tf_indices, evidence, neg_pool)


def load_splits(splits_dir: str, cell_type: str, tier: str) -> Dict[str, Dict[str, torch.Tensor]]:
    return torch.load(Path(splits_dir) / f"{cell_type}_{tier}_splits.pt")
