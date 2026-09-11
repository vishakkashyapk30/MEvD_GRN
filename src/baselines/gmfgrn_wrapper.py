"""GMFGRN baseline — the direct architectural competitor (GNN + matrix
factorization, RNA-only; plan 9.3).

Two entry points:
  run_gmfgrn()   -> uses the OFFICIAL repo if cloned to src/baselines/gmfgrn_repo/
                    (https://github.com/Lishuoyy/GMFGRN); raises with instructions
                    otherwise.
  run_gmf_gae()  -> a self-contained, faithful stand-in from the same family:
                    a Graph Auto-Encoder (Kipf & Welling) over a co-expression
                    kNN graph, scoring TF->target by embedding inner product.
                    Always runnable (no external deps beyond torch/sklearn).
Results are labeled distinctly so the official method is never misrepresented.
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path
from typing import Dict

import numpy as np
import torch
import torch.nn.functional as F

from src.baselines.common import (EmbeddingScorer, load_expression_df,
                                   universe_gene_names)
from src.models.gnn import GNNBackbone


# --------------------------------------------------------------------------- official
def run_gmfgrn(rna_path: str, gene_index: Dict[str, int], pcfg: dict, device: str = "cuda:0"):
    """Try the official GMFGRN repo. Clone it first:
        git clone https://github.com/Lishuoyy/GMFGRN src/baselines/gmfgrn_repo
    """
    repo = Path(__file__).parent / "gmfgrn_repo"
    if not repo.exists():
        raise RuntimeError(
            "Official GMFGRN not found. Clone it:\n"
            "  git clone https://github.com/Lishuoyy/GMFGRN src/baselines/gmfgrn_repo\n"
            "then adapt this function to the repo's train/predict entry points, or use "
            "run_gmf_gae() for the self-contained graph-MF baseline.")
    sys.path.insert(0, str(repo))
    importlib.import_module("GMFGRN")  # entry module — adapt to repo API
    raise NotImplementedError(
        "GMFGRN repo found; wire its training call here to return an EmbeddingScorer "
        "or DenseMatrixScorer aligned to gene_index.")


# --------------------------------------------------------------------------- self-contained GAE
class _GAE(torch.nn.Module):
    def __init__(self, in_dim, hidden=128, emb=64, dropout=0.2):
        super().__init__()
        # Single relation (the co-expression kNN graph); GNNBackbone always
        # expects a LIST of per-relation edge_index tensors, one per relation.
        self.enc = GNNBackbone(in_dim, hidden, num_layers=2, dropout=dropout, n_relations=1)
        self.proj = torch.nn.Linear(hidden, emb)

    def encode(self, x, edge_index):
        return self.proj(self.enc(x, [edge_index]))


def run_gmf_gae(rna_path: str, gene_index: Dict[str, int], pcfg: dict,
                device: str = "cuda:0", k_neighbors: int = 15, init_dim: int = 64,
                emb_dim: int = 64, epochs: int = 100, seed: int = 42) -> EmbeddingScorer:
    """Unsupervised graph auto-encoder over a co-expression kNN graph."""
    from sklearn.decomposition import TruncatedSVD
    from sklearn.neighbors import NearestNeighbors

    torch.manual_seed(seed)
    dev = device if (str(device).startswith("cuda") and torch.cuda.is_available()) else "cpu"
    genes = universe_gene_names(gene_index)
    expr = load_expression_df(rna_path, gene_index, pcfg)[genes].to_numpy(np.float32)  # cells×genes
    gene_profiles = expr.T                                        # genes × cells
    n_genes = gene_profiles.shape[0]

    # gene node features via truncated SVD of gene profiles
    d = min(init_dim, max(2, min(gene_profiles.shape) - 1))
    X = TruncatedSVD(n_components=d, random_state=seed).fit_transform(gene_profiles)
    X = torch.tensor(X, dtype=torch.float32, device=dev)

    # co-expression kNN graph (cosine)
    k = min(k_neighbors, n_genes - 1)
    nn = NearestNeighbors(n_neighbors=k + 1, metric="cosine").fit(gene_profiles)
    _, idx = nn.kneighbors(gene_profiles)
    src = np.repeat(np.arange(n_genes), k)
    dst = idx[:, 1:].reshape(-1)                                  # drop self
    edge_index = torch.tensor(np.stack([src, dst]), dtype=torch.long, device=dev)
    print(f"[gmf-gae] {n_genes} genes, {edge_index.shape[1]} kNN edges, "
          f"init_dim={d}, emb={emb_dim} on {dev}", flush=True)

    model = _GAE(X.shape[1], hidden=128, emb=emb_dim).to(dev)
    opt = torch.optim.Adam(model.parameters(), lr=1e-2, weight_decay=1e-5)
    n_edges = edge_index.shape[1]
    for ep in range(1, epochs + 1):
        model.train()
        z = model.encode(X, edge_index)
        pos = (z[edge_index[0]] * z[edge_index[1]]).sum(-1)
        neg_dst = torch.randint(0, n_genes, (n_edges,), device=dev)
        neg = (z[edge_index[0]] * z[neg_dst]).sum(-1)
        logits = torch.cat([pos, neg])
        labels = torch.cat([torch.ones_like(pos), torch.zeros_like(neg)])
        loss = F.binary_cross_entropy_with_logits(logits, labels)
        opt.zero_grad(); loss.backward(); opt.step()
        if ep % 25 == 0 or ep == 1:
            print(f"[gmf-gae] epoch {ep:3d} loss={loss.item():.4f}", flush=True)

    model.eval()
    with torch.no_grad():
        Z = model.encode(X, edge_index).cpu().numpy()
    return EmbeddingScorer(Z)
