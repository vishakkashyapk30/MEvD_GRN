"""scMultiomeGRN baseline adapted to MEvD-GRN's directed TF->gene benchmark.

The implementation follows Xu et al. (NAR 2025, gkaf138) and the authors'
Zenodo source (10.5281/zenodo.14848389): modality-specific edge-aware
attention, cross-modal attention, a graph convolution, and a symmetric MLP
link decoder.

The published method is TF->TF only and materializes an N x N score matrix.
MEvD-GRN instead evaluates TF->all-gene edges (roughly 20k nodes), so this
adapter makes two explicit scalability/generalization changes:
  * all benchmark genes are nodes, with the already aligned RNA/ATAC features;
  * only requested edge pairs are decoded instead of materializing N x N.

One model is fitted per evidence tier using only that tier's training graph.
This is necessary because training on (say) localization and evaluating on the
nested dual-evidence split can otherwise leak nominal test positives.
"""
from __future__ import annotations

import copy
import random
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.utils import softmax

from src.baselines.common import TieredScorer, load_expression_df


def _standardize(x: torch.Tensor) -> torch.Tensor:
    mean = x.mean(0, keepdim=True)
    std = x.std(0, keepdim=True, unbiased=False).clamp_min(1e-6)
    return torch.nan_to_num((x - mean) / std)


class _EdgeAwareConv(nn.Module):
    """Official GraFrankConv computation, expressed without dense adjacency."""

    def __init__(self, in_dim: int, out_dim: int, edge_dim: int = 36):
        super().__init__()
        self.self_linear = nn.Linear(in_dim, out_dim)
        self.message_linear = nn.Linear(in_dim + edge_dim, out_dim)
        self.attn = nn.Linear(out_dim, 1)
        self.attn_i = nn.Linear(out_dim, 1)
        self.lin_l = nn.Linear(out_dim, out_dim)
        self.lin_r = nn.Linear(out_dim, out_dim)
        # Two valid 3x3 convolutions + 2x2 pooling map 16x16 -> 6x6.
        self.edge_conv = nn.Sequential(
            nn.Conv2d(1, 1, 3),
            nn.Conv2d(1, 1, 3),
            nn.MaxPool2d(2),
            nn.Dropout(0.1),
        )

    def forward(self, x: torch.Tensor, edge_index: torch.Tensor,
                edge_attr: torch.Tensor) -> torch.Tensor:
        src, dst = edge_index
        self_emb = self.self_linear(x)
        edge_emb = self.edge_conv(edge_attr.unsqueeze(1)).flatten(1)
        message = self.message_linear(torch.cat([x[src], edge_emb], dim=-1))
        alpha = F.leaky_relu(
            self.attn(message) + self.attn_i(self_emb)[dst], negative_slope=0.1
        )
        alpha = softmax(alpha, dst, num_nodes=x.shape[0])
        alpha = F.dropout(alpha, p=0.1, training=self.training)
        aggregate = torch.zeros(
            (x.shape[0], message.shape[1]), dtype=x.dtype, device=x.device
        )
        aggregate.index_add_(0, dst, message * alpha)
        return self.lin_l(aggregate) + self.lin_r(self_emb)


class _ModalityStack(nn.Module):
    def __init__(self, in_dim: int, hidden_dim: int, n_layers: int = 2):
        super().__init__()
        dims = [in_dim] + [hidden_dim] * n_layers
        self.layers = nn.ModuleList(
            [_EdgeAwareConv(dims[i], dims[i + 1]) for i in range(n_layers)]
        )

    def forward(self, x, edge_index, edge_attr):
        for i, layer in enumerate(self.layers):
            x = layer(x, edge_index, edge_attr)
            if i + 1 < len(self.layers):
                x = F.dropout(F.relu(x), p=0.1, training=self.training)
        return x


class _CrossModalAttention(nn.Module):
    def __init__(self, hidden_dim: int):
        super().__init__()
        self.transform = nn.Linear(hidden_dim, hidden_dim)
        self.attention = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, modalities):
        stacked = torch.stack([x.relu() for x in modalities], dim=1)
        weights = torch.softmax(self.attention(stacked).squeeze(-1), dim=1)
        return (weights.unsqueeze(-1) * self.transform(stacked)).sum(dim=1)


class _SymmetricPairDecoder(nn.Module):
    """Official concatenation MLP decoder without allocating an N x N tensor."""

    def __init__(self, node_dim: int):
        super().__init__()
        self.mlp = nn.Sequential(
            nn.Linear(node_dim * 2, node_dim),
            nn.Linear(node_dim, 1),
        )

    def forward(self, z_src, z_dst):
        uv = self.mlp(torch.cat([z_src, z_dst], dim=-1)).squeeze(-1)
        vu = self.mlp(torch.cat([z_dst, z_src], dim=-1)).squeeze(-1)
        return torch.sigmoid((uv + vu) * 0.5)


class ScMultiomeGRN(nn.Module):
    def __init__(self, rna_dim: int, atac_dim: int, fusion_dim: int = 512,
                 hidden_dim: int = 256, n_layers: int = 2):
        super().__init__()
        self.rna_stack = _ModalityStack(rna_dim, fusion_dim, n_layers)
        self.atac_stack = _ModalityStack(atac_dim, fusion_dim, n_layers)
        self.fusion = _CrossModalAttention(fusion_dim)
        self.graph_linear = nn.Linear(fusion_dim, hidden_dim)
        self.decoder = _SymmetricPairDecoder(hidden_dim)

    def encode(self, rna_x, atac_x, edge_index, edge_attr, adj_norm):
        rna_z = self.rna_stack(rna_x, edge_index, edge_attr)
        atac_z = self.atac_stack(atac_x, edge_index, edge_attr)
        fused = self.fusion([rna_z, atac_z])
        return torch.sparse.mm(adj_norm, self.graph_linear(fused))

    def decode(self, z, src, dst):
        return self.decoder(z[src], z[dst])


def _make_training_graph(train_pos: torch.Tensor, n_nodes: int, max_neighbors: int,
                         seed: int) -> torch.Tensor:
    """Build the initial network from training positives only, then symmetrize."""
    rng = np.random.default_rng(seed)
    edges = train_pos.cpu().numpy()
    kept = []
    for src in np.unique(edges[0]):
        ids = np.flatnonzero(edges[0] == src)
        if len(ids) > max_neighbors:
            ids = rng.choice(ids, size=max_neighbors, replace=False)
        kept.append(edges[:, ids])
    directed = np.concatenate(kept, axis=1) if kept else np.zeros((2, 0), dtype=np.int64)
    reverse = directed[::-1]
    self_edges = np.stack([np.arange(n_nodes), np.arange(n_nodes)])
    graph = np.concatenate([directed, reverse, self_edges], axis=1)
    graph = np.unique(graph, axis=1)
    return torch.as_tensor(graph, dtype=torch.long)


def _normalized_adjacency(edge_index: torch.Tensor, n_nodes: int,
                          device: str) -> torch.Tensor:
    src, dst = edge_index.to(device)
    degree = torch.bincount(dst, minlength=n_nodes).float().clamp_min(1.0)
    values = degree[src].rsqrt() * degree[dst].rsqrt()
    # Sparse matrix row=destination, col=source performs neighbor aggregation.
    return torch.sparse_coo_tensor(
        torch.stack([dst, src]), values, (n_nodes, n_nodes), device=device
    ).coalesce()


def _expression_bins(expression: np.ndarray, bins: int = 16) -> np.ndarray:
    x = np.log10(np.maximum(expression, 0.0) + 1e-2)
    lo = x.min(axis=0, keepdims=True)
    span = x.max(axis=0, keepdims=True) - lo
    scaled = np.floor((x - lo) * bins / np.maximum(span, 1e-12))
    return np.clip(scaled, 0, bins - 1).astype(np.uint8)


def _edge_histograms(expression: np.ndarray, edge_index: torch.Tensor,
                     bins: int = 16) -> torch.Tensor:
    """Paper equation (2), computed efficiently from pre-binned expression."""
    digitized = _expression_bins(expression, bins)
    edges = edge_index.cpu().numpy()
    attrs = np.empty((edges.shape[1], bins, bins), dtype=np.float32)
    cache = {}
    n_cells = max(expression.shape[0], 1)
    for k, (src, dst) in enumerate(edges.T):
        key = (int(src), int(dst))
        reverse = (int(dst), int(src))
        if reverse in cache:
            attrs[k] = cache[reverse].T
            cache[key] = attrs[k]
            continue
        flat = digitized[:, src].astype(np.int16) * bins + digitized[:, dst]
        hist = np.bincount(flat, minlength=bins * bins).reshape(bins, bins)
        attr = (np.log10(hist / n_cells + 1e-4) + 4.0) / 4.0
        attrs[k] = attr.astype(np.float32)
        cache[key] = attrs[k]
    return torch.from_numpy(attrs)


def _edge_arrays(split_part: dict):
    pos = split_part["pos"].long()
    neg = split_part["neg"].long()
    edges = torch.cat([pos, neg], dim=1)
    labels = torch.cat([
        torch.ones(pos.shape[1], dtype=torch.float32),
        torch.zeros(neg.shape[1], dtype=torch.float32),
    ])
    return edges, labels


def _sample_epoch_edges(edges: torch.Tensor, labels: torch.Tensor, cap: int,
                        rng: np.random.Generator):
    if cap <= 0 or edges.shape[1] <= cap:
        order = rng.permutation(edges.shape[1])
    else:
        pos = np.flatnonzero(labels.numpy() == 1)
        neg = np.flatnonzero(labels.numpy() == 0)
        n_pos = min(len(pos), max(1, round(cap * len(pos) / (len(pos) + len(neg)))))
        n_neg = min(len(neg), cap - n_pos)
        order = np.concatenate([
            rng.choice(pos, n_pos, replace=False),
            rng.choice(neg, n_neg, replace=False),
        ])
        rng.shuffle(order)
    return edges[:, order], labels[order]


@torch.no_grad()
def _score_edges(model, z, edges, batch_size):
    scores = []
    for start in range(0, edges.shape[1], batch_size):
        e = edges[:, start:start + batch_size].to(z.device)
        scores.append(model.decode(z, e[0], e[1]).cpu())
    return torch.cat(scores)


def _validation_aupr(model, rna_x, atac_x, graph, edge_attr, adj_norm,
                     val_edges, val_labels, batch_size):
    from sklearn.metrics import average_precision_score

    model.eval()
    with torch.no_grad():
        z = model.encode(rna_x, atac_x, graph, edge_attr, adj_norm)
        scores = _score_edges(model, z, val_edges, batch_size).numpy()
    return float(average_precision_score(val_labels.numpy(), scores))


class _ModelScorer:
    def __init__(self, model, z, batch_size=65536):
        self.model = model
        self.z = z
        self.batch_size = batch_size

    def score(self, tf_idx: np.ndarray, tg_idx: np.ndarray) -> np.ndarray:
        edges = torch.as_tensor(np.stack([tf_idx, tg_idx]), dtype=torch.long)
        self.model.eval()
        with torch.no_grad():
            return _score_edges(self.model, self.z, edges, self.batch_size).numpy()


def _fit_tier(rna_x: torch.Tensor, atac_x: torch.Tensor, expression: np.ndarray,
              split: dict, device: str, seed: int, max_epochs: int, patience: int,
              graph_max_neighbors: int, train_edge_cap: int, batch_size: int,
              learning_rate: float, ckpt_path: Optional[str] = None) -> _ModelScorer:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if str(device).startswith("cuda") and torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        dev = device
    else:
        dev = "cpu"

    n_nodes = rna_x.shape[0]
    graph_cpu = _make_training_graph(
        split["train"]["pos"], n_nodes, graph_max_neighbors, seed
    )
    print(
        f"[scMultiomeGRN] initial train-only graph: {n_nodes} nodes, "
        f"{graph_cpu.shape[1]} directed/self edges",
        flush=True,
    )
    hist_cpu = _edge_histograms(expression, graph_cpu)
    graph = graph_cpu.to(dev)
    edge_attr = hist_cpu.to(dev)
    adj_norm = _normalized_adjacency(graph_cpu, n_nodes, dev)
    rna_x = rna_x.to(dev)
    atac_x = atac_x.to(dev)

    model = ScMultiomeGRN(rna_x.shape[1], atac_x.shape[1]).to(dev)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    train_edges, train_labels = _edge_arrays(split["train"])
    val_edges, val_labels = _edge_arrays(split["val"])
    pos_weight = float((train_labels == 0).sum() / max(int((train_labels == 1).sum()), 1))
    rng = np.random.default_rng(seed)
    best_aupr = -np.inf
    best_state: Optional[dict] = None
    stale = 0

    for epoch in range(1, max_epochs + 1):
        sampled_edges, sampled_labels = _sample_epoch_edges(
            train_edges, train_labels, train_edge_cap, rng
        )
        model.train()
        optimizer.zero_grad(set_to_none=True)
        z = model.encode(rna_x, atac_x, graph, edge_attr, adj_norm)
        total = sampled_edges.shape[1]
        epoch_loss = 0.0
        for start in range(0, total, batch_size):
            end = min(start + batch_size, total)
            e = sampled_edges[:, start:end].to(dev)
            y = sampled_labels[start:end].to(dev)
            pred = model.decode(z, e[0], e[1])
            per_item = F.binary_cross_entropy(pred, y, reduction="none")
            weights = torch.where(y > 0, pos_weight, 1.0)
            loss = (per_item * weights).mean() * ((end - start) / total)
            loss.backward(retain_graph=end < total)
            epoch_loss += float(loss.detach())
        optimizer.step()

        if epoch == 1 or epoch % 5 == 0 or epoch == max_epochs:
            val_aupr = _validation_aupr(
                model, rna_x, atac_x, graph, edge_attr, adj_norm,
                val_edges, val_labels, batch_size,
            )
            print(
                f"[scMultiomeGRN] epoch={epoch:4d} loss={epoch_loss:.5f} "
                f"val_AUPR={val_aupr:.5f}",
                flush=True,
            )
            if val_aupr > best_aupr + 1e-6:
                best_aupr = val_aupr
                best_state = copy.deepcopy(model.state_dict())
                stale = 0
                # Persist best-so-far to disk on every improvement. Training
                # runs here can take hours (max_epochs=2000) and only return
                # a usable model at the very end otherwise -- if the job is
                # killed (e.g. a SLURM walltime or an OOM) mid-run, this is
                # the difference between losing everything and keeping the
                # best result seen up to that point.
                if ckpt_path is not None:
                    torch.save({"epoch": epoch, "best_val_aupr": best_aupr,
                               "model": best_state}, ckpt_path)
                    print(f"[scMultiomeGRN] checkpoint saved: epoch={epoch} "
                          f"val_AUPR={best_aupr:.5f} -> {ckpt_path}", flush=True)
            else:
                stale += 1
                if stale >= patience:
                    print(f"[scMultiomeGRN] early stop at epoch {epoch}", flush=True)
                    break

    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        z = model.encode(rna_x, atac_x, graph, edge_attr, adj_norm).detach()
    # Edge attributes and graph are no longer needed after embeddings are fixed.
    return _ModelScorer(model, z, batch_size=max(batch_size, 65536))


def run_scmultiomegrn(
    rna_path: str,
    data,
    splits_per_tier: Dict[str, dict],
    pcfg: dict,
    device: str = "cuda:0",
    seed: int = 42,
    max_epochs: int = 2000,
    patience: int = 100,
    graph_max_neighbors: int = 500,
    train_edge_cap: int = 200_000,
    batch_size: int = 8192,
    learning_rate: float = 1e-5,
    main_curriculum_tiers: Optional[list] = None,
    cell_type: str = "K562",
) -> TieredScorer:
    """Train leakage-safe scMultiomeGRN models on the shared benchmark splits.

    Only fits one model per tier in `main_curriculum_tiers` (default:
    localization, perturbation) -- matching MEvD-GRN's own protocol of
    treating tiers outside that set (default: dual_evidence) as held-out
    zero-shot inference rather than a fine-tuning target, since dual_evidence
    is heavily nested inside the other tiers (see configs/default.yaml). The
    returned TieredScorer falls back to the LAST main-tier model (the most
    "refined" one, e.g. perturbation) for any tier it never explicitly trained.
    """
    expression = load_expression_df(rna_path, data.gene_index, pcfg).to_numpy(np.float32)
    signature = data.rna_signature.float()
    if signature.ndim != 2 or signature.shape[1] <= 1 or not torch.isfinite(signature).all():
        signature = data.rna_features.float()
    rna_x = _standardize(torch.cat([signature, data.rna_features.float()], dim=1))
    # openness is now a (n_genes, 4) regulatory-potential descriptor (see
    # plan.md Section 2), not a single scalar -- no reshape needed.
    atac_x = _standardize(torch.cat([
        data.atac_features.float(),
        data.openness.float(),
    ], dim=1))

    main_tiers = list(main_curriculum_tiers) if main_curriculum_tiers else list(splits_per_tier)
    trainable = [t for t in splits_per_tier if t in main_tiers]
    skipped = [t for t in splits_per_tier if t not in main_tiers]
    if skipped:
        print(f"[scMultiomeGRN] treating {skipped} as held-out zero-shot inference "
              f"(no model trained; scored with the {trainable[-1] if trainable else '?'} "
              "model)", flush=True)

    ckpt_dir = Path("results/checkpoints") / f"scmultiomegrn_{cell_type}"
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    scorers = {}
    for offset, tier in enumerate(trainable):
        print(f"[scMultiomeGRN] fitting tier={tier}", flush=True)
        scorers[tier] = _fit_tier(
            rna_x, atac_x, expression, splits_per_tier[tier], device, seed + offset,
            max_epochs, patience, graph_max_neighbors, train_edge_cap,
            batch_size, learning_rate, ckpt_path=str(ckpt_dir / f"{tier}_best.pt"),
        )
    fallback = trainable[-1] if trainable else None
    return TieredScorer(scorers, fallback_tier=fallback)
