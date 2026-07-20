"""Graph construction (plan Part 4): gene universe, prior graph, evidence edges,
negative pool, and edge-level train/val/test splits.

Reference-network files are 3-column TSV WITH a header:  Source\\tTarget\\tRelationship
(the 'Relationship' column is Activation/Repression; v1 treats edges as binary).
Gene symbols are upper-cased to match the canonicalised single-cell symbols.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Set, Tuple

import numpy as np
import pandas as pd
import torch

Edge = Tuple[int, int]


def canon(name: str) -> str:
    return str(name).strip().upper()


# ----------------------------------------------------------------------------- networks
def load_network_pairs(path: str) -> List[Tuple[str, str]]:
    """Read a reference network -> list of (TF, Target) upper-cased pairs.

    Robust to presence/absence of the header and to 2- or 3-column files.
    """
    df = pd.read_csv(path, sep="\t", dtype=str, header=0)
    cols = [c.lower() for c in df.columns]
    if "source" in cols and "target" in cols:
        s = df.iloc[:, cols.index("source")]
        t = df.iloc[:, cols.index("target")]
    else:                                   # no usable header -> reread headerless
        df = pd.read_csv(path, sep="\t", dtype=str, header=None)
        s, t = df.iloc[:, 0], df.iloc[:, 1]
    pairs = [(canon(a), canon(b)) for a, b in zip(s, t)
             if isinstance(a, str) and isinstance(b, str)]
    return pairs


def build_gene_universe(rna_names: List[str],
                        network_paths: Dict[str, str]) -> Tuple[List[str], Set[str]]:
    """Universe = RNA genes ∩ (all genes appearing in any reference network).
    Returns (sorted universe list, set of TF names = any 'Source' gene)."""
    rna_set = set(rna_names)
    net_genes: Set[str] = set()
    tf_names: Set[str] = set()
    for path in network_paths.values():
        if not Path(path).exists():
            continue
        for tf, tg in load_network_pairs(path):
            net_genes.add(tf)
            net_genes.add(tg)
            tf_names.add(tf)
    universe = sorted(rna_set & net_genes)
    return universe, tf_names


def build_evidence_edges(network_paths: Dict[str, str],
                         gene_index: Dict[str, int]) -> Dict[str, torch.Tensor]:
    """tier -> LongTensor(2, E) of in-universe positive edges. Reports coverage."""
    out: Dict[str, torch.Tensor] = {}
    for tier, path in network_paths.items():
        if not Path(path).exists():
            print(f"[edges] WARNING missing {tier}: {path}", flush=True)
            continue
        pairs = load_network_pairs(path)
        src, dst = [], []
        for tf, tg in pairs:
            i, j = gene_index.get(tf), gene_index.get(tg)
            if i is not None and j is not None:
                src.append(i)
                dst.append(j)
        edges = torch.tensor([src, dst], dtype=torch.long) if src else torch.zeros((2, 0), dtype=torch.long)
        out[tier] = _unique_edges(edges)
        cov = 100 * out[tier].shape[1] / max(len(pairs), 1)
        print(f"[edges] {tier}: {out[tier].shape[1]}/{len(pairs)} in-universe ({cov:.1f}%)", flush=True)
    return out


def _unique_edges(edges: torch.Tensor) -> torch.Tensor:
    if edges.shape[1] == 0:
        return edges
    return torch.unique(edges, dim=1)


def edge_set(edges: torch.Tensor) -> Set[Edge]:
    return set(zip(edges[0].tolist(), edges[1].tolist()))


# ----------------------------------------------------------------------------- nesting
def verify_nesting(evidence: Dict[str, torch.Tensor],
                   hierarchy: List[str]) -> Dict[str, float]:
    """Check E_dual ⊆ E_pert ⊆/and E_loc (subset fractions). Curriculum assumes this."""
    report: Dict[str, float] = {}
    sets = {t: edge_set(e) for t, e in evidence.items()}
    for i in range(1, len(hierarchy)):
        hi = hierarchy[i]
        if hi not in sets or not sets[hi]:
            continue
        for j in range(i):
            lo = hierarchy[j]
            if lo not in sets:
                continue
            inter = len(sets[hi] & sets[lo])
            frac = inter / max(len(sets[hi]), 1)
            report[f"{hi}⊆{lo}"] = frac
            flag = "OK" if frac > 0.95 else "NOT nested"
            print(f"[nesting] {frac*100:5.1f}% of '{hi}' edges are in '{lo}'  [{flag}]", flush=True)
    return report


def reconstruct_dual(evidence: Dict[str, torch.Tensor]) -> torch.Tensor:
    """dual = localization ∩ perturbation (guarantees nesting)."""
    loc = edge_set(evidence.get("localization", torch.zeros((2, 0), dtype=torch.long)))
    pert = edge_set(evidence.get("perturbation", torch.zeros((2, 0), dtype=torch.long)))
    inter = sorted(loc & pert)
    if not inter:
        return torch.zeros((2, 0), dtype=torch.long)
    return torch.tensor(list(zip(*inter)), dtype=torch.long)


# ----------------------------------------------------------------------------- co-expression graph
def _topk_neighbors(signatures: np.ndarray, k: int, chunk: int = 512) -> List[List[int]]:
    """For each gene, the k most co-expressed OTHER genes by signature dot product
    (s_i . s_j ~ correlation). Chunked matmul so we never build the full NxN matrix."""
    n = signatures.shape[0]
    S = signatures.astype(np.float32)
    k = int(min(k, max(n - 1, 1)))
    neigh: List[List[int]] = []
    for start in range(0, n, chunk):
        block = S[start:start + chunk] @ S.T                 # (b, n) similarity
        for r in range(block.shape[0]):
            gi = start + r
            block[r, gi] = -np.inf                            # exclude self
            top = np.argpartition(block[r], -k)[-k:]
            neigh.append([int(j) for j in top])
    return neigh


def build_coexpression_graph(signatures: np.ndarray, k: int = 20) -> torch.Tensor:
    """Gene-gene co-expression kNN graph (UNDIRECTED). Nodes that co-vary across
    cells share regulatory context; message passing over this graph lets the GNN
    aggregate co-regulation modules. Returns a symmetric LongTensor(2, E).
    """
    n = signatures.shape[0]
    if n == 0 or signatures.shape[1] == 0 or not np.any(signatures):
        print("[coexpr] no signatures; empty co-expression graph", flush=True)
        return torch.zeros((2, 0), dtype=torch.long)
    neigh = _topk_neighbors(signatures, k)
    edges: Set[Edge] = set()
    for i, nbrs in enumerate(neigh):
        for j in nbrs:
            edges.add((i, j))
            edges.add((j, i))                                 # symmetric
    if not edges:
        return torch.zeros((2, 0), dtype=torch.long)
    src, dst = zip(*sorted(edges))
    g = torch.tensor([list(src), list(dst)], dtype=torch.long)
    print(f"[coexpr] co-expression kNN graph: k={k} -> {g.shape[1]} edges over {n} genes",
          flush=True)
    return g


# ----------------------------------------------------------------------------- TF candidate graph
def build_prior_graph(gene_index: Dict[str, int], tf_indices: List[int],
                      atac_features: np.ndarray, rna_features: np.ndarray,
                      top_k: int, signatures: np.ndarray | None = None,
                      openness: np.ndarray | None = None,
                      evidence: Dict[str, torch.Tensor] | None = None,
                      exclude_positives: bool = True) -> torch.Tensor:
    """TF-specific MESSAGE-PASSING candidate graph (directed TF -> target).

    For each TF, candidate targets are the top_k genes that are BOTH (a) accessible
    (positive `openness`, i.e. the locus has accessible peaks near its TSS) AND
    (b) most co-expressed with that TF (signature dot product). This yields a
    DIFFERENT neighbourhood per TF, unlike the old global accessibility ranking that
    gave every TF the same targets. Falls back to the global accessibility ranking
    when signatures/openness are unavailable (legacy-graph ablation).

    Used ONLY for GNN message passing — the decoder can score ANY (TF, gene) pair.
    Known positive edges are EXCLUDED so val/test labels can't leak via messages.
    Returns LongTensor(2, E).
    """
    n_genes = len(gene_index)
    pos_keys: Set[Edge] = set()
    if exclude_positives and evidence:
        for e in evidence.values():
            pos_keys |= edge_set(e)

    edges: Set[Edge] = set()
    have_coexpr = signatures is not None and signatures.size and np.any(signatures)

    if have_coexpr:
        S = signatures.astype(np.float32)
        accessible = (openness > 0) if (openness is not None and np.any(openness)) \
            else np.ones(n_genes, dtype=bool)
        acc_idx = np.where(accessible)[0]
        for tf in tf_indices:
            sim = S[acc_idx] @ S[tf]                          # co-expression of TF with accessible genes
            order = acc_idx[np.argsort(-sim)]
            cnt = 0
            for tgt in order:
                tgt = int(tgt)
                if tgt == tf or (tf, tgt) in pos_keys:
                    continue
                edges.add((tf, tgt))
                cnt += 1
                if cnt >= top_k:
                    break
        mode = "co-expression x accessibility (TF-specific)"
    else:                                                    # legacy global accessibility ranking
        score = atac_features[:, 0].astype(np.float64)
        if not np.any(score > 0):
            score = rna_features[:, 0].astype(np.float64)
        order_desc = np.argsort(-score)
        for tf in tf_indices:
            cnt = 0
            for tgt in order_desc:
                tgt = int(tgt)
                if tgt == tf or (tf, tgt) in pos_keys:
                    continue
                edges.add((tf, tgt))
                cnt += 1
                if cnt >= top_k:
                    break
        mode = "global accessibility (legacy)"

    src, dst = zip(*sorted(edges)) if edges else ([], [])
    prior = torch.tensor([list(src), list(dst)], dtype=torch.long)
    print(f"[prior] {len(tf_indices)} TFs, top_k={top_k}, mode={mode} -> {prior.shape[1]} "
          f"message-passing edges over {n_genes} genes (positives excluded={exclude_positives})",
          flush=True)
    return prior


# ----------------------------------------------------------------------------- negatives & splits
def create_negative_pool(all_positive: Dict[str, torch.Tensor], tf_indices: List[int],
                         n_genes: int, cap: int = 1_000_000, seed: int = 42) -> torch.Tensor:
    """Universal negative pool: random TF×gene pairs that are positive in NO tier.

    Sampling from the FULL bipartite candidate space (not the accessibility-ranked
    prior) is essential: restricting negatives to the prior's top-K most-accessible
    genes lets the model separate positives/negatives by an expression-level
    artifact, badly inflating AUPR. Field-standard practice is random negatives.
    """
    pos: Set[int] = set()
    for e in all_positive.values():
        for s, d in zip(e[0].tolist(), e[1].tolist()):
            pos.add(s * n_genes + d)
    tf_arr = np.asarray(tf_indices, dtype=np.int64)
    full_space = len(tf_arr) * n_genes
    target = min(cap, full_space - len(pos))
    rng = np.random.default_rng(seed)
    found: Set[int] = set()
    attempts, max_attempts = 0, 60
    while len(found) < target and attempts < max_attempts:
        batch = max(target * 2, 100_000)
        s = tf_arr[rng.integers(0, len(tf_arr), size=batch)]
        d = rng.integers(0, n_genes, size=batch)
        keys = s * n_genes + d
        for k, si, di in zip(keys.tolist(), s.tolist(), d.tolist()):
            if si == di or k in pos or k in found:
                continue
            found.add(k)
            if len(found) >= target:
                break
        attempts += 1
    src = [k // n_genes for k in found]
    dst = [k % n_genes for k in found]
    print(f"[neg] negative pool: {len(found)} random TF×gene non-positive edges "
          f"(space={full_space}, positives={len(pos)})", flush=True)
    if not src:
        return torch.zeros((2, 0), dtype=torch.long)
    return torch.tensor([src, dst], dtype=torch.long)


def _sample_cols(edges: torch.Tensor, n: int, generator: torch.Generator) -> torch.Tensor:
    total = edges.shape[1]
    if total == 0 or n <= 0:
        return torch.zeros((2, 0), dtype=torch.long)
    n = min(n, total)
    idx = torch.randperm(total, generator=generator)[:n]
    return edges[:, idx]


def create_edge_splits(pos_edges: torch.Tensor, neg_pool: torch.Tensor,
                       train_ratio: float, val_ratio: float,
                       neg_train_ratio: int, neg_eval_ratio: int,
                       seed: int) -> Dict[str, Dict[str, torch.Tensor]]:
    """Edge-level split (plan 4.4). Negatives drawn disjointly per split."""
    g = torch.Generator().manual_seed(seed)
    n_pos = pos_edges.shape[1]
    perm = torch.randperm(n_pos, generator=g)
    n_train = int(n_pos * train_ratio)
    n_val = int(n_pos * val_ratio)
    train_pos = pos_edges[:, perm[:n_train]]
    val_pos = pos_edges[:, perm[n_train:n_train + n_val]]
    test_pos = pos_edges[:, perm[n_train + n_val:]]

    # Disjoint negative draws. Partition the pool by split RATIO first so val/test
    # always receive negatives even when the pool is small vs. a dense tier (e.g.
    # localization), then cap each split to its desired neg:pos ratio.
    n_neg_total = neg_pool.shape[1]
    neg_shuf = neg_pool[:, torch.randperm(n_neg_total, generator=g)]
    c_tr = int(n_neg_total * train_ratio)
    c_va = int(n_neg_total * val_ratio)
    chunk_tr = neg_shuf[:, :c_tr]
    chunk_va = neg_shuf[:, c_tr:c_tr + c_va]
    chunk_te = neg_shuf[:, c_tr + c_va:]
    train_neg = chunk_tr[:, :min(neg_train_ratio * train_pos.shape[1], chunk_tr.shape[1])]
    val_neg = chunk_va[:, :min(neg_eval_ratio * val_pos.shape[1], chunk_va.shape[1])]
    test_neg = chunk_te[:, :min(neg_eval_ratio * test_pos.shape[1], chunk_te.shape[1])]
    return {
        "train": {"pos": train_pos, "neg": train_neg},
        "val": {"pos": val_pos, "neg": val_neg},
        "test": {"pos": test_pos, "neg": test_neg},
    }
