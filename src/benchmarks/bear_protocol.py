"""Fair supervised protocol for scoring MeVD-GRN inside BEAR-GRN's tables.

BEAR-GRN's methods are unsupervised: each emits one ranked TF->target edge list
per dataset, and the list is scored against a ground truth (GT). A supervised
model can only be put in the same table if no label it was trained on belongs
to the edges it is scored on. We use TF-disjoint cross-fitting:

  1. The TFs of the scored GT are partitioned into K folds (degree-stratified,
     seeded). For fold k, every label of a fold-k TF is removed from EVERY
     training label source (all SC-MO-GRN-DB tiers, all GTs), so a held-out TF
     contributes no positive and no negative anywhere in training.
  2. A model is trained on the remaining TFs (validation = a TF-disjoint slice
     of the training TFs, never the held-out fold).
  3. The model scores every candidate edge of the fold-k TFs.
  4. The K held-out score blocks are concatenated into ONE edge list covering
     all GT TFs. That list is then scored by BEAR-GRN's metric code exactly as
     any unsupervised method's output is (src/benchmarks/bear_metrics.py).

Training negatives are degree-matched (Stock et al. 2025): a negative
(tf, tgt) is drawn with tf ~ out-degree and tgt ~ in-degree of the training
positives, so a model cannot separate positives from negatives by target
hubness alone. Uniform negatives are available as an ablation.

Everything here is pure numpy on integer (tf_idx, gene_idx) pairs.
"""
from __future__ import annotations

from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


def edge_keys(edges: np.ndarray, n_genes: int) -> np.ndarray:
    """(2, E) or (E, 2) int array -> int64 keys tf * n_genes + gene."""
    e = np.asarray(edges, dtype=np.int64)
    if e.ndim == 2 and e.shape[0] != 2 and e.shape[1] == 2:
        e = e.T
    if e.size == 0:
        return np.zeros(0, dtype=np.int64)
    return e[0] * np.int64(n_genes) + e[1]


def keys_to_edges(keys: np.ndarray, n_genes: int) -> np.ndarray:
    keys = np.asarray(keys, dtype=np.int64)
    return np.stack([keys // n_genes, keys % n_genes]).astype(np.int64)


def tf_disjoint_folds(tfs: Sequence[int], out_degree: Dict[int, int], k: int,
                      seed: int) -> List[np.ndarray]:
    """Degree-stratified TF partition: sort TFs by GT out-degree, cut into
    blocks of size k, shuffle within each block and deal one TF per fold. Every
    fold then spans the whole degree range (hub TFs are not all in one fold)."""
    rng = np.random.default_rng(seed)
    tfs = np.asarray(sorted(set(int(t) for t in tfs)), dtype=np.int64)
    k = max(1, min(k, len(tfs)))
    order = tfs[np.argsort([-out_degree.get(int(t), 0) for t in tfs], kind="stable")]
    folds: List[List[int]] = [[] for _ in range(k)]
    for s in range(0, len(order), k):
        block = order[s:s + k].copy()
        rng.shuffle(block)
        slots = rng.permutation(k)[:len(block)]
        for t, f in zip(block, slots):
            folds[f].append(int(t))
    return [np.asarray(sorted(f), dtype=np.int64) for f in folds]


def drop_tfs(edges: np.ndarray, tfs: Iterable[int]) -> np.ndarray:
    """Remove every edge whose source TF is in `tfs`. edges: (2, E)."""
    e = np.asarray(edges, dtype=np.int64).reshape(2, -1)
    bad = np.isin(e[0], np.asarray(list(tfs), dtype=np.int64))
    return e[:, ~bad]


def keep_tfs(edges: np.ndarray, tfs: Iterable[int]) -> np.ndarray:
    e = np.asarray(edges, dtype=np.int64).reshape(2, -1)
    return e[:, np.isin(e[0], np.asarray(list(tfs), dtype=np.int64))]


def degree_matched_negatives(pos: np.ndarray, n: int, n_genes: int,
                             forbid_keys: np.ndarray, rng: np.random.Generator,
                             target_pool: np.ndarray | None = None,
                             smoothing: float = 1.0, max_rounds: int = 50) -> np.ndarray:
    """Draw `n` unique negatives (2, n) with source ~ out-degree(pos) and
    target ~ in-degree(pos) + `smoothing` (restricted to `target_pool` if
    given). `forbid_keys` (int64 keys) = every known positive of any label
    source plus any pair that must not be trained on."""
    pos = np.asarray(pos, dtype=np.int64).reshape(2, -1)
    if pos.shape[1] == 0 or n <= 0:
        return np.zeros((2, 0), dtype=np.int64)
    src_vals, src_cnt = np.unique(pos[0], return_counts=True)
    p_src = src_cnt / src_cnt.sum()
    tgt_space = (np.arange(n_genes, dtype=np.int64) if target_pool is None
                 else np.asarray(target_pool, dtype=np.int64))
    indeg = np.bincount(pos[1], minlength=n_genes).astype(np.float64)[tgt_space] + smoothing
    p_tgt = indeg / indeg.sum()
    forbid = np.unique(np.asarray(forbid_keys, dtype=np.int64))
    got = np.zeros(0, dtype=np.int64)
    for _ in range(max_rounds):
        m = int((n - got.size) * 1.3) + 64
        s = rng.choice(src_vals, size=m, p=p_src)
        t = rng.choice(tgt_space, size=m, p=p_tgt)
        keys = s * np.int64(n_genes) + t
        keys = keys[(s != t)]
        keys = keys[~np.isin(keys, forbid)]
        got = np.unique(np.concatenate([got, keys]))
        if got.size >= n:
            break
    if got.size > n:
        got = rng.choice(got, size=n, replace=False)
    return keys_to_edges(np.sort(got), n_genes)


def uniform_negatives(src_tfs: np.ndarray, n: int, n_genes: int, forbid_keys: np.ndarray,
                      rng: np.random.Generator, target_pool: np.ndarray | None = None,
                      max_rounds: int = 50) -> np.ndarray:
    """Ablation: uniform random TF x gene negatives over the training TFs."""
    src_tfs = np.asarray(src_tfs, dtype=np.int64)
    tgt_space = (np.arange(n_genes, dtype=np.int64) if target_pool is None
                 else np.asarray(target_pool, dtype=np.int64))
    forbid = np.unique(np.asarray(forbid_keys, dtype=np.int64))
    got = np.zeros(0, dtype=np.int64)
    for _ in range(max_rounds):
        m = int((n - got.size) * 1.3) + 64
        s = rng.choice(src_tfs, size=m)
        t = rng.choice(tgt_space, size=m)
        keys = (s * np.int64(n_genes) + t)[s != t]
        keys = keys[~np.isin(keys, forbid)]
        got = np.unique(np.concatenate([got, keys]))
        if got.size >= n:
            break
    if got.size > n:
        got = rng.choice(got, size=n, replace=False)
    return keys_to_edges(np.sort(got), n_genes)


def assert_tf_disjoint(train_label_edges: Dict[str, np.ndarray], held_tfs: Iterable[int],
                       what: str = "") -> None:
    """Hard audit: no training label (pos or neg, any source) may have a
    held-out TF as its source."""
    held = np.asarray(list(held_tfs), dtype=np.int64)
    for name, e in train_label_edges.items():
        e = np.asarray(e, dtype=np.int64).reshape(2, -1)
        n_bad = int(np.isin(e[0], held).sum())
        assert n_bad == 0, f"{what}: {n_bad} '{name}' training labels come from held-out TFs"


def candidate_block(tfs: np.ndarray, targets: np.ndarray) -> np.ndarray:
    """All (tf, target) pairs, tf in tfs, target in targets, tf != target -> (2, E)."""
    tfs = np.asarray(tfs, dtype=np.int64)
    targets = np.asarray(targets, dtype=np.int64)
    s = np.repeat(tfs, targets.size)
    t = np.tile(targets, tfs.size)
    keep = s != t
    return np.stack([s[keep], t[keep]])
