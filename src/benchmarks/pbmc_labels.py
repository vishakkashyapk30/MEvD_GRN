"""Training labels, splits and negatives for the PBMC10k / LINGER-Cistrome benchmark.

Pre-registered design (docs/experiments/pbmc10k_linger_benchmark.md, section 4):

* Label sources (none is the evaluation ChIP set, none is motif-derived):
    collectri   CollecTRI (OmniPath `datasets=collectri`), literature-curated.
                Multi-member complexes (AP-1/NF-kB heteromers) are dropped;
                single-member "COMPLEX:<uniprot>" rows map to that gene.
    dorothea_ab DoRothEA confidence levels A and B only (C-E are dominated by
                ChIP-only / motif-only / inferred evidence).
* Every TF in the LINGER Cistrome evaluation set is removed as a REGULATOR
  (TF-disjoint w.r.t. the evaluation). It may still appear as a target.
* Split regimes over the remaining label set (70/15/15, seeded):
    tf      TF-disjoint     (headline; unit = regulator)
    target  target-disjoint (unit = target gene; train pairs only touch train targets)
    random  edge split      (EpiAwareNet-style, comparability only)
    target_all  (amendment A8, 2026-10-08) target-disjoint over the WHOLE gene
            universe: the labelled-target partition is exactly `target`'s (same
            positives in every part), and the never-labelled genes are
            partitioned too, so train negatives may use never-labelled genes of
            the train partition and val/test negatives those of their own part.
            In `target`, train negatives only ever use the ~3.6k train targets,
            so ~20k never-labelled genes are never seen as negatives at all.
* Negatives:
    train   degree-matched pool: each positive (t, g) spawns `pool_mult * neg_ratio`
            negatives (t, g') with the same TF (TF out-degree preserved) and g'
            drawn with probability proportional to g's in-degree among the train
            positives; a `uniform_frac` share draws g' uniformly so genes with no
            known regulator are still seen as negatives. The trainer samples
            uniformly from this pool each epoch.
    val     uniform same-TF negatives (fixed sample, `val_neg_ratio`:1) --
            mirrors the headline per-TF ranking over all genes.
    test    uniform (as val) AND degree-matched 1:1 (Stock et al. 2025) sets.
  Negatives never contain a known positive of the label source (any split),
  never a self pair, and for the target regime only use targets of that part.
"""
from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

SOURCES = ("collectri", "dorothea_ab")
REGIMES = ("tf", "target", "random", "target_all")


def canon(x: str) -> str:
    return str(x).strip().upper()


# ----------------------------------------------------------------------------- loading
def load_label_source(name: str, path: str) -> pd.DataFrame:
    """Return a de-duplicated DataFrame[tf, target] of UPPERCASE symbols."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    if name == "collectri":
        multi = df["source"].str.startswith("COMPLEX:") & df["source"].str.contains("_")
        df = df[~multi]
    elif name == "dorothea_ab":
        lev = df["dorothea_level"].fillna("")
        df = df[lev.str.contains("A") | lev.str.contains("B")]
    else:
        raise ValueError(f"unknown label source {name!r}; expected one of {SOURCES}")
    out = pd.DataFrame({"tf": df["source_genesymbol"].map(canon),
                        "target": df["target_genesymbol"].map(canon)})
    out = out[(out.tf != out.target) & ~out.tf.str.contains("_") & ~out.target.str.contains("_")]
    return out.drop_duplicates().reset_index(drop=True)


def restrict(df: pd.DataFrame, genes: Iterable[str], exclude_tfs: Iterable[str]) -> Tuple[pd.DataFrame, dict]:
    genes = set(genes)
    ex = {canon(t) for t in exclude_tfs}
    n0 = len(df)
    in_uni = df[df.tf.isin(genes) & df.target.isin(genes)]
    dropped_eval = in_uni[in_uni.tf.isin(ex)]
    kept = in_uni[~in_uni.tf.isin(ex)].reset_index(drop=True)
    info = {"edges_raw": n0, "edges_in_universe": len(in_uni),
            "edges_removed_eval_tfs": len(dropped_eval),
            "eval_tfs_present_as_regulator": sorted(set(dropped_eval.tf)),
            "edges": len(kept), "tfs": int(kept.tf.nunique()), "targets": int(kept.target.nunique())}
    return kept, info


# ----------------------------------------------------------------------------- splits
def _partition(units: List[str], ratios, rng) -> Dict[str, set]:
    units = sorted(units)
    perm = rng.permutation(len(units))
    n_tr = int(round(ratios[0] * len(units)))
    n_va = int(round(ratios[1] * len(units)))
    idx = {"train": perm[:n_tr], "val": perm[n_tr:n_tr + n_va], "test": perm[n_tr + n_va:]}
    return {p: {units[i] for i in ix} for p, ix in idx.items()}


def split_labels(df: pd.DataFrame, regime: str, seed: int,
                 ratios=(0.70, 0.15, 0.15)) -> Dict[str, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    if regime == "tf":
        part = _partition(df.tf.unique().tolist(), ratios, rng)
        return {p: df[df.tf.isin(s)].reset_index(drop=True) for p, s in part.items()}
    if regime in ("target", "target_all"):        # identical labelled-target partition
        part = _partition(df.target.unique().tolist(), ratios, rng)
        return {p: df[df.target.isin(s)].reset_index(drop=True) for p, s in part.items()}
    if regime == "random":
        perm = rng.permutation(len(df))
        n_tr = int(round(ratios[0] * len(df)))
        n_va = int(round(ratios[1] * len(df)))
        ix = {"train": perm[:n_tr], "val": perm[n_tr:n_tr + n_va], "test": perm[n_tr + n_va:]}
        return {p: df.iloc[np.sort(i)].reset_index(drop=True) for p, i in ix.items()}
    raise ValueError(f"regime must be one of {REGIMES}")


def to_index(df: pd.DataFrame, gene_index: Dict[str, int]) -> np.ndarray:
    if len(df) == 0:
        return np.zeros((0, 2), dtype=np.int64)
    return np.stack([df.tf.map(gene_index).to_numpy(), df.target.map(gene_index).to_numpy()],
                    axis=1).astype(np.int64)


# ----------------------------------------------------------------------------- negatives
def sample_same_tf_negatives(pos: np.ndarray, n_genes: int, forbidden: np.ndarray,
                             ratio: int, rng: np.random.Generator,
                             target_probs: Optional[np.ndarray] = None,
                             allowed_targets: Optional[np.ndarray] = None,
                             max_rounds: int = 50) -> np.ndarray:
    """For each positive (t, g) draw `ratio` negatives (t, g'): g' ~ target_probs
    (uniform if None) over `allowed_targets` (all genes if None), rejecting self
    pairs and any key in `forbidden` (sorted int64 keys t*n_genes+g). Returns (M, 2)."""
    if len(pos) == 0 or ratio <= 0:
        return np.zeros((0, 2), dtype=np.int64)
    allowed = np.arange(n_genes) if allowed_targets is None else np.asarray(allowed_targets)
    p = None
    if target_probs is not None:
        p = np.asarray(target_probs, dtype=np.float64)[allowed]
        p = p / p.sum()
    tfs = np.repeat(pos[:, 0], ratio)
    tg = np.full(tfs.shape[0], -1, dtype=np.int64)
    todo = np.arange(tfs.shape[0])
    for _ in range(max_rounds):
        if todo.size == 0:
            break
        draw = allowed[rng.choice(allowed.size, size=todo.size, p=p)]
        key = tfs[todo] * n_genes + draw
        bad = (draw == tfs[todo]) | _isin_sorted(key, forbidden)
        tg[todo[~bad]] = draw[~bad]
        todo = todo[bad]
    ok = tg >= 0
    return np.stack([tfs[ok], tg[ok]], axis=1)


def _isin_sorted(keys: np.ndarray, sorted_keys: np.ndarray) -> np.ndarray:
    if sorted_keys.size == 0:
        return np.zeros(keys.shape, dtype=bool)
    i = np.searchsorted(sorted_keys, keys)
    i[i >= sorted_keys.size] = sorted_keys.size - 1
    return sorted_keys[i] == keys


def keys_of(pairs: np.ndarray, n_genes: int) -> np.ndarray:
    return np.unique(pairs[:, 0].astype(np.int64) * n_genes + pairs[:, 1].astype(np.int64))


def degree_matched_probs(train_pos: np.ndarray, n_genes: int, uniform_frac: float) -> np.ndarray:
    indeg = np.bincount(train_pos[:, 1], minlength=n_genes).astype(np.float64)
    p_deg = indeg / max(indeg.sum(), 1.0)
    return (1.0 - uniform_frac) * p_deg + uniform_frac / n_genes


def unlabelled_gene_partition(P: Dict[str, np.ndarray], n_genes: int, seed: int,
                              ratios=(0.70, 0.15, 0.15)) -> Dict[str, np.ndarray]:
    """target_all: each part's labelled targets plus a seeded 70/15/15 share of
    the genes that are a target in no part (sorted index arrays)."""
    labelled = np.unique(np.concatenate([P[p][:, 1] for p in P]))
    rest = np.setdiff1d(np.arange(n_genes), labelled)
    perm = np.random.default_rng(seed + 2000).permutation(rest.size)
    n_tr = int(round(ratios[0] * rest.size))
    n_va = int(round(ratios[1] * rest.size))
    extra = {"train": rest[perm[:n_tr]], "val": rest[perm[n_tr:n_tr + n_va]],
             "test": rest[perm[n_tr + n_va:]]}
    return {p: np.union1d(np.unique(P[p][:, 1]), extra[p]) for p in P}


def build_split_arrays(parts: Dict[str, pd.DataFrame], all_pos: pd.DataFrame, regime: str,
                       gene_index: Dict[str, int], seed: int, neg_ratio: int = 5,
                       pool_mult: int = 4, uniform_frac: float = 0.2,
                       val_neg_ratio: int = 20) -> Dict[str, Dict[str, np.ndarray]]:
    """Index arrays (M, 2) per part:
       train: pos, neg (degree-matched POOL for the trainer)
       val:   pos, neg (uniform same-TF)
       test:  pos, neg (uniform same-TF), neg_dm (degree-matched 1:1)."""
    n = len(gene_index)
    rng = np.random.default_rng(seed + 1000)
    forb = keys_of(to_index(all_pos, gene_index), n)
    P = {p: to_index(parts[p], gene_index) for p in ("train", "val", "test")}
    allowed = {p: None for p in P}
    if regime == "target":
        allowed = {p: np.unique(P[p][:, 1]) for p in P}
    elif regime == "target_all":
        allowed = unlabelled_gene_partition(P, n, seed)
    probs_tr = degree_matched_probs(P["train"], n, uniform_frac)
    out = {"train": {"pos": P["train"],
                     "neg": sample_same_tf_negatives(P["train"], n, forb, neg_ratio * pool_mult, rng,
                                                     probs_tr, allowed["train"]),
                     # uniform same-TF pool, only for the *_uniformneg sensitivity variant
                     "neg_uniform": sample_same_tf_negatives(P["train"], n, forb, neg_ratio * pool_mult,
                                                             rng, None, allowed["train"])}}
    for p in ("val", "test"):
        out[p] = {"pos": P[p],
                  "neg": sample_same_tf_negatives(P[p], n, forb, val_neg_ratio, rng, None, allowed[p])}
    # degree-matched test negatives use the TEST part's own in-degree (Stock et al.)
    probs_te = degree_matched_probs(P["test"], n, 0.0) if len(P["test"]) else None
    out["test"]["neg_dm"] = sample_same_tf_negatives(P["test"], n, forb, 1, rng, probs_te, allowed["test"])
    # the three parts' negatives must be disjoint from each other's positives by construction
    return out


def check_disjoint(arrays: Dict[str, Dict[str, np.ndarray]], regime: str, n: int) -> dict:
    """Hard assertions for the leakage rules; returns simple counts for the log."""
    tr = arrays["train"]
    tr_keys = keys_of(tr["pos"], n)
    for k in ("neg", "neg_uniform"):
        if k in tr and len(tr[k]):
            tr_keys = np.union1d(tr_keys, keys_of(tr[k], n))
    rep = {}
    for p in ("val", "test"):
        for k in ("pos", "neg", "neg_dm"):
            if k not in arrays[p] or len(arrays[p][k]) == 0:
                continue
            ov = np.intersect1d(keys_of(arrays[p][k], n), tr_keys).size
            if k == "pos":
                assert ov == 0, f"{p}.{k}: {ov} pairs also in train"
            rep[f"{p}_{k}_overlap_train"] = int(ov)
        if regime == "tf":
            assert not (set(arrays[p]["pos"][:, 0]) & set(tr["pos"][:, 0])), f"{p} TFs overlap train TFs"
        if regime in ("target", "target_all"):
            assert not (set(arrays[p]["pos"][:, 1]) & set(tr["pos"][:, 1])), f"{p} targets overlap train"
            for k in ("neg", "neg_uniform"):
                assert not (set(tr[k][:, 1]) & set(arrays[p]["pos"][:, 1])), f"train {k} touch {p} targets"
        if regime == "target_all":            # no held-out-part gene is ever a train target
            held = set(arrays[p]["pos"][:, 1]) | set(arrays[p]["neg"][:, 1])
            for k in ("pos", "neg", "neg_uniform"):
                assert not (set(tr[k][:, 1]) & held), f"train {k} touch a {p}-partition gene"
    return rep
