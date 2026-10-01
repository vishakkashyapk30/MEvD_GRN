#!/usr/bin/env python
"""BEAR-GRN step 3: baselines emitted in the same edge-list format as MeVD-GRN
(results/<ds>/<baseline>/seed<s>/grn.tsv.gz), so scripts/24 scores them alike.

  indegree   supervised hub prior under the SAME TF-disjoint cross-fitting as
             scripts/22 (same fold function and seed): score(tf, g) = in-degree
             of g among the training TFs' labels. If MeVD-GRN does not beat it,
             its gain is target hubness, not data (s6.7).
  pearson    |Pearson r| of TF and target over cells (CP10k + log1p RNA), unsupervised.
  coverage   score 1 for every (TF node, RNA gene) pair: no ranking at all. Measures
             how much of BEAR's AUPRC is earned just by emitting the measured genes
             (unpredicted universe pairs score 0), i.e. the dense-output floor.
  grnboost2  arboreto GRNBoost2 importances (dask-free per-target path,
             src/benchmarks/scmgrn_features.py), TF nodes as regulators, unsupervised.

The 9 BEAR methods are not re-run: their released outputs are re-scored (scripts/24).

  python scripts/23_bear_baselines.py --root $BEAR_ROOT --dataset K562 --baselines indegree pearson --seeds 42 43
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from src.benchmarks import bear_protocol as bp             # noqa: E402
from src.data.preprocessing import _normalize_log, extract_symbol, load_features_by_cells  # noqa: E402
from src.utils.io import load_config, load_json            # noqa: E402


def write(out_dir: Path, names, src, dst, score, fold=None, meta=None):
    out_dir.mkdir(parents=True, exist_ok=True)
    nm = np.asarray(names, dtype=object)
    df = pd.DataFrame({"Source": nm[src], "Target": nm[dst], "Score": score})
    if fold is not None:
        df["fold"] = fold
    df.to_csv(out_dir / "grn.tsv.gz", sep="\t", index=False, float_format="%.6g")
    with open(out_dir / "meta.json", "w") as f:
        json.dump(meta or {}, f, indent=1)
    print(f"[done] {out_dir}/grn.tsv.gz ({len(df):,} edges)", flush=True)


def expr_matrix(root: Path, ds_cfg: dict, gene_index: dict) -> np.ndarray:
    """cells x universe-genes, CP10k + log1p (same normalisation as scripts/21)."""
    X, raw = load_features_by_cells(str(root / ds_cfg["rna"]))
    X = _normalize_log(X, {"normalize_total": True, "target_sum": 1e4, "log1p": True})
    names = [extract_symbol(g) for g in raw]
    col = {}
    for j, g in enumerate(names):
        col.setdefault(g, j)
    idx = np.array([col[g] for g in sorted(gene_index, key=gene_index.get)])
    return np.asarray(X[:, idx].todense(), dtype=np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.environ.get("BEAR_ROOT", "data/bear"))
    ap.add_argument("--config", default=str(REPO / "configs/bear/datasets.yaml"))
    ap.add_argument("--model_config", default=str(REPO / "configs/bear/mevd_base.yaml"),
                    help="for bear.n_folds (must match scripts/22)")
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--baselines", nargs="+", default=["indegree", "pearson"])
    ap.add_argument("--seeds", nargs="+", type=int, default=[42, 43, 44, 45, 46])
    ap.add_argument("--n_jobs", type=int, default=8)
    args = ap.parse_args()
    root = Path(args.root)
    ds_cfg = load_config(args.config)["datasets"][args.dataset]
    K = int(load_config(args.model_config)["bear"]["n_folds"])
    pdir = root / "processed" / args.dataset
    gene_index = load_json(pdir / "gene_index.json")
    names = [None] * len(gene_index)
    for g, i in gene_index.items():
        names[i] = g
    n = len(names)
    tf_nodes = np.asarray(load_json(pdir / "tf_indices.json"), dtype=np.int64)
    labels = [torch.load(p).numpy().reshape(2, -1) for p in sorted((pdir / "labels").glob("*.pt"))]
    lab = bp.keys_to_edges(np.unique(np.concatenate([bp.edge_keys(e, n) for e in labels])), n)
    outdeg = {int(k): int(v) for k, v in zip(*np.unique(lab[0], return_counts=True))}
    res = root / "results" / args.dataset

    if "indegree" in args.baselines:
        for s in args.seeds:
            folds = bp.tf_disjoint_folds(tf_nodes, outdeg, K, s)
            S, D, SC, F = [], [], [], []
            for k, held in enumerate(folds):
                train = bp.drop_tfs(lab, held)
                bp.assert_tf_disjoint({"labels": train}, held, what=f"indegree fold {k}")
                indeg = np.bincount(train[1], minlength=n).astype(np.float32)
                pairs = bp.candidate_block(held, np.arange(n))
                S.append(pairs[0]); D.append(pairs[1]); SC.append(indeg[pairs[1]])
                F.append(np.full(pairs.shape[1], k, dtype=np.int16))
            write(res / "baseline_indegree" / f"seed{s}", names, np.concatenate(S), np.concatenate(D),
                  np.concatenate(SC), np.concatenate(F),
                  {"baseline": "indegree", "seed": s, "folds": [f.tolist() for f in folds]})

    if "coverage" in args.baselines:
        pairs = bp.candidate_block(tf_nodes, np.arange(n))
        write(res / "baseline_coverage" / "seed0", names, pairs[0], pairs[1],
              np.ones(pairs.shape[1], dtype=np.float32), meta={"baseline": "coverage"})
    if "pearson" in args.baselines or "grnboost2" in args.baselines:
        X = expr_matrix(root, ds_cfg, gene_index)
    if "pearson" in args.baselines:
        Z = (X - X.mean(0)) / (X.std(0) + 1e-8)
        R = (Z[:, tf_nodes].T @ Z) / X.shape[0]            # TF x genes
        R[np.arange(len(tf_nodes)), tf_nodes] = 0.0
        pairs = bp.candidate_block(tf_nodes, np.arange(n))
        pos = {int(t): i for i, t in enumerate(tf_nodes)}
        sc = np.abs(R[[pos[int(t)] for t in pairs[0]], pairs[1]])
        write(res / "baseline_pearson" / "seed0", names, pairs[0], pairs[1], sc,
              meta={"baseline": "pearson_abs"})
    if "grnboost2" in args.baselines:
        from src.benchmarks.scmgrn_features import grnboost2_features
        t0 = time.time()
        expr = pd.DataFrame(X.T, index=names)
        net = grnboost2_features(expr, [names[i] for i in tf_nodes], seed=0, n_jobs=args.n_jobs)
        st = net.stack()
        st = st[st > 0]
        src = np.array([gene_index[a] for a in st.index.get_level_values(0)])
        dst = np.array([gene_index[b] for b in st.index.get_level_values(1)])
        keep = src != dst
        write(res / "baseline_grnboost2" / "seed0", names, src[keep], dst[keep],
              st.to_numpy()[keep], meta={"baseline": "grnboost2", "seconds": time.time() - t0})


if __name__ == "__main__":
    main()
