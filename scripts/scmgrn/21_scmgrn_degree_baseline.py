#!/usr/bin/env python
"""Trivial TF-degree baselines on the scMultiomeGRN benchmark (no features, no cell data, no neural network).

The benchmark's folds split LINKS, not TFs, so every TF appears in training. A pair can therefore be scored from how many
TRAINING links its two TFs already have. For each cell type and each of the 10 official folds:

  deg_u = number of the fold's TRAINING positive links of TF u (test and validation positives are never used)
  degree_product : score(u, v) = (deg_u + 1) * (deg_v + 1)
  degree_logsum  : score(u, v) = log1p(deg_u) + log1p(deg_v)          (same ranking as the product; kept as a sanity check)
  degree_logreg  : logistic regression on [log1p deg_u + log1p deg_v, |log1p deg_u - log1p deg_v|, log1p deg_u * log1p deg_v]
                   fitted on the fold's training positives and training negatives only

Scored with the same code as every other method (src/benchmarks/scmgrn_protocol.py: the official metric_fn on the held-out
test pairs, and on the full matrix). Writes results/<ds>/<ct>/<variant>/seed0/foldNN.json like the other methods.

Usage: python scripts/scmgrn/21_scmgrn_degree_baseline.py --root ~/scmgrn
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.benchmarks import scmgrn_protocol as proto   # noqa: E402
from src.utils.io import load_config, save_json       # noqa: E402


def pair_feats(ld: np.ndarray, u: np.ndarray, v: np.ndarray) -> np.ndarray:
    a, b = ld[u], ld[v]
    return np.stack([a + b, np.abs(a - b), a * b], axis=1)


def run(root: Path, ds: str, cfg_path: str):
    cfg = load_config(cfg_path)
    s = cfg["split"]
    for ct in cfg["cell_types"]:
        gdir = root / cfg["graph_dir"] / ct
        sdir = proto.split_dir(str(root / cfg["splits_dir"] / ct), int(s["seed"]), int(s["n_splits"]),
                               float(s["train_split_rate"]), bool(s["rest_all_train"]))
        if not (gdir / "x_graph.txt").exists():
            print(f"[skip] {ds}/{ct}: not built"); continue
        adj, _ = proto.load_adj(str(gdir / "x_graph.txt"))
        n = adj.shape[0]
        for fold in range(1, int(s["n_splits"]) + 1):
            f = proto.load_fold(sdir, fold, n)
            tr_pos = f["train"]["pos"]                              # both directions of every training link
            deg = np.bincount(tr_pos[:, 0], minlength=n).astype(float)   # one count per training link of the TF
            ld = np.log1p(deg)
            uu, vv = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
            S = {"degree_product": (deg[uu] + 1.0) * (deg[vv] + 1.0),
                 "degree_logsum": ld[uu] + ld[vv]}
            X = np.concatenate([pair_feats(ld, f["train"]["pos"][:, 0], f["train"]["pos"][:, 1]),
                                pair_feats(ld, f["train"]["neg"][:, 0], f["train"]["neg"][:, 1])])
            y = np.concatenate([np.ones(len(f["train"]["pos"])), np.zeros(len(f["train"]["neg"]))])
            lr = LogisticRegression(max_iter=1000).fit(X, y)
            S["degree_logreg"] = lr.predict_proba(pair_feats(ld, uu.ravel(), vv.ravel()))[:, 1].reshape(n, n)
            for name, M in S.items():
                M = (M + M.T) / 2.0
                m = proto.evaluate_score_matrix(M.astype(np.float64), adj, f["test_eval"])
                out = root / cfg["results_dir"] / ct / name / "seed0"
                out.mkdir(parents=True, exist_ok=True)
                save_json({"method": name, "cell_type": ct, "dataset": cfg["dataset"], "fold": fold, "seed": 0,
                           "note": "uses the fold's training links only", "metrics": m}, out / f"fold{fold:02d}.json")
        print(f"[done] {ds}/{ct}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    root = Path(a.root).expanduser()
    for ds in ("lung", "hema"):
        run(root, ds, str(REPO / f"configs/scmgrn/{ds}.yaml"))


if __name__ == "__main__":
    main()
