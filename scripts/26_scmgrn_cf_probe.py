#!/usr/bin/env python
"""Probe: how far does pure matrix completion of the training links go (no cell data, no Geneformer)?

The folds split LINKS, not TFs, so every test TF has ~90% of its links in training, and the motif-based answer key
is family-structured (TFs with similar motifs have near-identical link rows). Item-based collaborative filtering:
  C = cosine similarity between TFs' training-link rows (diagonal 0)
  cf_lookalike(u, v) = sum_w C[u, w] * A[w, v]  (+ transpose): are u's look-alike TFs linked to v?
  cf_lookalike_norm: the same divided by sum_w C[u, w] (fraction of u's look-alikes linked to v) -> less degree-driven.
Writes results/<ds>/<ct>/<method>/seed0/foldNN.json (official test + degree-matched test), like scripts/22.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from src.benchmarks import scmgrn_protocol as proto   # noqa: E402
from src.utils.io import load_config, save_json       # noqa: E402


def scores(A):
    Z = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-12)
    C = Z @ Z.T
    np.fill_diagonal(C, 0)
    S = C @ A
    Sn = S / (C.sum(1, keepdims=True) + 1e-12)
    return {"cf_lookalike": S + S.T, "cf_lookalike_norm": Sn + Sn.T}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True); a = ap.parse_args()
    root = Path(a.root).expanduser()
    for ds in ("lung", "hema"):
        cfg = load_config(str(REPO / f"configs/scmgrn/{ds}.yaml")); s = cfg["split"]
        for ct in cfg["cell_types"]:
            adj, _ = proto.load_adj(str(root / cfg["graph_dir"] / ct / "x_graph.txt")); n = adj.shape[0]
            sdir = proto.split_dir(str(root / cfg["splits_dir"] / ct), int(s["seed"]), int(s["n_splits"]),
                                   float(s["train_split_rate"]), bool(s["rest_all_train"]))
            for fold in range(1, int(s["n_splits"]) + 1):
                f = proto.load_fold(sdir, fold, n)
                A = np.zeros((n, n)); e = f["train"]["pos"]; A[e[:, 0], e[:, 1]] = 1; A = np.maximum(A, A.T)
                np.fill_diagonal(A, 0)
                dm = proto.degree_matched_test(adj, f, seed=fold)
                for name, S in scores(A).items():
                    m = proto.evaluate_score_matrix(S, adj, f["test_eval"])
                    m["test_degmatched"] = proto.evaluate_degree_matched(S, dm)
                    out = root / cfg["results_dir"] / ct / name / "seed0"; out.mkdir(parents=True, exist_ok=True)
                    save_json({"method": name, "cell_type": ct, "dataset": ds, "fold": fold, "seed": 0,
                               "note": "training links only", "metrics": m}, out / f"fold{fold:02d}.json")
            print("[done]", ds, ct, flush=True)


if __name__ == "__main__":
    main()
