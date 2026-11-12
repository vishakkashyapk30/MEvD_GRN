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
import argparse, json, sys
from pathlib import Path
import numpy as np
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.benchmarks import scmgrn_protocol as proto   # noqa: E402
from src.utils.io import load_config, save_json       # noqa: E402


def family_scores(A, F, k=5):
    """F: TF x 768 Geneformer embeddings. 'Is v already linked (in training) to TFs that look like u in Geneformer?'"""
    Z = F - F.mean(0, keepdims=True)
    Z = Z / (np.linalg.norm(Z, axis=1, keepdims=True) + 1e-12)
    G = Z @ Z.T
    np.fill_diagonal(G, -np.inf)
    K = np.zeros_like(G)
    top = np.argsort(-G, axis=1)[:, :k]
    np.put_along_axis(K, top, 1.0 / k, axis=1)
    S = K @ A                                         # fraction of u's k Geneformer look-alikes linked to v
    return {"gf_family_top5": S + S.T}


def scores(A):
    Z = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-12)
    C = Z @ Z.T
    np.fill_diagonal(C, 0)
    S = C @ A
    Sn = S / (C.sum(1, keepdims=True) + 1e-12)
    return {"cf_lookalike": S + S.T, "cf_lookalike_norm": Sn + Sn.T}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--root", required=True)
    ap.add_argument("--split_kind", default="official", choices=proto.SPLIT_KINDS); a = ap.parse_args()
    root = Path(a.root).expanduser()
    for ds in ("lung", "hema"):
        cfg = load_config(str(REPO / f"configs/scmgrn/{ds}.yaml")); s = cfg["split"]
        for ct in cfg["cell_types"]:
            adj, nodes = proto.load_adj(str(root / cfg["graph_dir"] / ct / "x_graph.txt")); n = adj.shape[0]
            sdir = proto.resolve_split_dir(str(root / cfg["splits_dir"] / ct), s, a.split_kind)
            pdir = root / cfg["processed_dir"] / ct
            gi = json.load(open(pdir / "gene_index.json"))
            F = np.load(pdir / "fm_gene_embeddings.npy")[[gi[x.upper()] for x in nodes]].astype(np.float64)
            for fold in range(1, int(s["n_splits"]) + 1):
                f = proto.load_fold(sdir, fold, n)
                A = np.zeros((n, n)); e = f["train"]["pos"]; A[e[:, 0], e[:, 1]] = 1; A = np.maximum(A, A.T)
                np.fill_diagonal(A, 0)
                for name, S in {**scores(A), **family_scores(A, F)}.items():
                    m = proto.evaluate_all(S, adj, f, sdir, fold)
                    out = root / proto.results_subdir(cfg["results_dir"], a.split_kind) / ct / name / "seed0"; out.mkdir(parents=True, exist_ok=True)
                    save_json({"method": name, "cell_type": ct, "dataset": ds, "fold": fold, "seed": 0,
                               "note": "training links only", "metrics": m}, out / f"fold{fold:02d}.json")
            print("[done]", ds, ct, flush=True)


if __name__ == "__main__":
    main()
