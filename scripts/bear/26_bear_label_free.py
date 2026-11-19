#!/usr/bin/env python
"""BEAR-GRN: label-free, non-graph-neural models (2026-11-19, advisor item 5).

No label of any kind is used: one model per dataset scores every (TF node, gene) pair, like the released BEAR methods.
  lf_ggm_ridge   Gaussian graphical model by neighbourhood regression (Meinshausen & Buehlmann, Ann. Stat. 2006),
                 ridge form: every gene's standardised expression regressed on all TF nodes,
                 beta = (X_T' X_T + lam I)^-1 X_T' z_g, score(t, g) = |beta_t|. When g is itself a TF it is left out of
                 its own regressors, exactly, via the block-inverse identity beta_{-g} = -Theta[-g, g] / Theta[g, g]
                 with Theta = (X_T' X_T + lam I)^-1. lam = 0.1 x number of cells.
  lf_regdiffusion RegDiffusion (Zhu & Slonim, J. Comput. Biol. 2024): diffusion model whose denoiser contains a
                 learnable gene x gene adjacency; score(t, g) = |A[t, g]| (rows = regulators).
Expression = CP10k + log1p over the dataset's gene universe (scripts/23 expr_matrix).
Outputs: results/<ds>/<model>/seed0/grn.tsv.gz.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from src.utils.io import load_config, load_json            # noqa: E402

_spec = importlib.util.spec_from_file_location("b23", REPO / "scripts/bear/23_bear_baselines.py")
b23 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(b23)


def write(out_dir: Path, names, tf_nodes, S: np.ndarray, meta: dict):
    out_dir.mkdir(parents=True, exist_ok=True)
    nm = np.asarray(names, dtype=object)
    n = S.shape[1]
    src = np.repeat(tf_nodes, n)
    dst = np.tile(np.arange(n), len(tf_nodes))
    keep = src != dst
    df = pd.DataFrame({"Source": nm[src[keep]], "Target": nm[dst[keep]], "Score": S.reshape(-1)[keep]})
    df.to_csv(out_dir / "grn.tsv.gz", sep="\t", index=False, float_format="%.4g")
    with open(out_dir / "meta.json", "w") as f:
        json.dump(meta, f, indent=1)
    print(f"[done] {out_dir} ({len(df):,} edges)", flush=True)


def ggm_ridge(X: np.ndarray, tf_nodes: np.ndarray, lam_frac: float = 0.1) -> np.ndarray:
    Z = (X - X.mean(0)) / (X.std(0) + 1e-8)
    Xt = Z[:, tf_nodes]
    lam = lam_frac * Z.shape[0]
    Theta = np.linalg.inv(Xt.T @ Xt + lam * np.eye(len(tf_nodes)))
    B = Theta @ (Xt.T @ Z)                                   # TF x genes
    for j, g in enumerate(tf_nodes):                          # TF targets: leave the target out of its regressors
        col = -Theta[:, j] / Theta[j, j]
        col[j] = 0.0
        B[:, g] = col
    return np.abs(B).astype(np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--config", default=str(REPO / "configs/bear/datasets.yaml"))
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--models", nargs="+", default=["ggm_ridge", "regdiffusion"])
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--n_steps", type=int, default=1000)
    args = ap.parse_args()
    root = Path(args.root)
    ds_cfg = load_config(args.config)["datasets"][args.dataset]
    pdir = root / "processed" / args.dataset
    gene_index = load_json(pdir / "gene_index.json")
    names = [None] * len(gene_index)
    for g, i in gene_index.items():
        names[i] = g
    tf_nodes = np.asarray(load_json(pdir / "tf_indices.json"), dtype=np.int64)
    X = b23.expr_matrix(root, ds_cfg, gene_index)             # cells x genes
    res = root / "results" / args.dataset
    print(f"[{args.dataset}] {X.shape[0]} cells x {X.shape[1]} genes, {len(tf_nodes)} TF nodes", flush=True)
    if "ggm_ridge" in args.models:
        t0 = time.time()
        S = ggm_ridge(X, tf_nodes)
        write(res / "lf_ggm_ridge" / "seed0", names, tf_nodes, S,
              {"model": "ggm_ridge", "lam": "0.1 x cells", "seconds": time.time() - t0, "labels_used": False})
    if "regdiffusion" in args.models:
        import regdiffusion as rd
        import torch
        t0 = time.time()
        torch.manual_seed(0)
        tr = rd.RegDiffusionTrainer(X.astype(np.float32), device=args.device, n_steps=args.n_steps,
                                    memory_efficient=True)
        tr.train()
        A = np.abs(np.asarray(tr.get_adj(), dtype=np.float32))
        write(res / "lf_regdiffusion" / "seed0", names, tf_nodes, A[tf_nodes],
              {"model": "regdiffusion", "version": rd.__version__ if hasattr(rd, "__version__") else "",
               "n_steps": args.n_steps, "seconds": time.time() - t0, "labels_used": False})


if __name__ == "__main__":
    main()
