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
  lf_deepsem     DeepSEM (Shu et al., Nat. Comput. Sci. 2021): beta-VAE whose encoder/decoder is a structural equation
                 model with a sparse gene x gene weight matrix W (L1 penalty); score(t, g) = |W[g, t]|. The authors'
                 model code (github.com/HantaoShu/DeepSEM, --deepsem_dir) is used unchanged with the default
                 hyper-parameters of their main.py (120 epochs, MLP and W optimised alternately, alpha 100, beta 1,
                 lr 1e-4); the final W is kept (their ground-truth file only prints progress, so none is given:
                 the TF set is the dataset's TF nodes). Genes = TF nodes + the 2,000 most variable genes, as in
                 the DeepSEM paper; other genes are not output (score 0).
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


def deepsem(X: np.ndarray, tf_nodes: np.ndarray, deepsem_dir: str, n_hvg: int = 2000, epochs: int = 120,
            seed: int = 0):
    import torch
    from torch.utils.data import DataLoader, TensorDataset
    spec = importlib.util.spec_from_file_location(               # authors' model file, unchanged (loaded by path:
        "deepsem_model", Path(deepsem_dir).expanduser().resolve() / "src" / "Model.py")   # this repo has a src too)
    dm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dm)
    VAE_EAD = dm.VAE_EAD
    torch.manual_seed(seed)
    np.random.seed(seed)
    var = X.var(0)
    genes = np.union1d(tf_nodes, np.argsort(-var)[:n_hvg])
    D = X[:, genes]
    D = (D - D.mean(0)) / (D.std(0) + 1e-8)
    ng = len(genes)
    A0 = np.ones([ng, ng]) / (ng - 1) + (np.random.rand(ng * ng) * 0.0002).reshape([ng, ng])
    np.fill_diagonal(A0, 0)
    is_tf = np.isin(genes, tf_nodes)

    class Opt:                                                  # main.py defaults
        alpha, beta, lr, lr_step_size, gamma, n_hidden, K, K1, K2, batch_size = 100, 1, 1e-4, 0.99, 0.95, 128, 1, 1, 2, 64
    opt = Opt()
    vae = VAE_EAD(A0, 1, opt.n_hidden, opt.K).float().cuda()
    o1 = torch.optim.RMSprop(vae.parameters(), lr=opt.lr)
    o2 = torch.optim.RMSprop([vae.adj_A], lr=opt.lr * 0.2)
    sch = torch.optim.lr_scheduler.StepLR(o1, step_size=opt.lr_step_size, gamma=opt.gamma)
    dl = DataLoader(TensorDataset(torch.FloatTensor(D)), batch_size=opt.batch_size, shuffle=True)
    vae.train()
    for ep in range(epochs + 1):
        w_step = ep % (opt.K1 + opt.K2) >= opt.K1
        vae.adj_A.requires_grad = w_step
        for (x,) in dl:
            o1.zero_grad()
            loss = vae(x.cuda(), dropout_mask=None, temperature=max(0.95 ** ep, 0.5), opt=opt)[0]
            loss = loss + opt.alpha * torch.mean(torch.abs(vae.adj_A))
            loss.backward()
            (o2 if w_step else o1).step()
        sch.step()
    W = np.abs(vae.adj_A.detach().cpu().numpy())               # W[receiver, sender]
    S = np.zeros((len(tf_nodes), X.shape[1]), dtype=np.float32)
    pos = {int(g): i for i, g in enumerate(genes)}
    for r, t in enumerate(tf_nodes):
        S[r, genes] = W[:, pos[int(t)]]
    return S, int(ng)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--config", default=str(REPO / "configs/bear/datasets.yaml"))
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--models", nargs="+", default=["ggm_ridge", "regdiffusion"])
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--n_steps", type=int, default=1000)
    ap.add_argument("--deepsem_dir", default="~/.cache/ext/DeepSEM")
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
    if "deepsem" in args.models:
        t0 = time.time()
        S, ng = deepsem(X, tf_nodes, args.deepsem_dir)
        write(res / "lf_deepsem" / "seed0", names, tf_nodes, S,
              {"model": "deepsem", "genes_modelled": ng, "seconds": time.time() - t0, "labels_used": False})
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
