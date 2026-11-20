#!/usr/bin/env python
"""LINGER (LingerGRN 1.110, method 'LINGER') re-run on one BEAR-GRN dataset under an ATAC condition (2026-11-19,
advisor item 3: how much does ATAC add to LINGER?). Runs in the separate `linger` env; imports nothing from this repo.

Inputs are BEAR-GRN's own filtered matrices (INPUT.DATA: genes x cells RNA, peaks x cells ATAC), the ones every
released BEAR method was run on. Steps follow LINGER's PBMC tutorial (as in src/benchmarks/pbmc_linger_run.py) with
one cell population: pseudo-bulk -> preprocess with the bulk reference -> training -> TF-RE binding -> cis -> trans.

ATAC conditions (--mode):
  full         unchanged
  peakshuffle  every peak takes another peak's cell profile (one fixed permutation): accessibility values and the
               cell structure are kept, the genomic location -> motif / nearby-gene mapping is broken
  random       every ATAC value permuted across all peaks and cells: no ATAC information left, same value
               distribution (a literal constant matrix breaks LINGER's ATAC PCA step)
Seeds (--seed, 2026-11-20): LINGER fixes its own training seed at 42 and pseudo-bulk draws with python's `random`.
Seed 0 keeps LINGER's 42; seed s > 0 sets python/numpy/torch seeds to s, the k-means seed to s, the ATAC permutation
seed to s and LINGER's training seed to 42 + s.
Pseudo-labels (--labels): `auto` keeps the one population when it already gives >= 50 metacells (floor(sqrt(N)) + 1,
e.g. iPS with 5,711 cells: LINGER exactly as published) and otherwise splits the cells by k-means on RNA PCA.
Output: <export_dir>/grn.tsv.gz with Source (TF), Target (gene), Score (cell-population trans-regulatory score).
"""
from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def build_adata(rna_csv: str, atac_csv: str, mode: str, seed: int = 0, max_cells: int = 0):
    import anndata
    R = pd.read_csv(rna_csv, index_col=0)
    a_cols = pd.read_csv(atac_csv, index_col=0, nrows=0).columns
    cells = [c for c in R.columns if c in set(a_cols)]
    if max_cells and len(cells) > max_cells:       # fixed cell subsample (seed 0), the same in every condition and seed
        keep = set(np.random.default_rng(0).choice(len(cells), max_cells, replace=False).tolist())
        cells = [c for i, c in enumerate(cells) if i in keep]
    A = pd.read_csv(atac_csv, index_col=0, usecols=[0] + [1 + list(a_cols).index(c) for c in cells])
    R, A = R[cells], A[cells]
    M = A.to_numpy(dtype=np.float32)
    rng = np.random.default_rng(seed)
    if mode == "peakshuffle":
        M = M[rng.permutation(M.shape[0])]
    elif mode == "random":
        M = rng.permutation(M.reshape(-1)).reshape(M.shape)
    elif mode != "full":
        raise ValueError(mode)
    peaks = [p.replace("-", ":", 1) for p in A.index.astype(str)]       # chr1-10112-10254 -> chr1:10112-10254
    a_rna = anndata.AnnData(X=sp.csc_matrix(R.to_numpy(dtype=np.float32).T))
    a_rna.var_names = R.index.astype(str)
    a_rna.var["gene_ids"] = a_rna.var_names
    a_atac = anndata.AnnData(X=sp.csc_matrix(M.T))
    a_atac.var_names = peaks
    a_atac.var["gene_ids"] = peaks
    for a in (a_rna, a_atac):
        a.obs_names = cells
        a.obs["barcode"] = cells
        a.obs["sample"] = "1"
        a.obs["label"] = "all"
    a_rna.var_names_make_unique()
    return a_rna, a_atac


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna", required=True)
    ap.add_argument("--atac", required=True)
    ap.add_argument("--grndir", required=True, help="unpacked LINGER data_bulk/ (trailing slash)")
    ap.add_argument("--outdir", required=True, help="LINGER working dir (trailing slash)")
    ap.add_argument("--export_dir", required=True)
    ap.add_argument("--mode", default="full", choices=["full", "peakshuffle", "random"])
    ap.add_argument("--genome", default="hg38")
    ap.add_argument("--species", default="Human")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--labels", default="auto", choices=["auto", "kmeans", "none"])
    ap.add_argument("--max_cells", type=int, default=0, help="random cell subsample (seed 0); 0 = all cells")
    args = ap.parse_args()
    import random
    import torch
    random.seed(args.seed); np.random.seed(args.seed); torch.manual_seed(args.seed)
    import scanpy as sc
    import LingerGRN.LINGER_tr as LINGER_tr
    import LingerGRN.LL_net as LL_net
    from LingerGRN.preprocess import preprocess
    from LingerGRN.pseudo_bulk import pseudo_bulk
    # LINGER hard-codes 8 joblib workers in preprocess (and sizes others from node RAM), which OOM-kills a 27 GB job.
    # Cap every joblib pool at LINGER_MAX_JOBS without editing the package: module-level imports and the
    # function-local "from joblib import Parallel" both resolve to the patched class. Results are unchanged
    # (the pools only map independent chromosomes).
    import joblib
    import LingerGRN.preprocess as _lp
    cap = int(os.environ.get("LINGER_MAX_JOBS", "2"))

    class _Capped(joblib.Parallel):
        def __init__(self, n_jobs=None, *a, **k):
            super().__init__(*a, n_jobs=min(cap, n_jobs or cap), **k)
    joblib.Parallel = _Capped
    for m in (_lp, LL_net, LINGER_tr):
        if hasattr(m, "Parallel"):
            m.Parallel = _Capped
    print(f"[linger {args.mode}] joblib pools capped at {cap} workers", flush=True)
    if args.seed:                       # LINGER calls set_global_seed() (default 42) at the start of training
        _orig_seed = LINGER_tr.set_global_seed
        LINGER_tr.set_global_seed = lambda seed=None: _orig_seed(42 + args.seed)
        print(f"[linger {args.mode}] seed {args.seed}: LINGER training seed {42 + args.seed}", flush=True)
    t0 = time.time()
    out = args.outdir
    os.makedirs(out, exist_ok=True)
    a_rna, a_atac = build_adata(args.rna, args.atac, args.mode, seed=args.seed, max_cells=args.max_cells)
    sc.pp.filter_cells(a_rna, min_genes=200)
    sc.pp.filter_genes(a_rna, min_cells=3)
    sc.pp.filter_cells(a_atac, min_genes=200)
    sc.pp.filter_genes(a_atac, min_cells=3)
    sel = sorted(set(a_rna.obs["barcode"]) & set(a_atac.obs["barcode"]))
    a_rna = a_rna[pd.Index(a_rna.obs["barcode"]).get_indexer(sel)].copy()
    a_atac = a_atac[pd.Index(a_atac.obs["barcode"]).get_indexer(sel)].copy()
    print(f"[linger {args.mode}] cells {a_rna.shape[0]}, genes {a_rna.shape[1]}, peaks {a_atac.shape[1]}", flush=True)
    # LINGER draws floor(sqrt(N)) + 1 metacells per cell-type label, and its training needs >= 50 of them (SHAP
    # background of 50 drawn without replacement). BEAR-GRN datasets are one population (474 cells -> 22 metacells),
    # so cells are split into expression clusters as pseudo-labels: k-means on RNA PCA (label-free, fixed seed),
    # the smallest k from 2 giving >= 60 metacells. The same labels are used in every ATAC condition.
    n_one = int(np.floor(np.sqrt(a_rna.shape[0]))) + 1
    if args.labels == "none" or (args.labels == "auto" and n_one >= 50):
        k, n_meta = 1, n_one
        a_rna.obs["label"] = a_atac.obs["label"] = "all"
    else:
        from sklearn.cluster import KMeans
        tmp = a_rna.copy()
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
        sc.pp.pca(tmp, n_comps=15, random_state=args.seed)
        for k in range(2, 21):
            lab = KMeans(k, n_init=10, random_state=args.seed).fit_predict(tmp.obsm["X_pca"])
            n_meta = sum(int(np.floor(np.sqrt(c))) + 1 for c in np.bincount(lab) if c >= 10)
            if n_meta >= 60:
                break
        a_rna.obs["label"] = a_atac.obs["label"] = [f"c{x}" for x in lab]
    print(f"[linger {args.mode}] pseudo-labels: k={k} -> {n_meta} metacells", flush=True)
    TG, RE = pseudo_bulk(a_rna.copy(), a_atac.copy(), False)
    RE[RE > 100] = 100
    TG, RE = TG.fillna(0), RE.fillna(0)
    a_atac.write(out + "adata_ATAC.h5ad")
    a_rna.write(out + "adata_RNA.h5ad")
    TG.to_csv(out + "TG_pseudobulk.tsv")
    RE.to_csv(out + "RE_pseudobulk.tsv")
    os.chdir(out)
    os.makedirs("data", exist_ok=True)
    pd.DataFrame(a_atac.var["gene_ids"]).to_csv("data/Peaks.txt", header=None, index=None)
    pd.DataFrame(a_atac.var["gene_ids"]).to_csv(out + "Peaks.txt", header=None, index=None)
    preprocess(TG, RE, args.grndir, args.genome, "LINGER", out)
    print(f"[linger {args.mode}] preprocess done {time.time() - t0:.0f}s", flush=True)
    LINGER_tr.training(args.grndir, "LINGER", out, "ReLU", args.species)
    print(f"[linger {args.mode}] training done {time.time() - t0:.0f}s", flush=True)
    import anndata
    a_rna, a_atac = anndata.read_h5ad(out + "adata_RNA.h5ad"), anndata.read_h5ad(out + "adata_ATAC.h5ad")
    LL_net.TF_RE_binding(args.grndir, a_rna, a_atac, args.genome, "LINGER", out)
    LL_net.cis_reg(args.grndir, a_rna, a_atac, args.genome, "LINGER", out)
    LL_net.trans_reg(args.grndir, "LINGER", out, args.genome)
    S = pd.read_csv(out + "cell_population_trans_regulatory.txt", sep="\t", index_col=0)     # TG x TF
    long = S.stack().reset_index()
    long.columns = ["Target", "Source", "Score"]
    ed = Path(args.export_dir)
    ed.mkdir(parents=True, exist_ok=True)
    long[["Source", "Target", "Score"]].to_csv(ed / "grn.tsv.gz", sep="\t", index=False, float_format="%.5g")
    json.dump({"method": "LINGER re-run", "mode": args.mode, "seed": args.seed, "labels": args.labels, "max_cells": args.max_cells,
               "k_pseudo_labels": int(k), "metacells": int(n_meta), "cells": int(a_rna.shape[0]),
               "tg": int(S.shape[0]), "tf": int(S.shape[1]), "seconds": time.time() - t0},
              open(ed / "meta.json", "w"), indent=1)
    print(f"[linger {args.mode}] exported {len(long):,} edges, {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
