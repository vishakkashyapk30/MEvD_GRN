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


def build_adata(rna_csv: str, atac_csv: str, mode: str, seed: int = 0):
    import anndata
    R = pd.read_csv(rna_csv, index_col=0)
    A = pd.read_csv(atac_csv, index_col=0)
    cells = [c for c in R.columns if c in set(A.columns)]
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
    args = ap.parse_args()
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
    t0 = time.time()
    out = args.outdir
    os.makedirs(out, exist_ok=True)
    a_rna, a_atac = build_adata(args.rna, args.atac, args.mode)
    sc.pp.filter_cells(a_rna, min_genes=200)
    sc.pp.filter_genes(a_rna, min_cells=3)
    sc.pp.filter_cells(a_atac, min_genes=200)
    sc.pp.filter_genes(a_atac, min_cells=3)
    sel = sorted(set(a_rna.obs["barcode"]) & set(a_atac.obs["barcode"]))
    a_rna = a_rna[pd.Index(a_rna.obs["barcode"]).get_indexer(sel)].copy()
    a_atac = a_atac[pd.Index(a_atac.obs["barcode"]).get_indexer(sel)].copy()
    print(f"[linger {args.mode}] cells {a_rna.shape[0]}, genes {a_rna.shape[1]}, peaks {a_atac.shape[1]}", flush=True)
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
    json.dump({"method": "LINGER re-run", "mode": args.mode, "cells": int(a_rna.shape[0]),
               "tg": int(S.shape[0]), "tf": int(S.shape[1]), "seconds": time.time() - t0},
              open(ed / "meta.json", "w"), indent=1)
    print(f"[linger {args.mode}] exported {len(long):,} edges, {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
