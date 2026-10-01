#!/usr/bin/env python
"""Re-run LINGER (LingerGRN==1.110, method='LINGER') on the benchmark's exact cells,
following docs/PBMC.md of github.com/Durenlab/LINGER step by step, then export
what the PBMC benchmark needs. Run inside the separate `linger` conda env
(python 3.10, LingerGRN==1.110, bedtools); it imports nothing from this repo.

Only deviation from the tutorial (runbook amendment A1b): cells are the
benchmark's cells.tsv (label file + >=200 genes RNA + >=200 peaks ATAC), i.e.
LingerGRN 1.110's later-added mt<5% filter in get_adata is NOT applied (it
would keep 514/9,543 cells; the published LINGER/KEGNI cell counts are
unfiltered). The adata objects are built exactly as get_adata does otherwise.

Exports (into --root):
  linger/tg_<ct>.txt        LINGER's cell-type-specific TG set (= candidate space `linger_tg`)
  results/<ct>/unsupervised/none/linger/eval_scores.npz + metrics.json(stub)
                            trans scores for the eval TFs over the benchmark gene universe
                            (genes outside LINGER's TG set get min-1, i.e. rank last)

Usage (in the linger env):
  python src/benchmarks/pbmc_linger_run.py --root $ROOT --grndir $ROOT/linger/data_bulk/ \
      --outdir $ROOT/linger/out/ --steps prep,train,celltype,export
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

CTS = {"classical_monocyte": "classical monocytes", "naive_cd4_t": "naive CD4 T cells",
       "naive_b": "naive B cells", "mdc": "myeloid DC"}
EVAL_TFS = ["CTCF", "ETS1", "FOXP3", "IRF1", "IRF4", "MYC", "REST", "RUNX1", "SPI1", "STAT1"]


def build_adata(root: Path):
    import anndata
    import h5py
    with h5py.File(root / "raw/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5") as f:
        m = f["matrix"]
        M = sp.csc_matrix((m["data"][:].astype(np.float32), m["indices"][:], m["indptr"][:]),
                          shape=tuple(m["shape"][:]))
        names = np.array([x.decode() for x in m["features"]["name"][:]])
        inter = np.array([x.decode() for x in m["features"]["interval"][:]])
        ft = np.array([x.decode() for x in m["features"]["feature_type"][:]])
        bc = np.array([x.decode() for x in m["barcodes"][:]])
    # LINGER's features[1] = name for genes; for peaks the 10x features.tsv name is chr:start-end
    feat_name = np.where(ft == "Peaks", inter, names)
    cells = pd.read_csv(root / "cells/cells.tsv", sep="\t")
    keep = pd.Index(bc).get_indexer(cells["barcode"])
    X = M.T.tocsr()[keep]
    adata = anndata.AnnData(X=sp.csc_matrix(X))
    adata.var["gene_ids"] = feat_name
    adata.var_names = feat_name
    adata.obs["barcode"] = cells["barcode"].values
    adata.obs_names = cells["barcode"].values
    adata.obs["sample"] = [b.split("-")[1] for b in cells["barcode"].values]
    adata.obs["label"] = cells["label"].values
    a_rna = adata[:, np.where(ft == "Gene Expression")[0]].copy()
    a_atac = adata[:, np.where(ft == "Peaks")[0]].copy()
    a_rna.var.index = a_rna.var["gene_ids"].values
    a_rna.var_names_make_unique()
    a_rna.var["gene_ids"] = a_rna.var.index
    return a_rna, a_atac


def step_prep(root, outdir, grndir):
    import scanpy as sc
    from LingerGRN.preprocess import preprocess
    from LingerGRN.pseudo_bulk import pseudo_bulk
    a_rna, a_atac = build_adata(root)
    sc.pp.filter_cells(a_rna, min_genes=200)
    sc.pp.filter_genes(a_rna, min_cells=3)
    sc.pp.filter_cells(a_atac, min_genes=200)
    sc.pp.filter_genes(a_atac, min_cells=3)
    sel = list(set(a_rna.obs["barcode"].values) & set(a_atac.obs["barcode"].values))
    a_rna = a_rna[pd.Index(a_rna.obs["barcode"]).get_indexer(sel)]
    a_atac = a_atac[pd.Index(a_atac.obs["barcode"]).get_indexer(sel)]
    print(f"[linger] cells {a_rna.shape[0]}, genes {a_rna.shape[1]}, peaks {a_atac.shape[1]}", flush=True)
    samplelist = list(set(a_atac.obs["sample"].values))
    TG, RE = pd.DataFrame([]), pd.DataFrame([])
    single = (a_rna.obs["sample"].unique().shape[0] * a_rna.obs["sample"].unique().shape[0] > 100)
    for s in samplelist:
        t, r = pseudo_bulk(a_rna[a_rna.obs["sample"] == s], a_atac[a_atac.obs["sample"] == s], single)
        TG = pd.concat([TG, t], axis=1)
        RE = pd.concat([RE, r], axis=1)
        RE[RE > 100] = 100
    os.makedirs(outdir, exist_ok=True)
    a_atac.write(outdir + "adata_ATAC.h5ad")
    a_rna.write(outdir + "adata_RNA.h5ad")
    TG, RE = TG.fillna(0), RE.fillna(0)
    pd.DataFrame(a_atac.var["gene_ids"]).to_csv(outdir + "Peaks.txt", header=None, index=None)
    TG.to_csv(outdir + "TG_pseudobulk.tsv")
    RE.to_csv(outdir + "RE_pseudobulk.tsv")
    cwd = os.getcwd()
    os.chdir(outdir)          # the tutorial runs preprocess from the dir holding data/Peaks.txt
    os.makedirs("data", exist_ok=True)
    pd.DataFrame(a_atac.var["gene_ids"]).to_csv("data/Peaks.txt", header=None, index=None)
    preprocess(TG, RE, grndir, "hg38", "LINGER", outdir)
    os.chdir(cwd)


def load_ad(outdir):
    import anndata
    return anndata.read_h5ad(outdir + "adata_RNA.h5ad"), anndata.read_h5ad(outdir + "adata_ATAC.h5ad")


def step_train(outdir, grndir):
    import LingerGRN.LINGER_tr as LINGER_tr
    os.chdir(outdir)
    LINGER_tr.training(grndir, "LINGER", outdir, "ReLU", "Human")


def step_celltype(outdir, grndir):
    import LingerGRN.LL_net as LL_net
    os.chdir(outdir)
    a_rna, a_atac = load_ad(outdir)
    LL_net.TF_RE_binding(grndir, a_rna, a_atac, "hg38", "LINGER", outdir)
    LL_net.cis_reg(grndir, a_rna, a_atac, "hg38", "LINGER", outdir)
    LL_net.trans_reg(grndir, "LINGER", outdir, "hg38")
    for ct in CTS.values():
        LL_net.cell_type_specific_TF_RE_binding(grndir, a_rna, a_atac, "hg38", ct, outdir, "LINGER")
        LL_net.cell_type_specific_cis_reg(grndir, a_rna, a_atac, "hg38", ct, outdir, "LINGER")
        LL_net.cell_type_specific_trans_reg(grndir, a_rna, ct, outdir)


def step_export(root: Path, outdir):
    genes = [l.strip() for l in open(root / "cells/genes.txt") if l.strip()]
    gi = {g: i for i, g in enumerate(genes)}
    (root / "linger").mkdir(parents=True, exist_ok=True)
    for ct, lb in CTS.items():
        S = pd.read_csv(outdir + f"cell_type_specific_trans_regulatory_{lb}.txt", sep="\t", index_col=0)
        S.index = [str(g).upper() for g in S.index]
        S.columns = [str(c).upper() for c in S.columns]
        (root / "linger" / f"tg_{ct}.txt").write_text("\n".join(S.index) + "\n")
        tfs = [t for t in EVAL_TFS if t in S.columns]
        out = np.zeros((len(tfs), len(genes)), dtype=np.float32)
        for k, t in enumerate(tfs):
            col = S[t]
            col = col[~col.index.duplicated()]
            out[k] = float(np.nanmin(col.values)) - 1.0
            idx = [gi[g] for g in col.index if g in gi]
            out[k, idx] = col[[g for g in col.index if g in gi]].values
        d = root / "results" / ct / "unsupervised" / "none" / "linger"
        d.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(d / "eval_scores.npz", scores=out, tfs=np.array(tfs))
        json.dump({"method": "linger", "cell_type": ct, "source": "unsupervised", "regime": "none",
                   "n_tg": int(S.shape[0]), "tg_not_in_universe": int(sum(g not in gi for g in S.index)),
                   "tfs": tfs, "note": "scored by scripts/29 (pbmc_eval)"},
                  open(d / "metrics.json", "w"), indent=2)
        print(f"[linger export] {ct}: {S.shape[0]} TGs x {S.shape[1]} TFs; eval TFs {tfs}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--grndir", required=True, help="unpacked data_bulk/ (trailing slash)")
    ap.add_argument("--outdir", required=True, help="LINGER output dir (trailing slash)")
    ap.add_argument("--steps", default="prep,train,celltype,export")
    a = ap.parse_args()
    root = Path(a.root).resolve()
    outdir = str(Path(a.outdir).resolve()) + "/"
    grndir = str(Path(a.grndir).resolve()) + "/"
    for s in a.steps.split(","):
        print(f"===== LINGER step {s} =====", flush=True)
        {"prep": lambda: step_prep(root, outdir, grndir), "train": lambda: step_train(outdir, grndir),
         "celltype": lambda: step_celltype(outdir, grndir), "export": lambda: step_export(root, outdir)}[s]()


if __name__ == "__main__":
    main()
