#!/usr/bin/env python
"""Build scMultiomeGRN's SECONDARY benchmark input (Buenrostro et al. 2018 human adult
hematopoiesis, hg19; Supp. Fig. S1: HSC / CMP / GMP) in the same per-cell-type layout the
lung pipeline uses:  <out>/<HSC|CMP|GMP>/{atac,scrna}/{matrix.mtx,barcodes.tsv,peaks.tsv|genes.tsv}.
After this, scripts/16 / 17 / 18 / 19 run unchanged with configs/scmgrn/hema.yaml.

  rna   scRNA from Buenrostro 2018 Data S2 (mmc4; despite its name a tar.gz holding
        supplementary_code/input/scRNAseq_data.mat: counts 12558 genes x 14432 cells, labelled
        cd34/cmp0/cmp1/cmp2/gmp/hsc/mono). HSC = 'hsc'; GMP = 'gmp'; CMP = cmp0 + cmp1 + cmp2
        (ASSUMPTION: the paper does not say which scRNA labels it used; scATAC and scRNA are
        separate experiments, as in the paper).
  atac  scATAC from the pinellolab scATAC-benchmarking inputs: metadata.tsv (cell -> label),
        combined.sorted.merged.bed (peaks) and the deduplicated single-cell BAMs (Dropbox
        tarball, README of the repo), counted per peak with `bedtools coverage -counts`.
        HSC 347 / CMP 502 / GMP 402 cells, matching the paper. Needs bedtools; run it as a batch
        job (slurm/scmgrn_hema_prep.sh), not on a login node.

Usage:
  python scripts/14b_scmgrn_hema_prepare.py rna  --mmc4 mmc4.tar.gz --out sep_data/hema
  python scripts/14b_scmgrn_hema_prepare.py atac --bam_dir sc-bams_nodup --metadata metadata.tsv \
         --peaks combined.sorted.merged.bed --out sep_data/hema --n_jobs 10
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tarfile
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

MEMBER = "supplementary_code/input/scRNAseq_data.mat"
RNA_LABELS = {"HSC": ["hsc"], "CMP": ["cmp0", "cmp1", "cmp2"], "GMP": ["gmp"]}
ATAC_LABELS = {"HSC": "HSC", "CMP": "CMP", "GMP": "GMP"}


def write_mtx(d: Path, M: sp.spmatrix, barcodes, feats, feat_file: str, var_features=None):
    d.mkdir(parents=True, exist_ok=True)
    sio.mmwrite(str(d / "matrix.mtx"), sp.coo_matrix(M), field="integer" if M.dtype.kind in "iu" else "real")
    (d / "barcodes.tsv").write_text("".join(f"{b}\n" for b in barcodes))
    (d / feat_file).write_text("".join(f"{f}\n" for f in feats))
    if var_features:
        (d / "var_features.tsv").write_text("".join(f"{f}\n" for f in var_features))


def cmd_rna(a):
    with tempfile.TemporaryDirectory() as tmp:
        with tarfile.open(a.mmc4) as t:
            m = t.getmember(MEMBER)
            t.extract(m, tmp, filter="data")
        mat = sio.loadmat(str(Path(tmp) / MEMBER))
    counts = sp.csc_matrix(mat["counts"])                                  # genes x cells
    genes = [str(np.ravel(x)[0]) for x in mat["geneLabels"].ravel()]
    labels = np.array([str(np.ravel(x)[0]) for x in mat["cellLabels"].ravel()])
    keep_g = ~pd.Series(genes).duplicated().to_numpy()
    for ct, labs in RNA_LABELS.items():
        idx = np.where(np.isin(labels, labs))[0]
        M = counts[:, idx][np.where(keep_g)[0], :]
        names = [f"{ct}_rna_{i}" for i in range(len(idx))]
        write_mtx(Path(a.out) / ct / "scrna", M, names, [g for g, k in zip(genes, keep_g) if k], "genes.tsv")
        print(f"[rna] {ct}: {len(idx)} cells x {int(keep_g.sum())} genes (labels {labs})", flush=True)


def _count_one(args):
    bam, bed = args
    r = subprocess.run(["bedtools", "coverage", "-counts", "-a", bed, "-b", bam],
                       capture_output=True, text=True, check=True)
    return np.fromiter((int(l.rsplit("\t", 1)[1]) for l in r.stdout.splitlines()), dtype=np.int32)


def cmd_atac(a):
    meta = pd.read_csv(a.metadata, sep="\t")
    cell_col, lab_col = meta.columns[0], meta.columns[1]
    bed = pd.read_csv(a.peaks, sep="\t", header=None, usecols=[0, 1, 2])
    peaks = [f"{c}:{s}-{e}" for c, s, e in bed.itertuples(index=False)]
    for ct, lab in ATAC_LABELS.items():
        cells = meta.loc[meta[lab_col] == lab, cell_col].astype(str).tolist()
        bams = [str(Path(a.bam_dir) / f"{c}.bam") for c in cells]
        missing = [b for b in bams if not Path(b).exists()]
        if missing:
            sys.exit(f"{ct}: {len(missing)} BAMs missing, e.g. {missing[:2]}")
        with ThreadPoolExecutor(max_workers=a.n_jobs) as ex:
            cols = list(ex.map(_count_one, [(b, a.peaks) for b in bams]))
        M = sp.csc_matrix(np.stack(cols, axis=1))                          # peaks x cells
        write_mtx(Path(a.out) / ct / "atac", M, cells, peaks, "peaks.tsv")
        print(f"[atac] {ct}: {len(cells)} cells x {len(peaks)} peaks, nnz={M.nnz}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("rna"); r.add_argument("--mmc4", required=True); r.add_argument("--out", required=True)
    t = sub.add_parser("atac")
    t.add_argument("--bam_dir", required=True); t.add_argument("--metadata", required=True)
    t.add_argument("--peaks", required=True); t.add_argument("--out", required=True)
    t.add_argument("--n_jobs", type=int, default=8)
    a = ap.parse_args()
    {"rna": cmd_rna, "atac": cmd_atac}[a.cmd](a)


if __name__ == "__main__":
    main()
