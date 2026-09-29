"""scMultiomeGRN's node features, re-implemented for the reproduction anchor.

* ATAC: MAESTRO `scatac-genescore --model Enhanced --genedistance 1000
  --species GRCh38` (MAESTRO v1.2.1, the version in the official docker image
  `winterdongqing/maestro:1.2.1`). `RP_AddExonRemovePromoter` is re-implemented
  here in vectorized form with identical scoring rules:
    phase 1 (all transcripts): peak centre inside the gene area
      [txStart-2kb | txStart, txEnd | txEnd+2kb] ->
        in an exon         : 1 / (total exon length in kb)
        in TSS +/- 2 kb    : 2^(-|centre - TSS| / decay)
        otherwise (intron) : no score
      Peaks scored in phase 1 for ANY gene are "in body" and excluded below.
    phase 2 (remaining peaks): |centre - TSS| <= 15*decay -> 2^(-|d|/decay)
  Scores x peak-count matrix, then per symbol the transcript with the largest
  total score is kept. The official pipeline always passes species GRCh38,
  including for the hg19 lung data; we keep that for fidelity.
* RNA: GRNBoost2 TF->gene importances (arboreto, seed 0) over variable
  genes U TFs, computed per target with arboreto's own
  `infer_partial_network` (the exact routine `grnboost2()` maps over), in a
  process pool instead of a dask cluster (dask-version independent).
"""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp

from src.benchmarks.scmgrn_groundtruth import clean_peak_name


# ----------------------------------------------------------------------------- MAESTRO RP
def load_refgenes(path: str) -> pd.DataFrame:
    bed = pd.read_csv(path, sep="\t", header=0, index_col=False, dtype={"chrom": str})
    bed["chrom"] = [c if str(c).startswith("chr") else f"chr{c}" for c in bed["chrom"]]
    plus = bed["strand"] == "+"
    bed["tss"] = np.where(plus, bed["txStart"], bed["txEnd"]).astype(np.int64)
    bed["gstart"] = np.where(plus, bed["txStart"] - 2000, bed["txStart"]).astype(np.int64)
    bed["gend"] = np.where(~plus, bed["txEnd"] + 2000, bed["txEnd"]).astype(np.int64)
    ex_s = [np.array([int(v) for v in s.strip(",").split(",")]) for s in bed["exonStarts"]]
    ex_e = [np.array([int(v) for v in s.strip(",").split(",")]) for s in bed["exonEnds"]]
    bed["exs"], bed["exe"] = ex_s, ex_e
    bed["length"] = [float(np.sum((e - s) / 1000.0)) for s, e in zip(ex_s, ex_e)]
    bed["uid"] = [f"{n}@{s}@{e}" for n, s, e in zip(bed["name2"], bed["gstart"], bed["gend"])]
    bed = bed.drop_duplicates(subset="uid", keep="first").reset_index(drop=True)
    return bed


def maestro_enhanced_rp(X_cells_by_peaks: sp.spmatrix, peak_names: Sequence[str],
                        refgenes: pd.DataFrame, symbols: Sequence[str],
                        decay: float = 1000.0) -> pd.DataFrame:
    """RP score (symbols x cells) for the requested gene symbols."""
    parsed = [clean_peak_name(p) for p in peak_names]
    chrom = np.array([p[0] for p in parsed])
    centre = np.array([(p[1] + p[2]) / 2.0 for p in parsed])
    n_peaks = len(parsed)
    gene_distance = 15 * decay
    by_chr: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for c in np.unique(chrom):
        idx = np.where(chrom == c)[0]
        o = np.argsort(centre[idx], kind="mergesort")
        by_chr[c] = (centre[idx][o], idx[o])

    want = set(symbols)
    inbody = np.zeros(n_peaks, dtype=bool)
    rows, cols, vals = [], [], []            # phase-1 scores, wanted transcripts only
    for t, g in enumerate(refgenes.itertuples(index=False)):
        pc = by_chr.get(g.chrom)
        if pc is None:
            continue
        cs, ids = pc
        lo = np.searchsorted(cs, g.gstart, side="left")
        hi = np.searchsorted(cs, g.gend, side="right")
        if hi <= lo:
            continue
        c, pid = cs[lo:hi], ids[lo:hi]
        in_exon = np.zeros(len(c), dtype=bool)
        for s, e in zip(g.exs, g.exe):
            in_exon |= (c >= s) & (c <= e)
        in_prom = (~in_exon) & (c >= g.tss - 2000) & (c <= g.tss + 2000)
        hit = in_exon | in_prom
        inbody[pid[hit]] = True
        if g.name2 in want:
            sc = np.where(in_exon, 1.0 / g.length if g.length > 0 else 0.0,
                          2.0 ** (-np.abs(c - g.tss) / decay))
            rows += [t] * int(hit.sum())
            cols += pid[hit].tolist()
            vals += sc[hit].tolist()
    # phase 2: out-of-body peaks within 15*decay of the TSS
    for t, g in enumerate(refgenes.itertuples(index=False)):
        if g.name2 not in want or g.chrom not in by_chr:
            continue
        cs, ids = by_chr[g.chrom]
        lo = np.searchsorted(cs, g.tss - gene_distance, side="left")
        hi = np.searchsorted(cs, g.tss + gene_distance, side="right")
        c, pid = cs[lo:hi], ids[lo:hi]
        ok = ~inbody[pid]
        rows += [t] * int(ok.sum())
        cols += pid[ok].tolist()
        vals += (2.0 ** (-np.abs(c[ok] - g.tss) / decay)).tolist()
    W = sp.csr_matrix((vals, (rows, cols)), shape=(len(refgenes), n_peaks))
    W.sum_duplicates()                       # (never happens: phases are disjoint)
    counts = sp.csc_matrix(X_cells_by_peaks).T.tocsr().astype(np.float64)  # peaks x cells
    G = (W @ counts).tocsr()                  # transcripts x cells
    tot = np.asarray(G.sum(axis=1)).ravel()
    best: Dict[str, Tuple[float, int]] = {}
    for t, name in enumerate(refgenes["name2"]):
        if name in want and (name not in best or tot[t] > best[name][0]):
            best[name] = (tot[t], t)
    keep = sorted(best)
    M = G[[best[s][1] for s in keep], :].toarray() if keep else np.zeros((0, counts.shape[1]))
    return pd.DataFrame(M, index=keep)


# ----------------------------------------------------------------------------- GRNBoost2
_GB_STATE: dict = {}


def _gb_init(tf_matrix, tf_names, seed):
    _GB_STATE.update(tf_matrix=tf_matrix, tf_names=tf_names, seed=seed)


def _gb_one(args):
    from arboreto.core import EARLY_STOP_WINDOW_LENGTH, SGBM_KWARGS, infer_partial_network
    target, expr = args
    df = infer_partial_network("GBM", SGBM_KWARGS, _GB_STATE["tf_matrix"], _GB_STATE["tf_names"],
                               target, expr, include_meta=False,
                               early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
                               seed=_GB_STATE["seed"])
    return df


def grnboost2_features(expr_genes_by_cells: pd.DataFrame, tf_nodes: Sequence[str],
                       seed: int = 0, n_jobs: int = 8) -> pd.DataFrame:
    """extract_scrna_feature2(method='grnboost2'): TF x (TFs + other genes)
    importance matrix (zeros where GRNBoost2 reports no link)."""
    data = expr_genes_by_cells
    genes = sorted(set(tf_nodes) & set(data.index))
    other = sorted(set(data.index) - set(genes))
    all_gene = genes + other
    X = data.loc[all_gene].T                              # cells x genes
    tf_matrix = X[genes].to_numpy(dtype=np.float64)
    jobs = [(g, X[g].to_numpy(dtype=np.float64)) for g in all_gene]
    parts = []
    with ProcessPoolExecutor(max_workers=max(1, n_jobs), initializer=_gb_init,
                             initargs=(tf_matrix, genes, seed)) as ex:
        for df in ex.map(_gb_one, jobs, chunksize=8):
            parts.append(df)
    assoc = pd.concat(parts, ignore_index=True)
    net = pd.DataFrame(0.0, index=genes, columns=all_gene)
    for tf, tg, imp in zip(assoc["TF"], assoc["target"], assoc["importance"]):
        net.at[tf, tg] = imp
    return net
