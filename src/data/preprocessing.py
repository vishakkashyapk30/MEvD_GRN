"""scRNA / scATAC preprocessing for SC-MO-GRN-DB (plan Part 3).

Reality check vs. the original plan (which assumed .h5ad, cells x genes):
  * SC-MO-GRN-DB single-cell files are TAB-DELIMITED TEXT, oriented
    FEATURES x CELLS:
        scRNA.Expression_Matrix.*.txt   row0="Symbol\\t<cell> ..."   floats
        scATAC.Peak_Matrix.*.txt        row0="peak_coord\\t<cell> ..." ints, peaks="chr1:start-end"
  * RNA gene symbols are mixed-case (Esrrb); reference-network symbols are
    UPPERCASE (ESRRB). We canonicalise everything to UPPERCASE for matching.
  * RNA values are already normalised floats; we re-apply CP10K + log1p so the
    transform is well-defined regardless of upstream scaling.

The model operates on GENE NODES, so every modality is reduced to a per-gene
(mean, variance) feature pair across cells.
"""
from __future__ import annotations

import glob
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp


# ----------------------------------------------------------------------------- helpers
def canon(name: str) -> str:
    """Canonical gene symbol: strip + uppercase (matches network files)."""
    return str(name).strip().upper()


def extract_symbol(label: str) -> str:
    """Pull a gene symbol out of heterogeneous SC-MO-GRN-DB row labels.

    Seen formats: "Esrrb" (plain) and
    "chr3:108107280-108146146|-|ENSMUSG00000000001.4|Gnai3" (coord|strand|ensembl|symbol).
    For the pipe form the symbol is the last field.
    """
    s = str(label).strip()
    if "|" in s:
        s = s.split("|")[-1]
    return canon(s)


def _resolve_glob(pattern: str) -> str:
    matches = sorted(glob.glob(pattern, recursive=True))
    if not matches:
        raise FileNotFoundError(f"No file matches pattern: {pattern}")
    return matches[0]


def _col_mean_var(X) -> Tuple[np.ndarray, np.ndarray]:
    """Per-column (per-feature) mean and variance, sparse- or dense-safe."""
    n = X.shape[0]
    if sp.issparse(X):
        mean = np.asarray(X.mean(axis=0)).ravel()
        mean_sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
        var = np.maximum(mean_sq - mean ** 2, 0.0)
    else:
        X = np.asarray(X)
        mean = X.mean(axis=0)
        var = X.var(axis=0)
    return mean.astype(np.float32), var.astype(np.float32)


def _col_detection_rate(X) -> np.ndarray:
    """Per-gene fraction of cells with non-zero signal (regulatory-activity proxy)."""
    n = max(X.shape[0], 1)
    nnz = np.asarray((X > 0).sum(axis=0)).ravel()
    return (nnz / n).astype(np.float32)


def coexpression_signatures(X, d: int) -> np.ndarray:
    """Per-gene co-expression signature s_g in R^d such that s_i . s_j approximates
    the Pearson correlation between genes i and j across cells.

    Built from a truncated SVD of the column-standardized, mean-centered expression
    matrix W (W = (X - col_mean) / (col_std * sqrt(C))), for which W^T W is exactly
    the gene-gene correlation matrix. With W = U S V^T, we set s_g = (V S)_g so that
    s_i . s_j = (V S^2 V^T)_ij approximates (W^T W)_ij = corr(i, j).

    X: cells x genes (sparse or dense), log-normalized. Centering is applied
    implicitly via a LinearOperator so the sparse matrix is never densified.
    Returns (n_genes, d') float32; d' <= d (clamped by matrix rank).
    """
    from scipy.sparse.linalg import LinearOperator, svds

    C, G = X.shape
    if sp.issparse(X):
        Xc = X.tocsr().astype(np.float64)
        mean = np.asarray(Xc.mean(axis=0)).ravel()
        msq = np.asarray(Xc.multiply(Xc).mean(axis=0)).ravel()
        var = np.maximum(msq - mean ** 2, 0.0)
    else:
        Xc = np.asarray(X, dtype=np.float64)
        mean = Xc.mean(axis=0)
        var = Xc.var(axis=0)
    std = np.sqrt(var)
    std[std == 0] = 1.0
    scale = 1.0 / (std * np.sqrt(max(C, 1)))            # per-gene column scaling

    k = int(min(d, min(C, G) - 1))
    if k < 1:                                            # degenerate (tiny matrix)
        return np.zeros((G, 1), dtype=np.float32)

    def matvec(x):                                       # W @ x, x in R^G -> R^C
        y = np.asarray(x).ravel() * scale
        return Xc.dot(y) - float(mean.dot(y))

    def rmatvec(u):                                      # W^T @ u, u in R^C -> R^G
        u = np.asarray(u).ravel()
        return (Xc.T.dot(u) - mean * float(u.sum())) * scale

    W = LinearOperator((C, G), matvec=matvec, rmatvec=rmatvec, dtype=np.float64)
    try:
        _, s, vt = svds(W, k=k)
    except Exception as e:                               # pragma: no cover
        print(f"[scRNA] WARNING: co-expression SVD failed ({e}); signatures = 0.", flush=True)
        return np.zeros((G, k), dtype=np.float32)
    sig = (vt.T * s).astype(np.float32)                  # (G, k)
    print(f"[scRNA] co-expression signatures: {sig.shape[0]} genes x {sig.shape[1]} dims "
          f"(s_i.s_j ~ corr)", flush=True)
    return sig


def _sniff_delimiter(path: str) -> str:
    """Detect tab vs comma from the first line (SC-MO-GRN-DB uses both)."""
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        first = f.readline()
    return "," if first.count(",") > first.count("\t") else "\t"


def load_10x_h5(path: str, feature_type: str,
                keep_prefix: Optional[str] = None) -> Tuple[sp.csr_matrix, List[str]]:
    """Load one feature type from a CellRanger 10x .h5 (features×barcodes CSC).

    Returns (X cells×features CSR float32, labels). For Peaks, labels come from
    'interval'; for genes from 'name'. `keep_prefix` (e.g. 'hg38.') filters a
    barnyard matrix to one species and strips the prefix.
    """
    import h5py
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        n_feat, n_bc = int(m["shape"][0]), int(m["shape"][1])
        M = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]),
                          shape=(n_feat, n_bc))                    # features × barcodes
        ftype = m["features"]["feature_type"][:]
        label_key = "interval" if feature_type == "Peaks" else "name"
        labels_all = [x.decode() if isinstance(x, bytes) else str(x)
                      for x in m["features"][label_key][:]]
    mask = ftype == feature_type.encode()
    rows = np.where(mask)[0]
    labels = [labels_all[i] for i in rows]
    if keep_prefix:
        keep = [k for k, lab in enumerate(labels) if lab.startswith(keep_prefix)]
        rows = rows[keep]
        labels = [labels[k][len(keep_prefix):] for k in keep]
    X = M[rows, :].T.tocsr().astype(np.float32)                   # cells × features
    return X, labels


def load_features_by_cells(path: str) -> Tuple[sp.csr_matrix, List[str]]:
    """Load a FEATURES x CELLS matrix (tab- or comma-delimited); return
    (X_cells_by_feat CSR float32, feature_names). Row labels = features."""
    sep = _sniff_delimiter(path)
    df = pd.read_csv(path, sep=sep, index_col=0, engine="c")
    feature_names = [str(f) for f in df.index.tolist()]
    X_feat_by_cell = df.to_numpy(dtype=np.float32)          # features x cells
    X = sp.csr_matrix(X_feat_by_cell.T)                      # cells x features
    return X, feature_names


# ----------------------------------------------------------------------------- scRNA
def _normalize_log(X: sp.csr_matrix, cfg: dict) -> sp.csr_matrix:
    """CP-target_sum library normalisation + log1p (both optional via cfg)."""
    X = X.astype(np.float32).tocsr()
    X.data = np.nan_to_num(X.data, nan=0.0, posinf=0.0, neginf=0.0)
    X.data[X.data < 0] = 0.0
    if cfg.get("normalize_total", True):
        target = float(cfg.get("target_sum", 1e4))
        lib = np.asarray(X.sum(axis=1)).ravel()
        lib[lib == 0] = 1.0
        scale = target / lib
        X = sp.diags(scale) @ X
    if cfg.get("log1p", True):
        X = X.tocsr()
        X.data = np.log1p(X.data)
    return X.tocsr()


def preprocess_scrna(rna_path: str, cfg: dict
                     ) -> Tuple[np.ndarray, List[str], np.ndarray, np.ndarray]:
    """Returns (gene_features [n_genes,3], gene_names UPPERCASE, hvg_mask [n_genes],
    signatures [n_genes, d]).

    gene_features columns = [mean, variance, detection_rate]; detection_rate is the
    fraction of cells expressing the gene (a regulatory-activity proxy). signatures
    are the co-expression embeddings (see `coexpression_signatures`).
    """
    X, raw_names = load_features_by_cells(rna_path)          # cells x genes
    gene_names = [extract_symbol(g) for g in raw_names]

    # QC: cells with >= min_genes detected, genes in >= min_cells.
    min_genes = int(cfg.get("min_genes_per_cell", 200))
    min_cells = int(cfg.get("min_cells_per_gene", 3))
    genes_per_cell = np.asarray((X > 0).sum(axis=1)).ravel()
    keep_cells = genes_per_cell >= min_genes
    if keep_cells.sum() == 0:                                # tiny/odd matrices: keep all
        keep_cells[:] = True
    X = X[keep_cells]
    cells_per_gene = np.asarray((X > 0).sum(axis=0)).ravel()
    keep_genes = cells_per_gene >= min_cells
    X = X[:, keep_genes]
    gene_names = [g for g, k in zip(gene_names, keep_genes) if k]

    X = _normalize_log(X, cfg)
    mean, var = _col_mean_var(X)
    detection = _col_detection_rate(X)
    gene_features = np.stack([mean, var, detection], axis=1).astype(np.float32)  # (n_genes, 3)

    # Co-expression signatures (s_i . s_j ~ corr(i, j)) — preserves the pairwise
    # co-expression signal that mean/var alone discard.
    sig_dim = int(cfg.get("signature_dim", 50))
    signatures = coexpression_signatures(X, sig_dim)          # (n_genes, d)

    # HVG mask (normalized dispersion); kept for optional restriction/ablation.
    hvg_mask = _hvg_mask(mean, var, int(cfg.get("n_hvg", 3000)))

    # Collapse duplicate symbols (rare) by keeping the higher-variance row.
    gene_features, gene_names, hvg_mask, signatures = _dedup_genes(
        gene_features, gene_names, hvg_mask, signatures)
    print(f"[scRNA] {X.shape[0]} cells x {len(gene_names)} genes after QC; "
          f"{int(hvg_mask.sum())} HVGs", flush=True)
    return gene_features, gene_names, hvg_mask, signatures


def _hvg_mask(mean: np.ndarray, var: np.ndarray, n_hvg: int) -> np.ndarray:
    """Rough normalized-dispersion HVG flag on log-normalized stats."""
    disp = np.divide(var, mean, out=np.zeros_like(var), where=mean > 0)
    n_hvg = min(n_hvg, len(disp))
    if n_hvg <= 0:
        return np.zeros(len(disp), dtype=bool)
    thresh_idx = np.argsort(disp)[::-1][:n_hvg]
    mask = np.zeros(len(disp), dtype=bool)
    mask[thresh_idx] = True
    return mask


def _dedup_genes(feats: np.ndarray, names: List[str], hvg: np.ndarray,
                 signatures: Optional[np.ndarray] = None):
    seen: Dict[str, int] = {}
    keep = []
    for i, g in enumerate(names):
        if g not in seen:
            seen[g] = i
            keep.append(i)
        else:
            j = seen[g]
            if feats[i, 1] > feats[j, 1]:          # keep higher-variance copy
                keep[keep.index(j)] = i
                seen[g] = i
    keep = sorted(keep)
    kept_feats = feats[keep]
    kept_names = [names[i] for i in keep]
    kept_hvg = hvg[keep]
    if signatures is None:
        return kept_feats, kept_names, kept_hvg
    return kept_feats, kept_names, kept_hvg, signatures[keep]


# ----------------------------------------------------------------------------- scATAC
_PEAK_RE = re.compile(r"^(chr[\w]+)[:\-_](\d+)[\-_](\d+)$")


def _parse_peaks(peak_names: List[str]) -> pd.DataFrame:
    chrom, start, end, idx = [], [], [], []
    for i, p in enumerate(peak_names):
        m = _PEAK_RE.match(str(p).strip())
        if not m:
            continue
        chrom.append(m.group(1))
        start.append(int(m.group(2)))
        end.append(int(m.group(3)))
        idx.append(i)
    df = pd.DataFrame({"chrom": chrom, "start": start, "end": end, "peak_idx": idx})
    df["mid"] = (df["start"] + df["end"]) // 2
    return df


def _load_tss(genome: str, gene_names: List[str], gtf_path: Optional[str]) -> pd.DataFrame:
    """TSS per gene: DataFrame[gene(UPPER), chrom, tss]. Tries GTF, then BioMart."""
    want = set(gene_names)
    if gtf_path and Path(gtf_path).exists():
        return _tss_from_gtf(gtf_path, want)
    try:                                                     # online fallback (login node only)
        return _tss_from_biomart(genome, want)
    except Exception as e:                                   # pragma: no cover
        print(f"[scATAC] WARNING: no TSS source ({e}); ATAC features will be zero.", flush=True)
        return pd.DataFrame(columns=["gene", "chrom", "tss"])


def _tss_from_gtf(gtf_path: str, want: set) -> pd.DataFrame:
    """Vectorized GENCODE-GTF parse -> DataFrame[gene(UPPER), chrom, tss].
    Reads .gtf or .gtf.gz. TSS = start(+ strand) / end(- strand)."""
    cols = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attr"]
    parts = []
    for chunk in pd.read_csv(gtf_path, sep="\t", comment="#", header=None, names=cols,
                             usecols=["seqname", "feature", "start", "end", "strand", "attr"],
                             dtype={"seqname": str, "feature": str, "start": np.int64,
                                    "end": np.int64, "strand": str, "attr": str},
                             chunksize=500_000):
        g = chunk[chunk["feature"] == "gene"].copy()
        if g.empty:
            continue
        g["gene"] = g["attr"].str.extract(r'gene_name "([^"]+)"', expand=False).map(
            lambda x: canon(x) if isinstance(x, str) else x)
        g = g[g["gene"].isin(want)]
        if g.empty:
            continue
        g["chrom"] = g["seqname"].map(lambda c: c if c.startswith("chr") else "chr" + c)
        g["tss"] = np.where(g["strand"] == "+", g["start"], g["end"]).astype(np.int64)
        parts.append(g[["gene", "chrom", "tss"]])
    if not parts:
        return pd.DataFrame(columns=["gene", "chrom", "tss"])
    out = pd.concat(parts, ignore_index=True).drop_duplicates("gene")
    print(f"[scATAC] TSS from GTF: {len(out)} genes matched", flush=True)
    return out


def _tss_from_biomart(genome: str, want: set) -> pd.DataFrame:
    from pybiomart import Server
    server = Server(host="http://www.ensembl.org")
    ds_name = "mmusculus_gene_ensembl" if genome.startswith("mm") else "hsapiens_gene_ensembl"
    ds = server.marts["ENSEMBL_MART_ENSEMBL"].datasets[ds_name]
    res = ds.query(attributes=["external_gene_name", "chromosome_name",
                               "transcription_start_site"])
    res.columns = ["gene", "chrom", "tss"]
    res["gene"] = res["gene"].map(canon)
    res = res[res["gene"].isin(want)].dropna()
    res["chrom"] = res["chrom"].astype(str).map(lambda c: c if c.startswith("chr") else "chr" + c)
    res["tss"] = res["tss"].astype(int)
    print(f"[scATAC] TSS from BioMart: {len(res)} genes matched", flush=True)
    return res.drop_duplicates("gene")


def _peak_gene_incidence(peaks: pd.DataFrame, tss: pd.DataFrame,
                         gene_index: Dict[str, int], n_peaks: int,
                         window_bp: int) -> sp.csr_matrix:
    """Sparse (n_peaks x n_genes) 0/1 incidence: peak within +/-window of gene TSS."""
    n_genes = len(gene_index)
    rows, cols = [], []
    for chrom, gss in tss.groupby("chrom"):
        pk = peaks[peaks["chrom"] == chrom]
        if pk.empty:
            continue
        mids = pk["mid"].to_numpy()
        order = np.argsort(mids)
        mids_sorted = mids[order]
        pk_idx_sorted = pk["peak_idx"].to_numpy()[order]
        for _, gr in gss.iterrows():
            g = gr["gene"]
            gi = gene_index.get(g)
            if gi is None:
                continue
            lo = np.searchsorted(mids_sorted, gr["tss"] - window_bp, side="left")
            hi = np.searchsorted(mids_sorted, gr["tss"] + window_bp, side="right")
            for pi in pk_idx_sorted[lo:hi]:
                rows.append(int(pi))
                cols.append(gi)
    if not rows:
        return sp.csr_matrix((n_peaks, n_genes), dtype=np.float32)
    data = np.ones(len(rows), dtype=np.float32)
    return sp.csr_matrix((data, (rows, cols)), shape=(n_peaks, n_genes))


def preprocess_scatac(atac_path: str, gene_names: List[str], gene_index: Dict[str, int],
                      cfg: dict) -> Tuple[np.ndarray, np.ndarray]:
    """Gene-activity-score features + locus openness per gene.

    Returns (atac_features [n_genes, 2] = [mean, var], openness [n_genes]).
    `openness` is a normalized count of accessible peaks within +/-window of the
    gene's TSS (a motif-free accessibility proxy for how regulatable the locus is):
    openness_g = log1p(n_peaks_near_TSS(g)) / max_g. Genes with no mapped peaks get
    (0, 0) activity and 0 openness.
    """
    n_genes = len(gene_names)
    atac = cfg.get("atac", {}) if isinstance(cfg.get("atac"), dict) else {}
    genome = cfg.get("genome", "hg38")
    window_bp = int(cfg.get("atac_window_bp", 100_000))
    gtf_path = atac.get("gtf_path")
    zero = np.zeros((n_genes, 2), dtype=np.float32), np.zeros(n_genes, dtype=np.float32)

    if not atac_path:
        print("[scATAC] no ATAC path; returning zero features.", flush=True)
        return zero

    if str(atac_path).endswith(".h5"):                          # 10x multiome h5 (Peaks)
        prefix = "hg38." if not str(genome).startswith("mm") else "mm10."
        X, peak_names = load_10x_h5(atac_path, "Peaks", keep_prefix=prefix)
        print(f"[scATAC] 10x h5: {X.shape[0]} barcodes × {X.shape[1]} {prefix} peaks", flush=True)
    else:
        X, peak_names = load_features_by_cells(atac_path)      # cells x peaks (text)
    # QC on cells
    min_peaks = int(cfg.get("min_peaks_per_cell", 100))
    peaks_per_cell = np.asarray((X > 0).sum(axis=1)).ravel()
    keep = peaks_per_cell >= min_peaks
    if keep.sum() == 0:
        keep[:] = True
    X = X[keep]
    X = _normalize_log(X, cfg)                                 # CP10K + log1p

    peaks = _parse_peaks(peak_names)
    tss = _load_tss(genome, gene_names, gtf_path)
    if tss.empty or peaks.empty:
        return zero

    incidence = _peak_gene_incidence(peaks, tss, gene_index, len(peak_names), window_bp)
    # Gene activity (cells x genes) = X(cells x peaks) @ incidence(peaks x genes)
    gene_act = X @ incidence                                   # sparse cells x genes
    mean, var = _col_mean_var(gene_act)
    feats = np.stack([mean, var], axis=1).astype(np.float32)

    # Openness = normalized log peak-count near TSS (accessibility / regulatory potential).
    peak_count = np.asarray(incidence.sum(axis=0)).ravel().astype(np.float64)  # (n_genes,)
    openness = np.log1p(peak_count)
    mx = float(openness.max())
    openness = (openness / mx if mx > 0 else openness).astype(np.float32)

    nz = int((feats[:, 0] > 0).sum())
    print(f"[scATAC] {X.shape[0]} cells; {nz}/{n_genes} genes with non-zero activity "
          f"({100*nz/max(n_genes,1):.1f}%); mean peaks/gene={peak_count.mean():.2f}", flush=True)
    return feats, openness
