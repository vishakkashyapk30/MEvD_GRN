"""Re-implementation of scMultiomeGRN's ground-truth construction.

Xu et al., NAR 2025 (gkaf138), code: Zenodo 10.5281/zenodo.14848389,
`src/utils.py` functions ATAC_filter_from_MTX, extract_chromosome
(GetSequence.pl), fimo, get_TFBS_region, get_TFBS_from_promoter,
get_contained_region + build_Graph, build_cross_graph. The official code
has no license, so this is a line-by-line *re-implementation* (not a copy)
that reproduces the same numerical semantics, noted inline:

  1. peaks detected in > `threshold` (0.1) of cells, `rate > threshold`
  2. sequences via Bio::DB::Fasta subseq(start=>end) = 1-based inclusive
     on the BED start/end as written; header ">chr-start-end-"; chrY skipped
  3. `fimo --oc <dir> --thresh 1e-4 --no-qvalue <motif> <fasta>` per motif
  4. keep hits with p <= 1e-6; TF = motif_id.split("_")[0] (merges redundant
     motifs); abs coords = peak_start + fimo_start/stop - 1
  5. promoters = rows of the GENCODE protein-coding promoter table whose gene
     is a motif TF (TSS +/- 2 kb)
  6. edge TF_a <-> TF_b (symmetric) if a TFBS of TF_b lies fully inside TF_a's
     promoter; nodes = promoter TFs with >= 1 contained TFBS; sorted
  (bedops -e -20% before the containment loop is implied by containment.)
"""
from __future__ import annotations

import os
import shutil
import subprocess
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp


# ----------------------------------------------------------------------------- peaks
def clean_peak_name(name: str) -> Tuple[str, int, int]:
    """Their clean_chromosome_name(...).split('_'): tab/':'/'-' -> '_'."""
    s = str(name).replace("\t", "_").replace(":", "_").replace("-", "_") \
        .replace("'", "").replace('"', "")
    chrom, start, end = s.split("_")[:3]
    return chrom, int(start), int(end)


def filter_peaks(X_cells_by_peaks: sp.spmatrix, peak_names: List[str],
                 threshold: float = 0.1) -> pd.DataFrame:
    """Peaks with detection rate strictly > threshold, sort-bed ordered."""
    X = sp.csc_matrix(X_cells_by_peaks)
    n_cells = X.shape[0]
    nnz = np.diff((X > 0).astype(np.int8).tocsc().indptr)
    rate = nnz / max(n_cells, 1)
    keep = np.where(rate > threshold)[0]
    rows = [clean_peak_name(peak_names[i]) + (float(rate[i]),) for i in keep]
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "rate"])
    # sort-bed: lexicographic chromosome, then numeric start, end
    return df.sort_values(["chrom", "start", "end"], kind="mergesort").reset_index(drop=True)


# ----------------------------------------------------------------------------- FASTA
class IndexedFasta:
    """Minimal samtools-faidx-compatible random access (no pyfaidx needed).
    Builds `<fasta>.fai` on first use if absent (fixed line width required,
    true for UCSC hg19/hg38 .fa)."""

    def __init__(self, path: str):
        self.path = str(path)
        fai = self.path + ".fai"
        if not os.path.exists(fai):
            self._build_fai(fai)
        self.index: Dict[str, Tuple[int, int, int, int]] = {}
        with open(fai) as f:
            for line in f:
                name, length, offset, lb, lw = line.rstrip("\n").split("\t")[:5]
                self.index[name] = (int(length), int(offset), int(lb), int(lw))
        self._fh = open(self.path, "rb")

    def _build_fai(self, fai: str) -> None:
        entries = []
        with open(self.path, "rb") as f:
            name, length, offset, lb, lw = None, 0, 0, 0, 0
            pos = 0
            for raw in f:
                if raw.startswith(b">"):
                    if name is not None:
                        entries.append((name, length, offset, lb, lw))
                    name = raw[1:].split()[0].decode()
                    length, lb, lw = 0, 0, 0
                    offset = pos + len(raw)
                else:
                    if lb == 0:
                        lb, lw = len(raw.rstrip(b"\r\n")), len(raw)
                    length += len(raw.rstrip(b"\r\n"))
                pos += len(raw)
            if name is not None:
                entries.append((name, length, offset, lb, lw))
        with open(fai, "w") as out:
            for e in entries:
                out.write("\t".join(map(str, e)) + "\n")

    def fetch(self, chrom: str, start0: int, end0: int) -> Optional[str]:
        """0-based half-open [start0, end0)."""
        if chrom not in self.index:
            return None
        length, offset, lb, lw = self.index[chrom]
        start0, end0 = max(0, start0), min(length, end0)
        if end0 <= start0:
            return ""
        b0 = offset + (start0 // lb) * lw + start0 % lb
        b1 = offset + (end0 // lb) * lw + end0 % lb
        self._fh.seek(b0)
        return self._fh.read(b1 - b0).replace(b"\n", b"").replace(b"\r", b"").decode().upper()


def write_peak_fasta(peaks: pd.DataFrame, genome_fasta: str, out_fasta: str) -> int:
    """GetSequence.pl semantics: subseq(start=>end), 1-based inclusive, i.e.
    0-based [start-1, end); header '>chr-start-end-'; chrY skipped."""
    fa = IndexedFasta(genome_fasta)
    n = 0
    with open(out_fasta, "w") as out:
        for chrom, start, end in zip(peaks["chrom"], peaks["start"], peaks["end"]):
            if chrom == "chrY":
                continue
            seq = fa.fetch(chrom, int(start) - 1, int(end))
            if seq is None:
                print(f"[fasta] WARNING {chrom} not in genome; skipped", flush=True)
                continue
            out.write(f">{chrom}-{start}-{end}-\n{seq}\n")
            n += 1
    return n


# ----------------------------------------------------------------------------- motifs / FIMO
def unpack_motifs(motif_zip_or_dir: str, dest: str) -> str:
    """Return a directory containing the HOCOMOCO v11 .meme files."""
    p = Path(motif_zip_or_dir)
    if p.is_dir():
        return str(p)
    d = Path(dest)
    if not any(d.rglob("*.meme")):
        d.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(p) as z:
            z.extractall(d)
    return str(next(d.rglob("*.meme")).parent)


def motif_files(motif_dir: str) -> List[str]:
    return sorted(os.path.join(motif_dir, f) for f in os.listdir(motif_dir) if f.endswith(".meme"))


def motif_tf_names(motif_dir: str) -> set:
    """get_TFBS_from_promoter: file.split('_')[0] for .meme files."""
    return {f.split("_")[0] for f in os.listdir(motif_dir) if f.endswith(".meme")}


def _fimo_out_name(motif_file: str) -> str:
    tmp = os.path.basename(motif_file).split(".")
    return f"{tmp[0].split('_HUMAN')[0]}.{'.'.join(tmp[2:-1])}"


def run_fimo(fasta: str, motif_dir: str, out_dir: str, fimo_bin: str = "fimo",
             thresh: float = 1e-4, n_jobs: int = 8,
             only_motifs: Optional[Iterable[str]] = None) -> str:
    """One FIMO call per motif file, same flags as the official pipeline.
    Resumable: motifs whose output already exists are skipped."""
    out = Path(out_dir)
    (out / "tmp").mkdir(parents=True, exist_ok=True)
    files = motif_files(motif_dir)
    if only_motifs is not None:
        keep = set(only_motifs)
        files = [f for f in files if os.path.basename(f).split("_")[0] in keep]

    def one(mf: str) -> Tuple[str, int]:
        name = _fimo_out_name(mf)
        dst = out / f"{name}.txt"
        if dst.exists():
            return name, 0
        tmp_dir = out / "tmp" / name
        cmd = [fimo_bin, "--oc", str(tmp_dir), "--thresh", str(thresh), "--no-qvalue",
               "--verbosity", "1", mf, fasta]
        r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"fimo failed on {mf}: {r.stderr[-500:]}")
        shutil.move(str(tmp_dir / "fimo.tsv"), str(dst))
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return name, 1

    with ThreadPoolExecutor(max_workers=max(1, n_jobs)) as ex:
        done = sum(n for _, n in ex.map(one, files))
    shutil.rmtree(out / "tmp", ignore_errors=True)
    print(f"[fimo] {len(files)} motifs ({done} newly scanned) -> {out}", flush=True)
    return str(out)


def parse_tfbs(fimo_dir: str, p_thresh: float = 1e-6) -> pd.DataFrame:
    """get_TFBS_region: rows with p <= p_thresh -> (chrom, abs_start, abs_end, tf)."""
    rows = []
    for fn in sorted(f for f in os.listdir(fimo_dir) if f.endswith(".txt")):
        with open(os.path.join(fimo_dir, fn)) as f:
            head = f.readline().strip()
            if not head:
                continue
            assert head.startswith("motif_id"), fn
            for line in f:
                line = line.strip()
                if not line:
                    break
                t = line.split("\t")
                if float(t[7]) > p_thresh:
                    continue
                tf = t[0].split("_")[0]
                chro, start, _end = t[2].split("-")[:3]
                rows.append((chro, int(start) + int(t[3]) - 1, int(start) + int(t[4]) - 1, tf))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "tf"])
    return df.sort_values(["chrom", "start", "end"], kind="mergesort").reset_index(drop=True)


# ----------------------------------------------------------------------------- graph
def load_tf_promoters(promoter_file: str, tf_names: set) -> pd.DataFrame:
    df = pd.read_csv(promoter_file, sep="\t", header=0, dtype={"gene_name": str})
    df.columns = ["chrom", "start", "end", "strand", "gene"][: df.shape[1]] + list(df.columns[5:])
    df = df[df["gene"].isin(tf_names)].copy()
    return df.sort_values(["chrom", "start", "end"], kind="mergesort").reset_index(drop=True)


def build_tf_graph(tfbs: pd.DataFrame, promoters: pd.DataFrame) -> pd.DataFrame:
    """build_Graph: symmetric TF x TF 0/1 adjacency over sorted node set."""
    pairs = set()
    tf_set = set()
    by_chr = {c: g for c, g in tfbs.groupby("chrom")}
    for chrom, pg in promoters.groupby("chrom"):
        q = by_chr.get(chrom)
        if q is None or q.empty:
            continue
        qs, qe, qt = q["start"].to_numpy(), q["end"].to_numpy(), q["tf"].to_numpy()
        order = np.argsort(qs, kind="mergesort")
        qs, qe, qt = qs[order], qe[order], qt[order]
        for start, end, tf in zip(pg["start"].to_numpy(), pg["end"].to_numpy(), pg["gene"].to_numpy()):
            lo = np.searchsorted(qs, start, side="left")
            hi = np.searchsorted(qs, end, side="right")
            hit = qt[lo:hi][qe[lo:hi] <= end]
            if hit.size:
                tf_set.add(tf)
                for other in set(hit.tolist()):
                    pairs.add((tf, other))
                    pairs.add((other, tf))
    nodes = sorted(tf_set)
    idx = {t: i for i, t in enumerate(nodes)}
    A = np.zeros((len(nodes), len(nodes)), dtype=np.int64)
    for a, b in pairs:
        if a in idx and b in idx:
            A[idx[a], idx[b]] = 1
    return pd.DataFrame(A, index=nodes, columns=nodes)


def write_their_adj(adj: pd.DataFrame, path: str) -> None:
    """Same TSV layout as build_Graph / build_cross_graph (header 'TF')."""
    out = adj.copy()
    out.index.name = "TF"
    out.to_csv(path, sep="\t")


def cross_nodes(adj: pd.DataFrame, *feature_indices: Iterable[str]) -> List[str]:
    """build_cross_graph: keep node order of the graph, intersect with features."""
    keep = set(adj.index)
    for fi in feature_indices:
        keep &= set(fi)
    return [n for n in adj.index if n in keep]


def write_promoter_tss_gtf(promoter_file: str, out_gtf: str) -> str:
    """Minimal GTF (gene lines) whose start=end=TSS (= promoter midpoint), so
    src.data.preprocessing._tss_from_gtf gives exactly the TSS used by the
    ground truth -- no separate GENCODE GTF download needed."""
    df = pd.read_csv(promoter_file, sep="\t", header=0)
    df.columns = ["chrom", "start", "end", "strand", "gene"][: df.shape[1]] + list(df.columns[5:])
    tss = ((df["start"] + df["end"]) // 2).astype(int)
    with open(out_gtf, "w") as f:
        for c, t, s, g in zip(df["chrom"], tss, df["strand"], df["gene"]):
            f.write(f"{c}\tpromoter_table\tgene\t{t}\t{t}\t.\t{s}\t.\tgene_name \"{g}\";\n")
    return out_gtf
