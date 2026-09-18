"""Real TF motif scanning (plan.md Section 2, "still open" item).

Builds a genuinely motif-informed TF -> target-gene relation, replacing the
accessibility-only heuristic in `graph_builder.build_prior_graph` (which
only asks "is this locus open," never "does this specific TF's DNA binding
pattern actually appear here"). This closes a fidelity gap where MEvD-GRN's
own model never did what the scMultiomeGRN baseline wrapper already does
(FIMO-style motif scanning against the paper's spec).

Pipeline:
  1. A genome FASTA (hg38) + a JASPAR CORE PFM flat file are the only two
     external resources needed -- both standard, small enough to keep in
     `data/raw/genome/` and `data/raw/motifs/` respectively.
  2. Each TF with both (a) a JASPAR motif and (b) a reference-network entry
     gets its PFM converted to a log-odds PWM (`pfm_to_pwm`).
  3. For each candidate (TF, target gene) pair, the target's top-K
     accessible peaks (already computed by `preprocess_scatac`, see
     `gene_top_peak_coords`) are fetched as sequence and scanned for the
     TF's PWM (`best_pwm_score`). A hit above a fixed fraction of the PWM's
     max possible score creates an edge.

No external scanning binary (FIMO/MOODS) is required -- this uses a plain
NumPy sliding-window log-odds scan, which is exact (not an approximation)
for a single best-hit-per-sequence query, at the cost of not producing
calibrated p-values the way FIMO does. That tradeoff is acceptable here
since the resulting edges are used for a fixed top-K message-passing graph,
not for a p-value-thresholded genome-wide report.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import torch

_BASES = "ACGT"
_BASE_IDX = {b: i for i, b in enumerate(_BASES)}
_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


# --------------------------------------------------------------------------- genome FASTA
class FastaIndex:
    """Minimal random-access FASTA reader (one file, indexed once on load).

    Avoids a pysam/samtools dependency: builds an in-memory
    {chrom: byte_offset} index by scanning the file once, then seeks +
    reads fixed-width-wrapped sequence lines directly. Fine for a single
    hg38-scale genome (~3GB) read a handful of times per gene.
    """

    def __init__(self, fasta_path: str):
        self.path = Path(fasta_path)
        self._offsets: Dict[str, int] = {}
        self._line_len: Dict[str, int] = {}
        self._seq_len: Dict[str, int] = {}
        self._fh = None  # persistent handle -- fetch() is called millions of times
        self._index()

    def _index(self) -> None:
        with open(self.path, "rb") as f:
            chrom = None
            while True:
                pos = f.tell()
                line = f.readline()
                if not line:
                    break
                if line.startswith(b">"):
                    chrom = line[1:].split()[0].decode()
                    self._offsets[chrom] = f.tell()
                    self._seq_len[chrom] = 0
                elif chrom is not None:
                    if chrom not in self._line_len:
                        self._line_len[chrom] = len(line.rstrip(b"\n"))
                    self._seq_len[chrom] += len(line.rstrip(b"\n"))
        print(f"[motif_scan] indexed {len(self._offsets)} sequences from {self.path}", flush=True)

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """0-based half-open [start, end), like BED/UCSC coordinates."""
        if chrom not in self._offsets or chrom not in self._line_len:
            return ""
        line_len = self._line_len[chrom]
        if line_len == 0:
            return ""
        start = max(0, start)
        end = min(end, self._seq_len[chrom])
        if end <= start:
            return ""
        start_line, start_col = divmod(start, line_len)
        n_lines_needed = (end - start) // line_len + 2
        byte_start = self._offsets[chrom] + start_line * (line_len + 1)
        if self._fh is None:
            self._fh = open(self.path, "rb")
        self._fh.seek(byte_start)
        raw = self._fh.read(n_lines_needed * (line_len + 1))
        seq = raw.replace(b"\n", b"").decode()
        return seq[start_col:start_col + (end - start)].upper()


# --------------------------------------------------------------------------- JASPAR PFMs
def load_jaspar_pfms(path: str) -> Dict[str, np.ndarray]:
    """Parse a JASPAR CORE flat file (`pfm_meme.txt`-style or the plain
    ">MOTIF_ID  TF_NAME" + 4-row-counts format used by the JASPAR bundle
    download) into {TF_SYMBOL: pfm (4, motif_len) float32 counts}.

    Symbols are canon-uppercased for direct matching against SC-MO-GRN-DB's
    gene/TF symbols. A TF name shared by multiple JASPAR matrices keeps the
    first one encountered.
    """
    from src.data.preprocessing import canon

    pfms: Dict[str, np.ndarray] = {}
    name = None
    rows: List[List[float]] = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None and len(rows) == 4:
                    pfms.setdefault(name, np.array(rows, dtype=np.float32))
                parts = line[1:].split()
                name = canon(parts[1]) if len(parts) > 1 else canon(parts[0])
                rows = []
            elif line.strip():
                # Rows may be "A [ 1 2 3 ... ]" (JASPAR) or plain space-separated counts.
                cleaned = line.replace("[", " ").replace("]", " ")
                toks = cleaned.split()
                toks = [t for t in toks if t not in _BASES]
                rows.append([float(t) for t in toks])
        if name is not None and len(rows) == 4:
            pfms.setdefault(name, np.array(rows, dtype=np.float32))
    print(f"[motif_scan] loaded {len(pfms)} PFMs from {path}", flush=True)
    return pfms


def pfm_to_pwm(pfm: np.ndarray, pseudocount: float = 0.8,
              background: Tuple[float, float, float, float] = (0.25, 0.25, 0.25, 0.25)
              ) -> np.ndarray:
    """(4, L) count matrix -> (4, L) log-odds PWM, log2(freq / background)."""
    counts = pfm + pseudocount
    freqs = counts / counts.sum(axis=0, keepdims=True)
    bg = np.asarray(background, dtype=np.float32).reshape(4, 1)
    return np.log2(freqs / bg).astype(np.float32)


def _one_hot(seq: str) -> np.ndarray:
    """(4, len(seq)) one-hot; unrecognized bases (N, etc.) get an all-zero column."""
    out = np.zeros((4, len(seq)), dtype=np.float32)
    for i, ch in enumerate(seq):
        j = _BASE_IDX.get(ch)
        if j is not None:
            out[j, i] = 1.0
    return out


def best_pwm_score(seq: str, pwm: np.ndarray) -> float:
    """Best log-odds score of `pwm` (and its reverse complement) anywhere
    in `seq`; -inf if `seq` is shorter than the motif."""
    L = pwm.shape[1]
    if len(seq) < L:
        return float("-inf")
    onehot = _one_hot(seq)
    best = float("-inf")
    for pwm_variant in (pwm, pwm[:, ::-1][::-1, :]):  # forward, reverse-complement
        for start in range(onehot.shape[1] - L + 1):
            score = float((onehot[:, start:start + L] * pwm_variant).sum())
            if score > best:
                best = score
    return best


def pwm_max_score(pwm: np.ndarray) -> float:
    return float(pwm.max(axis=0).sum())


# --------------------------------------------------------------------------- graph builder
def build_motif_graph(gene_index: Dict[str, int], tf_indices: List[int],
                      tf_names_by_index: Dict[int, str],
                      gene_top_peaks: Dict[int, List[Tuple[str, int, int]]],
                      fasta_path: str, jaspar_path: str,
                      score_frac: float = 0.75,
                      evidence: Optional[Dict[str, torch.Tensor]] = None,
                      exclude_positives: bool = True,
                      candidate_edges: Optional[torch.Tensor] = None) -> torch.Tensor:
    """Directed TF -> target-gene edges backed by an ACTUAL DNA motif hit,
    not an accessibility proxy. For each TF with a JASPAR motif, scans that
    TF's PWM against its candidate target genes' top-K accessible peak
    sequences (from `gene_top_peak_coords`); an edge is added if the best
    hit reaches `score_frac` of the PWM's max possible score.

    `candidate_edges`: restrict scanning to (tf, gene) pairs already in this
    LongTensor(2, E) -- normally `build_prior_graph`'s co-expression x
    accessibility candidates. Without this, the scan is a full TF x
    all-genes-with-a-peak Cartesian product (~119 real-motif TFs x ~20,000
    genes here), which is ~20x more (TF, gene) pairs than the intended
    candidate set and takes hours; passing the existing candidate edges
    (already a deliberately narrowed shortlist) is both the intended scope
    -- this graph is meant to REFINE the candidate set with a real sequence
    check, not scan the whole genome per TF -- and the only way this
    finishes in a reasonable time on CPU.

    Known positive edges are excluded (same label-leakage-safety contract
    as `build_prior_graph`). Returns LongTensor(2, E); empty if the genome
    FASTA / JASPAR file aren't available (graceful degradation, matching
    the rest of the ATAC pipeline's missing-resource handling).
    """
    if not fasta_path or not Path(fasta_path).exists():
        print(f"[motif_scan] WARNING: no genome FASTA at {fasta_path!r}; "
              "motif graph will be empty.", flush=True)
        return torch.zeros((2, 0), dtype=torch.long)
    if not jaspar_path or not Path(jaspar_path).exists():
        print(f"[motif_scan] WARNING: no JASPAR motif file at {jaspar_path!r}; "
              "motif graph will be empty.", flush=True)
        return torch.zeros((2, 0), dtype=torch.long)

    fasta = FastaIndex(fasta_path)
    pfms = load_jaspar_pfms(jaspar_path)

    pos_keys: Set[Tuple[int, int]] = set()
    if exclude_positives and evidence:
        for e in evidence.values():
            pos_keys |= set(zip(e[0].tolist(), e[1].tolist()))

    # Restrict to (tf, gene) pairs already shortlisted by build_prior_graph,
    # not the full TF x all-genes Cartesian product -- see docstring.
    genes_by_tf: Dict[int, List[int]]
    if candidate_edges is not None and candidate_edges.numel():
        genes_by_tf = {}
        for tf_i, gi_i in zip(candidate_edges[0].tolist(), candidate_edges[1].tolist()):
            genes_by_tf.setdefault(tf_i, []).append(gi_i)
    else:
        print("[motif_scan] WARNING: no candidate_edges given -- falling back to a full "
              "TF x all-genes scan, which is much slower.", flush=True)
        genes_by_tf = {tf: list(gene_top_peaks) for tf in tf_indices}

    seq_cache: Dict[int, List[str]] = {}  # gene -> fetched peak sequences, shared across TFs

    def peak_seqs(gi: int) -> List[str]:
        if gi not in seq_cache:
            seq_cache[gi] = [fasta.fetch(c, s, e) for c, s, e in gene_top_peaks.get(gi, [])]
        return seq_cache[gi]

    edges: Set[Tuple[int, int]] = set()
    n_tf_with_motif = 0
    n_pairs_scanned = 0
    for tf in tf_indices:
        tf_name = tf_names_by_index.get(tf)
        pfm = pfms.get(tf_name) if tf_name else None
        if pfm is None:
            continue
        n_tf_with_motif += 1
        pwm = pfm_to_pwm(pfm)
        thresh = score_frac * pwm_max_score(pwm)
        for gi in genes_by_tf.get(tf, []):
            if gi == tf or (tf, gi) in pos_keys:
                continue
            n_pairs_scanned += 1
            hit = False
            for seq in peak_seqs(gi):
                if seq and best_pwm_score(seq, pwm) >= thresh:
                    hit = True
                    break
            if hit:
                edges.add((tf, gi))
    print(f"[motif_scan] {n_tf_with_motif}/{len(tf_indices)} TFs had a JASPAR motif, "
          f"{n_pairs_scanned} candidate (TF, gene) pairs scanned -> "
          f"{len(edges)} motif-backed edges over {len(gene_index)} genes "
          f"(positives excluded={exclude_positives})", flush=True)
    if not edges:
        return torch.zeros((2, 0), dtype=torch.long)
    src, dst = zip(*sorted(edges))
    return torch.tensor([list(src), list(dst)], dtype=torch.long)
