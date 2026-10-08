"""TF-motif x target-accessibility pair features for MeVD-GRN under BEAR-GRN
(pre-registered model M3, docs/experiments/bear_grn_benchmark.md s12.6).

For each (TF node, gene): the best JASPAR 2024 PWM log-odds score of the TF,
relative to the PWM's maximum, over the gene's top-5 RP-weighted accessible
peaks (from scripts/21, `gene_top_peaks.json`); the number of those peaks with
a hit >= 0.8 x max; and a has-motif indicator. Label-free: only sequences,
PWMs and the dataset's own peaks are used. ChIP / KO ground truths are not
motif-derived, and every BEAR method uses TF motifs, so this is not circular.

Pieces:
  * pack_genome(): stream a UCSC *.fa.gz once into a 2-bit-packed store
    (~775 MB for hg38; N runs kept separately), so no 3 GB FASTA is ever on disk.
  * Genome2bit.fetch(): codes 0-3 (ACGT), 4 = N.
  * load_jaspar(): all matrices (redundant set), keyed by upper-cased TF name;
    heterodimers "FOS::JUN" are assigned to both partners.
  * scan_best(): vectorised torch conv1d scan (both strands, CPU or GPU),
    best relative score per (peak, motif), position-validity masked.
  * pair_features(): -> float16 array (n_tf_nodes, n_genes, 3).
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import torch

_LUT = np.full(256, 4, dtype=np.uint8)
for _c, _v in zip(b"ACGTacgt", [0, 1, 2, 3, 0, 1, 2, 3]):
    _LUT[_c] = _v
MAIN_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]}


# --------------------------------------------------------------------------- genome store
def _pack_chrom(seq: bytes, out: Path, name: str) -> int:
    codes = _LUT[np.frombuffer(seq, dtype=np.uint8)]
    n = codes.size
    isn = codes == 4
    edges = np.flatnonzero(np.diff(np.concatenate([[0], isn.view(np.int8), [0]])))
    nruns = edges.reshape(-1, 2) if edges.size else np.zeros((0, 2), dtype=np.int64)
    c = np.where(isn, 0, codes).astype(np.uint8)
    pad = (-n) % 4
    if pad:
        c = np.concatenate([c, np.zeros(pad, dtype=np.uint8)])
    c = c.reshape(-1, 4)
    packed = (c[:, 0] | (c[:, 1] << 2) | (c[:, 2] << 4) | (c[:, 3] << 6)).astype(np.uint8)
    np.savez(out / f"{name}.npz", packed=packed, n=np.int64(n), nruns=nruns.astype(np.int64))
    return n


def pack_genome(fa_gz: str, out_dir: str, chroms=MAIN_CHROMS) -> Dict[str, int]:
    """Stream a gzipped FASTA once; write <out_dir>/<chrom>.npz (2-bit + N runs)."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    sizes, name, buf = {}, None, []
    with gzip.open(fa_gz, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if name in chroms and buf:
                    sizes[name] = _pack_chrom(b"".join(buf), out, name)
                    print(f"[genome] {name}: {sizes[name]:,} bp", flush=True)
                name = line[1:].split()[0].decode()
                buf = []
            elif name in chroms:
                buf.append(line.rstrip(b"\n\r"))
        if name in chroms and buf:
            sizes[name] = _pack_chrom(b"".join(buf), out, name)
    (out / "sizes.json").write_text(json.dumps(sizes))
    return sizes


class Genome2bit:
    def __init__(self, store_dir: str):
        self.dir = Path(store_dir)
        self._cache: Dict[str, Tuple[np.ndarray, int, np.ndarray]] = {}

    def _load(self, chrom: str):
        if chrom not in self._cache:
            p = self.dir / f"{chrom}.npz"
            if not p.exists():
                return None
            z = np.load(p)
            self._cache[chrom] = (z["packed"], int(z["n"]), z["nruns"])
        return self._cache[chrom]

    def fetch(self, chrom: str, start: int, end: int) -> np.ndarray:
        """0-based half-open [start, end) -> uint8 codes (0-3, 4 = N/outside)."""
        rec = self._load(chrom)
        L = max(0, end - start)
        if rec is None or L == 0:
            return np.full(L, 4, dtype=np.uint8)
        packed, n, nruns = rec
        s, e = max(0, start), min(n, end)
        out = np.full(L, 4, dtype=np.uint8)
        if e <= s:
            return out
        b0, b1 = s // 4, (e + 3) // 4
        by = packed[b0:b1]
        codes = np.stack([by & 3, (by >> 2) & 3, (by >> 4) & 3, (by >> 6) & 3], axis=1).reshape(-1)
        codes = codes[s - b0 * 4: s - b0 * 4 + (e - s)]
        out[s - start: e - start] = codes
        if nruns.size:
            i0 = np.searchsorted(nruns[:, 1], s, side="right")
            for a, b in nruns[i0:]:
                if a >= e:
                    break
                a2, b2 = max(a, s), min(b, e)
                if b2 > a2:
                    out[a2 - start: b2 - start] = 4
        return out


# --------------------------------------------------------------------------- motifs
def load_jaspar(path: str) -> List[Tuple[str, str, np.ndarray]]:
    """-> [(matrix_id, TF name, pfm (4, L))]; every matrix of the file."""
    out, mid, name, rows = [], None, None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if mid and len(rows) == 4:
                    out.append((mid, name, np.array(rows, dtype=np.float32)))
                parts = line[1:].split()
                mid, name, rows = parts[0], (parts[1] if len(parts) > 1 else parts[0]), []
            elif line:
                toks = line.replace("[", " ").replace("]", " ").split()
                rows.append([float(t) for t in toks if t not in ("A", "C", "G", "T")])
        if mid and len(rows) == 4:
            out.append((mid, name, np.array(rows, dtype=np.float32)))
    return out


def pfm_to_pwm(pfm: np.ndarray, pseudocount: float = 0.8) -> np.ndarray:
    c = pfm + pseudocount
    return np.log2(c / c.sum(0, keepdims=True) / 0.25).astype(np.float32)


def motifs_for_tfs(jaspar: List[Tuple[str, str, np.ndarray]], tf_names: Sequence[str]
                   ) -> Dict[str, List[np.ndarray]]:
    """TF name (upper) -> list of PWMs; 'A::B' matrices count for A and for B."""
    want = {t.upper() for t in tf_names}
    by_tf: Dict[str, List[np.ndarray]] = {}
    for _mid, name, pfm in jaspar:
        for part in name.upper().split("::"):
            part = part.strip()
            if part in want:
                by_tf.setdefault(part, []).append(pfm_to_pwm(pfm))
    return by_tf


@torch.no_grad()
def scan_best(seqs: List[np.ndarray], pwms: List[np.ndarray], device: str = "cpu",
              peak_batch: int = 64, motif_batch: int = 256) -> np.ndarray:
    """Best (forward or reverse-complement) log-odds score / PWM max, per
    (sequence, PWM). Only windows fully inside the sequence count; -inf-like
    (= -1) where the sequence is shorter than the motif. -> (n_seqs, n_pwms) float32."""
    n, m = len(seqs), len(pwms)
    out = np.full((n, m), -1.0, dtype=np.float32)
    if n == 0 or m == 0:
        return out
    Lm = max(p.shape[1] for p in pwms)
    W = np.zeros((2 * m, 4, Lm), dtype=np.float32)
    mlen = np.zeros(2 * m, dtype=np.int64)
    pmax = np.zeros(2 * m, dtype=np.float32)
    for j, p in enumerate(pwms):
        L = p.shape[1]
        W[j, :, :L] = p
        W[m + j, :, :L] = p[::-1, ::-1]                    # reverse complement (ACGT -> TGCA)
        mlen[j] = mlen[m + j] = L
        pmax[j] = pmax[m + j] = p.max(0).sum()
    Wt = torch.from_numpy(W).to(device)
    mlen_t = torch.from_numpy(mlen).to(device)
    pmax_t = torch.from_numpy(pmax).to(device)
    order = np.argsort([len(s) for s in seqs])
    eye = torch.eye(5, device=device)[:, :4]                # code 4 (N) -> all-zero column
    for b0 in range(0, n, peak_batch):
        idx = order[b0:b0 + peak_batch]
        lens = np.array([len(seqs[i]) for i in idx])
        Lb = int(max(lens.max(), Lm))
        codes = np.full((len(idx), Lb), 4, dtype=np.int64)
        for r, i in enumerate(idx):
            codes[r, :len(seqs[i])] = seqs[i]
        x = eye[torch.from_numpy(codes).to(device)].permute(0, 2, 1).contiguous()   # B,4,Lb
        lens_t = torch.from_numpy(lens).to(device)
        T = Lb - Lm + 1
        pos = torch.arange(T, device=device)
        res = torch.empty((len(idx), m), device=device)
        for m0 in range(0, m, motif_batch):
            m1 = min(m, m0 + motif_batch)
            sel = torch.cat([torch.arange(m0, m1), torch.arange(m + m0, m + m1)]).to(device)
            y = torch.nn.functional.conv1d(x, Wt[sel])      # B, 2k, T  (windows start at 0..T-1)
            valid = (pos[None, None, :] + mlen_t[sel][None, :, None]) <= lens_t[:, None, None]
            y = torch.where(valid, y, torch.full_like(y, -1e9))
            best = y.max(dim=2).values / pmax_t[sel][None, :]
            k = m1 - m0
            res[:, m0:m1] = torch.maximum(best[:, :k], best[:, k:])
        res = torch.where(res < -1e6, torch.full_like(res, -1.0), res)
        out[idx] = res.cpu().numpy()
    return out


def pair_features(gene_top_peaks: Dict[str, list], tf_nodes: Sequence[int], names: Sequence[str],
                  genome: Genome2bit, jaspar_path: str, n_genes: int, top_k: int = 5,
                  max_len: int = 2000, hit_frac: float = 0.8, device: str = "cpu") -> Tuple[np.ndarray, dict]:
    """-> (F float16 (n_tf_nodes, n_genes, 3) = [best rel score, #hit peaks / top_k, has_motif], info)."""
    tf_names = [names[t] for t in tf_nodes]
    by_tf = motifs_for_tfs(load_jaspar(jaspar_path), tf_names)
    pwms, owner = [], []
    for r, t in enumerate(tf_names):
        for p in by_tf.get(t.upper(), []):
            pwms.append(p)
            owner.append(r)
    owner = np.asarray(owner, dtype=np.int64)
    # unique peaks over the genes' top-k lists (centre-cropped to max_len)
    peak_id: Dict[Tuple[str, int, int], int] = {}
    gene_peaks: Dict[int, List[int]] = {}
    for g, pk in gene_top_peaks.items():
        ids = []
        for c, s, e in pk[:top_k]:
            s, e = int(s), int(e)
            if e - s > max_len:
                mid = (s + e) // 2
                s, e = mid - max_len // 2, mid + max_len // 2
            key = (str(c), s, e)
            if key not in peak_id:
                peak_id[key] = len(peak_id)
            ids.append(peak_id[key])
        gene_peaks[int(g)] = ids
    keys = sorted(peak_id, key=peak_id.get)
    seqs = [genome.fetch(c, s, e) for c, s, e in keys]
    n_known = sum(int((q < 4).mean() > 0.5) for q in seqs)
    S = scan_best(seqs, pwms, device=device)                # n_peaks x n_pwms
    n_tf = len(tf_nodes)
    # per peak x TF: max over the TF's PWMs
    P = np.full((len(seqs), n_tf), -1.0, dtype=np.float32)
    for r in np.unique(owner):
        P[:, r] = S[:, owner == r].max(axis=1)
    F = np.zeros((n_tf, n_genes, 3), dtype=np.float16)
    has = np.zeros(n_tf, dtype=np.float32)
    has[np.unique(owner)] = 1.0
    F[:, :, 2] = has[:, None]
    for g, ids in gene_peaks.items():
        if not ids or g >= n_genes:
            continue
        sub = P[ids]                                         # k x n_tf
        F[:, g, 0] = np.clip(sub.max(axis=0), 0.0, 1.0)
        F[:, g, 1] = (sub >= hit_frac).sum(axis=0) / float(top_k)
    F[has == 0, :, :2] = 0.0
    info = {"n_tf_nodes": n_tf, "n_tf_with_motif": int(has.sum()), "n_pwms": len(pwms),
            "n_peaks_scanned": len(seqs), "n_peaks_mostly_ACGT": int(n_known),
            "top_k": top_k, "max_len": max_len, "hit_frac": hit_frac,
            "jaspar": Path(jaspar_path).name}
    return F, info
