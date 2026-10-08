#!/usr/bin/env python
"""Fetch BEAR-GRN's released data (Zenodo 10.5281/zenodo.20704929) member by
member, using HTTP range requests into the remote zips, so only the needed
files are downloaded (the full record is ~21.7 GB).

Run on a machine with internet (Ada LOGIN node; compute nodes are offline).

  python scripts/20_bear_download.py --root $BEAR_ROOT --datasets K562 Macrophage_S1 \
      --what input gt inferred
  python scripts/20_bear_download.py --root ~/.cache/bear/data --datasets K562 --what gt input
  python scripts/20_bear_download.py --list INFERRED.GRNS.zip

--what:
  input      INPUT.DATA/<ds>/*            (preprocessed RNA + ATAC used by every method)
  gt         every GT variant file for <ds> (GROUND.TRUTHS, .KO, .UNION, .INETRSECTION,
             CORE_GROUND_TRUTH, CELL.TYPE.EXCLUSIVE.GROUND.TRUTH)
  inferred   INFERRED.GRNS/<ds>/<method>/*  (released outputs of the 9 methods)
  stability  STABILITY_GRNS/<ds>/**         (released subsample networks; large)
  stabinput  INPUT.DATA.STABILITY/<ds>/**   (subsampled inputs for our stability runs)
Files keep their in-zip paths under --root. Existing files of the right size are skipped.
"""
from __future__ import annotations

import argparse
import io
import os
import shutil
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

try:  # INFERRED.GRNS.zip (and others) use Deflate64 (method 9), unsupported by stdlib zipfile
    import zipfile_deflate64  # noqa: F401  (pip install zipfile-deflate64; patches zipfile)
except ImportError:
    zipfile_deflate64 = None

RECORD = "20704929"
URL = "https://zenodo.org/api/records/{rec}/files/{name}/content"

# dataset name (INPUT.DATA / INFERRED.GRNS dir) -> substrings identifying its GT files
DATASETS = {
    "K562": {"inferred": "K562", "gt": ["_K562.tsv"]},
    "Macrophage_S1": {"inferred": "Macrophage_S1", "gt": ["_Buffer1.tsv"]},
    "Macrophage_S2": {"inferred": "Macrophage_S2", "gt": ["_Buffer2.tsv", "_Buffer1.tsv"]},
    "iPS": {"inferred": "iPS", "gt": ["_iPS.tsv"]},
    "mESC_E7.5_rep1": {"inferred": "E7.5_rep1", "gt": ["_E7.5_rep1.tsv"]},
    "mESC_E7.5_rep2": {"inferred": "E7.5_rep2", "gt": ["_E7.5_rep2.tsv"]},
    "mESC_E8.5_rep1": {"inferred": "E8.5_rep1", "gt": ["_E8.5_rep1.tsv"]},
    "mESC_E8.5_rep2": {"inferred": "E8.5_rep2", "gt": ["_E8.5_rep2.tsv"]},
    "Naive_mESC": {"inferred": "Naive_mESC", "gt": ["_Naive_mESC.tsv"]},
}
GT_ZIPS = ["GROUND.TRUTHS.zip", "GROUND.TRUTHS.KO.zip", "GROUND.TRUTHS.UNION.zip",
           "GROUND.TRUTHS.INETRSECTION.zip", "CORE_GROUND_TRUTH.zip",
           "CELL.TYPE.EXCLUSIVE.GROUND.TRUTH.zip"]


class HTTPFile(io.RawIOBase):
    """Seekable read-only file over HTTP range requests."""

    def __init__(self, url: str):
        self.url, self.pos, self.fetched = url, 0, 0
        for k in range(10):                     # Zenodo returns sporadic 5xx / timeouts
            try:
                with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=120) as r:
                    self.size = int(r.headers["Content-Length"])
                break
            except Exception:
                if k == 9:
                    raise
                time.sleep(10 * (k + 1))

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off)
        return self.pos

    def read(self, n=-1):
        if n < 0:
            n = self.size - self.pos
        if n == 0 or self.pos >= self.size:
            return b""
        end = min(self.pos + n, self.size) - 1
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
        for k in range(8):
            try:
                with urllib.request.urlopen(req, timeout=300) as r:
                    data = r.read()
                break
            except Exception:
                if k == 7:
                    raise
                time.sleep(5 * (k + 1))
        self.pos += len(data)
        self.fetched += len(data)
        return data

    def readinto(self, b):
        d = self.read(len(b))
        b[:len(d)] = d
        return len(d)


def open_zip(name: str):
    f = HTTPFile(URL.format(rec=RECORD, name=name))
    return f, zipfile.ZipFile(io.BufferedReader(f, buffer_size=1 << 20))


def _fetch_ranges(url: str, ranges, dst: str, conns: int, chunk: int = 4 << 20) -> int:
    """Write the byte ranges [a, b) of `url` into the (sparse) file `dst` at the
    same offsets, with `conns` parallel range requests (Zenodo throttles a single
    connection to ~10-50 KB/s, so parallelism is what makes this usable)."""
    from concurrent.futures import ThreadPoolExecutor

    jobs = []
    for a, b in ranges:
        for s in range(a, b, chunk):
            jobs.append((s, min(s + chunk, b)))

    def get(job):
        a, b = job
        req = urllib.request.Request(url, headers={"Range": f"bytes={a}-{b - 1}"})
        for k in range(12):
            try:
                with urllib.request.urlopen(req, timeout=300) as r:
                    data = r.read()
                if len(data) != b - a:
                    raise IOError(f"short read {len(data)}/{b - a}")
                break
            except Exception:
                if k == 11:
                    raise
                time.sleep(min(60, 5 * (k + 1)))
        with open(dst, "r+b") as f:
            f.seek(a)
            f.write(data)
        return len(data)

    with ThreadPoolExecutor(conns) as ex:
        return sum(ex.map(get, jobs))


NO_EXTRACT = False
TOUCHED = set()      # sparse zips written by this run (only these are cleaned up)


def extract(zname: str, want, root: Path, conns: int = 16) -> int:
    """Extract the members selected by `want` from a remote zip, downloading only
    their bytes (+ the central directory) into a sparse local copy of the zip,
    then unpacking with zipfile (Deflate: stdlib; Deflate64: zipfile-deflate64)
    or Info-ZIP `unzip` (handles Deflate64 natively)."""
    import struct
    import subprocess

    url = URL.format(rec=RECORD, name=zname)
    f, z = open_zip(zname)                      # reads only the central directory
    todo = [i for i in z.infolist() if not i.is_dir() and want(i.filename)
            and not ((root / i.filename).exists() and (root / i.filename).stat().st_size == i.file_size)]
    for i in z.infolist():
        if not i.is_dir() and want(i.filename) and i not in todo:
            print(f"[skip] {root / i.filename}", flush=True)
    if not todo:
        return 0
    size = f.size
    sparse = root / "zips" / (zname + ".sparse")
    if not NO_EXTRACT:
        TOUCHED.add(sparse)
    sparse.parent.mkdir(parents=True, exist_ok=True)
    if not sparse.exists():
        with open(sparse, "wb") as g:
            g.truncate(size)
    # central directory + EOCD (zip64 records included): from the first CD entry to EOF
    cd_start = min(z.start_dir, size)
    _fetch_ranges(url, [(cd_start, size)], str(sparse), conns)
    # each member: local header (30 B + name + extra) + data (+ 16..24 B descriptor)
    hdr = [(i.header_offset, min(i.header_offset + 30, size)) for i in todo]
    _fetch_ranges(url, hdr, str(sparse), conns, chunk=30)
    ranges = []
    with open(sparse, "rb") as g:
        for i in todo:
            g.seek(i.header_offset)
            h = g.read(30)
            n_name, n_extra = struct.unpack("<HH", h[26:30])
            a = i.header_offset
            ranges.append((a, min(size, a + 30 + n_name + n_extra + i.compress_size + 24)))
    t0 = time.time()
    got = _fetch_ranges(url, ranges, str(sparse), conns)
    print(f"[fetch] {zname}: {len(todo)} members, {got/1e6:.1f} MB in {time.time()-t0:.0f}s", flush=True)
    if NO_EXTRACT:      # members stay compressed inside the sparse zip (read with bear_data.read_csv_sparse)
        print(f"[keep] {len(todo)} members compressed in {sparse}", flush=True)
        return len(todo)
    n = 0
    for i in todo:
        out = root / i.filename
        out.parent.mkdir(parents=True, exist_ok=True)
        if i.compress_type == 9 and zipfile_deflate64 is None:
            subprocess.run(["unzip", "-o", "-q", str(sparse), i.filename, "-d", str(root)], check=True)
        else:
            with zipfile.ZipFile(sparse) as zz, zz.open(i.filename) as src, \
                    open(str(out) + ".part", "wb") as dst:
                shutil.copyfileobj(src, dst, length=1 << 22)
            os.replace(str(out) + ".part", out)
        assert out.stat().st_size == i.file_size, f"size mismatch for {out}"
        n += 1
        print(f"[get] {i.filename} ({i.file_size/1e6:.1f} MB)", flush=True)
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.environ.get("BEAR_ROOT", "data/bear"))
    ap.add_argument("--datasets", nargs="*", default=["K562"])
    ap.add_argument("--what", nargs="*", default=["input", "gt"])
    ap.add_argument("--methods", nargs="*", default=None, help="subset of INFERRED.GRNS methods")
    ap.add_argument("--list", default=None, help="print a zip's member list and exit")
    ap.add_argument("--conns", type=int, default=16, help="parallel range requests")
    ap.add_argument("--no_extract", action="store_true",
                    help="do not decompress; keep members inside <root>/zips/<zip>.sparse (implies --keep_sparse)")
    ap.add_argument("--keep_sparse", action="store_true",
                    help="keep <root>/zips/*.sparse (holes take no disk; reuse speeds reruns)")
    args = ap.parse_args()
    if args.list:
        f, z = open_zip(args.list)
        for i in z.infolist():
            print(f"{i.file_size:>14d} {i.compress_size:>14d} {i.filename}")
        return
    root = Path(args.root)
    root.mkdir(parents=True, exist_ok=True)
    global NO_EXTRACT
    NO_EXTRACT = bool(args.no_extract)
    if NO_EXTRACT:
        args.keep_sparse = True
    for ds in args.datasets:
        if ds not in DATASETS:
            sys.exit(f"unknown dataset {ds}; choose from {list(DATASETS)}")
        meta = DATASETS[ds]
        inf = meta["inferred"]
        for what in args.what:
            if what == "input":
                extract("INPUT.DATA.zip", lambda p: p.startswith(f"INPUT.DATA/{ds}/"), root, args.conns)
            elif what == "stabinput":
                extract("INPUT.DATA.STABILITY.zip", lambda p: f"/{ds}/" in p or f"/{inf}/" in p, root, args.conns)
            elif what == "gt":
                for zn in GT_ZIPS:
                    extract(zn, lambda p: any(p.endswith(s) for s in meta["gt"]), root, args.conns)
            elif what == "inferred":
                def want(p, inf=inf):
                    parts = p.split("/")
                    return (len(parts) >= 4 and parts[1] == inf and
                            (args.methods is None or parts[2] in args.methods))
                extract("INFERRED.GRNS.zip", want, root, args.conns)
            elif what == "stability":
                extract("STABILITY_GRNS.zip", lambda p: f"/{inf}/" in p or f"/{ds}/" in p, root, args.conns)
            else:
                sys.exit(f"unknown --what {what}")
    if not args.keep_sparse:
        for sp in TOUCHED:
            if sp.exists():
                sp.unlink()


if __name__ == "__main__":
    main()
