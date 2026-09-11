#!/usr/bin/env python
"""Parallel ranged-HTTP downloader for scmogrndb.psu.edu (single connections are
server-throttled to ~40-70KB/s; many concurrent range requests bypass that
per-connection cap). Usage:

  python scripts/_parallel_download.py <url> <outpath> [n_workers]
"""
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

N_WORKERS_DEFAULT = 16


def head_size(url: str) -> int:
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req) as r:
        return int(r.headers["Content-Length"])


def fetch_range(url: str, start: int, end: int) -> bytes:
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()


def download(url: str, out_path: str, n_workers: int = N_WORKERS_DEFAULT) -> None:
    size = head_size(url)
    chunk = max(size // n_workers, 1)
    ranges = [(i, min(i + chunk - 1, size - 1)) for i in range(0, size, chunk)]
    buf = [None] * len(ranges)
    done = 0
    with ThreadPoolExecutor(max_workers=n_workers) as ex:
        futs = {ex.submit(fetch_range, url, s, e): i for i, (s, e) in enumerate(ranges)}
        for fut in as_completed(futs):
            i = futs[fut]
            buf[i] = fut.result()
            done += 1
            print(f"\r[{out_path}] {done}/{len(ranges)} chunks ({size/1e6:.1f} MB total)",
                  end="", flush=True)
    print()
    with open(out_path, "wb") as f:
        for b in buf:
            f.write(b)


if __name__ == "__main__":
    url, out = sys.argv[1], sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else N_WORKERS_DEFAULT
    download(url, out, workers)
