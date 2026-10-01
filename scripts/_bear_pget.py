#!/usr/bin/env python3
"""Parallel HTTP range downloader (stdlib only, Python >= 3.6) for Zenodo, which
throttles single connections to ~50-200 KB/s. Resumable: finished chunks are
kept in <out>.parts/ and skipped on rerun.

  python3 scripts/_bear_pget.py URL OUT [--conns 16] [--chunk_mb 16]
"""
import argparse
import os
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed


def head_size(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=120) as r:
        return int(r.headers["Content-Length"])


def get_range(url, a, b, path):
    if os.path.exists(path) and os.path.getsize(path) == b - a + 1:
        return 0
    for k in range(12):
        try:
            req = urllib.request.Request(url, headers={"Range": "bytes=%d-%d" % (a, b)})
            with urllib.request.urlopen(req, timeout=300) as r:
                data = r.read()
            if len(data) != b - a + 1:
                raise IOError("short read %d/%d" % (len(data), b - a + 1))
            with open(path + ".tmp", "wb") as f:
                f.write(data)
            os.replace(path + ".tmp", path)
            return len(data)
        except Exception as e:  # noqa: BLE001
            if k == 11:
                raise
            time.sleep(min(60, 5 * (k + 1)))
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("url")
    ap.add_argument("out")
    ap.add_argument("--conns", type=int, default=16)
    ap.add_argument("--chunk_mb", type=int, default=16)
    a = ap.parse_args()
    size = head_size(a.url)
    if os.path.exists(a.out) and os.path.getsize(a.out) == size:
        print("[skip] %s (%d bytes)" % (a.out, size), flush=True)
        return
    pdir = a.out + ".parts"
    os.makedirs(pdir, exist_ok=True)
    cs = a.chunk_mb << 20
    ranges = [(s, min(s + cs, size) - 1) for s in range(0, size, cs)]
    t0, done = time.time(), 0
    with ThreadPoolExecutor(a.conns) as ex:
        futs = {ex.submit(get_range, a.url, s, e, os.path.join(pdir, "%06d" % i)): i
                for i, (s, e) in enumerate(ranges)}
        for n, f in enumerate(as_completed(futs), 1):
            done += f.result()
            if n % 10 == 0 or n == len(ranges):
                el = time.time() - t0
                print("  %s: %d/%d chunks, %.1f MB/s" % (os.path.basename(a.out), n, len(ranges),
                                                       done / 1e6 / max(el, 1e-6)), flush=True)
    with open(a.out + ".tmp", "wb") as out:
        for i in range(len(ranges)):
            with open(os.path.join(pdir, "%06d" % i), "rb") as f:
                out.write(f.read())
    assert os.path.getsize(a.out + ".tmp") == size, "size mismatch"
    os.replace(a.out + ".tmp", a.out)
    for i in range(len(ranges)):
        os.remove(os.path.join(pdir, "%06d" % i))
    os.rmdir(pdir)
    print("[done] %s %d bytes in %.0fs" % (a.out, size, time.time() - t0), flush=True)


if __name__ == "__main__":
    sys.exit(main())
