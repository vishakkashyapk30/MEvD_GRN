#!/usr/bin/env python
"""Recompute every metric (scmgrn_protocol.evaluate_all) for fold results that saved their score matrix
(scores_foldNN.npy next to foldNN.json), and write the new metric sets into the JSON. Existing keys are replaced
by the recomputed values; nothing else in the file changes. Used when a new evaluation set is added.

Usage: python scripts/31_scmgrn_backfill_metrics.py --root ~/scmgrn
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from src.benchmarks import scmgrn_protocol as proto   # noqa: E402
from src.utils.io import load_config, save_json       # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    root = Path(a.root).expanduser()
    cfgs = {ds: load_config(str(REPO / f"configs/scmgrn/{ds}.yaml")) for ds in ("lung", "hema")}
    cache, n_done = {}, 0
    for sp in sorted(root.glob("results/*/*/*/seed*/scores_fold*.npy")):
        res_ds, ct = sp.parts[-5], sp.parts[-4]
        ds, _, kind = res_ds.partition("_")
        kind = kind or "official"
        cfg = cfgs[ds]
        fold = int(sp.stem.replace("scores_fold", ""))
        jp = sp.with_name(f"fold{fold:02d}.json")
        if not jp.exists():
            continue
        key = (ds, ct, kind)
        if key not in cache:
            adj, _ = proto.load_adj(str(root / cfg["graph_dir"] / ct / "x_graph.txt"))
            cache[key] = (adj, proto.resolve_split_dir(str(root / cfg["splits_dir"] / ct), cfg["split"], kind))
        adj, sdir = cache[key]
        f = proto.load_fold(sdir, fold, adj.shape[0])
        d = json.load(open(jp))
        d["metrics"].update(proto.evaluate_all(np.load(sp), adj, f, sdir, fold))
        save_json(d, jp)
        n_done += 1
    print(f"[backfill] updated {n_done} fold files")


if __name__ == "__main__":
    main()
