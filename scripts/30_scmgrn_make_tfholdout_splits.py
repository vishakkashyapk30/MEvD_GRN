#!/usr/bin/env python
"""Write the held-out-TF folds (src/benchmarks/scmgrn_protocol.make_tf_holdout_splits) for every built cell type.

Usage: python scripts/30_scmgrn_make_tfholdout_splits.py --root ~/scmgrn
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from src.benchmarks import scmgrn_protocol as proto   # noqa: E402
from src.utils.io import load_config                  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    root = Path(a.root).expanduser()
    report = {}
    for ds in ("lung", "hema"):
        cfg = load_config(str(REPO / f"configs/scmgrn/{ds}.yaml"))
        s = cfg["split"]
        for ct in cfg["cell_types"]:
            adj_file = root / cfg["graph_dir"] / ct / "x_graph.txt"
            if not adj_file.exists():
                print(f"[skip] {ds}/{ct}: not built"); continue
            sdir = proto.make_tf_holdout_splits(str(adj_file), str(root / cfg["splits_dir"] / ct),
                                                seed=int(s["seed"]), n_splits=int(s["n_splits"]))
            adj, _ = proto.load_adj(str(adj_file))
            n_c2 = n_c3 = 0
            for k in range(1, int(s["n_splits"]) + 1):
                f = proto.load_fold(sdir, k, adj.shape[0])
                meta = json.load(open(Path(sdir) / f"split-{k:03d}_tfs.json"))
                hid = set(meta["test_tfs"])
                tp = f["test"]["pos"]
                tp = tp[tp[:, 0] < tp[:, 1]]
                c3 = sum(1 for u, v in tp if u in hid and v in hid)
                n_c3 += c3; n_c2 += len(tp) - c3
            report[f"{ds}/{ct}"] = {"n_tfs": int(adj.shape[0]), "test_links_c2": n_c2, "test_links_c3": n_c3}
            print(f"{ds}/{ct}: {report[f'{ds}/{ct}']}", flush=True)
    json.dump(report, open(root / "splits" / "tfholdout_report.json", "w"), indent=1)


if __name__ == "__main__":
    main()
