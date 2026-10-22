#!/usr/bin/env python
"""Label-free TF-TF prior network from STRING (human, v12.0), one file per cell type.

For the scMultiomeGRN ablation 'string_prior': the official model normally message-passes over the fold's
training links (part of the motif answer key). This replaces that graph with STRING protein-protein
associations among the same TFs, which never uses the answer key.

  STRING proteins are mapped to TF symbols through protein.info `preferred_name` (uppercase match);
  an edge is kept when combined_score >= --min_score (default 700, STRING's 'high confidence').
  Output: graph/<ds>/<ct>/x_string_graph.txt, same format as x_graph.txt (tab-separated, TF names as
  index and header, 0/1, symmetric, zero diagonal), plus a coverage / overlap report.

Usage: python scripts/24_build_string_tf_graph.py --root ~/scmgrn --links <links.txt.gz> --info <info.txt.gz>
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from src.benchmarks import scmgrn_protocol as proto   # noqa: E402
from src.utils.io import load_config                  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--links", required=True)
    ap.add_argument("--info", required=True)
    ap.add_argument("--min_score", type=int, default=700)
    a = ap.parse_args()
    root = Path(a.root).expanduser()

    info = pd.read_csv(a.info, sep="\t", usecols=[0, 1], dtype=str)
    info.columns = ["protein", "name"]
    prot2name = dict(zip(info["protein"], info["name"].str.upper()))
    links = pd.read_csv(a.links, sep=" ", dtype={"protein1": str, "protein2": str, "combined_score": int})
    links = links[links["combined_score"] >= a.min_score]
    links["a"] = links["protein1"].map(prot2name)
    links["b"] = links["protein2"].map(prot2name)
    links = links.dropna(subset=["a", "b"])
    print(f"STRING links with score >= {a.min_score}: {len(links):,}", flush=True)

    report = {}
    for ds in ("lung", "hema"):
        cfg = load_config(str(REPO / f"configs/scmgrn/{ds}.yaml"))
        for ct in cfg["cell_types"]:
            gdir = root / cfg["graph_dir"] / ct
            if not (gdir / "x_graph.txt").exists():
                continue
            adj, nodes = proto.load_adj(str(gdir / "x_graph.txt"))
            idx = {n.upper(): i for i, n in enumerate(nodes)}
            sub = links[links["a"].isin(idx) & links["b"].isin(idx)]
            S = np.zeros((len(nodes), len(nodes)), dtype=int)
            S[sub["a"].map(idx).to_numpy(), sub["b"].map(idx).to_numpy()] = 1
            S = np.maximum(S, S.T)
            np.fill_diagonal(S, 0)
            pd.DataFrame(S, index=nodes, columns=nodes).to_csv(gdir / "x_string_graph.txt", sep="\t")
            iu = np.triu_indices(len(nodes), 1)
            s_e, a_e = S[iu].astype(bool), adj[iu].astype(bool)
            report[f"{ds}/{ct}"] = {
                "n_tfs": len(nodes),
                "tfs_with_string_edge": int((S.sum(1) > 0).sum()),
                "string_edges": int(s_e.sum()),
                "answer_key_edges": int(a_e.sum()),
                "overlap_edges": int((s_e & a_e).sum()),
                "frac_answer_key_in_string": float((s_e & a_e).sum() / max(a_e.sum(), 1)),
                "frac_string_in_answer_key": float((s_e & a_e).sum() / max(s_e.sum(), 1)),
            }
            print(f"{ds}/{ct}: {report[f'{ds}/{ct}']}", flush=True)
    out = root / "graph" / "string_report.json"
    json.dump(report, open(out, "w"), indent=1)
    print("report ->", out)


if __name__ == "__main__":
    main()
