#!/usr/bin/env python
"""BEAR pooled AUROC / AUPRC (verified port) for the MID-GRN outputs, cached next to the existing scores:
$BEAR_ROOT/results/<ds>/scores/<variant>__seed<s>__<gt>.json (same format as scripts/bear/24).
mid_full: 20 draws (as pre-registered); other variants: 5 draws (README D3)."""
from __future__ import annotations

import argparse
import json

import numpy as np

import mdh_lib as L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datasets", nargs="*", default=L.DATASETS)
    ap.add_argument("--variants", nargs="*", default=["mid_full", "mid_dataonly", "mid_bce", "mid2_full", "mid2_marginal"])
    args = ap.parse_args()
    for ds in args.datasets:
        gts = [g for g in ("ChIP", "KO", "Union", "Intersection", "Core") if L.gt_full(ds, g) is not None]
        for v in args.variants:
            for f in sorted((L.ROOT / "results" / ds / v).glob("seed*/grn.tsv.gz")):
                seed = f.parent.name
                e = None
                for gt in gts:
                    out = L.ROOT / "results" / ds / "scores" / f"{v}__{seed}__{gt}.json"
                    if out.exists():
                        continue
                    e = L.read_ours(f) if e is None else e
                    r = L.pooled_score(e, ds, gt, n_rep=20 if v in ("mid_full", "mid2_full") else 5)
                    r.update({"dataset": ds, "gt": gt, "method": f"{v}/{seed}", "kind": "ours", "file": str(f)})
                    json.dump(r, open(out, "w"), indent=1)
                    print(ds, v, seed, gt, round(r["AUROC"], 4), round(r["AUPRC"], 4), flush=True)


if __name__ == "__main__":
    main()
