#!/usr/bin/env python
"""Build the PBMC10k / LINGER-Cistrome benchmark inputs (CPU; label-free features).

Steps (--steps, comma-separated; default cells,processed,splits):
  cells      LINGER cell QC (label file, mt<5%, >=200 genes RNA, >=200 peaks ATAC,
             intersection) + gene universe (filter_genes(min_cells=3) on RNA, + eval TFs)
             -> <root>/cells/cells.tsv, genes.txt
  processed  one MeVD-GRN processed dir per evaluated cell type (label-free graphs)
             -> <root>/processed/<ct>/
  fm         Geneformer V2-104M embeddings for the universe (HF cache needed)
  splits     per label source: TF-disjoint / target-disjoint / random splits with
             degree-matched + uniform negatives, eval TFs removed as regulators
             -> <root>/splits/<source>/{tf,target,random}.npz, meta.json

Usage:
  python scripts/26_pbmc_build.py --config configs/pbmc/benchmark.yaml --root $PBMC_ROOT
  python scripts/26_pbmc_build.py ... --steps processed --cell_types naive_b --max_cells 300  # smoke
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from src.benchmarks import pbmc_data as pdx      # noqa: E402
from src.benchmarks import pbmc_eval as pe       # noqa: E402
from src.benchmarks import pbmc_labels as pl     # noqa: E402
from src.utils.io import load_config, save_json  # noqa: E402


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


def eval_tfs(cfg) -> list:
    return sorted({d["tf"].upper() for d in cfg["groundtruth"]["datasets"]})


def load_all(cfg, root):
    rna, genes, atac, peaks, bc = pdx.load_multiome(str(root / cfg["h5"]))
    rna, genes = pdx.dedup_columns(rna, genes)
    return rna, genes, atac, peaks, bc


def step_cells(cfg, root, rna, genes, atac, bc):
    q = cfg["qc"]
    lab = pd.read_csv(root / cfg["cell_labels"], sep="\t")
    lab = dict(zip(lab["barcode_use"], lab["label"]))
    in_lab = np.array([b in lab for b in bc])
    mt = np.array([g.startswith("MT-") for g in genes])
    tot = np.asarray(rna.sum(axis=1)).ravel()
    pct_mt = 100 * np.asarray(rna[:, mt].sum(axis=1)).ravel() / np.maximum(tot, 1)
    ng = np.asarray((rna > 0).sum(axis=1)).ravel()
    npk = np.asarray((atac > 0).sum(axis=1)).ravel()
    mt_ok = np.ones_like(in_lab) if q.get("max_pct_mt") is None else (pct_mt < q["max_pct_mt"])
    rna_ok = in_lab & mt_ok & (ng >= q["rna_min_genes"])
    atac_ok = in_lab & (npk >= q["atac_min_peaks"])
    keep = rna_ok & atac_ok
    # genes: LINGER filters genes on the RNA object after its own cell filter (pre-intersection)
    det = np.asarray((rna[rna_ok] > 0).sum(axis=0)).ravel()
    uni = sorted({g for g, d in zip(genes, det) if d >= q["min_cells_per_gene"]} |
                 {t for t in eval_tfs(cfg) if t in set(genes)})
    out = root / "cells"
    out.mkdir(parents=True, exist_ok=True)
    cells = pd.DataFrame({"barcode": [b for b, k in zip(bc, keep) if k],
                          "label": [lab[b] for b, k in zip(bc, keep) if k]})
    cells.to_csv(out / "cells.tsv", sep="\t", index=False)
    (out / "genes.txt").write_text("\n".join(uni) + "\n")
    counts = cells["label"].value_counts().to_dict()
    summ = {"labelled": int(in_lab.sum()), "pass_rna": int(rna_ok.sum()), "pass_atac": int(atac_ok.sum()),
            "kept": int(keep.sum()), "genes": len(uni), "per_label": counts,
            "missing_eval_tfs": [t for t in eval_tfs(cfg) if t not in set(uni)]}
    save_json(summ, out / "summary.json")
    log("cells:", json.dumps({k: v for k, v in summ.items() if k != "per_label"}))
    log("evaluated cell types:", {ct: counts.get(lb, 0) for ct, lb in cfg["cell_types"].items()})


def regulator_names(cfg, root, universe):
    regs = set(eval_tfs(cfg))
    for name, f in cfg["label_files"].items():
        df = pl.load_label_source(name, str(root / f))
        regs |= set(df.tf)
    return sorted(regs & set(universe))


def step_processed(cfg, mcfg, root, rna, genes, atac, peaks, bc, cts, max_cells):
    cells = pd.read_csv(root / "cells" / "cells.tsv", sep="\t")
    universe = [l.strip() for l in open(root / "cells" / "genes.txt") if l.strip()]
    regs = regulator_names(cfg, root, universe)
    bidx = {b: i for i, b in enumerate(bc)}
    for ct in cts:
        lb = cfg["cell_types"][ct]
        rows = np.array([bidx[b] for b in cells.loc[cells.label == lb, "barcode"]])
        if max_cells:
            rows = rows[:max_cells]
        mask = np.zeros(len(bc), dtype=bool)
        mask[rows] = True
        log(f"processed {ct}: {mask.sum()} cells, {len(universe)} genes, {len(regs)} regulators")
        summ = pdx.build_processed(root / "processed" / ct, rna, genes, atac, peaks, mask, universe, regs,
                                   str(root / cfg["gtf"]), mcfg["data"], mcfg.get("atac", {}) or {}, ct)
        log(f"processed {ct} ->", json.dumps(summ))


def step_fm(root, cts):
    spec = importlib.util.spec_from_file_location("fm11", REPO / "scripts" / "11_extract_fm_embeddings.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    universe = [l.strip() for l in open(root / "cells" / "genes.txt") if l.strip()]
    emb = mod.extract_embeddings(universe)
    for ct in cts:
        p = root / "processed" / ct
        gi = json.load(open(p / "gene_index.json"))
        assert [g for g, _ in sorted(gi.items(), key=lambda x: x[1])] == universe
        np.save(p / "fm_gene_embeddings.npy", emb.astype(np.float32))
    log(f"FM embeddings {emb.shape}, matched rows {(np.abs(emb).sum(1) > 0).sum()}")


def step_splits(cfg, root, sources):
    universe = [l.strip() for l in open(root / "cells" / "genes.txt") if l.strip()]
    gi = {g: i for i, g in enumerate(universe)}
    ev = eval_tfs(cfg)
    s = cfg["split"]
    for src in sources:
        df = pl.load_label_source(src, str(root / cfg["label_files"][src]))
        kept, info = pl.restrict(df, universe, ev)
        all_pos = df[df.tf.isin(gi) & df.target.isin(gi)]
        out = root / "splits" / src
        out.mkdir(parents=True, exist_ok=True)
        meta = {"source": src, "n_genes": len(universe), "eval_tfs": ev, "filter": info, "regimes": {}}
        for reg in pl.REGIMES:
            parts = pl.split_labels(kept, reg, int(s["seed"]), tuple(s["ratios"]))
            arr = pl.build_split_arrays(parts, all_pos, reg, gi, int(s["seed"]))
            rep = pl.check_disjoint(arr, reg, len(gi))
            evi = {gi[t] for t in ev if t in gi}
            for p in arr:
                assert not (set(arr[p]["pos"][:, 0]) & evi), f"eval TF is a regulator in {reg}.{p}"
            np.savez_compressed(out / f"{reg}.npz", **{f"{p}__{k}": v for p in arr for k, v in arr[p].items()})
            meta["regimes"][reg] = {"sizes": {p: {k: int(len(v)) for k, v in arr[p].items()} for p in arr},
                                    "n_tfs": {p: int(len(np.unique(arr[p]["pos"][:, 0]))) for p in arr},
                                    "checks": rep}
            log(f"splits {src}/{reg}:", json.dumps(meta["regimes"][reg]["sizes"]))
        save_json(meta, out / "meta.json")
        log(f"splits {src}: {json.dumps(info)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(REPO / "configs/pbmc/benchmark.yaml"))
    ap.add_argument("--model_config", default=str(REPO / "configs/pbmc/mevd_base.yaml"),
                    help="only its data/atac blocks are used (graph degrees, RP decay)")
    ap.add_argument("--root", default=os.environ.get("PBMC_ROOT", "data/pbmc"))
    ap.add_argument("--steps", default="cells,processed,splits")
    ap.add_argument("--cell_types", default=None, help="comma list; default all in config")
    ap.add_argument("--sources", default=",".join(pl.SOURCES))
    ap.add_argument("--max_cells", type=int, default=0, help="smoke tests only")
    args = ap.parse_args()
    cfg = load_config(args.config)
    mcfg = load_config(args.model_config)
    root = Path(args.root)
    steps = args.steps.split(",")
    cts = args.cell_types.split(",") if args.cell_types else list(cfg["cell_types"])
    data = None
    if "cells" in steps or "processed" in steps:
        t = time.time()
        data = load_all(cfg, root)
        log(f"loaded h5 in {time.time() - t:.0f}s: RNA {data[0].shape}, ATAC {data[2].shape}")
    if "cells" in steps:
        step_cells(cfg, root, data[0], data[1], data[2], data[4])
    if "processed" in steps:
        step_processed(cfg, mcfg, root, *data, cts, args.max_cells)
    if "fm" in steps:
        step_fm(root, cts)
    if "splits" in steps:
        step_splits(cfg, root, args.sources.split(","))


if __name__ == "__main__":
    main()
