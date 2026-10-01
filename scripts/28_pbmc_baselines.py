#!/usr/bin/env python
"""Baselines on the PBMC10k / LINGER-Cistrome benchmark, scored by the same code
(src/benchmarks/pbmc_eval.py) over the same gene universe and cells as MeVD-GRN.

Methods (--methods): degree, geneid (label-trained, per source/regime; same
train positives + degree-matched negative pool as MeVD-GRN), pearson, pearson_signed,
grnboost2 (label-free; CPU-heavy, run separately with --methods grnboost2).

Outputs: <root>/results/<ct>/<source>/<regime>/<method>/metrics.json (+ eval_scores.npz)
         label-free methods use source=unsupervised, regime=none.

Usage:
  python scripts/28_pbmc_baselines.py --root $PBMC_ROOT --cell_type naive_b \
      --methods degree,geneid,pearson,pearson_signed --sources collectri --regimes tf
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from src.benchmarks import pbmc_baselines as pb   # noqa: E402
from src.benchmarks import pbmc_data as pdx       # noqa: E402
from src.benchmarks import pbmc_eval as pe        # noqa: E402
from src.benchmarks import pbmc_labels as pl      # noqa: E402
from src.data.preprocessing import _normalize_log  # noqa: E402
from src.utils.io import load_config, load_json, save_json  # noqa: E402

LABEL_FREE = ("pearson", "pearson_signed", "grnboost2")


def save(root, ct, source, regime, method, S_eval, eval_tfs, chip, internal, extra=None):
    out = root / "results" / ct / source / regime / method
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "eval_scores.npz", scores=S_eval.astype(np.float32), tfs=np.array(eval_tfs))
    rec = {"method": method, "cell_type": ct, "source": source, "regime": regime,
           "internal_test": internal, "chip": chip, **(extra or {})}
    save_json(rec, out / "metrics.json")
    s = chip["summary"]
    print(f"[{ct} {source}/{regime} {method}] Cistrome n={s['n']} mean AUROC={s['auroc_mean']:.4f} "
          f"AUPR ratio={s['aupr_ratio_mean']:.3f}" +
          (f" | internal per-TF AUROC={internal.get('per_tf_auroc_mean', float('nan')):.4f} "
           f"dm-AUROC={internal.get('dm_auroc', float('nan')):.4f}" if internal else ""), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(REPO / "configs/pbmc/benchmark.yaml"))
    ap.add_argument("--root", default=os.environ.get("PBMC_ROOT", "data/pbmc"))
    ap.add_argument("--cell_type", required=True)
    ap.add_argument("--methods", default="degree,geneid,pearson,pearson_signed")
    ap.add_argument("--sources", default="collectri")
    ap.add_argument("--regimes", default="tf")
    ap.add_argument("--n_jobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    ap.add_argument("--skip_existing", action="store_true")
    args = ap.parse_args()
    dcfg = load_config(args.config)
    root = Path(args.root)
    ct = args.cell_type
    methods = args.methods.split(",")
    gene_index = load_json(root / "processed" / ct / "gene_index.json")
    genes = [g for g, _ in sorted(gene_index.items(), key=lambda x: x[1])]
    n = len(genes)
    gt = pe.load_groundtruth(dcfg, root)
    eval_tfs = sorted({d["tf"] for d in gt["datasets"] if d["tf"] in gene_index})
    eval_idx = np.array([gene_index[t] for t in eval_tfs])

    X = None
    if any(m in LABEL_FREE for m in methods):          # this cell type's log-normalised RNA
        rna, hgenes, _atac, _peaks, bc = pdx.load_multiome(str(root / dcfg["h5"]))
        rna, hgenes = pdx.dedup_columns(rna, hgenes)
        cells = pd.read_csv(root / "cells" / "cells.tsv", sep="\t")
        keep = set(cells.loc[cells.label == dcfg["cell_types"][ct], "barcode"])
        rows = np.array([i for i, b in enumerate(bc) if b in keep])
        col = {g: i for i, g in enumerate(hgenes)}
        X = _normalize_log(rna[rows][:, [col[g] for g in genes]], {"normalize_total": True,
                                                                   "target_sum": 1e4, "log1p": True})
        print(f"[{ct}] expression {X.shape}", flush=True)

    for m in methods:
        if m in LABEL_FREE:
            out = root / "results" / ct / "unsupervised" / "none" / m / "metrics.json"
            if args.skip_existing and out.exists():
                continue
            t = time.time()
            if m == "grnboost2":
                regs = sorted(set(load_json(root / "processed" / ct / "tf_indices.json")))
                S = pb.grnboost2_scores(X, genes, [genes[i] for i in regs], eval_tfs, n_jobs=args.n_jobs)
            else:
                S = pb.pearson_scores(X, eval_idx, signed=(m == "pearson_signed"))
            chip = pe.evaluate_chip({t_: S[k] for k, t_ in enumerate(eval_tfs)}, genes, ct, gt, dcfg, root)
            save(root, ct, "unsupervised", "none", m, S, eval_tfs, chip, {}, {"seconds": time.time() - t})
            continue
        for src in args.sources.split(","):
            for reg in args.regimes.split(","):
                out = root / "results" / ct / src / reg / m / "metrics.json"
                if args.skip_existing and out.exists():
                    continue
                z = np.load(root / "splits" / src / f"{reg}.npz")
                arrays = {p: {k.split("__")[1]: z[k] for k in z.files if k.startswith(p + "__")}
                          for p in ("train", "val", "test")}
                tr = arrays["train"]
                test_regs = np.unique(arrays["test"]["pos"][:, 0])
                regs_all = np.concatenate([eval_idx, test_regs])
                if m == "degree":
                    S = pb.degree_scores(tr["pos"], n, regs_all)
                elif m == "geneid":
                    S = pb.geneid_scores(tr["pos"], tr["neg"], n, regs_all)
                else:
                    raise SystemExit(f"unknown method {m}")
                S_eval, S_test = S[:len(eval_idx)], S[len(eval_idx):]
                internal = pe.internal_test_metrics(S_test, test_regs, arrays, n, reg)
                chip = pe.evaluate_chip({t_: S_eval[k] for k, t_ in enumerate(eval_tfs)}, genes, ct, gt, dcfg, root)
                save(root, ct, src, reg, m, S_eval, eval_tfs, chip, internal)


if __name__ == "__main__":
    main()
