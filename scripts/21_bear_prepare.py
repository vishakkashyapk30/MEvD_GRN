#!/usr/bin/env python
"""BEAR-GRN step 1: label-free processed dir per dataset (+ optional FM embeddings).

  python scripts/21_bear_prepare.py --root $BEAR_ROOT --datasets K562 Macrophage_S1 [--fm]

Writes $BEAR_ROOT/processed/<dataset>/ (see src/benchmarks/bear_data.py).
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from src.benchmarks.bear_data import prepare_dataset   # noqa: E402
from src.utils.io import load_config                    # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.environ.get("BEAR_ROOT", "data/bear"))
    ap.add_argument("--config", default=str(REPO / "configs/bear/datasets.yaml"))
    ap.add_argument("--model_config", default=str(REPO / "configs/bear/mevd_base.yaml"),
                    help="only its data/atac blocks are used (graph degrees, QC switches)")
    ap.add_argument("--datasets", nargs="+", required=True)
    ap.add_argument("--fm", action="store_true", help="also extract Geneformer embeddings (scripts/11)")
    ap.add_argument("--gtf", default=None, help="override the GTF path")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--motif", action="store_true",
                    help="build motif_pair_feats.npz (s12.6 M3; needs --genome_store and --jaspar)")
    ap.add_argument("--genome_store", default=None, help="2-bit store dir from bear_motif.pack_genome")
    ap.add_argument("--jaspar", default=None, help="JASPAR 2024 CORE vertebrates (redundant) PFM file")
    ap.add_argument("--motif_device", default="cpu")
    ap.add_argument("--regime", default="L1", choices=["L1", "L2"],
                    help="L2 = compendium labels only -> processed/<ds>__L2 (s6.2)")
    args = ap.parse_args()
    root = Path(args.root)
    dcfg = load_config(args.config)
    mcfg = load_config(args.model_config)
    pcfg = dict(mcfg["data"])
    pcfg["atac"] = dict(mcfg.get("atac", {}))
    for ds in args.datasets:
        c = dcfg["datasets"][ds]
        out = root / "processed" / (ds + ("__L2" if args.regime == "L2" else ""))
        comp = dcfg["compendium"][c["species"]] if args.regime == "L2" else None
        if (out / "summary.json").exists() and not args.overwrite:
            print(f"[skip] {out} exists", flush=True)
        else:
            gtf = args.gtf or str(root / dcfg["gtf"][c["species"]])
            genome = "hg38" if c["species"] == "human" else "mm10"
            t0 = time.time()
            print(f"\n##### {ds} ({c['species']}, {genome}) #####", flush=True)
            s = prepare_dataset(c, root, out, pcfg, gtf, genome, compendium=comp)
            print(f"[done] {ds} in {time.time()-t0:.0f}s: {s}", flush=True)
        if args.motif and not (out / "motif_pair_feats.npz").exists():
            import json as _json
            import numpy as _np
            from src.benchmarks import bear_motif
            from src.utils.io import load_json
            t0 = time.time()
            gi = load_json(out / "gene_index.json")
            names = [None] * len(gi)
            for g, i in gi.items():
                names[i] = g
            tf_nodes = _np.asarray(load_json(out / "tf_indices.json"), dtype=_np.int64)
            F, info = bear_motif.pair_features(load_json(out / "gene_top_peaks.json"), tf_nodes, names,
                                               bear_motif.Genome2bit(args.genome_store), args.jaspar,
                                               len(names), device=args.motif_device)
            _np.savez_compressed(out / "motif_pair_feats.npz", F=F, tf_nodes=tf_nodes)
            info["seconds"] = round(time.time() - t0, 1)
            (out / "motif_info.json").write_text(_json.dumps(info, indent=1))
            print(f"[motif] {ds}: {info}", flush=True)
        if args.fm and not (out / "fm_gene_embeddings.npy").exists():
            subprocess.run([sys.executable, str(REPO / "scripts/11_extract_fm_embeddings.py"),
                            "--cell_type", ds, "--processed_dir", str(out)], check=True)


if __name__ == "__main__":
    main()
