#!/usr/bin/env python
"""Reproduction anchor: run the OFFICIAL scMultiomeGRN code (Zenodo
10.5281/zenodo.14848389, unmodified, `Trainer.run_single_fold`) on our
regenerated ground truth + official-format features + the same split files
MeVD-GRN uses. If these numbers land near the paper's, the regenerated
benchmark is faithful and MeVD-GRN's numbers are directly comparable.

The official code is executed in a subprocess (both repos have a top-level
`src` package) with cwd = a per-fold work dir and PYTHONPATH = the official
repo. It needs `lightning` (or `pytorch_lightning`) in the environment.
Metrics are read from its CSV logs (`test/best_test/*`, `test/best_all/*`)
AND recomputed from its saved best_all score matrix with our
src/benchmarks/scmgrn_protocol.py, so both methods go through one metric code.

Usage:
  python scripts/scmgrn/18_scmgrn_official_repro.py --config configs/scmgrn/lung.yaml \
      --root $SCMGRN_ROOT --cell_type Lymphoid_cells --folds 1-10
"""
from __future__ import annotations

import argparse
import glob
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from src.benchmarks import scmgrn_protocol as proto   # noqa: E402
from src.utils.io import load_config, save_json       # noqa: E402

DRIVER = r"""
import inspect, json, sys, torch
# compatibility shim only: torch>=2.7 removed ReduceLROnPlateau(verbose=...), a
# print-only flag the official code passes (it targets torch 2.4). No-op on older torch.
_R = torch.optim.lr_scheduler.ReduceLROnPlateau
if "verbose" not in inspect.signature(_R.__init__).parameters:
    class _RP(_R):
        def __init__(self, *a, verbose=None, **k):
            super().__init__(*a, **k)
    torch.optim.lr_scheduler.ReduceLROnPlateau = _RP
from src.trainer import Trainer
kw = json.loads(sys.argv[1])
parser = Trainer.add_argparse_args(**kw)
config = parser.parse_args([])
config.accelerator = "gpu" if torch.cuda.is_available() else "cpu"
print("LOGDIR=" + Trainer.run_single_fold(**vars(config)), flush=True)
"""


def last_metrics(csv_path: str) -> dict:
    df = pd.read_csv(csv_path)
    out = {}
    for c in df.columns:
        v = df[c].dropna()
        if len(v) and c.startswith("test/"):
            out[c] = float(v.iloc[-1])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--root", default=os.environ.get("SCMGRN_ROOT", "data/scmgrn"))
    ap.add_argument("--cell_type", required=True)
    ap.add_argument("--folds", default="1-10")
    ap.add_argument("--model_seed", type=int, default=666, help="official default")
    ap.add_argument("--max_epochs", type=int, default=2000)
    ap.add_argument("--patience", type=int, default=100)
    ap.add_argument("--work_dir", default=None, help="scratch for lightning logs")
    ap.add_argument("--skip_existing", action="store_true")
    ap.add_argument("--official_repo", default=None, help="code dir relative to --root (default: config official_repo); "
                    "ablations use the patched copy from scripts/scmgrn/23_patch_official_ablation.py")
    ap.add_argument("--ablation", default=None, help="tag: results go to official_scmultiomegrn_<tag>/")
    ap.add_argument("--method_name", default=None, help="results directory name (overrides the ablation-based name; "
                    "used for other models of the official repo, e.g. --extra_kw '{\"model_cls\": \"DeepTFni\"}')")
    ap.add_argument("--extra_kw", default="{}", help="JSON of extra official kwargs (scripts/23 switches); "
                    "'{graph_dir}' in strings -> this cell type's graph dir; gf_from_processed=true writes the "
                    "Geneformer TSV (x_node.txt order) and sets gf_feature_file")
    ap.add_argument("--save_scores", action="store_true", help="keep the best_all score matrix (scores_foldNN.npy)")
    ap.add_argument("--split_kind", default="official", choices=proto.SPLIT_KINDS)
    args = ap.parse_args()

    dcfg = load_config(args.config)
    root = Path(args.root)
    ct = args.cell_type
    official = str(root / (args.official_repo or dcfg["official_repo"]))
    gdir = root / dcfg["graph_dir"] / ct
    s = dcfg["split"]
    sdir = proto.resolve_split_dir(str(root / dcfg["splits_dir"] / ct), s, args.split_kind)
    adj, _ = proto.load_adj(str(gdir / "x_graph.txt"))
    folds = (list(range(int(args.folds.split("-")[0]), int(args.folds.split("-")[1]) + 1))
             if "-" in args.folds else [int(x) for x in args.folds.split(",")])
    import json
    method = args.method_name or ("official_scmultiomegrn" + (f"_{args.ablation}" if args.ablation else ""))
    out_dir = root / proto.results_subdir(dcfg["results_dir"], args.split_kind) / ct / method / f"seed{args.model_seed}"
    out_dir.mkdir(parents=True, exist_ok=True)
    work = Path(args.work_dir or (out_dir / "work"))
    extra = {k: (v.replace("{graph_dir}", str(gdir)) if isinstance(v, str) else v)
             for k, v in json.loads(args.extra_kw).items()}
    if extra.pop("gf_from_processed", False):
        pdir = root / dcfg["processed_dir"] / ct
        gi = json.load(open(pdir / "gene_index.json"))
        names = pd.read_csv(gdir / "x_node.txt", names=["name"])["name"].tolist()
        E = np.load(pdir / "fm_gene_embeddings.npy")
        work.mkdir(parents=True, exist_ok=True)
        gf = work / "gf_feature.tsv"
        pd.DataFrame(E[[gi[n.upper()] for n in names]], index=names).to_csv(gf, sep="\t")
        extra["gf_feature_file"] = str(gf)
    if extra:
        print(f"[official {ct}] ablation={args.ablation} extra kwargs: {extra}", flush=True)

    for fold in folds:
        out = out_dir / f"fold{fold:02d}.json"
        if args.skip_existing and out.exists():
            continue
        files = proto.fold_files(sdir, fold)
        wd = work / f"fold{fold:02d}"
        wd.mkdir(parents=True, exist_ok=True)
        kw = dict(name=ct, adj_file=str(gdir / "x_graph.txt"), node_file=str(gdir / "x_node.txt"),
                  scrna_file=str(gdir / "tmp" / "scrna.csv"), comment="mevd_repro" + (f"_{args.ablation}" if args.ablation else "") + ("" if args.split_kind == "official" else f"_{args.split_kind}"),
                  atac_feature_file=str(gdir / "x_atac_rp_score_feature.txt"),
                  scrna_feature_file=str(gdir / "x_scRNA_grnboost2_feature.txt"),
                  model_cls="scMultiomeGRN", mask_rate=0.0, data_dir=str(wd / "data"),
                  train_edge_file=files["train"], val_edge_file=files["val"],
                  test_edge_file=files["test"], split_id=fold, max_epochs=args.max_epochs,
                  patience=args.patience, model_seed=args.model_seed)
        kw.update(extra)                       # extra kwargs may override defaults (e.g. model_cls)
        t0 = time.time()
        r = subprocess.run([sys.executable, "-u", "-c", DRIVER, json.dumps(kw)], cwd=str(wd),
                           env={**os.environ, "PYTHONPATH": os.pathsep.join(filter(None, [official, os.environ.get("PYTHONPATH")]))}, capture_output=True, text=True)
        secs = time.time() - t0
        (wd / "stdout.log").write_text(r.stdout[-200000:])
        (wd / "stderr.log").write_text(r.stderr[-200000:])
        if r.returncode != 0:
            raise SystemExit(f"official run failed (fold {fold}); see {wd/'stderr.log'}\n{r.stderr[-3000:]}")
        log_dir = [l for l in r.stdout.splitlines() if l.startswith("LOGDIR=")][-1][len("LOGDIR="):]
        log_dir = str(wd / log_dir) if not os.path.isabs(log_dir) else log_dir
        logged = last_metrics(os.path.join(log_dir, "metrics.csv"))
        S = np.loadtxt(os.path.join(log_dir, "best_all_score_matrix.txt"))
        f = proto.load_fold(sdir, fold, adj.shape[0])
        m = proto.evaluate_all(S, adj, f, sdir, fold)
        if args.save_scores:
            np.save(out_dir / f"scores_fold{fold:02d}.npy", S.astype(np.float32))
        for ck in glob.glob(os.path.join(log_dir, "checkpoints")):
            shutil.rmtree(ck, ignore_errors=True)
        rec = {"method": method, "model_cls": kw["model_cls"], "ablation": args.ablation, "split_kind": args.split_kind, "extra_kw": extra, "cell_type": ct, "dataset": dcfg["dataset"],
               "fold": fold, "seed": args.model_seed, "train_seconds": secs, "log_dir": log_dir,
               "metrics": m, "official_logged_metrics": logged}
        save_json(rec, out)
        print(f"[official {ct} fold {fold}] test AUROC={m['test']['auroc']:.4f} AUPR={m['test']['aupr']:.4f} "
              f"| all AUROC={m['all']['auroc']:.4f} AUPR={m['all']['aupr']:.4f} "
              f"(logged best_all AUROC={logged.get('test/best_all/auroc', float('nan')):.4f}) "
              f"({secs:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
