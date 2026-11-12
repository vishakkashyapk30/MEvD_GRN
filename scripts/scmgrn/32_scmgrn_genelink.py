#!/usr/bin/env python
"""GENELink (Chen & Liu, Brief Bioinform 2022) on the scMultiomeGRN benchmark, run with the copy shipped in the
official scMultiomeGRN repository (src/benchmark/models/genelink.py, its default hyper-parameters) on our fold files.

The official benchmark Evaluator imports DeepWalk / GENIE3 wrappers that need gensim and pyscenic, so this driver
imports only GENELink and the benchmark GraphDataModule and repeats Evaluator.run_single_fold's steps:
seed, train on the fold's training pairs, score every TF pair (sigmoid, symmetrised as in GENELink.evaluate).
The data module's pair histograms (edge_attr) are replaced by zeros: GENELink never reads them, and they are the
slow part of loading. Scores go through scmgrn_protocol.evaluate_all like every other method.

Usage: python scripts/scmgrn/32_scmgrn_genelink.py --config configs/scmgrn/lung.yaml --root ~/scmgrn --cell_type Myeloid_cells \
           --folds 1-10 --seed 8 [--split_kind tfholdout]
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from src.benchmarks import scmgrn_protocol as proto   # noqa: E402
from src.utils.io import load_config, save_json       # noqa: E402

DRIVER = r"""
import json, sys, random, numpy as np, torch
kw = json.loads(sys.argv[1])
from src.benchmark.dataset import GraphDataModule
GraphDataModule.compute_edge_attr = lambda self, scrna_file, adj_norm, bins=16: np.zeros(
    (adj_norm._nnz() if hasattr(adj_norm, "_nnz") else adj_norm.nnz, 16, 16))
from src.benchmark.models.genelink import GENELink
import argparse
p = GENELink.add_argparse_args(argparse.ArgumentParser())
cfg = vars(p.parse_args([]))
cfg.update(kw)
cfg["seed"] = kw["model_seed"]
random.seed(cfg["seed"]); np.random.seed(cfg["seed"]); torch.manual_seed(cfg["seed"]); torch.cuda.manual_seed_all(cfg["seed"])
data = GraphDataModule(**cfg)
data.prepare_data()
model = GENELink(**cfg)                  # cfg carries root_dir
model.train(**data.train_data.data)
m, feat, adj = model.model, model.data_feature, model.adj
n = adj.shape[0]
pairs = torch.from_numpy(np.stack(np.meshgrid(np.arange(n), np.arange(n))).T.reshape(-1, 2)).to(feat.device)
m.eval()
with torch.no_grad():
    s = torch.sigmoid(m(feat, adj, pairs)).reshape(n, n)
s = ((s + s.T) / 2).cpu().numpy()
np.save(kw["out_scores"], s)
print("SCORES_OK", flush=True)
"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--root", default=os.environ.get("SCMGRN_ROOT", "data/scmgrn"))
    ap.add_argument("--cell_type", required=True)
    ap.add_argument("--folds", default="1-10")
    ap.add_argument("--seed", type=int, default=8, help="GENELink's default seed")
    ap.add_argument("--official_repo", default=None)
    ap.add_argument("--split_kind", default="official", choices=proto.SPLIT_KINDS)
    ap.add_argument("--work_dir", default=None)
    ap.add_argument("--skip_existing", action="store_true")
    ap.add_argument("--save_scores", action="store_true")
    args = ap.parse_args()

    dcfg = load_config(args.config)
    root = Path(args.root)
    ct = args.cell_type
    official = str(root / (args.official_repo or dcfg["official_repo"]))
    gdir = root / dcfg["graph_dir"] / ct
    sdir = proto.resolve_split_dir(str(root / dcfg["splits_dir"] / ct), dcfg["split"], args.split_kind)
    adj, _ = proto.load_adj(str(gdir / "x_graph.txt"))
    folds = (list(range(int(args.folds.split("-")[0]), int(args.folds.split("-")[1]) + 1))
             if "-" in args.folds else [int(x) for x in args.folds.split(",")])
    out_dir = root / proto.results_subdir(dcfg["results_dir"], args.split_kind) / ct / "genelink" / f"seed{args.seed}"
    out_dir.mkdir(parents=True, exist_ok=True)
    work = Path(args.work_dir or (out_dir / "work"))

    for fold in folds:
        out = out_dir / f"fold{fold:02d}.json"
        if args.skip_existing and out.exists():
            continue
        files = proto.fold_files(sdir, fold)
        wd = work / f"fold{fold:02d}"
        wd.mkdir(parents=True, exist_ok=True)
        kw = dict(name=ct, adj_file=str(gdir / "x_graph.txt"), node_file=str(gdir / "x_node.txt"),
                  scrna_file=str(gdir / "tmp" / "scrna.csv"),
                  atac_feature_file=str(gdir / "x_atac_rp_score_feature.txt"),
                  scrna_feature_file=str(gdir / "x_scRNA_grnboost2_feature.txt"),
                  model_cls="GENELink", mask_rate=0.0, data_dir=str(wd / "data"),
                  train_edge_file=files["train"], val_edge_file=files["val"], test_edge_file=files["test"],
                  split_id=fold, model_seed=args.seed, root_dir=str(wd / "logs"),
                  out_scores=str(wd / "scores.npy"))
        t0 = time.time()
        r = subprocess.run([sys.executable, "-u", "-c", DRIVER, json.dumps(kw)], cwd=str(wd), capture_output=True,
                           text=True, env={**os.environ, "PYTHONPATH": os.pathsep.join(
                               filter(None, [official, os.environ.get("PYTHONPATH")]))})
        (wd / "stdout.log").write_text(r.stdout[-200000:])
        (wd / "stderr.log").write_text(r.stderr[-200000:])
        if r.returncode != 0 or "SCORES_OK" not in r.stdout:
            raise SystemExit(f"GENELink failed (fold {fold}); see {wd/'stderr.log'}\n{r.stderr[-3000:]}")
        S = np.load(wd / "scores.npy")
        f = proto.load_fold(sdir, fold, adj.shape[0])
        m = proto.evaluate_all(S, adj, f, sdir, fold)
        if args.save_scores:
            np.save(out_dir / f"scores_fold{fold:02d}.npy", S.astype(np.float32))
        save_json({"method": "genelink", "model": "GENELink (official scMultiomeGRN benchmark copy)",
                   "split_kind": args.split_kind, "cell_type": ct, "dataset": dcfg["dataset"], "fold": fold,
                   "seed": args.seed, "train_seconds": time.time() - t0, "metrics": m}, out)
        print(f"[genelink {ct} fold {fold}] test AUROC={m['test']['auroc']:.4f} AUPR={m['test']['aupr']:.4f} "
              f"({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
