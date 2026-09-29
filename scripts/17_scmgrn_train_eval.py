#!/usr/bin/env python
"""Train + evaluate MeVD-GRN on scMultiomeGRN's benchmark, one cell type x
one model seed x a set of folds, with the official protocol:

  * the fold's official split files (scripts/16), both directions of every pair
  * training negatives come ONLY from the fold's train negatives (the processed
    dir's global negative_pool is replaced per fold, so val/test negatives are
    never seen in training)
  * single evidence tier, one training stage, model selection on the fold's
    val split only (val AUPR, MeVD-GRN's native rule), <= 2000 epochs,
    patience 100 epochs (configs/scmgrn/mevd_base.yaml)
  * scores symmetrised, S = (sigmoid(i->j) + sigmoid(j->i)) / 2, because the
    ground truth is an undirected TF-TF graph
  * metrics = the official metric_fn on the "test" (held-out fold, 1:1) and
    "all" (full N x N matrix) sets -- see src/benchmarks/scmgrn_protocol.py

--obs_graph: information parity with scMultiomeGRN, whose encoder message-
passes over the fold's TRAINING adjacency. The fold's train positives are fed
as MeVD-GRN's 3rd relation (the slot the motif graph uses on SC-MO-GRN-DB;
that motif graph is disabled here because the ground truth is a motif scan).

Usage:
  python scripts/17_scmgrn_train_eval.py --config configs/scmgrn/lung.yaml \
     --model_config configs/scmgrn/mevd_fm_h384.yaml --obs_graph \
     --root $SCMGRN_ROOT --cell_type Lymphoid_cells --seed 42 --device cuda:0
"""
from __future__ import annotations

import argparse
import os
import random
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import torch

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from src.benchmarks import scmgrn_protocol as proto        # noqa: E402
from src.data.dataset import load_celltype_data            # noqa: E402
from src.models.mevd_grn import MEvDGRN                    # noqa: E402
from src.training.curriculum import stages_from_config     # noqa: E402
from src.training.trainer import MEvDTrainer               # noqa: E402
from src.utils.io import load_config, save_json            # noqa: E402


def T(a: np.ndarray) -> torch.Tensor:
    a = np.asarray(a, dtype=np.int64).reshape(-1, 2)
    return torch.from_numpy(a.T.copy())


def build_model(mcfg: dict, use_motif: bool) -> MEvDGRN:
    return MEvDGRN(rna_in_dim=mcfg["rna_in_dim"], atac_in_dim=mcfg["atac_in_dim"],
                   hidden_dim=mcfg["hidden_dim"], n_gnn_layers=mcfg["n_gnn_layers"],
                   dropout=mcfg["dropout"], integration=mcfg.get("integration", "role_aware"),
                   graph_mode=mcfg.get("graph_mode", "both"),
                   use_edge_mlp=bool(mcfg.get("use_edge_mlp", False)),
                   combine_mode=mcfg.get("combine_mode", "sum"),
                   use_fm=bool(mcfg.get("use_fm", False)), fm_in_dim=int(mcfg.get("fm_in_dim", 768)),
                   use_motif=use_motif)


@torch.no_grad()
def score_matrix(trainer: MEvDTrainer, n: int, chunk: int = 200_000) -> np.ndarray:
    trainer.model.eval()
    emb = trainer._encode()
    ii, jj = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    src = torch.from_numpy(ii.reshape(-1)).to(trainer.device)
    dst = torch.from_numpy(jj.reshape(-1)).to(trainer.device)
    out = []
    for s in range(0, src.shape[0], chunk):
        out.append(torch.sigmoid(trainer._decode(emb, src[s:s + chunk], dst[s:s + chunk])).cpu())
    P = torch.cat(out).numpy().reshape(n, n).astype(np.float64)
    return (P + P.T) / 2.0


def run_fold(args, dcfg, fold: int, out_path: Path) -> dict:
    seed = args.seed
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed); torch.cuda.manual_seed_all(seed)
    cfg = load_config(args.model_config)
    cfg["data"]["seed"] = seed
    tier = cfg["data"]["evidence_tiers"][0]
    root = Path(args.root)
    pdir = root / dcfg["processed_dir"] / args.cell_type
    gdir = root / dcfg["graph_dir"] / args.cell_type
    s = dcfg["split"]
    sdir = proto.split_dir(str(root / dcfg["splits_dir"] / args.cell_type), int(s["seed"]),
                           int(s["n_splits"]), float(s["train_split_rate"]), bool(s["rest_all_train"]))
    adj, nodes = proto.load_adj(str(gdir / "x_graph.txt"))
    n = adj.shape[0]

    data = load_celltype_data(str(pdir), [tier])
    if cfg["model"].get("use_fm") and data.fm_embeddings.shape[1] == 0:
        raise SystemExit(f"use_fm=true but {pdir}/fm_gene_embeddings.npy is missing "
                         "(run scripts/16 --steps fm on the login node)")
    f = proto.load_fold(sdir, fold, n)
    data.negative_pool = T(f["train"]["neg"])                     # fold-train negatives only
    data.motif_edges = T(f["train"]["pos"]) if args.obs_graph else torch.zeros((2, 0), dtype=torch.long)
    splits = {p: {"pos": T(f[p]["pos"]), "neg": T(f[p]["neg"])} for p in ("train", "val", "test")}
    for p in ("val", "test"):                                    # sanity: nothing held out is trained on
        held = set(map(tuple, f[p]["pos"].tolist())) | set(map(tuple, f[p]["neg"].tolist()))
        trained = set(map(tuple, f["train"]["pos"].tolist())) | set(map(tuple, f["train"]["neg"].tolist()))
        assert not (held & trained), f"fold {fold}: {p} pairs overlap the train split"

    ckpt = Path(tempfile.mkdtemp(prefix=f"ckpt_{args.cell_type}_f{fold}_", dir=args.tmp_dir))
    cfg["checkpoint_dir"] = str(ckpt)
    model = build_model(cfg["model"], use_motif=args.obs_graph)
    trainer = MEvDTrainer(model, data, cfg, args.device)
    stage = stages_from_config(cfg)[0]
    if args.max_epochs:
        stage.n_epochs = args.max_epochs
    t0 = time.time()
    trainer.restrict_negative_pool({tier: splits})   # no-op: pool is already fold-train only
    res = trainer.train_stage(stage, splits)
    secs = time.time() - t0
    S = score_matrix(trainer, n)
    m = proto.evaluate_score_matrix(S, adj, f["test_eval"])
    shutil.rmtree(ckpt, ignore_errors=True)
    rec = {"method": "MeVD-GRN", "variant": args.variant, "cell_type": args.cell_type,
           "dataset": dcfg["dataset"], "fold": fold, "seed": seed, "obs_graph": bool(args.obs_graph),
           "model_config": str(args.model_config), "n_params": model.count_parameters(),
           "n_nodes": n, "n_train_pos_directed": int(splits["train"]["pos"].shape[1]),
           "epochs_run": len(res["history"]), "best_val_aupr": res["best_val_aupr"],
           "train_seconds": secs, "metrics": m}
    save_json(rec, out_path)
    print(f"[{args.cell_type} fold {fold} seed {seed} {args.variant}] epochs={rec['epochs_run']} "
          f"val={res['best_val_aupr']:.4f} | test AUROC={m['test']['auroc']:.4f} "
          f"AUPR={m['test']['aupr']:.4f} acc={m['test']['acc']:.4f} | all AUROC={m['all']['auroc']:.4f} "
          f"AUPR={m['all']['aupr']:.4f} acc={m['all']['acc']:.4f} ({secs:.0f}s)", flush=True)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="dataset config, e.g. configs/scmgrn/lung.yaml")
    ap.add_argument("--model_config", default=str(REPO / "configs/scmgrn/mevd_fm_h384.yaml"))
    ap.add_argument("--root", default=os.environ.get("SCMGRN_ROOT", "data/scmgrn"))
    ap.add_argument("--cell_type", required=True)
    ap.add_argument("--folds", default="1-10", help="e.g. 1-10 or 1,3,5")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--obs_graph", action="store_true")
    ap.add_argument("--variant", default=None, help="name used in output paths")
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--max_epochs", type=int, default=None, help="smoke tests only")
    ap.add_argument("--tmp_dir", default=None)
    ap.add_argument("--skip_existing", action="store_true")
    args = ap.parse_args()
    if args.device.startswith("cuda") and not torch.cuda.is_available():
        print("[warn] CUDA unavailable -> CPU", flush=True)
        args.device = "cpu"
    if args.variant is None:
        args.variant = Path(args.model_config).stem.replace("mevd_", "") + ("_obs" if args.obs_graph else "_noobs")
    dcfg = load_config(args.config)
    if "-" in args.folds:
        a, b = args.folds.split("-")
        folds = list(range(int(a), int(b) + 1))
    else:
        folds = [int(x) for x in args.folds.split(",")]
    out_dir = Path(args.root) / dcfg["results_dir"] / args.cell_type / f"mevd_{args.variant}" / f"seed{args.seed}"
    out_dir.mkdir(parents=True, exist_ok=True)
    for fold in folds:
        out = out_dir / f"fold{fold:02d}.json"
        if args.skip_existing and out.exists():
            print(f"[skip] {out}", flush=True)
            continue
        run_fold(args, dcfg, fold, out)


if __name__ == "__main__":
    main()
