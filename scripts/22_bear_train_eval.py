#!/usr/bin/env python
"""BEAR-GRN step 2: MeVD-GRN, TF-disjoint cross-fitted, -> one GRN edge list
per (dataset, variant, seed) that BEAR's scoring code can read unchanged.

Protocol (pre-registered, docs/experiments/bear_grn_benchmark.md s6):
  * TF nodes are split into K degree-stratified folds (seeded). For fold k every
    label whose source is a fold-k TF is dropped from EVERY label tier; 15% of
    the remaining TFs are an inner validation set (early stopping only).
  * Training = one all_at_once stage on the union of the training TFs' labels;
    negatives = degree-matched (tf ~ out-degree, gene ~ in-degree) non-positive
    pairs of training TFs, resampled from a fixed per-fold pool every epoch.
  * The fold-k model scores every (fold-k TF, gene) pair. The K blocks are
    concatenated into <out>/grn.tsv.gz (Source, Target, Score, fold).
  * Hard audits: no training/val label or negative has a held-out source TF;
    restrict_negative_pool (leak fix) is called; graphs are label-free (scripts/21).

  python scripts/22_bear_train_eval.py --root $BEAR_ROOT --dataset K562 \
      --model_config configs/bear/mevd_fm_h384.yaml --seed 42 --device cuda:0
"""
from __future__ import annotations

import argparse
import json
import os
import random
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from src.benchmarks import bear_protocol as bp             # noqa: E402
from src.data.dataset import load_celltype_data            # noqa: E402
from src.models.mevd_grn import MEvDGRN                    # noqa: E402
from src.training.curriculum import stages_from_config     # noqa: E402
from src.training.trainer import MEvDTrainer               # noqa: E402
from src.utils.io import load_config, load_json, save_json  # noqa: E402


def build_model(m: dict) -> MEvDGRN:
    return MEvDGRN(rna_in_dim=m["rna_in_dim"], atac_in_dim=m["atac_in_dim"],
                   hidden_dim=m["hidden_dim"], n_gnn_layers=m["n_gnn_layers"],
                   dropout=m["dropout"], integration=m.get("integration", "role_aware"),
                   graph_mode=m.get("graph_mode", "both"),
                   use_edge_mlp=bool(m.get("use_edge_mlp", False)),
                   combine_mode=m.get("combine_mode", "sum"),
                   use_fm=bool(m.get("use_fm", False)), fm_in_dim=int(m.get("fm_in_dim", 768)),
                   use_motif=False)


def load_edges(p: Path) -> np.ndarray:
    return torch.load(p).numpy().astype(np.int64).reshape(2, -1)


def T(e: np.ndarray) -> torch.Tensor:
    return torch.from_numpy(np.ascontiguousarray(e, dtype=np.int64))


@torch.no_grad()
def score_block(trainer: MEvDTrainer, tfs: np.ndarray, n_genes: int, chunk: int = 262_144) -> np.ndarray:
    trainer.model.eval()
    emb = trainer._encode()
    pairs = bp.candidate_block(tfs, np.arange(n_genes))
    out = np.empty(pairs.shape[1], dtype=np.float32)
    for s in range(0, pairs.shape[1], chunk):
        src = torch.from_numpy(pairs[0, s:s + chunk]).to(trainer.device)
        dst = torch.from_numpy(pairs[1, s:s + chunk]).to(trainer.device)
        out[s:s + chunk] = torch.sigmoid(trainer._decode(emb, src, dst)).float().cpu().numpy()
    return pairs, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.environ.get("BEAR_ROOT", "data/bear"))
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--model_config", default=str(REPO / "configs/bear/mevd_fm_h384.yaml"))
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--variant", default=None)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--folds", default=None, help="subset, e.g. 0,1 (smoke tests); default all")
    ap.add_argument("--max_epochs", type=int, default=None, help="smoke tests only")
    ap.add_argument("--tmp_dir", default=None)
    ap.add_argument("--skip_existing", action="store_true")
    args = ap.parse_args()
    if args.device.startswith("cuda") and not torch.cuda.is_available():
        print("[warn] CUDA unavailable -> CPU", flush=True)
        args.device = "cpu"
    variant = args.variant or Path(args.model_config).stem.replace("mevd_", "")
    root = Path(args.root)
    cfg0 = load_config(args.model_config)
    regime = cfg0.get("bear", {}).get("regime", "L1")
    pdir = root / "processed" / (args.dataset + ("__L2" if regime == "L2" else ""))
    out_dir = root / "results" / args.dataset / f"mevd_{variant}" / f"seed{args.seed}"
    if args.skip_existing and (out_dir / "grn.tsv.gz").exists():
        print(f"[skip] {out_dir}/grn.tsv.gz", flush=True)
        return
    out_dir.mkdir(parents=True, exist_ok=True)

    cfg = load_config(args.model_config)
    bcfg = cfg["bear"]
    seed = args.seed
    cfg["data"]["seed"] = seed
    gene_index = load_json(pdir / "gene_index.json")
    names = [None] * len(gene_index)
    for g, i in gene_index.items():
        names[i] = g
    n_genes = len(names)
    tf_nodes = np.asarray(load_json(pdir / "tf_indices.json"), dtype=np.int64)
    labels = {p.stem: load_edges(p) for p in sorted((pdir / "labels").glob("*.pt"))}
    gts = {p.stem: load_edges(p) for p in sorted((pdir / "gts").glob("*.pt"))}
    all_lab = np.unique(np.concatenate([bp.edge_keys(e, n_genes) for e in labels.values()]))
    all_known = np.unique(np.concatenate([all_lab] + [bp.edge_keys(e, n_genes) for e in gts.values()]))
    lab_union = bp.keys_to_edges(all_lab, n_genes)
    outdeg = dict(zip(*np.unique(lab_union[0], return_counts=True)))
    outdeg = {int(k): int(v) for k, v in outdeg.items()}

    if regime == "L2":
        # one model scores every GT TF; it never sees a GT TF as a training source
        gt_tfs = np.unique(np.concatenate([e[0] for e in gts.values()]))
        folds = [np.intersect1d(tf_nodes, gt_tfs)]
    else:
        folds = bp.tf_disjoint_folds(tf_nodes, outdeg, int(bcfg["n_folds"]), seed)
    run_folds = range(len(folds)) if args.folds is None else [int(x) for x in args.folds.split(",")]
    print(f"[{args.dataset}] {n_genes} genes, {len(tf_nodes)} TF nodes, tiers={list(labels)} "
          f"({', '.join(f'{k}:{v.shape[1]}' for k, v in labels.items())}), "
          f"{len(folds)} TF-disjoint folds of sizes {[len(f) for f in folds]}", flush=True)

    rng_master = np.random.default_rng(seed)
    blocks, fold_log = [], []
    t_all = time.time()
    for k in run_folds:
        random.seed(seed + k); np.random.seed(seed + k)
        torch.manual_seed(seed * 100 + k); torch.cuda.manual_seed_all(seed * 100 + k)
        rng = np.random.default_rng(rng_master.integers(1 << 31))
        held = folds[k]
        train_tfs = np.setdiff1d(tf_nodes, held)
        inner = bp.tf_disjoint_folds(train_tfs, outdeg, max(2, int(round(1 / float(bcfg["val_tf_frac"])))),
                                     seed + 1000 + k)
        val_tfs = inner[0]
        fit_tfs = np.setdiff1d(train_tfs, val_tfs)

        train_pos = bp.keep_tfs(lab_union, fit_tfs)
        val_pos = bp.keep_tfs(lab_union, val_tfs)
        stage = stages_from_config(cfg)[0]
        if args.max_epochs:
            stage.n_epochs = args.max_epochs
        n_avail = len(fit_tfs) * (n_genes - 1) - train_pos.shape[1]
        n_pool = int(min(float(bcfg["neg_pool_mult"]) * stage.neg_ratio * train_pos.shape[1],
                         float(bcfg.get("neg_pool_max_frac", 0.5)) * n_avail))
        if bcfg["negatives"] == "degree_matched":
            pool = bp.degree_matched_negatives(train_pos, n_pool, n_genes, all_known, rng)
            vneg = bp.degree_matched_negatives(val_pos, int(bcfg["val_neg_ratio"]) * val_pos.shape[1],
                                               n_genes, all_known, rng)
        else:
            pool = bp.uniform_negatives(fit_tfs, n_pool, n_genes, all_known, rng)
            vneg = bp.uniform_negatives(val_tfs, int(bcfg["val_neg_ratio"]) * val_pos.shape[1],
                                        n_genes, all_known, rng)
        bp.assert_tf_disjoint({"train_pos": train_pos, "train_neg_pool": pool, "val_pos": val_pos,
                               "val_neg": vneg}, held, what=f"fold {k}")
        assert not np.isin(bp.edge_keys(pool, n_genes), all_known).any(), "negative pool holds a known positive"

        data = load_celltype_data(str(pdir), [])
        if cfg["model"].get("use_fm") and data.fm_embeddings.shape[1] == 0:
            raise SystemExit(f"use_fm=true but {pdir}/fm_gene_embeddings.npy is missing (scripts/21 --fm)")
        if bcfg.get("zero_atac"):
            data.atac_features = torch.zeros_like(data.atac_features)
            data.openness = torch.zeros_like(data.openness)
        data.evidence = {"all": T(train_pos)}
        data.negative_pool = T(pool)
        data.motif_edges = torch.zeros((2, 0), dtype=torch.long)
        splits = {"train": {"pos": T(train_pos), "neg": T(pool)},
                  "val": {"pos": T(val_pos), "neg": T(vneg)},
                  "test": {"pos": torch.zeros((2, 0), dtype=torch.long),
                           "neg": torch.zeros((2, 0), dtype=torch.long)}}
        if args.tmp_dir:
            Path(args.tmp_dir).mkdir(parents=True, exist_ok=True)
        tmp = Path(tempfile.mkdtemp(prefix=f"bear_{args.dataset}_s{seed}_f{k}_", dir=args.tmp_dir))
        cfg["checkpoint_dir"] = str(tmp)
        model = build_model(cfg["model"])
        trainer = MEvDTrainer(model, data, cfg, args.device)
        trainer.restrict_negative_pool({"all": splits})
        t0 = time.time()
        res = trainer.train_stage(stage, splits)
        pairs, sc = score_block(trainer, held, n_genes)
        shutil.rmtree(tmp, ignore_errors=True)
        blocks.append((pairs, sc, k))
        rec = {"fold": k, "n_held_tfs": int(len(held)), "n_fit_tfs": int(len(fit_tfs)),
               "n_val_tfs": int(len(val_tfs)), "n_train_pos": int(train_pos.shape[1]),
               "n_neg_pool": int(pool.shape[1]), "n_val_pos": int(val_pos.shape[1]),
               "n_val_neg": int(vneg.shape[1]), "best_val_aupr": float(res["best_val_aupr"]),
               "epochs_run": len(res["history"]), "seconds": time.time() - t0,
               "n_params": int(model.count_parameters())}
        fold_log.append(rec)
        print(f"[{args.dataset} seed {seed} {variant}] fold {k}: {rec}", flush=True)
        del trainer, model, data
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    src = np.concatenate([b[0][0] for b in blocks])
    dst = np.concatenate([b[0][1] for b in blocks])
    score = np.concatenate([b[1] for b in blocks])
    fold = np.concatenate([np.full(b[1].shape[0], b[2], dtype=np.int16) for b in blocks])
    nm = np.asarray(names, dtype=object)
    df = pd.DataFrame({"Source": nm[src], "Target": nm[dst], "Score": score, "fold": fold})
    df.to_csv(out_dir / "grn.tsv.gz", sep="\t", index=False, float_format="%.6g")
    meta = {"dataset": args.dataset, "variant": variant, "seed": seed,
            "model_config": str(args.model_config), "regime": regime,
            "protocol": "TF-disjoint cross-fit (s6.1)" if regime == "L1" else "compendium, all GT TFs held out (s6.2)",
            "negatives": bcfg["negatives"], "zero_atac": bool(bcfg.get("zero_atac", False)),
            "folds": [f.tolist() for f in folds], "fold_tf_names": [[names[i] for i in f] for f in folds],
            "fold_log": fold_log, "n_edges": int(len(df)), "seconds": time.time() - t_all,
            "complete": args.folds is None and args.max_epochs is None}
    with open(out_dir / "meta.json", "w") as f:
        json.dump(meta, f, indent=1)
    print(f"[done] {out_dir}/grn.tsv.gz: {len(df):,} edges in {time.time()-t_all:.0f}s", flush=True)


if __name__ == "__main__":
    main()
