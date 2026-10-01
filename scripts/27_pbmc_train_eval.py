#!/usr/bin/env python
"""Train + evaluate MeVD-GRN on the PBMC10k / LINGER-Cistrome benchmark:
one cell type x one label source x one split regime x one model seed x one variant.

Protocol (pre-registered, docs/experiments/pbmc10k_linger_benchmark.md s4):
  * labels = the build's split files (scripts/26); every Cistrome evaluation TF
    was removed as a regulator at build time (re-asserted here)
  * training negatives = the split's degree-matched train pool (or uniform, for
    the *_uniformneg variant); restrict_negative_pool stays on
  * single stage, model selection on val AUPR only; test + Cistrome are scored
    once, after training, from the restored best-val checkpoint
  * --rna_only: ATAC features + locus descriptor zeroed, use_atac=False, and the
    TF-candidate graph without the accessibility gate (the key ablation)

Outputs (under <root>/results/<cell_type>/<source>/<regime>/<variant>/seed<S>/):
  metrics.json      internal test metrics + Cistrome metrics + run info
  eval_scores.npz   float32 score rows for the evaluation TFs (re-scorable)

Usage:
  python scripts/27_pbmc_train_eval.py --config configs/pbmc/benchmark.yaml \
      --model_config configs/pbmc/mevd_fm_h384.yaml --root $PBMC_ROOT \
      --cell_type classical_monocyte --source collectri --regime tf --seed 42
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

from src.benchmarks import pbmc_eval as pe                  # noqa: E402
from src.benchmarks import pbmc_labels as pl                # noqa: E402
from src.data.dataset import load_celltype_data             # noqa: E402
from src.models.mevd_grn import MEvDGRN                     # noqa: E402
from src.training.curriculum import stages_from_config      # noqa: E402
from src.training.trainer import MEvDTrainer                # noqa: E402
from src.utils.io import load_config, load_json, save_json  # noqa: E402


def T(a: np.ndarray) -> torch.Tensor:
    a = np.asarray(a, dtype=np.int64).reshape(-1, 2)
    return torch.from_numpy(a.T.copy())


def build_model(m: dict, rna_only: bool) -> MEvDGRN:
    return MEvDGRN(rna_in_dim=m["rna_in_dim"], atac_in_dim=m["atac_in_dim"],
                   hidden_dim=m["hidden_dim"], n_gnn_layers=m["n_gnn_layers"], dropout=m["dropout"],
                   use_atac=not rna_only, integration=m.get("integration", "role_aware"),
                   graph_mode=m.get("graph_mode", "both"), use_edge_mlp=bool(m.get("use_edge_mlp", False)),
                   combine_mode=m.get("combine_mode", "sum"), use_fm=bool(m.get("use_fm", False)),
                   fm_in_dim=int(m.get("fm_in_dim", 768)), use_motif=False)


@torch.no_grad()
def score_rows(trainer: MEvDTrainer, regs: np.ndarray, n: int) -> np.ndarray:
    """sigmoid scores for every (reg, gene) pair: (len(regs), n) float32."""
    trainer.model.eval()
    emb = trainer._encode()
    tg = torch.arange(n, device=trainer.device)
    out = np.zeros((len(regs), n), dtype=np.float32)
    for k, r in enumerate(regs):
        src = torch.full((n,), int(r), dtype=torch.long, device=trainer.device)
        out[k] = torch.sigmoid(trainer._decode(emb, src, tg)).float().cpu().numpy()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(REPO / "configs/pbmc/benchmark.yaml"))
    ap.add_argument("--model_config", default=str(REPO / "configs/pbmc/mevd_fm_h384.yaml"))
    ap.add_argument("--root", default=os.environ.get("PBMC_ROOT", "data/pbmc"))
    ap.add_argument("--cell_type", required=True)
    ap.add_argument("--source", default="collectri", choices=pl.SOURCES)
    ap.add_argument("--regime", default="tf", choices=pl.REGIMES)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--rna_only", action="store_true")
    ap.add_argument("--variant", default=None)
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--max_epochs", type=int, default=None, help="smoke tests only")
    ap.add_argument("--tmp_dir", default=None)
    ap.add_argument("--skip_existing", action="store_true")
    args = ap.parse_args()
    if args.device.startswith("cuda") and not torch.cuda.is_available():
        print("[warn] CUDA unavailable -> CPU", flush=True)
        args.device = "cpu"
    if args.variant is None:
        args.variant = Path(args.model_config).stem.replace("mevd_", "") + ("_rnaonly" if args.rna_only else "")
    dcfg = load_config(args.config)
    root = Path(args.root)
    out_dir = (root / "results" / args.cell_type / args.source / args.regime / args.variant /
               f"seed{args.seed}")
    if args.skip_existing and (out_dir / "metrics.json").exists():
        print(f"[skip] {out_dir}", flush=True)
        return
    seed = args.seed
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed); torch.cuda.manual_seed_all(seed)

    cfg = load_config(args.model_config)
    cfg["data"]["seed"] = seed
    tier = cfg["data"]["evidence_tiers"][0]
    pdir = root / "processed" / args.cell_type
    sfile = root / "splits" / args.source / f"{args.regime}.npz"
    sp_ = np.load(sfile)
    arrays = {p: {k.split("__")[1]: sp_[k] for k in sp_.files if k.startswith(p + "__")}
              for p in ("train", "val", "test")}
    meta = load_json(root / "splits" / args.source / "meta.json")
    gene_index = load_json(pdir / "gene_index.json")
    n = len(gene_index)
    assert n == meta["n_genes"], "processed dir and split files use different gene universes"
    genes = [None] * n
    for g, i in gene_index.items():
        genes[i] = g
    gt = pe.load_groundtruth(dcfg, root)
    eval_tfs = sorted({d["tf"] for d in gt["datasets"]})
    eval_idx = np.array([gene_index[t] for t in eval_tfs if t in gene_index], dtype=np.int64)
    train_regs = set(arrays["train"]["pos"][:, 0].tolist())
    assert not (set(eval_idx.tolist()) & train_regs), "a Cistrome evaluation TF is a training regulator"

    data = load_celltype_data(str(pdir), [tier])
    if cfg["model"].get("use_fm") and data.fm_embeddings.shape[1] == 0:
        raise SystemExit(f"use_fm=true but {pdir}/fm_gene_embeddings.npy is missing (scripts/26 --steps fm)")
    if args.rna_only:
        data.atac_features = torch.zeros_like(data.atac_features)
        data.openness = torch.zeros_like(data.openness)
        data.tf_candidate_edges = torch.load(pdir / "tf_candidate_edges_rnaonly.pt")
    neg_mode = cfg.get("negatives", {}).get("mode", "degree_matched")
    if neg_mode == "uniform":
        pool = arrays["train"]["neg_uniform"]
    else:
        pool = arrays["train"]["neg"]
    data.negative_pool = T(pool)
    data.motif_edges = torch.zeros((2, 0), dtype=torch.long)
    splits = {p: {"pos": T(arrays[p]["pos"]), "neg": T(arrays[p]["neg"])} for p in ("train", "val", "test")}

    ckpt = Path(tempfile.mkdtemp(prefix=f"ckpt_{args.cell_type}_{args.regime}_", dir=args.tmp_dir))
    cfg["checkpoint_dir"] = str(ckpt)
    model = build_model(cfg["model"], args.rna_only)
    trainer = MEvDTrainer(model, data, cfg, args.device)
    stage = stages_from_config(cfg)[0]
    if args.max_epochs:
        stage.n_epochs = args.max_epochs
    removed = trainer.restrict_negative_pool({tier: splits})
    t0 = time.time()
    res = trainer.train_stage(stage, splits)
    secs = time.time() - t0

    # ---- internal test (label-set held-out part)
    test_regs = np.unique(arrays["test"]["pos"][:, 0])
    S_test = score_rows(trainer, test_regs, n)
    internal = pe.internal_test_metrics(S_test, test_regs, arrays, n, args.regime)
    # ---- Cistrome (LINGER protocol)
    S_eval = score_rows(trainer, eval_idx, n)
    chip = pe.evaluate_chip({genes[r]: S_eval[k] for k, r in enumerate(eval_idx)}, genes,
                            args.cell_type, gt, dcfg)
    shutil.rmtree(ckpt, ignore_errors=True)

    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_dir / "eval_scores.npz", scores=S_eval.astype(np.float32),
                        tfs=np.array([genes[r] for r in eval_idx]))
    rec = {"method": "MeVD-GRN", "variant": args.variant, "cell_type": args.cell_type,
           "source": args.source, "regime": args.regime, "seed": seed, "rna_only": bool(args.rna_only),
           "negatives": neg_mode, "model_config": str(args.model_config),
           "n_params": model.count_parameters(), "n_genes": n,
           "n_train_pos": int(len(arrays["train"]["pos"])), "neg_pool": int(len(pool)),
           "val_test_negatives_removed_from_pool": int(removed),
           "epochs_run": len(res["history"]), "best_val_aupr": res["best_val_aupr"],
           "train_seconds": secs, "internal_test": internal, "chip": chip}
    save_json(rec, out_dir / "metrics.json")
    s = chip.get("summary", {})
    print(f"[{args.cell_type} {args.source}/{args.regime} {args.variant} seed {seed}] "
          f"epochs={rec['epochs_run']} val_aupr={res['best_val_aupr']:.4f} | internal per-TF AUROC="
          f"{internal.get('per_tf_auroc_mean', float('nan')):.4f} dm-AUROC={internal.get('dm_auroc', float('nan')):.4f}"
          f" | Cistrome n={s.get('n', 0)} mean AUROC={s.get('auroc_mean', float('nan')):.4f} "
          f"AUPR ratio={s.get('aupr_ratio_mean', float('nan')):.3f} ({secs:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
