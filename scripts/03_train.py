#!/usr/bin/env python
"""Train MEvD-GRN with the multi-evidence curriculum (plan Step 11).

Usage:
  python scripts/03_train.py --config configs/k562.yaml --device cuda:0
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.data.dataset import load_celltype_data, load_splits
from src.models.mevd_grn import MEvDGRN
from src.training.curriculum import stages_from_config
from src.training.trainer import MEvDTrainer
from src.utils.io import ensure_dir, load_config, save_json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--cell_type", default=None)
    ap.add_argument("--device", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    cell_type = args.cell_type or cfg["cell_type"]
    device = args.device or cfg["hardware"]["device"]
    if device.startswith("cuda") and not torch.cuda.is_available():
        print("[warn] CUDA unavailable, falling back to CPU", flush=True)
        device = "cpu"
    cfg["checkpoint_dir"] = f"results/checkpoints/{cell_type}"

    tiers = cfg["data"]["evidence_tiers"]
    data = load_celltype_data(cfg["paths"]["processed_dir"], tiers)
    splits_per_tier = {t: load_splits("data/splits", cell_type, t)
                       for t in tiers
                       if (Path("data/splits") / f"{cell_type}_{t}_splits.pt").exists()}

    mcfg = cfg["model"]
    model = MEvDGRN(rna_in_dim=mcfg["rna_in_dim"], atac_in_dim=mcfg["atac_in_dim"],
                    hidden_dim=mcfg["hidden_dim"], n_gnn_layers=mcfg["n_gnn_layers"],
                    dropout=mcfg["dropout"], integration=mcfg.get("integration", "role_aware"))
    print(f"[model] MEvD-GRN parameters: {model.count_parameters():,}", flush=True)

    trainer = MEvDTrainer(model, data, cfg, device)
    stages = stages_from_config(cfg)

    t0 = time.time()
    curr_results = trainer.run_full_curriculum(splits_per_tier, stages)
    train_secs = time.time() - t0

    # final per-tier test evaluation
    print("\n=== Final test evaluation (all tiers) ===", flush=True)
    final = {}
    for tier, sp in splits_per_tier.items():
        m = trainer.evaluate_split(sp["test"])
        final[tier] = m
        print(f"  {tier:14s} AUPR={m['aupr']:.4f} AUROC={m['auroc']:.4f} "
              f"EP={m['early_precision']:.4f} EPR={m['epr']:.2f}", flush=True)

    ensure_dir("results")
    out = {"cell_type": cell_type, "train_seconds": train_secs,
           "n_params": model.count_parameters(),
           "curriculum": {k: v["best_val_aupr"] for k, v in curr_results.items()},
           "test": final}
    save_json(out, f"results/{cell_type}_results.json")
    torch.save({"model": model.state_dict()}, f"results/checkpoints/{cell_type}/final_model.pt")
    print(f"\n[done] trained in {train_secs/60:.1f} min -> results/{cell_type}_results.json", flush=True)


if __name__ == "__main__":
    main()
