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
from src.evaluation.benchmarker import eval_edges_for_tier
from src.models.mevd_grn import MEvDGRN
from src.training.curriculum import build_all_at_once_stage, stages_from_config
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
    model = MEvDGRN(
        rna_in_dim=mcfg["rna_in_dim"], atac_in_dim=mcfg["atac_in_dim"],
        hidden_dim=mcfg["hidden_dim"], n_gnn_layers=mcfg["n_gnn_layers"],
        dropout=mcfg["dropout"], integration=mcfg.get("integration", "role_aware"),
        graph_mode=mcfg.get("graph_mode", "both"),
        use_edge_mlp=bool(mcfg.get("use_edge_mlp", False)),
        combine_mode=mcfg.get("combine_mode", "sum"),
        use_fm=bool(mcfg.get("use_fm", False)), fm_in_dim=int(mcfg.get("fm_in_dim", 768)),
    )
    print(f"[model] MEvD-GRN parameters: {model.count_parameters():,} "
          f"(graph_mode={model.graph_mode}, edge_mlp={model.use_edge_mlp})", flush=True)

    trainer = MEvDTrainer(model, data, cfg, device)
    # Train only on the "main curriculum" tiers (default: localization,
    # perturbation). Tiers outside this set (default: dual_evidence) are
    # deliberately never trained on -- they are heavily nested inside the
    # trained tiers (see configs/default.yaml), so held-out performance on
    # them is reported as zero-shot INFERENCE, not fine-tuned performance.
    main_tiers = set(cfg["curriculum"].get("main_curriculum_tiers", tiers))
    held_out_tiers = [t for t in splits_per_tier if t not in main_tiers]

    # Default protocol is "sequential" (localization then perturbation).
    # "all_at_once" (joint training on the union of both tiers at once)
    # structurally avoids catastrophic forgetting, but on a full-budget K562
    # run it measurably UNDERPERFORMS sequential on perturbation and on the
    # headline zero-shot dual-evidence metric -- see the detailed comment on
    # configs/default.yaml's curriculum.protocol key for the re-measured
    # numbers and working theory (perturbation's positives get diluted by
    # the ~4.6x larger localization pool). Sequential's own mild, remaining
    # forgetting is an acceptable trade-off, not a bug, once the real
    # training bug (hard-negative leakage across tiers) is fixed.
    protocol = cfg["curriculum"].get("protocol", "sequential")
    t0 = time.time()
    if protocol == "all_at_once":
        stage, merged = build_all_at_once_stage(cfg, data, splits_per_tier)
        print(f"[curriculum] protocol=all_at_once; training jointly on tiers "
              f"{sorted(main_tiers)}; held out as zero-shot inference: {held_out_tiers}",
              flush=True)
        curr_results = {stage.name: trainer.train_stage(stage, merged)}
    else:
        stages = [s for s in stages_from_config(cfg) if s.evidence_tier in main_tiers]
        print(f"[curriculum] protocol=sequential; training on {[s.name for s in stages]}; "
              f"held out as zero-shot inference: {held_out_tiers}", flush=True)
        curr_results = trainer.run_full_curriculum(splits_per_tier, stages)
    train_secs = time.time() - t0

    # final per-tier evaluation: TEST-only for trained tiers, VAL+TEST for
    # held-out (never-trained) tiers -- see eval_edges_for_tier.
    print("\n=== Final evaluation (all tiers) ===", flush=True)
    final = {}
    for tier, sp in splits_per_tier.items():
        pos, neg = eval_edges_for_tier(sp, tier, held_out_tiers)
        m = trainer.evaluate_split({"pos": pos, "neg": neg})
        final[tier] = m
        tag = "inference" if tier in held_out_tiers else "test"
        print(f"  {tier:14s} [{tag:9s}] AUPR={m['aupr']:.4f} AUROC={m['auroc']:.4f} "
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
