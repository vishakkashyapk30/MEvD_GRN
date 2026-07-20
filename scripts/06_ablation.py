#!/usr/bin/env python
"""Ablation runner (plan Part 10, Step 12). One ablation per invocation.

  full_curriculum | loc_only | pert_only | dual_only | all_at_once
  rna_only | gated_fusion | concat_fusion | no_gnn | with_replay

The `gated_fusion` and `concat_fusion` ablations swap the default biologically
role-aware integration (ATAC-as-accessibility-gate + co-expression + openness) for
the legacy convex gated fusion / plain concatenation, isolating the value of the
role-aware integration. Each ablation changes exactly one component vs. the full
model and is evaluated on the dual-evidence test split. Usage:
  python scripts/06_ablation.py --config configs/k562.yaml --ablation rna_only
"""
from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.data.dataset import load_celltype_data, load_splits
from src.data.graph_builder import create_edge_splits
from src.models.mevd_grn import MEvDGRN
from src.training.curriculum import CurriculumStage, stages_from_config
from src.training.trainer import MEvDTrainer
from src.utils.io import ensure_dir, load_config, save_json

ABLATIONS = ["full_curriculum", "loc_only", "pert_only", "dual_only", "all_at_once",
             "rna_only", "gated_fusion", "concat_fusion", "no_gnn", "with_replay"]


def build_model(cfg, ablation):
    m = cfg["model"]
    integration = m.get("integration", "role_aware")
    if ablation == "concat_fusion":
        integration = "concat"
    elif ablation == "gated_fusion":
        integration = "gated"
    return MEvDGRN(
        rna_in_dim=m["rna_in_dim"], atac_in_dim=m["atac_in_dim"], hidden_dim=m["hidden_dim"],
        n_gnn_layers=m["n_gnn_layers"], dropout=m["dropout"],
        use_atac=(ablation != "rna_only"),
        use_gnn=(ablation != "no_gnn"),
        integration=integration,
    )


def select_stages(cfg, ablation):
    stages = stages_from_config(cfg)
    by_tier = {s.evidence_tier: s for s in stages}
    if ablation in ("full_curriculum", "rna_only", "gated_fusion", "concat_fusion",
                    "no_gnn", "with_replay", "all_at_once"):
        return stages
    if ablation == "loc_only":
        return [by_tier["localization"]]
    if ablation == "pert_only":
        return [by_tier["perturbation"]]
    if ablation == "dual_only":
        return [by_tier["dual_evidence"]]
    return stages


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--ablation", required=True, choices=ABLATIONS)
    ap.add_argument("--device", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    cell_type = cfg["cell_type"]
    device = args.device or cfg["hardware"]["device"]
    if device.startswith("cuda") and not torch.cuda.is_available():
        device = "cpu"
    cfg = copy.deepcopy(cfg)
    cfg["checkpoint_dir"] = f"results/checkpoints/{cell_type}_ablation_{args.ablation}"
    if args.ablation == "with_replay":
        cfg["curriculum"]["use_memory_replay"] = True

    tiers = cfg["data"]["evidence_tiers"]
    data = load_celltype_data(cfg["paths"]["processed_dir"], tiers)
    if args.ablation == "rna_only":                     # remove ALL ATAC-derived signal
        data.atac_features = torch.zeros_like(data.atac_features)
        data.openness = torch.zeros_like(data.openness)

    splits_per_tier = {t: load_splits("data/splits", cell_type, t) for t in tiers
                       if (Path("data/splits") / f"{cell_type}_{t}_splits.pt").exists()}

    model = build_model(cfg, args.ablation)
    trainer = MEvDTrainer(model, data, cfg, device)

    if args.ablation == "all_at_once":
        # Multi-task baseline: one stage on the UNION of all tiers' positives.
        merged = _merge_tiers(splits_per_tier, data, cfg)
        total_epochs = sum(int(s["n_epochs"]) for s in cfg["curriculum"]["stages"])
        stage = CurriculumStage("AllAtOnce", "all", total_epochs, 1e-3, 5, False)
        cfg["curriculum"]["use_hard_negatives"] = False   # no tier hierarchy here
        trainer.train_stage(stage, merged)
    else:
        trainer.run_full_curriculum(splits_per_tier, select_stages(cfg, args.ablation))

    # Evaluate on dual-evidence test (the paper's key comparison)
    key_tier = "dual_evidence" if "dual_evidence" in splits_per_tier else list(splits_per_tier)[-1]
    result = {t: trainer.evaluate_split(sp["test"]) for t, sp in splits_per_tier.items()}
    ensure_dir("results/ablations")
    save_json({"ablation": args.ablation, "cell_type": cell_type,
               "key_tier": key_tier, "results": result},
              f"results/ablations/{args.ablation}_{cell_type}.json")
    m = result[key_tier]
    print(f"[ablation {args.ablation}] {key_tier} AUPR={m['aupr']:.4f} "
          f"AUROC={m['auroc']:.4f} EP={m['early_precision']:.4f}", flush=True)


def _merge_tiers(splits_per_tier, data, cfg):
    """Union all tiers' positives into a single split (multi-task baseline)."""
    pos = torch.cat([data.evidence[t] for t in splits_per_tier], dim=1)
    pos = torch.unique(pos, dim=1)
    dcfg = cfg["data"]
    return create_edge_splits(pos, data.negative_pool, float(dcfg["train_ratio"]),
                              float(dcfg["val_ratio"]), int(dcfg["neg_train_ratio"]), 5,
                              int(dcfg["seed"]))


if __name__ == "__main__":
    main()
