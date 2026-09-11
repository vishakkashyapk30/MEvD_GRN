#!/usr/bin/env python
"""Ablation runner (plan Part 10, Step 12). One ablation per invocation.

  full_curriculum | loc_only | pert_only | dual_only | all_at_once
  rna_only | gated_fusion | concat_fusion | no_gnn | with_replay
  coexpr_only | tf_cand_only | edge_mlp | gated_relations
  with_fm | fm_only

The `gated_fusion` and `concat_fusion` ablations swap the default biologically
role-aware integration (ATAC-as-accessibility-gate + co-expression + openness) for
the legacy convex gated fusion / plain concatenation, isolating the value of the
role-aware integration. Each ablation changes exactly one component vs. the main
model.

Graph / edge ablations (new):
  - coexpr_only: GNN message-passes only over the co-expression kNN prior
  - tf_cand_only: GNN message-passes only over the TF-candidate prior
  - edge_mlp: role-aware decoder uses a nonlinear edge embedding (MPNN-style)
    instead of linear w_coexp / w_open scalars
  - gated_relations: the GNN learns a per-layer softmax weight over relations
    (coexpr vs. TF-candidate) instead of a fixed unweighted sum -- see
    src/models/gnn.py GNNBackbone combine_mode. The learned weights are an
    interpretability artifact (which relation the model actually leans on),
    not expected to change accuracy much with only 2 relations; more useful
    once a 3rd (motif-based) relation exists.

Foundation-model ablations (new, plan.md Section 5; run
scripts/11_extract_fm_embeddings.py first):
  - with_fm: full model + a pretrained Geneformer gene embedding fused
    additively into the RNA channel.
  - fm_only: FM embedding ALONE -- hand-crafted RNA features/co-expression
    signature zeroed, isolating what the pretrained embedding contributes
    by itself (mirrors the rna_only / no_gnn pattern).

Training scope per ablation, in terms of which tier(s) get a gradient step
(everything else is evaluated as held-out zero-shot inference, val+test
combined -- see benchmarker.eval_edges_for_tier):
  - full_curriculum, rna_only, gated_fusion, concat_fusion, no_gnn, all_at_once,
    coexpr_only, tf_cand_only, edge_mlp, gated_relations, with_fm, fm_only:
    the MAIN model's protocol -- `curriculum.main_curriculum_tiers` only
    (default localization+perturbation); dual_evidence held out.
  - loc_only / pert_only / dual_only: train on ONLY that one named tier
  - with_replay: the OLD 3-stage curriculum with memory replay enabled

Usage:
  python scripts/06_ablation.py --config configs/k562_fast.yaml --ablation coexpr_only
"""
from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.data.dataset import load_celltype_data, load_splits
from src.evaluation.benchmarker import eval_edges_for_tier
from src.models.mevd_grn import MEvDGRN
from src.training.curriculum import build_all_at_once_stage, stages_from_config
from src.training.trainer import MEvDTrainer
from src.utils.io import ensure_dir, load_config, save_json

ABLATIONS = [
    "full_curriculum", "loc_only", "pert_only", "dual_only", "all_at_once",
    "rna_only", "gated_fusion", "concat_fusion", "no_gnn", "with_replay",
    "coexpr_only", "tf_cand_only", "edge_mlp", "gated_relations",
    "with_fm", "fm_only",
]


def build_model(cfg, ablation):
    m = cfg["model"]
    integration = m.get("integration", "role_aware")
    if ablation == "concat_fusion":
        integration = "concat"
    elif ablation == "gated_fusion":
        integration = "gated"
    graph_mode = m.get("graph_mode", "both")
    if ablation == "coexpr_only":
        graph_mode = "coexpr"
    elif ablation == "tf_cand_only":
        graph_mode = "tf_candidate"
    use_edge_mlp = bool(m.get("use_edge_mlp", False)) or (ablation == "edge_mlp")
    combine_mode = m.get("combine_mode", "sum")
    if ablation == "gated_relations":
        combine_mode = "gated"
    use_fm = bool(m.get("use_fm", False)) or (ablation in ("with_fm", "fm_only"))
    return MEvDGRN(
        rna_in_dim=m["rna_in_dim"], atac_in_dim=m["atac_in_dim"], hidden_dim=m["hidden_dim"],
        n_gnn_layers=m["n_gnn_layers"], dropout=m["dropout"],
        use_atac=(ablation != "rna_only"),
        use_gnn=(ablation != "no_gnn"),
        integration=integration,
        graph_mode=graph_mode,
        use_edge_mlp=use_edge_mlp,
        combine_mode=combine_mode,
        use_fm=use_fm, fm_in_dim=int(m.get("fm_in_dim", 768)),
    )


def select_stages(cfg, ablation):
    stages = stages_from_config(cfg)
    by_tier = {s.evidence_tier: s for s in stages}
    main_tiers = set(cfg["curriculum"].get("main_curriculum_tiers", by_tier))
    if ablation == "with_replay":                # explicit 3-stage reference ablation
        return stages
    if ablation == "loc_only":
        return [by_tier["localization"]]
    if ablation == "pert_only":
        return [by_tier["perturbation"]]
    if ablation == "dual_only":                  # train ONLY on dual_evidence's own split
        return [by_tier["dual_evidence"]]
    # full_curriculum + architectural ablations: match the main model's protocol
    # (main_curriculum_tiers only; dual_evidence held out as zero-shot inference).
    return [s for s in stages if s.evidence_tier in main_tiers]


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
    elif args.ablation == "fm_only":     # isolate the FM embedding's own contribution:
        data.rna_signature = torch.zeros_like(data.rna_signature)  # remove hand-crafted
        data.rna_features = torch.zeros_like(data.rna_features)    # co-expression + stats

    splits_per_tier = {t: load_splits("data/splits", cell_type, t) for t in tiers
                       if (Path("data/splits") / f"{cell_type}_{t}_splits.pt").exists()}

    model = build_model(cfg, args.ablation)
    trainer = MEvDTrainer(model, data, cfg, device)

    main_tiers = set(cfg["curriculum"].get("main_curriculum_tiers", splits_per_tier))
    # Tiers this ablation actually trains a gradient step on -- everything else
    # is evaluated as held-out zero-shot inference (val+test combined).
    if args.ablation == "all_at_once":
        trained_tiers = set(main_tiers)
    else:
        trained_tiers = {s.evidence_tier for s in select_stages(cfg, args.ablation)}
    held_out_tiers = [t for t in splits_per_tier if t not in trained_tiers]

    if args.ablation == "all_at_once":
        # Multi-task baseline: one stage on the UNION of the main tiers'
        # positives (dual_evidence excluded -- same held-out-inference
        # treatment as the main model gets, for a fair comparison).
        stage, merged = build_all_at_once_stage(cfg, data, splits_per_tier)
        trainer.train_stage(stage, merged)
    else:
        trainer.run_full_curriculum(splits_per_tier, select_stages(cfg, args.ablation))

    # Evaluate on dual-evidence (the paper's key comparison); TEST-only for
    # trained tiers, VAL+TEST for tiers this ablation never trained on.
    key_tier = "dual_evidence" if "dual_evidence" in splits_per_tier else list(splits_per_tier)[-1]
    result = {}
    for t, sp in splits_per_tier.items():
        pos, neg = eval_edges_for_tier(sp, t, held_out_tiers)
        result[t] = trainer.evaluate_split({"pos": pos, "neg": neg})
    ensure_dir("results/ablations")
    relation_weights = {name: w.tolist() for name, w in model.get_relation_weights().items()}
    save_json({"ablation": args.ablation, "cell_type": cell_type,
               "key_tier": key_tier, "trained_tiers": sorted(trained_tiers),
               "results": result, "relation_weights": relation_weights},
              f"results/ablations/{args.ablation}_{cell_type}.json")
    m = result[key_tier]
    print(f"[ablation {args.ablation}] {key_tier} AUPR={m['aupr']:.4f} "
          f"AUROC={m['auroc']:.4f} EP={m['early_precision']:.4f}", flush=True)


if __name__ == "__main__":
    main()
