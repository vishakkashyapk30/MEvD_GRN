#!/usr/bin/env python
"""Run baselines and evaluate them on the SAME test splits as MEvD-GRN (plan Step 12).

Usage:
  python scripts/05_run_baselines.py --config configs/k562.yaml \
      --baselines grnboost2,regdiffusion,gmf_gae,scmultiomegrn --device cuda:0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.baselines import grnboost2_wrapper as gb2
from src.baselines import gmfgrn_wrapper as gmf
from src.baselines import regdiffusion_wrapper as rdw
from src.baselines import scmultiomegrn_wrapper as smg
from src.data.dataset import load_celltype_data, load_splits
from src.data.preprocessing import _resolve_glob
from src.evaluation.benchmarker import eval_edges_for_tier, evaluate_scorer_on_splits
from src.utils.io import ensure_dir, load_config, save_json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--baselines", default="grnboost2,regdiffusion,gmf_gae")
    ap.add_argument("--device", default=None)
    ap.add_argument("--scmultiomegrn-epochs", type=int, default=2000,
                    help="Paper default is 2000; use a smaller value for exploratory runs.")
    ap.add_argument("--scmultiomegrn-patience", type=int, default=100,
                    help="Early-stopping patience in 5-epoch validation checks.")
    ap.add_argument("--scmultiomegrn-edge-cap", type=int, default=200000,
                    help="Maximum labeled training edges sampled per epoch (0 = all).")
    ap.add_argument("--scmultiomegrn-graph-neighbors", type=int, default=500,
                    help="Maximum train-positive graph neighbors retained per TF.")
    args = ap.parse_args()

    cfg = load_config(args.config)
    cell_type = cfg["cell_type"]
    device = args.device or cfg["hardware"]["device"]
    tiers = cfg["data"]["evidence_tiers"]

    data = load_celltype_data(cfg["paths"]["processed_dir"], tiers)
    inv = {i: g for g, i in data.gene_index.items()}
    tf_names = [inv[int(i)] for i in data.tf_indices.tolist()]

    splits_per_tier = {t: load_splits("data/splits", cell_type, t) for t in tiers
                       if (Path("data/splits") / f"{cell_type}_{t}_splits.pt").exists()}

    main_tiers = list(cfg["curriculum"].get("main_curriculum_tiers", tiers))
    held_out_tiers = [t for t in splits_per_tier if t not in main_tiers]

    # target genes actually scored (limits GRNBoost2 to needed fits); use the
    # SAME eval-set composition as the metrics themselves (test-only for
    # trained tiers, val+test for held-out/inference tiers).
    eval_targets = set()
    for tier, sp in splits_per_tier.items():
        pos, neg = eval_edges_for_tier(sp, tier, held_out_tiers)
        eval_targets.update(pos[1].tolist())
        eval_targets.update(neg[1].tolist())
    eval_target_syms = [inv[int(i)] for i in eval_targets]

    pcfg = dict(cfg["data"]); pcfg["genome"] = cfg.get("genome"); pcfg["atac"] = cfg.get("atac", {})
    rna_path = _resolve_glob(cfg["paths"]["single_cell"]["rna_glob"])
    seed = int(cfg["data"]["seed"])
    workers = int(cfg["hardware"].get("n_dataloader_workers", 4))

    ensure_dir("results/baselines")
    for name in [b.strip() for b in args.baselines.split(",") if b.strip()]:
        print(f"\n########## Baseline: {name} ##########", flush=True)
        try:
            if name == "grnboost2":
                scorer = gb2.run_grnboost2(rna_path, data.gene_index, tf_names, pcfg, seed,
                                           workers, target_symbols=eval_target_syms)
            elif name == "regdiffusion":
                scorer = rdw.run_regdiffusion(rna_path, data.gene_index, pcfg, device,
                                              tf_names=tf_names)
            elif name in ("gmf_gae", "gmfgrn"):
                scorer = gmf.run_gmf_gae(rna_path, data.gene_index, pcfg, device, seed=seed)
            elif name in ("scmultiomegrn", "scmultiome_grn"):
                scorer = smg.run_scmultiomegrn(
                    rna_path, data, splits_per_tier, pcfg, device=device, seed=seed,
                    max_epochs=args.scmultiomegrn_epochs,
                    patience=args.scmultiomegrn_patience,
                    graph_max_neighbors=args.scmultiomegrn_graph_neighbors,
                    train_edge_cap=args.scmultiomegrn_edge_cap,
                    batch_size=int(cfg["training"]["batch_size"]),
                    main_curriculum_tiers=main_tiers,
                    cell_type=cell_type,
                )
            else:
                print(f"[skip] unknown baseline {name}", flush=True)
                continue
        except Exception as e:
            print(f"[error] {name} failed: {e}", flush=True)
            save_json({"baseline": name, "cell_type": cell_type, "error": str(e)},
                      f"results/baselines/{name}_{cell_type}.json")
            continue

        res = evaluate_scorer_on_splits(scorer, splits_per_tier, held_out_tiers)
        save_json({"baseline": name, "cell_type": cell_type, "results": res},
                  f"results/baselines/{name}_{cell_type}.json")
        for tier, m in res.items():
            print(f"  {name:12s} {tier:14s} AUPR={m['aupr']:.4f} AUROC={m['auroc']:.4f} "
                  f"EP={m['early_precision']:.4f}", flush=True)


if __name__ == "__main__":
    main()
