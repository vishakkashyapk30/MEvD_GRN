#!/usr/bin/env python
"""Evaluate a trained checkpoint, incl. cross-cell-type transfer (plan Step 12, 8.5).

Usage:
  python scripts/04_evaluate.py --config configs/k562.yaml \
      --checkpoint results/checkpoints/K562/final_model.pt --transfer_config configs/esc.yaml
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.data.dataset import load_celltype_data, load_splits
from src.models.mevd_grn import MEvDGRN
from src.training.trainer import MEvDTrainer
from src.utils.io import load_config, save_json


def build(cfg, device):
    mcfg = cfg["model"]
    # NOTE: must mirror every architecture flag scripts/03_train.py passes,
    # or a checkpoint trained with a non-default integration/graph_mode/
    # combine_mode/use_edge_mlp will silently load into the WRONG
    # architecture (mismatched or missing state_dict keys).
    return MEvDGRN(rna_in_dim=mcfg["rna_in_dim"], atac_in_dim=mcfg["atac_in_dim"],
                   hidden_dim=mcfg["hidden_dim"], n_gnn_layers=mcfg["n_gnn_layers"],
                   dropout=mcfg["dropout"], integration=mcfg.get("integration", "role_aware"),
                   graph_mode=mcfg.get("graph_mode", "both"),
                   use_edge_mlp=bool(mcfg.get("use_edge_mlp", False)),
                   combine_mode=mcfg.get("combine_mode", "sum"),
                   use_fm=bool(mcfg.get("use_fm", False)),
                   fm_in_dim=int(mcfg.get("fm_in_dim", 768))).to(device)


def eval_celltype(cfg, model, device):
    cell_type = cfg["cell_type"]
    tiers = cfg["data"]["evidence_tiers"]
    data = load_celltype_data(cfg["paths"]["processed_dir"], tiers)
    trainer = MEvDTrainer(model, data, cfg, device)
    res = {}
    for t in tiers:
        p = Path("data/splits") / f"{cell_type}_{t}_splits.pt"
        if not p.exists():
            continue
        res[t] = trainer.evaluate_split(load_splits("data/splits", cell_type, t)["test"])
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--transfer_config", default=None, help="evaluate transfer on another cell type")
    ap.add_argument("--device", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    device = args.device or cfg["hardware"]["device"]
    if device.startswith("cuda") and not torch.cuda.is_available():
        device = "cpu"

    model = build(cfg, device)
    ck = torch.load(args.checkpoint, map_location=device)
    model.load_state_dict(ck["model"])
    model.eval()

    print(f"=== In-domain eval ({cfg['cell_type']}) ===", flush=True)
    own = eval_celltype(cfg, model, device)
    for t, m in own.items():
        print(f"  {t:14s} AUPR={m['aupr']:.4f} AUROC={m['auroc']:.4f}", flush=True)
    save_json(own, f"results/{cfg['cell_type']}_eval.json")

    if args.transfer_config:
        tcfg = load_config(args.transfer_config)
        print(f"\n=== Transfer eval ({cfg['cell_type']} -> {tcfg['cell_type']}) ===", flush=True)
        # transfer needs a model whose dims match; reuse same hyperparams
        trans = eval_celltype(tcfg, model, device)
        for t, m in trans.items():
            print(f"  {t:14s} AUPR={m['aupr']:.4f} AUROC={m['auroc']:.4f}", flush=True)
        save_json(trans, f"results/{cfg['cell_type']}_to_{tcfg['cell_type']}_transfer.json")


if __name__ == "__main__":
    main()
