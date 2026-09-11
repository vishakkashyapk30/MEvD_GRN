#!/usr/bin/env python
"""Joint multi-cell-type training (plan.md Section 6, Workstream 5).

Trains ONE shared MEvD-GRN model on the LOCALIZATION tier of multiple cell
types at once (each epoch takes one full gradient pass per cell type, on
that cell type's own graph/features), then evaluates it zero-shot on a
held-out cell type it never saw at all.

Motivated by a clear, reproducible finding: a model trained on K562 alone
with a foundation-model gene embedding (Geneformer) transfers WORSE to new
cell types than one without it (see plan.md Section 6) -- consistent across
two independent target cell types (Macrophage, MCF7). The working
hypothesis is that the small network reading the embedding (`FMEncoder`)
over-specializes to K562-only patterns. Training on more than one cell type
at once is the cheapest way to test whether that's fixable.

Usage:
  python scripts/12_joint_train.py --configs configs/k562.yaml,configs/mcf7.yaml \
      --held_out_config configs/macrophage.yaml --device cuda --epochs 30
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
from src.training.curriculum import CurriculumStage
from src.training.trainer import MEvDTrainer
from src.utils.io import ensure_dir, load_config, save_json


def _build_trainer(cfg_path: str, model, device: str):
    cfg = load_config(cfg_path)
    cell_type = cfg["cell_type"]
    data = load_celltype_data(cfg["paths"]["processed_dir"], cfg["data"]["evidence_tiers"])
    splits = load_splits("data/splits", cell_type, "localization")
    trainer = MEvDTrainer(model, data, cfg, device)
    trainer.set_train_positives({"localization": splits})
    return trainer, splits, cell_type


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--configs", required=True, help="comma-separated training configs")
    ap.add_argument("--held_out_config", required=True,
                    help="zero-shot eval only, never trained on")
    ap.add_argument("--device", default=None)
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--neg_ratio", type=int, default=5)
    args = ap.parse_args()

    cfg_paths = args.configs.split(",")
    first_cfg = load_config(cfg_paths[0])
    device = args.device or first_cfg["hardware"]["device"]
    if device.startswith("cuda") and not torch.cuda.is_available():
        device = "cpu"

    mcfg = first_cfg["model"]
    model = MEvDGRN(
        rna_in_dim=mcfg["rna_in_dim"], atac_in_dim=mcfg["atac_in_dim"],
        hidden_dim=mcfg["hidden_dim"], n_gnn_layers=mcfg["n_gnn_layers"],
        dropout=mcfg["dropout"], integration=mcfg.get("integration", "role_aware"),
        graph_mode=mcfg.get("graph_mode", "both"),
        use_edge_mlp=bool(mcfg.get("use_edge_mlp", False)),
        combine_mode=mcfg.get("combine_mode", "sum"),
        use_fm=bool(mcfg.get("use_fm", False)), fm_in_dim=int(mcfg.get("fm_in_dim", 768)),
    ).to(device)
    print(f"[model] MEvD-GRN parameters: {model.count_parameters():,} (use_fm={model.use_fm})",
          flush=True)

    # Multiple MEvDTrainer facades sharing the SAME underlying model: each
    # cell type keeps its own graph/features/negative-pool/RNG, but every
    # gradient step updates the one shared set of parameters.
    trainers, names = [], []
    for p in cfg_paths:
        trainer, splits, cell_type = _build_trainer(p, model, device)
        trainers.append((trainer, splits))
        names.append(cell_type)
    print(f"[joint] training jointly on: {names}", flush=True)

    stage = CurriculumStage("Joint", "localization", args.epochs, args.lr, args.neg_ratio, False)
    tcfg = first_cfg["training"]
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr,
                            weight_decay=float(tcfg["weight_decay"]))
    clip = float(tcfg["clip_grad_norm"])
    bs = int(tcfg["batch_size"])

    t0 = time.time()
    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        for trainer, splits in trainers:
            edges, labels, sw = trainer._build_train_edges(stage, splits["train"]["pos"])
            perm = torch.randperm(edges.shape[1], generator=trainer.gen)
            edges, labels, sw = edges[:, perm], labels[perm], sw[perm]
            n_p = float((labels == 1).sum().item())
            n_n = float((labels == 0).sum().item())
            pos_weight = max(n_n / max(n_p, 1.0), 1.0)
            for start in range(0, edges.shape[1], bs):
                b = slice(start, start + bs)
                emb = trainer._encode()
                logits = trainer._decode(emb, edges[0, b], edges[1, b])
                loss = trainer._loss(logits, labels[b], stage, sw[b], pos_weight)
                opt.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), clip)
                opt.step()
                total_loss += float(loss.item())
        if epoch % 5 == 0 or epoch == args.epochs:
            msg = f"  [Joint] epoch {epoch:3d} loss={total_loss:.4f}"
            for (trainer, splits), name in zip(trainers, names):
                val = trainer.evaluate_split(splits["val"])
                msg += f" | {name}_val_aupr={val['aupr']:.4f}"
            print(msg, flush=True)

    train_secs = time.time() - t0
    print(f"\n=== Final evaluation ({train_secs / 60:.1f} min) ===", flush=True)
    final = {}
    for (trainer, splits), name in zip(trainers, names):
        m = trainer.evaluate_split(splits["test"])
        final[name] = m
        print(f"  {name:12s} [test, trained]       AUPR={m['aupr']:.4f} AUROC={m['auroc']:.4f}",
              flush=True)

    # Held-out cell type: zero-shot, never trained on -- use its FULL
    # evidence (train+val+test all held out), same convention as
    # eval_edges_for_tier uses for tiers never trained on within one cell type.
    held_cfg = load_config(args.held_out_config)
    held_data = load_celltype_data(held_cfg["paths"]["processed_dir"],
                                   held_cfg["data"]["evidence_tiers"])
    held_trainer = MEvDTrainer(model, held_data, held_cfg, device)
    held_splits = load_splits("data/splits", held_cfg["cell_type"], "localization")
    held_pos = torch.cat([held_splits["train"]["pos"], held_splits["val"]["pos"],
                          held_splits["test"]["pos"]], dim=1)
    held_neg = torch.cat([held_splits["train"]["neg"], held_splits["val"]["neg"],
                          held_splits["test"]["neg"]], dim=1)
    m = held_trainer.evaluate_split({"pos": held_pos, "neg": held_neg})
    final[held_cfg["cell_type"]] = m
    print(f"  {held_cfg['cell_type']:12s} [zero-shot, held out] AUPR={m['aupr']:.4f} "
          f"AUROC={m['auroc']:.4f}", flush=True)

    ensure_dir("results/joint_training")
    ensure_dir("results/checkpoints")
    # Tag includes use_fm -- otherwise a with-FM and a without-FM run on the
    # same cell types silently overwrite each other's results/checkpoint.
    tag = "_".join(names) + ("_fm" if model.use_fm else "_nofm")
    save_json({"trained_on": names, "held_out": held_cfg["cell_type"],
               "use_fm": model.use_fm, "results": final},
              f"results/joint_training/{tag}_holdout_{held_cfg['cell_type']}.json")
    torch.save({"model": model.state_dict()}, f"results/checkpoints/joint_{tag}_final.pt")
    print(f"\n[done] trained in {train_secs / 60:.1f} min -> "
          f"results/joint_training/{tag}_holdout_{held_cfg['cell_type']}.json", flush=True)


if __name__ == "__main__":
    main()
