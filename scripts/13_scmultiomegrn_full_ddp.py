#!/usr/bin/env python
"""Full scMultiomeGRN baseline run, data-parallel across all available GPUs.

Trains localization then perturbation to completion (matching MEvD-GRN's own
`curriculum.main_curriculum_tiers` protocol), scores dual_evidence zero-shot
from the perturbation model -- the SAME protocol `scripts/05_run_baselines.py`
uses, just distributed.

Why hand-rolled DDP instead of `nn.parallel.DistributedDataParallel`: the
wrapper's training loop calls `model.encode(...)` once per epoch and then
`model.decode(...)` + `loss.backward(retain_graph=...)` many times per epoch
on the SAME retained autograd graph before a single `optimizer.step()`. DDP's
gradient-bucket reducer expects one `forward()`/`backward()` pair per
iteration and resets its "buckets ready" state inside `forward()` -- calling
`encode`/`decode` directly bypasses that and risks silent hangs. Instead:
every rank runs the identical epoch loop (same seeded RNG on every rank, so
`_sample_epoch_edges` draws the SAME sampled edges everywhere), each rank
decodes+backprops only its shard of that batch (`sampled[rank::world_size]`,
weighted by the GLOBAL batch size so the math matches a single-GPU run
exactly), then gradients are summed across ranks with `dist.all_reduce`
before `optimizer.step()`. This also cuts the dominant per-epoch cost --
repeated `retain_graph=True` backward passes through the shared encode graph
-- roughly `world_size`-fold, since each rank now does ~1/world_size as many
backward calls.

Usage (single node, 4 GPUs):
  torchrun --standalone --nproc_per_node=4 scripts/13_scmultiomegrn_full_ddp.py \\
      --config configs/k562.yaml
"""
from __future__ import annotations

import argparse
import copy
import os
import random
import sys
from pathlib import Path
from typing import Dict, Optional

import numpy as np
import torch
import torch.distributed as dist
import torch.nn.functional as F

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.baselines.common import TieredScorer, load_expression_df
from src.baselines.scmultiomegrn_wrapper import (
    ScMultiomeGRN,
    _ModelScorer,
    _edge_arrays,
    _edge_histograms,
    _make_training_graph,
    _normalized_adjacency,
    _sample_epoch_edges,
    _standardize,
    _validation_aupr,
)
from src.data.dataset import load_celltype_data, load_splits
from src.data.preprocessing import _resolve_glob
from src.evaluation.benchmarker import eval_edges_for_tier, evaluate_scorer_on_splits
from src.utils.io import ensure_dir, load_config, save_json


def log0(rank: int, msg: str) -> None:
    if rank == 0:
        print(msg, flush=True)


def fit_tier_ddp(rank: int, world_size: int, rna_x: torch.Tensor, atac_x: torch.Tensor,
                  expression: np.ndarray, split: dict, seed: int, max_epochs: int,
                  patience: int, graph_max_neighbors: int, train_edge_cap: int,
                  batch_size: int, learning_rate: float,
                  ckpt_path: Optional[str] = None) -> Optional[_ModelScorer]:
    dev = f"cuda:{rank}"
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)          # identical CPU-side param init on every rank
    torch.cuda.manual_seed_all(seed)

    n_nodes = rna_x.shape[0]
    graph_cpu = _make_training_graph(split["train"]["pos"], n_nodes, graph_max_neighbors, seed)
    log0(rank, f"[scMultiomeGRN-DDP] train-only graph: {n_nodes} nodes, "
               f"{graph_cpu.shape[1]} directed/self edges, world_size={world_size}")
    hist_cpu = _edge_histograms(expression, graph_cpu)
    graph = graph_cpu.to(dev)
    edge_attr = hist_cpu.to(dev)
    adj_norm = _normalized_adjacency(graph_cpu, n_nodes, dev)
    rna_x_dev = rna_x.to(dev)
    atac_x_dev = atac_x.to(dev)

    model = ScMultiomeGRN(rna_x.shape[1], atac_x.shape[1]).to(dev)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

    train_edges, train_labels = _edge_arrays(split["train"])
    val_edges, val_labels = _edge_arrays(split["val"])
    pos_weight = float((train_labels == 0).sum() / max(int((train_labels == 1).sum()), 1))
    rng = np.random.default_rng(seed)  # same seed -> same draws on every rank, no broadcast needed

    best_aupr = -np.inf
    best_state: Optional[dict] = None
    stale = 0

    for epoch in range(1, max_epochs + 1):
        sampled_edges, sampled_labels = _sample_epoch_edges(train_edges, train_labels, train_edge_cap, rng)
        total = sampled_edges.shape[1]
        local_idx = np.arange(rank, total, world_size)
        local_edges = sampled_edges[:, local_idx]
        local_labels = sampled_labels[local_idx]

        model.train()
        optimizer.zero_grad(set_to_none=False)   # keep grads as zero-tensors, never None (all_reduce needs a tensor on every rank)
        z = model.encode(rna_x_dev, atac_x_dev, graph, edge_attr, adj_norm)
        local_total = local_edges.shape[1]
        epoch_loss_local = 0.0
        for start in range(0, local_total, batch_size):
            end = min(start + batch_size, local_total)
            e = local_edges[:, start:end].to(dev)
            y = local_labels[start:end].to(dev)
            pred = model.decode(z, e[0], e[1])
            per_item = F.binary_cross_entropy(pred, y, reduction="none")
            weights = torch.where(y > 0, pos_weight, 1.0)
            # Normalize by the GLOBAL batch size (not the local shard size) so
            # that summing gradients across ranks reproduces exactly the
            # single-GPU full-batch gradient over all `total` sampled edges.
            loss = (per_item * weights).sum() / total
            loss.backward(retain_graph=end < local_total)
            epoch_loss_local += float(loss.detach())

        for p in model.parameters():
            dist.all_reduce(p.grad, op=dist.ReduceOp.SUM)
        optimizer.step()

        loss_tensor = torch.tensor([epoch_loss_local], device=dev)
        dist.all_reduce(loss_tensor, op=dist.ReduceOp.SUM)
        epoch_loss = float(loss_tensor.item())

        stop_flag = torch.zeros(1, device=dev)
        if epoch == 1 or epoch % 5 == 0 or epoch == max_epochs:
            if rank == 0:
                val_aupr = _validation_aupr(
                    model, rna_x_dev, atac_x_dev, graph, edge_attr, adj_norm,
                    val_edges, val_labels, batch_size,
                )
                print(f"[scMultiomeGRN-DDP] epoch={epoch:4d} loss={epoch_loss:.5f} "
                      f"val_AUPR={val_aupr:.5f}", flush=True)
                if val_aupr > best_aupr + 1e-6:
                    best_aupr = val_aupr
                    best_state = copy.deepcopy(model.state_dict())
                    stale = 0
                    if ckpt_path is not None:
                        torch.save({"epoch": epoch, "best_val_aupr": best_aupr, "model": best_state}, ckpt_path)
                        print(f"[scMultiomeGRN-DDP] checkpoint saved: epoch={epoch} "
                              f"val_AUPR={best_aupr:.5f} -> {ckpt_path}", flush=True)
                else:
                    stale += 1
                    if stale >= patience:
                        print(f"[scMultiomeGRN-DDP] early stop at epoch {epoch}", flush=True)
                        stop_flag[0] = 1.0
            dist.broadcast(stop_flag, src=0)
            if stop_flag.item() > 0.5:
                break

    if rank == 0 and best_state is not None:
        model.load_state_dict(best_state)
    for p in model.state_dict().values():
        dist.broadcast(p, src=0)   # every rank ends on rank0's best checkpoint

    model.eval()
    with torch.no_grad():
        z = model.encode(rna_x_dev, atac_x_dev, graph, edge_attr, adj_norm).detach()
    if rank != 0:
        return None
    return _ModelScorer(model, z, batch_size=max(batch_size, 65536))


def run_scmultiomegrn_ddp(rank: int, world_size: int, rna_path: str, data, splits_per_tier: Dict[str, dict],
                           pcfg: dict, seed: int, max_epochs: int, patience: int,
                           graph_max_neighbors: int, train_edge_cap: int, batch_size: int,
                           learning_rate: float, main_curriculum_tiers, cell_type: str,
                           out_tag: str = "") -> Optional[TieredScorer]:
    expression = load_expression_df(rna_path, data.gene_index, pcfg).to_numpy(np.float32)
    signature = data.rna_signature.float()
    if signature.ndim != 2 or signature.shape[1] <= 1 or not torch.isfinite(signature).all():
        signature = data.rna_features.float()
    rna_x = _standardize(torch.cat([signature, data.rna_features.float()], dim=1))
    atac_x = _standardize(torch.cat([data.atac_features.float(), data.openness.float()], dim=1))

    main_tiers = list(main_curriculum_tiers) if main_curriculum_tiers else list(splits_per_tier)
    trainable = [t for t in splits_per_tier if t in main_tiers]
    skipped = [t for t in splits_per_tier if t not in main_tiers]
    if skipped:
        log0(rank, f"[scMultiomeGRN-DDP] treating {skipped} as held-out zero-shot inference "
                   f"(scored with the {trainable[-1] if trainable else '?'} model)")

    ckpt_dir = Path("results/checkpoints") / f"scmultiomegrn_{cell_type}_ddp{out_tag}"
    if rank == 0:
        ckpt_dir.mkdir(parents=True, exist_ok=True)
    dist.barrier()

    scorers = {}
    for offset, tier in enumerate(trainable):
        log0(rank, f"[scMultiomeGRN-DDP] fitting tier={tier}")
        scorer = fit_tier_ddp(
            rank, world_size, rna_x, atac_x, expression, splits_per_tier[tier], seed + offset,
            max_epochs, patience, graph_max_neighbors, train_edge_cap, batch_size, learning_rate,
            ckpt_path=str(ckpt_dir / f"{tier}_best.pt"),
        )
        if rank == 0:
            scorers[tier] = scorer

    if rank != 0:
        return None
    fallback = trainable[-1] if trainable else None
    return TieredScorer(scorers, fallback_tier=fallback)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--epochs", type=int, default=2000)
    ap.add_argument("--patience", type=int, default=100)
    ap.add_argument("--edge-cap", type=int, default=200_000)
    ap.add_argument("--graph-neighbors", type=int, default=500)
    ap.add_argument("--out-tag", default="", help="suffix for results/checkpoints paths, "
                    "e.g. '_smoketest' -- keeps a short exploratory run from clobbering "
                    "the real results/baselines/scmultiomegrn_<cell_type>.json")
    args = ap.parse_args()

    rank = int(os.environ["RANK"])
    world_size = int(os.environ["WORLD_SIZE"])
    local_rank = int(os.environ["LOCAL_RANK"])
    torch.cuda.set_device(local_rank)
    dist.init_process_group(backend="nccl")

    cfg = load_config(args.config)
    cell_type = cfg["cell_type"]
    tiers = cfg["data"]["evidence_tiers"]
    data = load_celltype_data(cfg["paths"]["processed_dir"], tiers)
    splits_per_tier = {t: load_splits("data/splits", cell_type, t) for t in tiers
                       if (Path("data/splits") / f"{cell_type}_{t}_splits.pt").exists()}
    main_tiers = list(cfg["curriculum"].get("main_curriculum_tiers", tiers))
    held_out_tiers = [t for t in splits_per_tier if t not in main_tiers]

    pcfg = dict(cfg["data"]); pcfg["genome"] = cfg.get("genome"); pcfg["atac"] = cfg.get("atac", {})
    rna_path = _resolve_glob(cfg["paths"]["single_cell"]["rna_glob"])
    seed = int(cfg["data"]["seed"])

    log0(rank, f"[scMultiomeGRN-DDP] world_size={world_size} tiers={tiers} main_tiers={main_tiers}")

    scorer = run_scmultiomegrn_ddp(
        rank, world_size, rna_path, data, splits_per_tier, pcfg, seed,
        max_epochs=args.epochs, patience=args.patience,
        graph_max_neighbors=args.graph_neighbors, train_edge_cap=args.edge_cap,
        batch_size=int(cfg["training"]["batch_size"]), learning_rate=1e-5,
        main_curriculum_tiers=main_tiers, cell_type=cell_type, out_tag=args.out_tag,
    )

    if rank == 0:
        ensure_dir("results/baselines")
        res = evaluate_scorer_on_splits(scorer, splits_per_tier, held_out_tiers)
        save_json({"baseline": "scmultiomegrn", "cell_type": cell_type,
                   "note": f"Full {len(main_tiers)}-tier DDP run (world_size={world_size}, "
                           f"epochs={args.epochs}, edge_cap={args.edge_cap}); "
                           f"{held_out_tiers} scored zero-shot from the last trained tier's model, "
                           "matching MEvD-GRN's own curriculum.main_curriculum_tiers protocol.",
                   "results": res}, f"results/baselines/scmultiomegrn_{cell_type}{args.out_tag}.json")
        for tier, m in res.items():
            print(f"  scmultiomegrn {tier:14s} AUPR={m['aupr']:.4f} AUROC={m['auroc']:.4f} "
                  f"EP={m['early_precision']:.4f}", flush=True)

    dist.barrier()
    dist.destroy_process_group()


if __name__ == "__main__":
    main()
