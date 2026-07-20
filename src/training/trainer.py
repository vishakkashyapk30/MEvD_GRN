"""Training loop with multi-evidence curriculum (plan Part 7, Step 8).

Strategy: the prior graph is small enough to encode in full at every gradient
step (plan 7.2 "full-graph GNN encoding"); only the decoder edges are
mini-batched. Negatives are resampled fresh each epoch.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import torch

from src.data.dataset import CellTypeData
from src.evaluation.metrics import compute_all_metrics
from src.training.curriculum import CurriculumStage
from src.training.losses import bce_weighted_loss, focal_loss
from src.training.sampler import get_hard_negatives, sample_negatives
from src.utils.io import ensure_dir


class MEvDTrainer:
    def __init__(self, model, data: CellTypeData, config: dict, device: str):
        self.model = model.to(device)
        self.data = data.to(device)
        self.cfg = config
        self.device = device
        self.tcfg = config["training"]
        self.hierarchy = config["data"]["tier_hierarchy"]
        self.ckpt_dir = ensure_dir(config.get("checkpoint_dir", "results/checkpoints"))
        self.gen = torch.Generator().manual_seed(int(config["data"]["seed"]))
        self.history: Dict[str, list] = {}

    # ------------------------------------------------------------------ helpers
    def _encode(self):
        return self.model.encode(self.data.rna_features, self.data.atac_features,
                                 self.data.coexpr_edges, self.data.tf_candidate_edges)

    def _decode(self, emb, tf_idx, target_idx) -> torch.Tensor:
        return self.model.decode(emb, tf_idx, target_idx,
                                 self.data.rna_signature, self.data.openness)

    def _build_train_edges(self, stage: CurriculumStage, train_pos: torch.Tensor,
                           replay: Optional[List[torch.Tensor]] = None):
        """Return (edges (2,M), labels (M,), sample_weight (M,)) for one epoch."""
        n_pos = train_pos.shape[1]
        n_neg = stage.neg_ratio * max(n_pos, 1)
        # split requested negatives between hard and random when hard negatives exist
        hard = torch.zeros((2, 0), dtype=torch.long)
        if self.cfg["curriculum"].get("use_hard_negatives", True):
            n_hard = n_neg // 2
            hard = get_hard_negatives(self.data.evidence, stage.evidence_tier,
                                      self.hierarchy, n_hard, self.gen)
        n_rand = n_neg - hard.shape[1]
        rand = sample_negatives(self.data.negative_pool, n_rand, self.gen)
        neg = torch.cat([hard, rand], dim=1) if hard.shape[1] else rand

        pos_list = [train_pos]
        w_list = [torch.ones(n_pos)]
        if replay:
            rw = float(self.cfg["curriculum"].get("replay_weight", 0.1))
            for r in replay:
                if r.shape[1]:
                    pos_list.append(r)
                    w_list.append(torch.full((r.shape[1],), rw))
        pos_all = torch.cat(pos_list, dim=1)
        edges = torch.cat([pos_all, neg], dim=1)
        labels = torch.cat([torch.ones(pos_all.shape[1]), torch.zeros(neg.shape[1])])
        sample_w = torch.cat(w_list + [torch.ones(neg.shape[1])])
        return edges.to(self.device), labels.to(self.device), sample_w.to(self.device)

    def _loss(self, logits, labels, stage, sample_w, pos_weight):
        if self.tcfg.get("loss_type") == "focal":
            return focal_loss(logits, labels, self.tcfg["focal_alpha"],
                              self.tcfg["focal_gamma"], sample_w)
        return bce_weighted_loss(logits, labels, pos_weight=pos_weight, sample_weight=sample_w)

    # ------------------------------------------------------------------ stage
    def train_stage(self, stage: CurriculumStage, splits: dict,
                    replay: Optional[List[torch.Tensor]] = None) -> dict:
        self.model.set_encoder_frozen(stage.freeze_encoder)
        params = [p for p in self.model.parameters() if p.requires_grad]
        opt = torch.optim.AdamW(params, lr=stage.learning_rate,
                                weight_decay=float(self.tcfg["weight_decay"]))
        sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=max(stage.n_epochs, 1))
        bs = int(self.tcfg["batch_size"])
        clip = float(self.tcfg["clip_grad_norm"])
        eval_every = int(self.tcfg["eval_every"])
        patience = int(self.tcfg["early_stopping_patience"])

        best_aupr, bad_checks = -1.0, 0
        best_path = self.ckpt_dir / f"best_model_{stage.name}.pt"
        hist = []
        for epoch in range(1, stage.n_epochs + 1):
            self.model.train()
            edges, labels, sw = self._build_train_edges(stage, splits["train"]["pos"], replay)
            perm = torch.randperm(edges.shape[1], generator=self.gen)
            edges, labels, sw = edges[:, perm], labels[perm], sw[perm]
            # actual neg:pos ratio (localization is dense within accessible space)
            n_p = float((labels == 1).sum().item())
            n_n = float((labels == 0).sum().item())
            pos_weight = max(n_n / max(n_p, 1.0), 1.0)
            if epoch == 1:
                print(f"  [{stage.name}] train pos={int(n_p)} neg={int(n_n)} "
                      f"pos_weight={pos_weight:.2f}", flush=True)

            total_loss, n_batches = 0.0, 0
            for start in range(0, edges.shape[1], bs):
                b = slice(start, start + bs)
                emb = self._encode()
                logits = self._decode(emb, edges[0, b], edges[1, b])
                loss = self._loss(logits, labels[b], stage, sw[b], pos_weight)
                opt.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(params, clip)
                opt.step()
                total_loss += float(loss.item())
                n_batches += 1
            sched.step()
            avg_loss = total_loss / max(n_batches, 1)

            rec = {"epoch": epoch, "train_loss": avg_loss}
            if epoch % eval_every == 0 or epoch == stage.n_epochs:
                val = self.evaluate_split(splits["val"])
                rec.update({f"val_{k}": v for k, v in val.items()})
                print(f"  [{stage.name}] epoch {epoch:3d} loss={avg_loss:.4f} "
                      f"val_aupr={val['aupr']:.4f} val_auroc={val['auroc']:.4f}", flush=True)
                if val["aupr"] > best_aupr:
                    best_aupr, bad_checks = val["aupr"], 0
                    self.save_checkpoint(best_path, best_aupr)
                else:
                    bad_checks += 1
                    if bad_checks >= patience:
                        print(f"  [{stage.name}] early stop at epoch {epoch}", flush=True)
                        hist.append(rec)
                        break
            hist.append(rec)

        if best_path.exists():               # restore best-in-stage weights
            self.load_checkpoint(best_path)
        self.history[stage.name] = hist
        return {"best_val_aupr": best_aupr, "history": hist}

    def run_full_curriculum(self, splits_per_tier: Dict[str, dict],
                            stages: List[CurriculumStage]) -> dict:
        use_replay = bool(self.cfg["curriculum"].get("use_memory_replay", False))
        results = {}
        seen_train_pos: List[torch.Tensor] = []
        for stage in stages:
            tier = stage.evidence_tier
            if tier not in splits_per_tier:
                print(f"  [skip] no splits for tier {tier}", flush=True)
                continue
            print(f"\n=== {stage.name} (tier={tier}, lr={stage.learning_rate}, "
                  f"freeze_enc={stage.freeze_encoder}) ===", flush=True)
            replay = seen_train_pos if (use_replay and seen_train_pos) else None
            results[stage.name] = self.train_stage(stage, splits_per_tier[tier], replay)
            seen_train_pos = seen_train_pos + [splits_per_tier[tier]["train"]["pos"]]
        return results

    # ------------------------------------------------------------------ eval
    @torch.no_grad()
    def evaluate_split(self, split: dict, n_total_candidates: Optional[int] = None) -> dict:
        self.model.eval()
        pos, neg = split["pos"].to(self.device), split["neg"].to(self.device)
        edges = torch.cat([pos, neg], dim=1)
        labels = torch.cat([torch.ones(pos.shape[1]), torch.zeros(neg.shape[1])]).numpy()
        emb = self._encode()
        logits = self._decode(emb, edges[0], edges[1])
        scores = torch.sigmoid(logits).cpu().numpy()
        if n_total_candidates is None:
            # EPR base rate over the full bipartite candidate space (BEELINE convention)
            n_total_candidates = int(self.data.tf_indices.numel() * (self.data.n_genes - 1))
        return compute_all_metrics(labels, scores, n_total_candidates)

    # ------------------------------------------------------------------ ckpt
    def save_checkpoint(self, path, best_val_aupr: float = -1.0) -> None:
        torch.save({"model": self.model.state_dict(), "best_val_aupr": best_val_aupr}, path)

    def load_checkpoint(self, path) -> None:
        ck = torch.load(path, map_location=self.device)
        self.model.load_state_dict(ck["model"])
