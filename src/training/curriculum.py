"""Multi-evidence curriculum stage manager (plan 6, 9.1-9.2)."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import torch


@dataclass
class CurriculumStage:
    name: str
    evidence_tier: str          # localization | perturbation | dual_evidence
    n_epochs: int
    learning_rate: float
    neg_ratio: int
    freeze_encoder: bool


def stages_from_config(cfg: dict) -> List[CurriculumStage]:
    out = []
    for s in cfg["curriculum"]["stages"]:
        out.append(CurriculumStage(
            name=s["name"], evidence_tier=s["evidence_tier"], n_epochs=int(s["n_epochs"]),
            learning_rate=float(s["lr"]), neg_ratio=int(s["neg_ratio"]),
            freeze_encoder=bool(s["freeze_encoder"]),
        ))
    return out


def build_all_at_once_stage(cfg: dict, splits_per_tier: Dict[str, dict]
                            ) -> Tuple[CurriculumStage, dict]:
    """Multi-task protocol: ONE stage trained on the UNION of the main
    curriculum tiers' positives, instead of sequential per-tier stages.

    Structurally avoids catastrophic forgetting -- there is no stage
    boundary to forget across. Whether this beats the sequential 2-stage
    curriculum is SIZE-DEPENDENT, not a fixed fact -- see results.md
    Section 3: at the small (307K-param, no-FM) size, joint training
    dilutes perturbation's positives with localization's larger pool and
    loses on perturbation/dual-evidence despite winning localization;
    at the larger FM-embedding size (h384/l2), that tradeoff disappears
    and all_at_once wins every tier. Check results.md for the current
    numbers at whatever size you're using -- don't assume either
    direction holds without checking. Mutates
    `cfg["curriculum"]["use_hard_negatives"]` to False, since tier-hierarchy
    hard negatives ("positive in a weaker tier, not in the current one")
    are undefined once tiers are merged into one.

    The merged train/val/test splits are built from the UNION of each
    tier's own already-persisted, leak-checked splits
    (`splits_per_tier[t][part]`) -- NOT by re-splitting the raw union of
    evidence from scratch. Re-splitting independently would only avoid
    leaking a tier's val/test positives into the merged train set by
    coincidence (matching seed, ratios, and sort order with the original
    per-tier split), not by construction; building from the existing
    per-tier train/val/test partitions guarantees it regardless of seed or
    ratio changes elsewhere. The assertion below checks this holds, rather
    than relying on it silently continuing to hold.
    """
    main_tiers = set(cfg["curriculum"].get("main_curriculum_tiers", splits_per_tier))
    merge_tiers = [t for t in splits_per_tier if t in main_tiers]

    def _union(part: str, key: str) -> torch.Tensor:
        pieces = [splits_per_tier[t][part][key] for t in merge_tiers]
        return torch.unique(torch.cat(pieces, dim=1), dim=1)

    merged = {part: {"pos": _union(part, "pos"), "neg": _union(part, "neg")}
             for part in ("train", "val", "test")}

    train_pos_set = set(map(tuple, merged["train"]["pos"].t().tolist()))
    for t in merge_tiers:
        for part in ("val", "test"):
            held_out = set(map(tuple, splits_per_tier[t][part]["pos"].t().tolist()))
            leaked = train_pos_set & held_out
            assert not leaked, (
                f"all_at_once leak: {len(leaked)} edge(s) from {t}.{part} "
                f"also appear in the merged all_at_once train split")

    total_epochs = sum(int(s["n_epochs"]) for s in cfg["curriculum"]["stages"]
                       if s["evidence_tier"] in main_tiers)
    cfg["curriculum"]["use_hard_negatives"] = False
    stage = CurriculumStage("AllAtOnce", "all", total_epochs, 1e-3, 5, False)
    return stage, merged


class CurriculumManager:
    def __init__(self, stages: List[CurriculumStage]):
        self.stages = stages
        self.idx = 0
        self.history: dict = {}

    @property
    def current_stage(self) -> CurriculumStage:
        return self.stages[self.idx]

    def advance(self) -> bool:
        if self.idx < len(self.stages) - 1:
            self.idx += 1
            return True
        return False

    def is_done(self) -> bool:
        return self.idx >= len(self.stages) - 1

    def get_encoder_freeze_state(self) -> bool:
        return self.current_stage.freeze_encoder
