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


def build_all_at_once_stage(cfg: dict, data, splits_per_tier: Dict[str, dict]
                            ) -> Tuple[CurriculumStage, dict]:
    """Multi-task protocol: ONE stage trained on the UNION of the main
    curriculum tiers' positives, instead of sequential per-tier stages.

    Structurally avoids catastrophic forgetting -- there is no stage
    boundary to forget across -- and already beats the sequential 2-stage
    curriculum on every metric in prior ablations (see results.md). Mutates
    `cfg["curriculum"]["use_hard_negatives"]` to False, since tier-hierarchy
    hard negatives ("positive in a weaker tier, not in the current one")
    are undefined once tiers are merged into one.
    """
    from src.data.graph_builder import create_edge_splits

    main_tiers = set(cfg["curriculum"].get("main_curriculum_tiers", splits_per_tier))
    merge_tiers = [t for t in splits_per_tier if t in main_tiers]
    pos = torch.cat([data.evidence[t] for t in merge_tiers], dim=1)
    pos = torch.unique(pos, dim=1)
    dcfg = cfg["data"]
    merged = create_edge_splits(pos, data.negative_pool, float(dcfg["train_ratio"]),
                                float(dcfg["val_ratio"]), int(dcfg["neg_train_ratio"]), 5,
                                int(dcfg["seed"]))
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
