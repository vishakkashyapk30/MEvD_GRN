"""Multi-evidence curriculum stage manager (plan 6, 9.1-9.2)."""
from __future__ import annotations

from dataclasses import dataclass
from typing import List


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
