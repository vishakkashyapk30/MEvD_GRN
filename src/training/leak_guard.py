"""Zero-leak guard (2026-10-08, docs/experiments/labelfree_graph_rerun.md s6).

Every MeVD-GRN training run must have both test->train leak fixes on:

1. Label-free message-passing graphs. The processed dir's summary.json must
   record that its TF-candidate graph was built without labels:
   `prior_exclude_positives: false` (scripts/02_preprocess.py), or the
   equivalent records written by the benchmark builders, `label_free: true`
   (PBMC, src/benchmarks/pbmc_data.py) and `label_free_graphs: true` (BEAR,
   src/benchmarks/bear_data.py). A missing record counts as leaky: every
   processed dir built before 2026-10-08 by 02_preprocess.py excluded all
   positives of all splits (label-dependent graph).
2. `training.exclude_eval_negatives: true` (the 2026-09-30 negative-pool fix;
   a missing key means true, the trainer's default).

`MEvDTrainer` calls `enforce` before any training edge is built, so no entry
point can train leaky by accident. The only way past it is an explicit
top-level `legacy_allow_leaks: true` in the config, meant ONLY for the
before/after legacy configs that reproduce pre-fix numbers; a loud warning is
printed and the status is recorded. `leak_status` is what scripts write into
every results JSON.
"""
from __future__ import annotations

from typing import Optional

LEGACY_KEY = "legacy_allow_leaks"


class LeakageError(RuntimeError):
    pass


def graph_label_free(provenance: Optional[dict]) -> bool:
    """True only if the processed dir explicitly records label-free graphs."""
    p = provenance or {}
    if "prior_exclude_positives" in p:
        return p["prior_exclude_positives"] is False
    return p.get("label_free") is True or p.get("label_free_graphs") is True


def leak_status(cfg: dict, data) -> dict:
    prov = getattr(data, "provenance", None) or {}
    st = {
        "graph_label_free": graph_label_free(prov),
        "eval_negatives_excluded": bool((cfg.get("training") or {}).get("exclude_eval_negatives", True)),
        LEGACY_KEY: bool(cfg.get(LEGACY_KEY, False)),
        "processed_dir": prov.get("_processed_dir"),
    }
    st["zero_leak"] = st["graph_label_free"] and st["eval_negatives_excluded"]
    return st


def enforce(cfg: dict, data) -> dict:
    """Raise LeakageError unless the run is leak-free or explicitly legacy."""
    st = leak_status(cfg, data)
    if st["zero_leak"]:
        return st
    problems = []
    if not st["graph_label_free"]:
        problems.append(f"the TF-candidate graph of {st['processed_dir'] or 'this data'} is not recorded "
                        "as label-free (summary.json needs prior_exclude_positives: false; build it with "
                        "scripts/02_preprocess.py and a *_labelfree config)")
    if not st["eval_negatives_excluded"]:
        problems.append("training.exclude_eval_negatives is false (val/test negatives stay in the "
                        "training negative pool)")
    if not st[LEGACY_KEY]:
        raise LeakageError("refusing to train with test->train leakage: " + "; ".join(problems) +
                           f". Set `{LEGACY_KEY}: true` only in a legacy before/after config.")
    bar = "!" * 100
    print(f"\n{bar}\n!!! LEGACY LEAKY RUN ({LEGACY_KEY}: true). Results are NOT valid evidence:\n!!! - "
          + "\n!!! - ".join(problems) + f"\n{bar}\n", flush=True)
    return st
