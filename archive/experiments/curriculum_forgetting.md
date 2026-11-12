# The curriculum's localization forgetting (K562, label-free graph)

**Status (2026-10-09): the variant runs are DEFERRED, not run; BEAR-GRN has
priority on the GPU (coordinator decision, 2026-10-08 22:20).** The code-reading
findings in section 1 stand on their own. The 13 queued runs are kept, ready to
resume, in `~/.cache/local_runs/queue_task2_deferred.txt` (format of
`~/.cache/local_runs/scripts/run_queue.sh`). They need the `full_curriculum`
label-free checkpoints of seeds 42-44, which exist in
`results/checkpoints/K562_ablation_full_curriculum_fm_h384l2_labelfree_seed{42,43,44}/`.

Setting: the label-free K562 setup of
[labelfree_graph_rerun.md](labelfree_graph_rerun.md) (FM + h384/l2, both leak
fixes on). `full_curriculum` = Stage 1 localization (30 epochs, lr 1e-3) then
Stage 2 perturbation (15 epochs, lr 3e-4). Dual evidence is never trained on;
it is scored zero-shot on val+test in every run below.

The problem (leakfix_rerun.md, 5 seeds, before the graph fix): the curriculum wins
perturbation (AUPR 0.817) and dual evidence (0.956) but loses localization
(0.745), against 0.954 for `all_at_once`.

## 1. What Stage 2 actually trains on (code reading + logs)

- **Memory replay is already on in `full_curriculum`.** `configs/default.yaml`
  sets `curriculum.use_memory_replay: true` and `replay_weight: 0.1`, and no K562
  config overrides them. So every `full_curriculum` number, including the paper's
  and leakfix_rerun.md's, already replays localization train positives in Stage 2:
  `round(0.1 x 160,934) = 16,093` per epoch, at loss weight 0.1. The Stage-2 log
  shows `train pos=177027` (160,934 perturbation + 16,093 replayed).
- **`--ablation with_replay` does not isolate replay.** `06_ablation.py` maps it
  to `select_stages(...) = all stages`, which adds `Stage3_DualEvidence` and
  **trains on dual evidence**. That breaks the held-out dual protocol. Since replay
  is already on, `with_replay` = `full_curriculum` + a dual fine-tuning stage. It
  was therefore **not run**. The replay knob is studied through `replay_weight`
  instead (sections 2-3).
- **Stage 2 uses no hard negatives at all, although `use_hard_negatives: true`.**
  - `_build_train_edges` passes `exclude = torch.cat(replay)` to
    `get_hard_negatives`.
  - With replay on, that exclusion is the *whole replay memory*, i.e. every
    localization train positive.
  - The Stage-2 hard-negative candidates are "localization train positives that
    are not perturbation edges", so all of them are excluded.
  - The Stage-2 log confirms it: `neg=700000`, exactly the size of the
    restricted random pool. With hard negatives the count would be
    5 x 160,934 = 804,670, half of them hard.
  - **Consequence for leakfix_rerun.md item 4:** its explanation of why the leak
    hit `all_at_once` harder ("the curriculum takes half its negatives from
    tier-hierarchy hard negatives") does not match the code. With replay on, both
    protocols draw 100% of their training negatives from the random pool.
  - The forgetting is therefore not "hard negatives teach loc-only edges as
    label 0". It is ordinary drift: Stage 2 sees perturbation positives against
    random negatives, plus a 10% replay at 0.1 weight. Localization-only edges
    are 98.2% of all localization positives (20,622 of 1,168,069 are also
    perturbation edges) and 98.3% of the localization test positives. They are
    almost absent from Stage 2.

## 2. Variants (3 seeds each, 42-44; dual held out in all)

Every variant changes only what happens after Stage 1 (or after Stage 2). So each
one is **resumed** from the matching seed's `full_curriculum` checkpoint
(`06_ablation.py --init_from ... --start_stage ...`, added 2026-10-08). Resuming
reloads the restored best-in-stage weights, so the starting model is identical. It
also skips the 38-minute Stage 1. The only other difference from a continuous run
is the RNG stream of negative sampling. A resumed default run (`resumeS2`)
measures that effect.

| Tag | Config | Resumed from | What changes |
|---|---|---|---|
| `resumeS2` (control) | `k562_fm_h384_l2_labelfree.yaml` | Stage-1 best | nothing: the default Stage 2, re-run |
| `nohardneg` (seed 42 only) | `..._labelfree_nohardneg.yaml` | Stage-1 best | `use_hard_negatives: false`. Expected to be a no-op (section 1); confirmation only |
| `replay03` | `..._labelfree_replay03.yaml` | Stage-1 best | `replay_weight: 0.3` (replayed fraction and loss weight 0.3) |
| `replay10` | `..._labelfree_replay10.yaml` | Stage-1 best | `replay_weight: 1.0` (as many localization positives as perturbation positives, full weight) |
| `locrefresh` | `..._labelfree_locrefresh.yaml` | Stage-2 best | adds `Stage3_LocRefresh`: localization, 3 epochs, lr 1e-4, with replay of the localization + perturbation train positives |

Checkpoints: `results/checkpoints/K562_ablation_full_curriculum_fm_h384l2_labelfree_seed{42,43,44}/`.
Results: `results/ablations/full_curriculum_fm_h384l2_labelfree_<variant>_seedN_K562.json`.

## 3. Results

Deferred (see Status). The baseline the variants would be compared against is the label-free
`full_curriculum` of [labelfree_graph_rerun.md](labelfree_graph_rerun.md) section 5.

Resume with:
`nohup ~/.cache/local_runs/scripts/run_queue.sh ~/.cache/local_runs/queue_task2_deferred.txt &`.
About 2.5 h on the RTX 4060: about 13 min per resumed Stage-2 run, 4 min per loc-refresh run.
