# K562 rerun with a label-free TF-candidate graph (second leak fixed)

**Status (2026-10-09 02:10 IST): seeds 42-44 done for both protocols (local
RTX 4060; Ada account locked). Seeds 45-46 are queued (second local queue) and
are added to section 5.2 when they finish.** Section 6 makes label-free graphs
the default and adds a hard guard against leaky training.

Context: [leakfix_rerun.md](leakfix_rerun.md) fixed the first leak
(val/test negatives in the random-negative pool). This document fixes the
second, which that doc lists as open (item 5), and reruns the K562 headline
with both fixes on.

## 1. The leak

- `scripts/02_preprocess.py` builds the TF-candidate message-passing graph
  (top-500 proximally accessible, most co-expressed genes per TF) with
  `evidence=<all tiers, all splits>` and `exclude_positives=True`
  (`data.prior_exclude_positives`, default `true`).
- So every known positive, including every val/test positive of every tier,
  is skipped when filling a TF's 500 slots. Absence from the graph depends on
  the test labels.
- Measured on the paper's graph (`data/processed/K562/tf_candidate_edges.pt`,
  112,500 edges): **0** val or test positives of any tier are in the graph,
  against **2.8-3.0%** of the val/test negatives. So "pair is an edge of the
  input graph" implies "label 0" with certainty.

| Split (K562) | positives in old graph | positives in label-free graph | negatives in old graph | negatives in label-free graph |
|---|---|---|---|---|
| localization val | 0 / 174,870 | 4,245 (2.43%) | 4,346 (2.90%) | 3,068 (2.05%) |
| localization test | 0 / 175,237 | 4,142 (2.36%) | 4,428 (2.95%) | 3,151 (2.10%) |
| perturbation val | 0 / 34,820 | 779 (2.24%) | 4,346 (2.90%) | 3,068 (2.05%) |
| perturbation test | 0 / 34,457 | 799 (2.32%) | 4,428 (2.95%) | 3,151 (2.10%) |
| dual val | 0 / 3,964 | 89 (2.25%) | 555 (2.80%) | 405 (2.04%) |
| dual test | 0 / 3,979 | 83 (2.09%) | 554 (2.78%) | 396 (1.99%) |

(Localization and perturbation share one negative set per split, since
negatives are drawn once from the shared pool.)

In the label-free graph positives and negatives are graph edges at nearly
the same rate (about 2.3% vs 2.1%), so membership is close to uninformative.
The old graph and the label-free graph share 79,571 of their 112,500 edges;
32,929 edges (29%) differ, and those 32,929 slots of the label-free graph are
all known positives (train, val or test).

## 2. What changed

- **New config** `configs/sweep/k562_fm_h384_l2_labelfree.yaml`: inherits
  `k562_fm_h384_l2.yaml` (FM + h384/l2, negative-sampling fix on) and sets
  `data.prior_exclude_positives: false`,
  `paths.processed_dir: data/processed/K562_labelfree`,
  `paths.reuse_processed_from: data/processed/K562`.
- **`scripts/02_preprocess.py`, graph-only rebuild mode** (used only when a
  config sets `paths.reuse_processed_from`):
  - copies every label-free artifact unchanged from the source dir
    (features, signatures, openness, gene index, TF indices, co-expression
    graph, evidence, negative pool, `fm_gene_embeddings.npy`);
  - rebuilds only the TF-candidate graph, with the config's
    `prior_exclude_positives`;
  - first rebuilds the graph with `exclude_positives=True` from the copied
    arrays and requires it to equal the source graph (it does, tensor-equal),
    which proves the arrays are the ones the paper's graph was built from;
  - regenerates the global splits in memory and requires every tensor to be
    edge-identical to `data/splits/K562_*_splits.pt`; writes nothing to
    `data/splits/`.
  - Why not rerun the raw pipeline: `coexpression_signatures` calls ARPACK
    `svds` without a fixed start vector, so a raw rerun would also perturb the
    signatures, the co-expression graph and the candidate ranking. The
    graph-only rebuild makes the TF-candidate graph the only difference.
  - `paths.splits_dir` (default `data/splits`) is now honoured by
    `02_preprocess.py`, so no config can overwrite the canonical splits by
    accident. Every existing config keeps the old behaviour.
- `configs/default.yaml` is unchanged (`prior_exclude_positives: true`), so
  `data/processed/K562`, `data/splits/K562_*` and every existing config
  reproduce the paper's inputs.

## 3. Sanity checks

| Check | Result |
|---|---|
| Graph rebuilt with `exclude_positives=True` from the copied arrays equals `data/processed/K562/tf_candidate_edges.pt` | **yes** (tensor-equal) |
| Regenerated splits equal `data/splits/K562_{localization,perturbation,dual_evidence}_splits.pt` | **yes**, all 18 tensors (3 tiers x train/val/test x pos/neg) |
| `K562_labelfree` files byte-identical (md5) to `K562`: `rna_features_aligned.npy`, `atac_features_aligned.npy`, `rna_signature.npy`, `openness.npy`, `gene_index.json`, `tf_indices.json`, `coexpr_edges.pt`, `negative_pool.pt`, `evidence_*.pt`, `motif_edges.pt` (empty in both), `fm_gene_embeddings.npy` | **yes** (only `tf_candidate_edges.pt` and `summary.json` differ) |
| FM embeddings aligned with the gene index | yes: `gene_index.json` and `fm_gene_embeddings.npy` are both byte-identical copies of the paper's pair |
| Label-free graph size | 112,500 edges (225 TFs x 500), same as before |
| `scripts/00_sanity_check.py` | **passes** (all checks) |
| Training negatives | `exclude_eval_negatives: true`: 300,000 val/test negatives removed (1,000,000 -> 700,000), as in the fixed rows of leakfix_rerun.md |

Build log: `python scripts/02_preprocess.py --config configs/sweep/k562_fm_h384_l2_labelfree.yaml`
(9 s).

## 4. Runs

`scripts/06_ablation.py --config configs/sweep/k562_fm_h384_l2_labelfree.yaml
--ablation {full_curriculum,all_at_once} --seed N --tag fm_h384l2_labelfree_seedN`,
one at a time on the local RTX 4060 Laptop (8 GB). About 75 s per Stage-1
epoch, so roughly 50-60 min per run. Results:
`results/ablations/{full_curriculum,all_at_once}_fm_h384l2_labelfree_seedN_K562.json`.

## 5. Results

### 5.1 Seed 42 (interim, 2026-10-08 22:10 IST)

Same seed, same splits, same negatives. The only difference is the TF-candidate graph.

| Run | loc AUPR | loc AUROC | pert AUPR | pert AUROC | dual AUPR | dual AUROC |
|---|---|---|---|---|---|---|
| full_curriculum, fixed (leakfix_rerun.md) | 0.7423 | 0.7589 | 0.8162 | 0.9601 | 0.9549 | 0.9899 |
| full_curriculum, label-free | 0.7307 | 0.7525 | 0.8147 | 0.9601 | 0.9551 | 0.9899 |
| all_at_once, fixed (leakfix_rerun.md) | 0.9540 | 0.9481 | 0.5788 | 0.8827 | 0.9063 | 0.9791 |
| all_at_once, label-free | 0.9525 | 0.9466 | 0.5869 | 0.8846 | 0.9007 | 0.9783 |

Run times on the RTX 4060: full_curriculum 49.6 min, all_at_once 62.7 min.

### 5.2 Seeds 42-44 (mean ± sample std, n = 3)

"Leak-1 fixed only" = the Ada runs of [leakfix_rerun.md](leakfix_rerun.md) (negative-pool fix, old
label-dependent graph), restricted to the same three seeds. "Both leaks fixed" = this rerun.

| Run | n | loc AUPR | loc AUROC | pert AUPR | pert AUROC | dual AUPR | dual AUROC |
|---|---|---|---|---|---|---|---|
| full_curriculum, leak-1 fixed only (Ada) | 3 | 0.7432±0.0032 | 0.7595±0.0040 | 0.8164±0.0003 | 0.9602±0.0001 | 0.9559±0.0009 | 0.9902±0.0002 |
| **full_curriculum, both leaks fixed** | 3 | 0.7302±0.0006 | 0.7519±0.0005 | **0.8145±0.0001** | **0.9601±0.0001** | **0.9541±0.0009** | **0.9898±0.0002** |
| all_at_once, leak-1 fixed only (Ada) | 3 | 0.9543±0.0004 | 0.9484±0.0003 | 0.5758±0.0061 | 0.8815±0.0018 | 0.9028±0.0037 | 0.9783±0.0009 |
| **all_at_once, both leaks fixed** | 3 | **0.9530±0.0006** | **0.9472±0.0008** | 0.5780±0.0117 | 0.8821±0.0035 | 0.9010±0.0017 | 0.9783±0.0002 |
| scMultiomeGRN re-run (`results/baselines/scmultiomegrn_K562.json`; never leaky) | 1 | 0.9059 | 0.9141 | 0.6550 | 0.9185 | 0.8466 | 0.9705 |
| *for reference: leak-1 fixed only, all 5 seeds (Ada)*: full_curriculum | 5 | 0.7450±0.0040 | 0.7602±0.0032 | 0.8165±0.0006 | 0.9602±0.0001 | 0.9562±0.0008 | 0.9902±0.0002 |
| *for reference: leak-1 fixed only, all 5 seeds (Ada)*: all_at_once | 5 | 0.9544±0.0007 | 0.9485±0.0007 | 0.5722±0.0114 | 0.8803±0.0039 | 0.9017±0.0040 | 0.9781±0.0009 |

Paired by seed, the label-free run minus the leak-1-only run (AUPR):

| Protocol | loc | pert | dual |
|---|---|---|---|
| full_curriculum (seeds 42 / 43 / 44) | -0.0116 / -0.0110 / -0.0164 | -0.0015 / -0.0024 / -0.0017 | +0.0002 / -0.0031 / -0.0024 |
| all_at_once (seeds 42 / 43 / 44) | -0.0015 / -0.0010 / -0.0014 | +0.0081 / -0.0040 / +0.0024 | -0.0056 / +0.0006 / -0.0006 |

Files: `results/ablations/{full_curriculum,all_at_once}_fm_h384l2_labelfree_seed{42,43,44}_K562.json`.
`leak_status` was backfilled into these six files on 2026-10-09 (section 6). The runs finished before
the guard existed; their config and processed dir are the verified label-free ones.

### 5.3 What it means

1. **The second leak was small.** Building the TF-candidate graph label-free moves perturbation and dual
   evidence by at most 0.003 AUPR in either protocol, inside or near the seed spread. The one consistent
   change is the curriculum's localization AUPR, -0.013 (all 3 seeds, 0.730 vs 0.743).
   - **Caveat:** the leak-1-only runs were on Ada (2080 Ti), the label-free runs on the laptop
     (RTX 4060), so these paired differences also contain cross-hardware non-determinism.
   - To attribute the -0.013 to the graph alone, run the leak-1-only setting on the same GPU:
     `configs/sweep/k562_fm_h384_l2_legacygraph.yaml` (`legacy_allow_leaks`, about 50 min per seed).
     **Not run here:** the GPU now belongs to the BEAR-GRN and PBMC queues.
2. **The curriculum still beats `all_at_once` on the zero-shot and perturbation tiers:**
   - perturbation AUPR 0.8145 vs 0.5780, AUROC 0.9601 vs 0.8821;
   - dual-evidence AUPR 0.9541 vs 0.9010, AUROC 0.9898 vs 0.9783.

   It still loses localization (0.730 vs 0.953).
3. **Against scMultiomeGRN** (pert AUPR 0.655, dual AUPR 0.847), with both leaks fixed:
   - **The curriculum** wins perturbation (AUPR 0.8145 vs 0.6550, AUROC 0.9601 vs 0.9185) and dual
     evidence (0.9541 vs 0.8466; AUROC 0.9898 vs 0.9705), i.e. 4 of 6 metrics. It loses localization
     (0.7302 vs 0.9059; AUROC 0.7519 vs 0.9141).
   - **`all_at_once`** wins localization and dual evidence but loses perturbation (AUPR 0.5780 vs 0.6550,
     AUROC 0.8821 vs 0.9185), also 4 of 6.
   - So the conclusion of leakfix_rerun.md is unchanged: the curriculum is the model that beats
     scMultiomeGRN on the tiers that test generalisation.
   - Caveats from `docs/mevd_vs_scmultiomegrn.md` §6 still apply: the scMultiomeGRN baseline is an
     adapted, single-seed run (R9); dual evidence also drove design choices (R4).
4. **Localization forgetting remains the curriculum's weakness**
   ([curriculum_forgetting.md](curriculum_forgetting.md); the variant runs there are deferred).

## 6. Zero-leak by default (2026-10-09)

Rule: every model we run has zero test->train leakage. A leaky configuration is never a default, and it
cannot run by accident.

### 6.1 Defaults and processed dirs
- `configs/default.yaml`: `data.prior_exclude_positives: false`. The old comment ("keep known edges OUT
  ... no label leakage") had it backwards: excluding positives *is* the leak. `training.exclude_eval_negatives`
  stays `true`. A new top-level `legacy_allow_leaks: false` is added (section 6.2).
- `configs/k562.yaml`, `macrophage.yaml` and `mcf7.yaml` now point to `data/processed/{K562,Macrophage,MCF7}_labelfree`.
  - Each has `reuse_processed_from: data/processed/<ct>`, so `02_preprocess.py --config configs/<ct>.yaml`
    does the graph-only rebuild. It never re-runs the raw pipeline, and never writes to `data/splits/`.
  - Every config inheriting from them is now label-free too: `k562_fm`, `k562_fast`, `macrophage_fm`,
    `mcf7_fm`, and all `configs/sweep/*` except the two legacy ones.
  - The old dirs `data/processed/{K562,Macrophage,MCF7}` are untouched. They are used only by the legacy
    configs.
- `slurm/train.sh` now reads `processed_dir` from the config.
- **New label-free dirs.** Built 2026-10-09 by the graph-only rebuild (section 2). Both checks passed:
  the source graph is reproduced exactly, and the regenerated splits are edge-identical to
  `data/splits/<ct>_localization_splits.pt`. Every other file is byte-identical to the source.

  | Dir | old graph edges | label-free edges | shared | positives in label-free graph | test positives in old / new graph | test negatives in old / new graph |
  |---|---|---|---|---|---|---|
  | K562_labelfree | 112,500 | 112,500 | 79,571 | 32,929 (29.3%) | 0 / 4,142 (2.36%, loc) | 4,428 (2.95%) / 3,151 (2.10%) |
  | Macrophage_labelfree | 10,867 | 11,000 | 6,943 | 4,057 (36.9%) | 0 / 596 (3.54%) | 1,628 (4.45%) / 1,009 (2.76%) |
  | MCF7_labelfree | 124,739 | 125,000 | 78,521 | 46,479 (37.2%) | 0 / 6,936 (3.82%) | 6,281 (4.19%) / 3,961 (2.64%) |

  The old Macrophage and MCF7 graphs had fewer than TFs x 500 edges because excluding the dense
  positives left some TFs short of candidates.
- **Stamped after verification.** These dirs, used by the queued PBMC runs and by BEAR, were *verified*
  label-free before being stamped with `prior_exclude_positives: false` (and a `label_free_verified`
  note):
  - `~/.cache/pbmc/root/processed/{classical_monocyte,naive_cd4_t,naive_b,mdc}`
  - `~/.cache/bear/data/processed/{K562,K562__L2}`
  - `data/processed/K562_labelfree`

  Verification (`~/.cache/local_runs/scripts/verify_labelfree_dirs.py`): the TF-candidate graph, the
  PBMC RNA-only candidate graph and the co-expression kNN were rebuilt from each dir's own saved arrays
  with no labels (`evidence=None, exclude_positives=False`). All are tensor-equal to the saved graphs,
  and all motif graphs are empty. The negative control, `data/processed/K562`, fails the same check (its
  graph differs) and was not stamped.
- **Builders now record the flag.**
  - `02_preprocess.py` writes `prior_exclude_positives` into `summary.json` (full pipeline and
    graph-only rebuild).
  - `src/benchmarks/pbmc_data.py`, `src/benchmarks/bear_data.py` and `scripts/16_scmgrn_build_benchmark.py`
    write `prior_exclude_positives: false`; their graphs are built with `evidence=None`.
  - `02_preprocess.py`'s fallback default for a missing config key is now `false`.

### 6.2 Hard guard (`src/training/leak_guard.py`)
- **What `MEvDTrainer` refuses.** It raises `LeakageError` before building any training edge (in
  `train_stage` and in `_build_train_edges`) unless both of the following hold:
  - The processed dir's `summary.json` records label-free graphs. The accepted records are
    `prior_exclude_positives: false`, or the benchmark builders' equivalents `label_free: true` (PBMC)
    and `label_free_graphs: true` (BEAR). An explicit `prior_exclude_positives: true` always wins. A
    missing record, or data with no provenance at all, counts as **leaky**.
  - `training.exclude_eval_negatives` is true.
- **Coverage.**
  - `load_celltype_data` attaches the summary as `CellTypeData.provenance`.
  - Every training entry point goes through these two methods:
    - 03 / 06 / 17 / 22 / 27 via `train_stage`;
    - 12 (joint training, its own loop) via `_build_train_edges`.
- **The only way past it** is a top-level `legacy_allow_leaks: true`. That prints a loud `!!!` banner
  and records the run as leaky. It is set only in:
  - `configs/sweep/k562_fm_h384_l2_legacyneg.yaml`: pre-2026-09-30 numbers; label-dependent graph +
    leaky negative pool; now explicitly points to `data/processed/K562`.
  - `configs/sweep/k562_fm_h384_l2_legacygraph.yaml` (new): leakfix_rerun.md's "fixed" rows; label-dependent
    graph, negative fix on.
  - To reproduce any other pre-2026-10-08 number, add `paths.processed_dir: data/processed/<ct>` and
    `legacy_allow_leaks: true` (plus `training.exclude_eval_negatives: false` for pre-09-30 numbers) to a
    copy of its config.
- **Results JSONs** written from now on carry
  `leak_status = {graph_label_free, eval_negatives_excluded, legacy_allow_leaks, processed_dir, zero_leak}`:
  - 03 `results/<ct>_results.json`
  - 06 `results/ablations/*.json`
  - 12 `results/joint_training/*.json` (per trained cell type, plus the held-out one)
  - 17 scMultiomeGRN-benchmark record
  - 22 BEAR `meta.json` (top level and per fold)
  - 27 PBMC `metrics.json`
  - 04 (`_leak_status` key; 07 skips `_` keys)
  - 08 transfer JSON

  For the eval-only scripts (04, 08), the status describes the data and config used at evaluation.

### 6.3 Verification (2026-10-09)
- **`scripts/00_sanity_check.py`** passes, 44 checks (26 before). The synthetic graph is now label-free
  with an explicit provenance. The new section 8b checks that the guard:
  - refuses a label-dependent graph;
  - refuses a missing record;
  - refuses unknown provenance;
  - refuses `exclude_eval_negatives: false`;
  - refuses via `_build_train_edges` too;
  - accepts the three label-free records;
  - lets an explicit `prior_exclude_positives: true` override `label_free: true`;
  - lets `legacy_allow_leaks` through and records the run as leaky.

  Section 9 checks that the guard reads `data/processed/{K562,Macrophage,MCF7}` as leaky and the
  `*_labelfree` dirs as label-free.
- **End-to-end CPU smokes of `06_ablation.py`** (tiny model, 1 epoch per stage, real K562 data):
  - label-free config: trains, `leak_status.zero_leak = true`;
  - legacy graph *without* the flag: `LeakageError` before the first epoch;
  - `k562_fm_h384_l2_legacygraph.yaml`: banner, trains, `zero_leak = false`.

  Smoke outputs were deleted.
- **End-to-end CPU smoke of `27_pbmc_train_eval.py`** (mDC, `target_all`, 1 epoch): passes the guard
  with the stamped dir, writes `leak_status`, and reproduces the pre-guard smoke exactly (Cistrome AUROC
  0.7715, AUPR ratio 2.496).
- **Static verdict for every config with a local processed dir.** `configs/{k562,k562_fast,k562_fm,macrophage,macrophage_fm,mcf7,mcf7_fm}.yaml`
  and every `configs/sweep/*.yaml` pass as zero-leak, including all queued `*_labelfree*` configs. The
  two legacy sweep configs pass only via `legacy_allow_leaks`. All 4 PBMC and both BEAR processed dirs
  pass.
- **Not runnable here** (dirs not on this laptop), with the expected behaviour:
  - Pre-2026-10-08 ESC builds (`data/processed/ESC`) and existing scMultiomeGRN-benchmark processed dirs
    (`configs/scmgrn/*`; their old summaries have no record) will be **refused** until rebuilt.
  - Rebuilding ESC with the new default gives a label-free graph. The scmgrn graphs are label-free by
    construction, so re-running `scripts/16 ... ` re-writes their summaries with the record.
