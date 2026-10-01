# K562 headline rerun after the negative-sampling leak fix

**Status (2026-10-02 04:00 IST): complete.** Ada job 467
(`slurm/leakfix_compare.sh`) ran 10 fixed runs (5 seeds × 2 protocols) and
2 legacy-pool sanity runs (seed 42).

**Sanity check passed:** with `exclude_eval_negatives: false`, the rerun
reproduces the paper's seed-42 numbers **exactly**, to 4 decimals on all 6
metrics (loc 0.9679/0.9611, pert 0.6984/0.9100, dual 0.9526/0.9883). Every
difference below therefore comes from the leak fix alone.

**What changed:** commit c26c03d (branch `leak-fix-and-benchmarks`). Before
the fix, training drew random negatives from the whole 1M-edge pool. That
pool included 100% of every tier's val/test negatives, so the model was
trained on the very negatives it was later scored on. The fixed runs
removed 300,000 val/test negatives from the training pool (1,000,000 →
700,000).

**Setting:** K562, Geneformer on, h384/l2 (`configs/sweep/k562_fm_h384_l2.yaml`),
run through `scripts/06_ablation.py`. Splits on disk are unchanged.
Dual evidence is never trained on: it is scored zero-shot on val+test.

## Results (mean ± std over seeds)

| Run | n | loc AUPR | loc AUROC | pert AUPR | pert AUROC | dual AUPR | dual AUROC |
|---|---|---|---|---|---|---|---|
| all_at_once, **pre-fix** (paper headline) | 5 | 0.9676±0.0002 | 0.9609±0.0002 | 0.6955±0.0017 | 0.9096±0.0004 | 0.9517±0.0007 | 0.9880±0.0002 |
| all_at_once, **fixed** | 5 | 0.9544±0.0007 | 0.9485±0.0007 | 0.5722±0.0114 | 0.8803±0.0039 | 0.9017±0.0040 | 0.9781±0.0009 |
| full_curriculum, **fixed** | 5 | 0.7450±0.0040 | 0.7602±0.0032 | **0.8165±0.0006** | **0.9602±0.0001** | **0.9562±0.0008** | **0.9902±0.0002** |
| scMultiomeGRN baseline (`results/baselines/scmultiomegrn_K562.json`; trains on train-split negatives only, so it was never leaky) | 1 | 0.9059 | 0.9141 | 0.6550 | 0.9185 | 0.8466 | 0.9705 |
| GRNBoost2 baseline | 1 | 0.5347 | 0.4987 | 0.2060 | 0.5504 | 0.1728 | 0.5245 |

## What it means

1. **The leak inflated the paper's headline.** For `all_at_once`, perturbation
   AUPR drops by 0.123 (0.696 → 0.572) and dual-evidence AUPR by 0.050
   (0.952 → 0.902). Both drops are far larger than the seed std.
2. **After the fix, the recommended model loses to scMultiomeGRN on
   perturbation** (AUPR 0.572 vs 0.655, AUROC 0.880 vs 0.919). It still wins
   localization and dual evidence, 4 of 6 metrics in all.
3. **After the fix, the curriculum beats `all_at_once` on the zero-shot
   headline.**
   - Dual evidence: AUPR 0.956 vs 0.902, AUROC 0.990 vs 0.978.
   - Perturbation: AUPR 0.817 vs 0.572.
   - Its weakness is forgetting localization (AUPR 0.745).
   - Against scMultiomeGRN it wins perturbation and dual evidence by wide
     margins, but loses localization.
   - This reverses the paper's current recommendation (`all_at_once`) and
     supports making the curriculum the headline model again (option (a) in
     the session notes).
4. **The legacy curriculum run settles the cause (seed 42, AUPR/AUROC):**

   | Run | loc | pert | dual |
   |---|---|---|---|
   | full_curriculum, legacy pool | 0.7509/0.7682 | 0.8329/0.9633 | 0.9619/0.9913 |
   | full_curriculum, fixed | 0.7423/0.7589 | 0.8162/0.9601 | 0.9549/0.9899 |
   | all_at_once, legacy pool | 0.9679/0.9611 | 0.6984/0.9100 | 0.9526/0.9883 |
   | all_at_once, fixed | 0.9540/0.9481 | 0.5788/0.8827 | 0.9063/0.9791 |

   - **The leak barely moved the curriculum** (pert −0.017, dual −0.007) but
     hit `all_at_once` hard (pert −0.120, dual −0.046).
   - A plausible mechanism: `build_all_at_once_stage` turns hard negatives
     off, so `all_at_once` draws *all* its negatives from the random pool,
     which is the leaky part. The curriculum takes half its negatives from
     tier-hierarchy hard negatives.
   - **Even before the fix, the curriculum already beat `all_at_once`** on
     perturbation and dual evidence at this size.
     `results/ablations/with_fm_h384l2_K562.json` (sequential + FM + h384/l2,
     pre-fix) gives pert 0.834 and dual 0.963, against all_at_once's 0.696 and
     0.952.
   - **The 2026-09-18 claim in results.md §3/§7, that `all_at_once` "wins all
     three tiers outright, including perturbation", is contradicted by
     results.md's own table (line 342).** `all_at_once` only wins
     localization.
   - So the paper's recommendation of `all_at_once` rested on a misread
     comparison, and the leak widened the gap. The curriculum should be the
     headline model. Its open weakness is localization forgetting
     (0.745 vs 0.954).
5. **Open checks before the paper changes:**
   - The second, still-unfixed leak: the TF-candidate graph is built with
     held-out positives excluded (`02_preprocess.py`; see
     `scmultiomegrn_generalization.md` §10).
   - The dual-evidence tier was also used for some design selection (see
     `docs/mevd_vs_scmultiomegrn.md` §6). So it is not a fully untouched test
     set.
