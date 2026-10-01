# K562 headline rerun after the negative-sampling leak fix

**Status: PARTIAL (2026-10-01 18:30 IST).** All_at_once has 4 of 5 seeds,
full_curriculum 3 of 5. The legacy-pool sanity runs (tasks 10–11) haven't
run yet. Ada job 467 is `slurm/leakfix_compare.sh`.

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
| all_at_once, **fixed** | 4 | 0.9546±0.0007 | 0.9487±0.0007 | 0.5703±0.0122 | 0.8796±0.0042 | 0.9012±0.0045 | 0.9779±0.0010 |
| full_curriculum, **fixed** | 3 | 0.7432±0.0032 | 0.7595±0.0040 | **0.8164±0.0003** | **0.9602±0.0001** | **0.9559±0.0009** | **0.9902±0.0002** |
| scMultiomeGRN baseline (`results/baselines/scmultiomegrn_K562.json`; trains on train-split negatives only, so it was never leaky) | 1 | 0.9059 | 0.9141 | 0.6550 | 0.9185 | 0.8466 | 0.9705 |
| GRNBoost2 baseline | 1 | 0.5347 | 0.4987 | 0.2060 | 0.5504 | 0.1728 | 0.5245 |

## What it means (provisional until all seeds and the legacy runs land)

1. **The leak inflated the paper's headline.** For `all_at_once`, perturbation
   AUPR drops by 0.125 (0.696 → 0.570) and dual-evidence AUPR by 0.051
   (0.952 → 0.901). Both drops are far larger than the seed std.
2. **After the fix, the recommended model loses to scMultiomeGRN on
   perturbation** (AUPR 0.570 vs 0.655, AUROC 0.880 vs 0.919). It still wins
   localization and dual evidence, 4 of 6 metrics in all.
3. **After the fix, the curriculum beats `all_at_once` on the zero-shot
   headline.**
   - Dual evidence: AUPR 0.956 vs 0.901, AUROC 0.990 vs 0.978.
   - Perturbation: AUPR 0.816 vs 0.570.
   - Its weakness is forgetting localization (AUPR 0.743).
   - Against scMultiomeGRN it wins perturbation and dual evidence by wide
     margins, but loses localization.
   - This reverses the paper's current recommendation (`all_at_once`) and
     supports making the curriculum the headline model again (option (a) in
     the session notes).
4. **Open checks before the paper changes:**
   - The legacy-pool full_curriculum run (task 11). It shows whether the
     curriculum's numbers were also inflated before the fix, which they
     likely were, and so whether the reversal comes from the fix rather than
     from model size.
   - The remaining seeds.
   - The second, still-unfixed leak: the TF-candidate graph is built with
     held-out positives excluded (`02_preprocess.py`; see
     `scmultiomegrn_generalization.md` §10).
   - The dual-evidence tier was also used for some design selection (see
     `docs/mevd_vs_scmultiomegrn.md` §6). So it is not a fully untouched test
     set.
