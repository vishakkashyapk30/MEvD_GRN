# Experiment MDH: testing the Marginal Dominance Hypothesis on BEAR-GRN

Branch `mdh`. The hypothesis and its predictions are in `HYPOTHESIS.md`. This file is the **pre-registration**,
written on 2026-11-14 before any run of the code in this folder. Results are appended below the line at the end, with
dates. Any change made after a result was seen is logged in "Deviations".

## Data and protocol (unchanged from `docs/benchmarks/bear_grn.md`)

- **Datasets.** The 9 BEAR-GRN datasets, using the local processed builds in `~/.cache/bear/data/processed/<ds>`
  (scripts/bear/21): measured genes, TF nodes, label tiers and GT edges indexed to the measured genes.
- **Ground truths.**
  - ChIP on all 9 datasets.
  - KO, Union and Intersection on K562 and the 5 mouse sets.
  - Core on the 4 human sets (secondary).
- **Folds.** TF-disjoint 5-fold cross-fitting: `bear_protocol.tf_disjoint_folds` with the same out-degrees, K and
  seed as scripts/22 and scripts/23, so fold *k* of seed *s* holds out exactly the same TFs as the in-degree
  baseline and MeVD-GRN. Checked against the stored `meta.json` folds before any scoring.
- **Inner validation.** The inner validation TFs are chosen exactly as in scripts/22.
- **Training labels.** Union of the dataset's label tiers (ChIP + KO where present) over the training TFs only,
  the same as the in-degree baseline and MeVD-GRN. A held-out TF contributes nothing to any training input, feature
  or statistic (asserted).
- **Seeds.** 42–46. The comparison with MeVD-GRN uses the seeds MeVD-GRN has locally (42–44; fewer on some mouse
  sets).

## Metrics

1. **Matched AUROC (primary; hypothesis-relevant).**
   - The probability that a positive pair scores above a negative pair **with the same TF** and **the same target
     in-degree stratum**.
   - Strata: the target's GT in-degree over all *other* GT TFs (leave-one-TF-out). Zero is its own stratum, then
     20 quantile bins of the non-zero values.
   - Universe: GT TFs that are TF nodes × measured genes, without self pairs. Pairs a method does not score get 0,
     as in BEAR.
   - Same TF removes the regulator marginal `α_i`. Same stratum removes the target marginal `β_j` (up to bin
     width). A pure marginal ranker scores ≈ 0.5.
   - Computed as Σ U / Σ n⁺n⁻ over all (TF, stratum) cells, ties counted ½.
2. **BEAR pooled AUROC and AUPRC (secondary; comparability with the paper).**
   - Scored with the verified port, `bear_metrics.score_roc_pr`, 20 draws, rng 43, exactly as
     `results/<ds>/scores/*.json`.

## Analyses

| ID | Prediction (HYPOTHESIS.md) | Analysis | Script |
|---|---|---|---|
| A1 | P1: GTs are near-additive in log-odds | Fit a two-way (β-model) logistic on each GT matrix and report the AUROC of `α_i + β_j`, of `α` alone and of `β` alone, next to the full-universe AUROC of the best released method, the in-degree baseline and MeVD-GRN | `01_gt_structure.py` |
| A2 | P2: published methods carry ≈ no pair-specific signal | Matched AUROC of the 9 released methods, Pearson, GRNBoost2, in-degree and MeVD-GRN on every dataset × GT | `02_audit.py` |
| A3 | P7: Intersection is anti-marginal | Spearman correlation over released methods of pooled ChIP AUROC vs Intersection AUROC (K562 + 5 mouse) | `02_audit.py` |
| A4 | P8: removing marginals in training recovers more pair-specific signal, and marginal + interaction beats the marginal baselines | MID-GRN (below) vs baselines | `03_midgrn.py`, `04_report.py` |

## MID-GRN (Marginal–Interaction Decomposition), frozen design

**Score.** For TF `i` and gene `j`, `s_ij = σ(w0 + w1·z_j + w2·z_j² + w3·α̂_i + w4·g_ij)`, where:
- `z_j` is the standardised `log1p` in-degree of `j` among the training TFs, i.e. the in-degree baseline. For a
  training TF's own row, the in-degree is computed leaving that TF out.
- `α̂_i` is the TF marginal: a ridge regression from TF node features to the logit of the TF's positive rate. The
  features are RNA (3), ATAC (3), 8 Geneformer PCs, has-motif, mean motif hit rate and co-expression signature
  norm. It is fit on the fit TFs.
- `g_ij` is the interaction head: an MLP (2 × 64, ReLU) on **pair** features only.
- `w` is fit by logistic regression on the inner validation TFs, on their natural distribution (every non-label
  pair is a negative).

**Pair features** (8; no node-level feature enters `g`):
1. Motif evidence of TF `i` in gene `j`'s top-5 linked peaks: best relative PWM score, and share of hit peaks
   (`motif_pair_feats.npz`).
2. Co-expression `s_i·s_j` over all 50 signature components (≈ Pearson r).
3. Background-removed co-expression: `s_i·s_j` without the 5 leading components, i.e. the dominant cell-state
   programs. This is the PSH analogue.
4. Geneformer cosine `cos(e_i, e_j)`.
5. Three **look-alike transfers**:
   - Formula: `LA_ij = Σ_k w_ik Y_kj − indeg_j / |T|`, over the 10 most similar training TFs `k ≠ i`, with
     `w = softmax(10·sim)`.
   - Similarities: Geneformer cosine, co-expression signature cosine, and motif-profile correlation across genes.
   - The subtraction centres out the target marginal. TFs without neighbours (no motif) get 0.

**Loss.** Matched-set conditional logit: each training positive `(i, j)` is contrasted with 8 negatives of the
**same TF** in the **same target in-degree stratum** (leave-own-TF-out in-degree among the training TFs), with a
softmax cross-entropy. `α_i` cancels exactly and `β_j` up to bin width. Training: Adam, lr 1e-3, batch 4096
positives, up to 8 epochs, early stopping on the inner-validation matched AUROC (patience 2).

**Variants.**

| Variant | Features | Loss | Role |
|---|---|---|---|
| `mid_full` | all 8 pair features | matched-set | headline |
| `mid_dataonly` | motif, co-expression ×2, Geneformer cosine (no look-alike, no labels inside `g`) | matched-set | Is the interaction signal from the cell data or from label transfer? |
| `mid_bce` | all 8 pair features **plus** `z_j`, `z_j²`, `α̂_i` as MLP inputs | ordinary BCE, uniform negatives 4:1 resampled each epoch; early stopping on inner-validation average precision over all validation pairs (the usual recipe); the MLP output is the score (no separate calibration) | The usual training recipe on the same information (tests P5/P8) |

For every variant the interaction-only score `g` is also written out, and scored on matched AUROC.

**Development.** Only on K562, seed 42. The design above was fixed before that run. If K562 forces a change (e.g. a
crash or a degenerate feature), it is logged under "Deviations", and the design is frozen again before the other 8
datasets and seeds 43–46 are run.

## Pre-registered success criteria

**Primary claim (model).**
- `mid_full` beats the TF-disjoint **in-degree baseline** and **MeVD-GRN (M3, `mevd_fm_h384_hubmotif`)** on
  matched AUROC for ChIP.
- "Significant" means the 95% CI of the across-dataset mean difference excludes 0. The CI comes from a cluster
  bootstrap over TFs within each dataset (2,000 resamples), with per-TF U statistics averaged over seeds.
- The same is reported for KO, Union and Intersection.
- Holm correction over the 2 comparators × 4 GTs.

**Secondary claim (BEAR's own metric).**
- `mid_full` pooled AUROC and AUPRC on ChIP above in-degree, MeVD-GRN and LINGER: the number of datasets out of 9,
  plus a two-sided Wilcoxon signed-rank test over datasets (seed means).
- The pre-registered BEAR win rule of `docs/benchmarks/bear_grn.md` §12.9 (above LINGER, in-degree and coverage) is
  counted over all dataset × GT × metric cells.

**Hypothesis claims.**
- **A1 supports P1** if the β-model AUROC ≥ the best released method's full-universe AUROC on ≥ 7/9 datasets (ChIP).
- **A2 supports P2** if every released method's matched AUROC is ≤ 0.55 on ChIP in ≥ 7/9 datasets.
- **A3 supports P7** if the Spearman ρ < 0.

**What would count against MDH.**
- Any released method ≥ 0.60 matched AUROC on most datasets.
- `mid_bce` ≥ `mid_full` on matched AUROC: in that case the conditional loss is not needed.

## Deviations

All three were made on 2026-11-14, after the K562 seed-42 development run of `mid_full` and before any other
dataset, seed or variant was run. The MID-GRN model, features and training are **unchanged** from the frozen design.
The changes concern evaluation validity only, and each was triggered by a failed negative control, not by
MID-GRN's own numbers.

**D1. Primary pair-specific metric: doubly-ranked AUROC (DR-AUROC), reported as excess over a marginal null.**

*Why the row-matched AUROC fails its controls:*
- The TF-disjoint in-degree baseline, which is target-only, scores 0.70 row-matched AUROC on K562 Intersection
  (0.63 on Core).
- With 19–37 GT TFs, a GT's leave-one-out in-degree is too noisy to remove the target marginal. ChIP hubness, which
  the training labels carry, leaks through.
- The row-matched AUROC is still reported, labelled as such.

*DR-AUROC, definition:*
1. Turn each TF row into percentiles against a reference gene set common to all rows (genes that are not an
   evaluated TF).
2. Within each target column, compare positive TFs' percentiles with negative TFs' percentiles. Ties count ½.
3. Any score that orders targets identically in every row ties exactly (every target marginal, whatever its
   source). So does any TF-only score.

*Further rules:*
- **Same-fold pairs only.** Only TF pairs held out in the same fold are compared. Across folds the in-degree
  baseline scored 0.45, because another fold's model was trained on TF *i*'s own labels. With both rules, in-degree
  scores exactly 0.500 on all 5 K562 GTs.
- **Marginal null.** 4 synthetic GTs are drawn from the GT's fitted two-way logistic model: same TF and target
  propensities, no pair structure. The same draws are used for every method of a dataset × GT.
- **Primary statistic:** DR excess = DR(real GT) − mean DR(null GTs), per method.
- **Inference:** a TF-cluster bootstrap (2,000 resamples, the same TF weights for real and null, averaged over
  seeds) gives CIs of method differences in DR excess.
- The success criteria above apply to DR excess in place of row-matched AUROC.

**D2. A1 oracle.**
- The in-sample β-model includes the scored pair's own label in its target effect (self-inclusion). It reaches 0.99
  row-matched AUROC on Naive_mESC KO.
- A1 now uses a leave-one-TF-out marginal oracle: logit(TF rate) + logit(target rate over the *other* GT TFs).
- The in-sample β-model AUROC is reported as an upper bound only.

**D3. Pooled BEAR scores.**
- `mid_full` is scored with 20 draws, as pre-registered.
- `mid_dataonly` and `mid_bce` use 5 draws, for time. The AUROC sd over draws is ~2e-4, and AUPRC is deterministic
  when the 10x negative cap is not binding.

## v2: post-hoc revision, frozen 2026-11-15 before any non-K562 v2 run (`06_midgrn_v2.py`)

**Status.** v1 stays the pre-registered model, and all of its results are reported. v2 is a **post-hoc** revision,
labelled as such everywhere.

**What I had seen before freezing v2:**
- v1 on all 9 datasets for seed 42: ChIP DR excess, and pooled ChIP on the 4 human sets.
- v2 test folds on K562 (the development dataset) only.

**Motivations:**
- MDH: target marginals are identifiable for unseen TFs, TF marginals are not.
- v1's training-side instability: calibration on 3–25 validation TFs, and its α weight flipping sign across folds.

**Design:**
1. **Target marginal from more labels.** Inputs are:
   - log in-degree rate among the training TFs, as v1;
   - **plus** log1p in-degree in the SC-MO-GRN-DB ChIP compendium (RN005/RN006, the L2 label source of
     `bear_grn.md` §6.2). Every source TF that appears in any GT or label file of the dataset is removed, so no
     evaluated TF ever contributes;
   - plus label-free target features: ATAC activity/openness (3), RNA mean/variance/detection (3), and whether the
     gene has a Geneformer token.
2. **No TF marginal.** α is dropped. On K562 it made the combiner rank whole TF rows by a noisy guess (ChIP AUROC
   0.558 with α vs 0.663 without).
3. **Out-of-fold combination.** 3 TF-disjoint inner folds over the training TFs give out-of-fold marginal
   features and `g` for every training TF.
   - The combiner is LightGBM: lr 0.05, 300 trees, 31 leaves, min_child 200, subsample 0.8.
   - It is monotone increasing in both in-degrees and in `g`.
   - It is fit on all training TFs, natural distribution.
4. **Interaction head.** `g` = v1's pair features, matched-set loss and MLP, with 4 fixed epochs. Held-out TFs get
   the mean of the 3 inner models.

**Outputs:**
- `mid2_full`: marginal + `g`, the headline.
- `mid2_marginal`: the same combiner without `g`. This is the strongest marginal-only baseline, which v2 must beat.
- `mid2_full_g`: `g` alone.

**Seeds:** 42–46 on all 9 datasets. Tests are as for v1, with `mid2_full` replacing `mid_full`. `mid2_marginal`
is added as a comparator in the primary family, so Holm runs over 3 comparators × 4 GTs.

---

## Results

Run 2026-11-14/15. Full tables: `results/REPORT.md`; raw: `results/a2_metrics.csv`, `results/pooled_bear.csv`, `results/a4_tests.json`.
Bug fixes after the runs, none of which changes a model or a score: `02_audit.py` skipped `mid2_*` folders (prefix filter), fixed and
re-run for `mid2_*` only; `05_report.py` used `df.gt` (resolves to `DataFrame.gt`), fixed.

**Pair-specific signal (DR-AUROC excess, mean over datasets, TF-cluster bootstrap, Holm).**
- v1 (pre-registered) beats in-degree on ChIP (+0.006, Holm p = 0.025) and Union (+0.009, p = 0.004); not on KO or Intersection.
  v1 does not beat MeVD-GRN on any GT and is significantly below it on KO and Union. The pre-registered A4 claim
  "MID-GRN > MeVD-GRN" fails.
- v2 (post-hoc) beats in-degree and its own marginal-only model on all 4 GTs (Δ +0.007 to +0.023, Holm p ≤ 0.02), and beats
  MeVD-GRN on ChIP (+0.006, p = 0.036). It ties MeVD-GRN on Union and Intersection and is below it on KO (−0.014).
- The pair-specific effects are small: DR excess ≈ 0.01–0.03 on a 0.5-centred scale.

**BEAR pooled AUROC/AUPRC (what reviewers compare).** v2 does not beat the strongest baseline, TF-disjoint in-degree, on the
headline GT. Detail:
- ChIP AUROC: mean Δ −0.009 vs in-degree. v2 is below it on Macrophage_S1 (−0.05), Macrophage_S2 (−0.04) and iPS (−0.01),
  about equal on mouse, and below LINGER on all 4 mESC time-course sets (0.60–0.61 vs 0.64–0.66).
- ChIP AUPRC equals in-degree. It beats LINGER on 9/9 datasets (Wilcoxon p = 0.004).
- Union: AUPRC beats MeVD-GRN and LINGER 6/6 (p = 0.031); AUROC beats MeVD-GRN 6/6 but is ≈ in-degree outside K562.
- Intersection AUROC beats in-degree 6/6 (p = 0.031).
- KO: everything is near or below 0.5 on mouse (in-degree 0.41); MeVD-GRN and the released methods are better.
- §12.9 win rule (above LINGER, in-degree and coverage): v2 19/52 cells, MeVD-GRN 16/52, v1 13/52.

**Verdict.** The MDH diagnosis holds: in-degree carries most of the pooled score, and the matched-set interaction head adds a
real but small pair-specific signal, which beats Pearson, GRNBoost2 and in-degree significantly on DR. The model built on it does
not win the BEAR pooled metrics: its gain is at most ≈ 0.01 AUROC except on K562 Union/Intersection, and it loses on the
small human datasets. v2 was designed after seeing v1 on all 9 datasets, so any further revision has no untouched BEAR data left
to be tested on.

A1 is evaluable only on K562 (released networks for the other 8 datasets are not local). A3: ρ = −0.16, p = 0.27, not significant.
