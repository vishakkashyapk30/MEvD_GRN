# Plan: advisor feedback of 19 Nov 2026 (deadline Sat 21 Nov, 01:00)

Status legend: DONE / RUNNING / NEXT / LATER. Times are targets, IST.

## Status, Fri 20 Nov 03:30
Deck: `docs/meetings/presentation_15th_nov_22nd_nov.md` (built by `docs/meetings/build_deck_11_22.py` from the
compiled tables; rebuilt on every new result). Brief: section 7 of `docs/meetings/brief_metrics_and_baselines.md`.

| Item | Status |
|---|---|
| 1. Citations | DONE: survey `docs/reference/lit_survey/degree_controlled_evaluation.md`; deck slide 2; brief 7.1 |
| 2. Hub-controlled tests on BEAR-GRN, every method | DONE: `docs/results/bear/matched.md` (released 8, baselines 4, MeVD-GRN + L2 + ablations, 16 simple, 3 label-free); KO also with the knockout-label in-degree as covariate (`matched_kocov.md`) |
| 3. ATAC value | scMultiomeGRN DONE (MeVD-GRN 3 seeds, scMultiomeGRN official). BEAR-GRN MeVD-GRN RUNNING: 32/54 units (3 workers, done ~07:00). LINGER re-run on Macrophage_S1: full DONE, randomised ~05:00, peaks shuffled ~06:40 (job 93928) |
| 4. Simple / Gaussian models | DONE: BEAR-GRN 8 models x 2 feature sets x 9 datasets, plus knockout-label-only training and the in-degree rule per label source; scMultiomeGRN parameter ladder (beta-model, latent-space, degree-corrected latent-space, LDA/QDA/GNB/GMM) |
| 5. Non-GNN models | DONE: ridge Gaussian graphical model, RegDiffusion, DeepSEM on 9 BEAR-GRN datasets; MLP and gradient boosting in item 4 |

Problems met on the way: LINGER OOM at 27 GB (it sizes worker pools from node RAM and hard-codes 8 in preprocess;
fixed by capping joblib pools) and LINGER needs >= 50 metacells (one 474-cell population gave 22; fixed with k-means
pseudo-labels, 62 metacells).

## Item 1. Citations for the degree-matched and caliper-matched tests
Goal: every degree-controlled test in the deck has a peer-reviewed citation, or is a direct application of one.
1. RUNNING (Thu 19:45 to ~23:00): background literature search, report in
   `docs/reference/lit_survey/degree_controlled_evaluation.md` (PPI, drug-target, knowledge graphs, recommender
   popularity bias, GNN link prediction, matched case-control / stratified AUC, GRN benchmarks).
2. Already found: Aiyappa et al., ICML 2025 (degree-corrected negatives); Janes & Pepe, Biometrika 2009
   (covariate-adjusted ROC); Yu et al., Bioinformatics 2010 (degree-preserving negatives).
3. NEXT (Fri): re-frame the deck: caliper test = covariate-adjusted AUROC with degree as the covariate; add the
   Aiyappa et al. degree-corrected test next to it; update the sources slide and the brief.

## Item 2. Degree-controlled tests on BEAR-GRN for every method
1. DONE: `src/benchmarks/bear_matched.py` (commit 7c161fd). Same pairs for every method; covariate = target
   in-degree among the training TFs of the fold (never the scored TF's own labels). Five numbers per method:
   in-space AUROC, degree-matched (same TF, 3 nearest in-degree), caliper (same TF, within 0.1 log-degree),
   covariate-adjusted AUROC (Janes & Pepe; covariates TF and in-degree group), degree-corrected (Aiyappa et al.).
2. RUNNING: pulling the released networks of the 8 other datasets from Ada (/share1) and the new MeVD-GRN seeds.
3. NEXT (Thu night): run on all 9 datasets x {ChIP, KO, Union, Intersection} x {9 released methods, 4 simple
   baselines, MeVD-GRN 5 seeds, L2 5 seeds, ablations}; table + slide.

## Item 3. How much ATAC adds (three conditions: full / ATAC removed / ATAC shuffled across genes)
"Removed" and "shuffled" cut every ATAC route: ATAC features, promoter openness, the peak-derived motif features and
the accessibility filter of the TF-candidate graph. Scored with the normal and the hub-controlled tests.

| Method | Benchmark | Full | ATAC removed | ATAC shuffled | Status |
|---|---|---|---|---|---|
| MeVD-GRN | BEAR-GRN (9 datasets, seeds 42-44) | have | job 93386 | job 93386 | RUNNING, ~7 h (54 units, 2 GPUs) |
| MeVD-GRN | scMultiomeGRN (12 cell types, seeds 42-44) | have | job 93391 | job 93391 | RUNNING |
| scMultiomeGRN | its own benchmark (seed 666) | have | have (RNA-only ablation) | job 93391 | RUNNING |
| LINGER | BEAR-GRN, one small dataset (Macrophage_S1) | released network + our re-run | ATAC constant | ATAC shuffled | NEXT (Fri morning, after GPUs free up) |
| released BEAR methods | BEAR-GRN | BEAR paper's own ATAC_Gene_Shuffle results | - | paper | cite |

Caveat to state: on the scMultiomeGRN benchmark the ground truth itself is built from ATAC peaks (circular), so the
clean measurement is BEAR-GRN. LINGER cannot run without ATAC at all, so "removed" = constant ATAC.

## Item 4. Simple models with few parameters, including Gaussian models
1. NEXT (Thu night, laptop CPU): `scripts/bear/25_bear_simple_models.py`, same TF-held-out folds as MeVD-GRN,
   20 interpretable pair features (in-degree, correlation, GRNBoost2, RNA/ATAC/openness summaries, motif scores,
   Geneformer cosine). Models: logistic regression (21 parameters), LDA, QDA, Gaussian naive Bayes, Gaussian-mixture
   classifier, Gaussian-kernel logistic regression, MLP (~700 parameters), gradient boosting; each on "cell data only"
   and "all features". 9 datasets, seed 42.
2. NEXT: scMultiomeGRN benchmark: parameter ladder 2-4 (link counts, done) -> one parameter per TF (beta-model) ->
   latent-space model (Hoff et al. 2002, k = 2, 4, 8) -> Gaussian models on pair features; official + held-out-TF
   splits, all tests.
3. Every model is scored with item 2's degree-controlled tests.

## Item 5. Non-GNN architectures (label-free and supervised), as extra experiments
1. NEXT (Thu night, laptop GPU): label-free on BEAR-GRN: RegDiffusion (diffusion model, installed), a Gaussian
   graphical model by ridge neighbourhood regression (Meinshausen & Buehlmann 2006), and DeepSEM (sparse
   structural-equation VAE) on TFs + highly variable genes if time allows. Label-free models score all TFs at once,
   which qualifies for BEAR-GRN as the released methods do.
2. Supervised non-GNN: the MLP and gradient boosting of item 4 (same TF-held-out protocol).

## Timeline
| When | Work |
|---|---|
| Thu 19:45-23:00 | item 2 run; item 4 BEAR models; RegDiffusion/GGM; literature agent |
| Thu night (Ada) | item 3 GPU runs (BEAR + scMultiomeGRN) |
| Fri morning | item 4 scMultiomeGRN ladder; LINGER Macrophage_S1 (3 conditions); DeepSEM |
| Fri afternoon | score everything with item 2's tests; ATAC tables |
| Fri evening | deck slides + brief update with citations (item 1) |
| Fri 23:00-Sat 01:00 | checks, final push |

Rules kept: zero test-to-train leakage; Ada batch jobs only; no force-push; commits carry only the user's name.
