# Degree-controlled evaluation, ATAC value, simple models and non-GNN models

Advisor meeting, follow-up on the feedback of 19 November 2026.

## 1. Summary

- **Citations (item 1).** Degree-matched test negatives are published (Aiyappa et al., ICML 2025; Cappelletti et al.,
  Bioinformatics Advances 2024; BioGraphBench-PPI, Frontiers in Bioinformatics 2026). The caliper test is the
  covariate-adjusted ROC of biostatistics (Janes & Pepe, Biometrika 2009) with degree as the covariate. No
  link-prediction paper uses a caliper on degree; the statistical method is established.
- **BEAR-GRN under degree control (item 2).** On ChIP, the best caliper-matched AUROC of any method is
  0.531 (MeVD-GRN); the best covariate-adjusted AUROC is 0.557 (MeVD-GRN). MeVD-GRN:
  0.531 caliper-matched,
  0.557 covariate-adjusted.
- **ATAC value (item 3), MeVD-GRN on BEAR-GRN, ChIP AUROC:** full 0.600,
  ATAC removed pending, ATAC shuffled
  pending.
- **Simple models (item 4) and non-GNN models (item 5):** slides 7-9.

## 2. Sources for the degree-controlled tests (item 1)

| Test | Closest published method | Published, or own |
|---|---|---|
| Degree-matched test (each test link paired with a non-link of similar endpoint degrees) | Degree-corrected negatives, Aiyappa et al., ICML 2025; degree-aware evaluation negatives, Cappelletti et al., Bioinformatics Advances 2024; "degree-matched" negatives preserving the endpoint log-degree-bin distribution, BioGraphBench-PPI, Frontiers in Bioinformatics 2026 | Published; ours matches each pair instead of the distribution |
| Caliper-matched test (non-link within 0.1 log-degree on both endpoints, else the link is dropped) | Matched case-control accuracy studies and the covariate-adjusted ROC (Janes & Pepe, Biometrics 2008, Biometrika 2009; Pepe, Fan & Seymour, Academic Radiology 2013); caliper matching (Rosenbaum & Rubin 1985; Austin 2011); dropped units (Stuart 2010) | Published statistical method, applied with degree as the covariate |
| AUROC within in-degree groups (BEAR-GRN) | Covariate-adjusted AUC = AUROC within covariate strata weighted by positives (Janes & Pepe 2009; Pepe et al. 2013); popularity-percentile evaluation (Bellogin et al., Information Retrieval Journal 2017) | Published statistical method |
| Degree-corrected test on BEAR-GRN (new) | Aiyappa et al., ICML 2025, directly | Published |

- No paper was found that applies a caliper on node degree to link-prediction negatives; the honest claim is the
  transfer of standard matching methods to a covariate (degree) that recent work shows confounds the evaluation.
- Full survey (33 sources, verification log): `docs/reference/lit_survey/degree_controlled_evaluation.md`.

## 3. BEAR-GRN under degree control: ChIP-seq ground truth (item 2)

- Covariate = target in-degree among the training TFs of the fold that holds the scored TF out (never the scored TF's
  own labels). Every method is scored on identical pairs; non-output pairs score 0.
- Degree-matched and caliper-matched pair each link (t, g) with a non-link (t, g') of the same TF, so the TF's own
  number of targets cannot help either. Covariate-adjusted = AUROC within (TF, in-degree group) cells, weighted by
  links. Degree-corrected = non-links drawn with the links' endpoint-degree distributions (Aiyappa et al.).
- The in-degree baseline is the covariate, so it scores about 0.5 on the caliper-matched and covariate-adjusted tests.

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | covariate-adjusted | degree-corrected |
|---|---|---|---|---|---|---|---|
| LINGER | 3 | 0.507 | pending | 0.532 | 0.514 | 0.517 | 0.544 |
| TRIPOD | 3 | 0.563 | pending | 0.519 | *0.517* | 0.520 | 0.559 |
| CellOracle | 3 | 0.515 | pending | 0.512 | 0.513 | 0.513 | 0.509 |
| SCENIC+ | 3 | 0.503 | pending | 0.504 | 0.504 | 0.503 | 0.504 |
| FigR | 3 | 0.498 | pending | 0.498 | 0.504 | 0.504 | 0.511 |
| GRaNIE | 3 | 0.500 | pending | 0.502 | 0.502 | 0.502 | 0.501 |
| Pando-GLM | 3 | 0.500 | pending | 0.506 | 0.502 | 0.503 | 0.503 |
| Pando-XGB | 3 | 0.501 | pending | 0.503 | 0.502 | 0.502 | 0.504 |
| In-degree (training TFs' links) | 2 | **0.609** | pending | 0.517 | 0.511 | *0.527* | *0.562* |
| Pearson | 2 | 0.520 | pending | 0.521 | 0.507 | 0.507 | 0.513 |
| GRNBoost2 | 2 | 0.563 | pending | 0.506 | 0.504 | 0.510 | 0.552 |
| All measured pairs | 2 | 0.500 | pending | 0.500 | 0.500 | 0.500 | 0.500 |
| MeVD-GRN | 2 | *0.600* | pending | *0.540* | **0.531** | **0.557** | **0.570** |
| MeVD-GRN, L2 (no cell-type labels) | 2 | 0.458 | pending | **0.553** | 0.496 | 0.498 | 0.486 |

## 4. BEAR-GRN under degree control: knockout ground truth (item 2)

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | covariate-adjusted | degree-corrected |
|---|---|---|---|---|---|---|---|
| LINGER | – | pending | pending | pending | pending | pending | pending |
| TRIPOD | – | pending | pending | pending | pending | pending | pending |
| CellOracle | – | pending | pending | pending | pending | pending | pending |
| SCENIC+ | – | pending | pending | pending | pending | pending | pending |
| FigR | – | pending | pending | pending | pending | pending | pending |
| GRaNIE | – | pending | pending | pending | pending | pending | pending |
| Pando-GLM | – | pending | pending | pending | pending | pending | pending |
| Pando-XGB | – | pending | pending | pending | pending | pending | pending |
| In-degree (training TFs' links) | – | pending | pending | pending | pending | pending | pending |
| Pearson | – | pending | pending | pending | pending | pending | pending |
| GRNBoost2 | – | pending | pending | pending | pending | pending | pending |
| All measured pairs | – | pending | pending | pending | pending | pending | pending |
| MeVD-GRN | – | pending | pending | pending | pending | pending | pending |
| MeVD-GRN, L2 (no cell-type labels) | – | pending | pending | pending | pending | pending | pending |

## 5. How much ATAC adds on BEAR-GRN (item 3)

- MeVD-GRN uses ATAC in four places: the ATAC features, promoter openness, the peak-derived motif features and the
  accessibility filter of the TF-candidate graph. "Removed" cuts all four; "shuffled" gives every gene another gene's
  ATAC values (same distribution, broken gene mapping). Seeds 42-44.
- LINGER is re-run on Macrophage_S1 with LINGER's own bulk reference: full data, ATAC peaks shuffled (genomic location
  broken, cell structure kept) and ATAC fully randomised (no ATAC information left).

ChIP-seq ground truth:

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | covariate-adjusted | degree-corrected |
|---|---|---|---|---|---|---|---|
| MeVD-GRN, full data | 2 | **0.600** | pending | **0.540** | **0.531** | **0.557** | **0.570** |
| MeVD-GRN, ATAC removed (all four routes) | – | pending | pending | pending | pending | pending | pending |
| MeVD-GRN, ATAC shuffled across genes | – | pending | pending | pending | pending | pending | pending |
| LINGER, released network (BEAR-GRN) | 3 | *0.507* | pending | *0.532* | *0.514* | *0.517* | *0.544* |
| LINGER re-run, full data | – | pending | pending | pending | pending | pending | pending |
| LINGER re-run, ATAC peaks shuffled | – | pending | pending | pending | pending | pending | pending |
| LINGER re-run, ATAC fully randomised | – | pending | pending | pending | pending | pending | pending |

## 6. How much ATAC adds on the scMultiomeGRN benchmark (item 3)

- The ground-truth network of this benchmark is itself built from the ATAC peaks (motif in an open promoter), so ATAC
  could only look more useful here than it is; BEAR-GRN (previous slide) is the clean test.

| Method | AUPR | AUROC | AUROC, degree-matched | AUROC, caliper-matched |
|---|---|---|---|---|
| MeVD-GRN, full data | **0.935** | **0.924** | **0.708** | **0.560** |
| MeVD-GRN, every ATAC input removed | pending | pending | pending | pending |
| MeVD-GRN, every ATAC input shuffled across genes | pending | pending | pending | pending |
| scMultiomeGRN, full data | *0.919* | *0.913* | *0.651* | *0.524* |
| scMultiomeGRN, ATAC removed (RNA only) | 0.916 | 0.910 | 0.650 | 0.511 |
| scMultiomeGRN, ATAC shuffled across TFs | pending | pending | pending | pending |

## 7. Simple and Gaussian models on BEAR-GRN, ChIP (item 4)

- Same TF-held-out folds as MeVD-GRN (seed 42), so no label of a scored TF is used. 20 interpretable features per
  (TF, gene) pair: target in-degree among training TFs (left out for the pair's own link), Pearson correlation,
  GRNBoost2 importance, RNA statistics of TF and gene, the gene's ATAC and promoter-openness summaries, 3 motif scores,
  Geneformer cosine. Training negatives are non-links of the same TF, so the models rank genes within a TF.
- Their scores are not calibrated across held-out TFs, so the per-TF and matched columns are the fair comparison.

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | covariate-adjusted | degree-corrected |
|---|---|---|---|---|---|---|---|
| In-degree | 2 | **0.609** | pending | 0.517 | 0.511 | 0.527 | *0.562* |
| MeVD-GRN | 2 | *0.600* | pending | **0.540** | **0.531** | **0.557** | **0.570** |
| logistic regression, all 20 features (21 parameters) | 2 | 0.411 | pending | 0.507 | 0.508 | 0.519 | 0.384 |
| logistic regression, cell data only (20 parameters) | 2 | 0.388 | pending | 0.514 | 0.513 | 0.532 | 0.377 |
| LDA, all 20 features (250 parameters) | 2 | 0.436 | pending | 0.498 | 0.502 | 0.514 | 0.407 |
| LDA, cell data only (228 parameters) | 2 | 0.416 | pending | 0.520 | 0.514 | 0.533 | 0.398 |
| QDA, all 20 features (460 parameters) | 2 | 0.428 | pending | 0.506 | 0.510 | 0.510 | 0.411 |
| QDA, cell data only (418 parameters) | 2 | 0.412 | pending | 0.509 | 0.512 | 0.515 | 0.406 |
| Gaussian naive Bayes, all 20 features (80 parameters) | 2 | 0.422 | pending | 0.521 | 0.519 | 0.528 | 0.435 |
| Gaussian naive Bayes, cell data only (76 parameters) | 2 | 0.401 | pending | *0.522* | *0.520* | 0.530 | 0.432 |
| Gaussian mixture classifier, all 20 features (1848 parameters) | 2 | 0.346 | pending | 0.477 | 0.480 | 0.484 | 0.443 |
| Gaussian mixture classifier, cell data only (1680 parameters) | 2 | 0.442 | pending | 0.488 | 0.497 | 0.501 | 0.500 |
| Gaussian-kernel logistic regression, all 20 features (301 parameters) | 2 | 0.468 | pending | 0.503 | 0.498 | 0.495 | 0.422 |
| Gaussian-kernel logistic regression, cell data only (301 parameters) | 2 | 0.469 | pending | 0.514 | 0.504 | 0.504 | 0.416 |
| MLP (32 units), all 20 features (705 parameters) | 2 | 0.503 | pending | 0.512 | 0.507 | 0.525 | 0.512 |
| MLP (32 units), cell data only (673 parameters) | 2 | 0.505 | pending | 0.505 | 0.508 | 0.531 | 0.502 |
| gradient boosting, all 20 features (200 trees) | 2 | 0.550 | pending | 0.507 | 0.507 | *0.535* | 0.547 |
| gradient boosting, cell data only (200 trees) | 2 | 0.545 | pending | 0.501 | 0.506 | 0.529 | 0.548 |

## 8. Parameter ladder on the scMultiomeGRN benchmark (item 4)

- Network models fitted on the training links only: the beta-model (one parameter per TF; Chatterjee, Diaconis & Sly
  2011), the latent-space model (Hoff, Raftery & Handcock 2002) and its degree-corrected form (Krivitsky et al. 2009).
  Gaussian classifiers use MAESTRO and GRNBoost2 pair features only (32 numbers, no training-graph input).

Official split:

| Method | AUPR | AUROC | AUROC, degree-matched | AUROC, caliper-matched |
|---|---|---|---|---|
| MeVD-GRN | **0.935** | **0.924** | **0.708** | 0.560 |
| Logistic regression on the two link counts (4 parameters) | 0.931 | 0.919 | 0.671 | 0.507 |
| beta-model (one parameter per TF) | *0.934* | *0.920* | 0.695 | 0.541 |
| Latent-space model, 2 dimensions | 0.571 | 0.628 | 0.590 | 0.548 |
| Latent-space model, 4 dimensions | 0.614 | 0.679 | 0.625 | **0.566** |
| Latent-space model, 8 dimensions | 0.672 | 0.734 | 0.656 | *0.565* |
| Degree-corrected latent-space model, 2 dimensions | 0.927 | 0.911 | 0.690 | 0.540 |
| Degree-corrected latent-space model, 8 dimensions | 0.918 | 0.902 | *0.701* | 0.561 |
| LDA on cell-data pair features (Gaussian) | 0.646 | 0.628 | 0.518 | 0.503 |
| QDA on cell-data pair features | 0.655 | 0.633 | 0.504 | 0.481 |
| Gaussian naive Bayes on cell-data pair features | 0.615 | 0.601 | 0.504 | 0.485 |
| Gaussian mixture classifier on cell-data pair features | 0.610 | 0.602 | 0.516 | 0.468 |

Held-out-TF split:

| Method | AUROC, pairs with one unseen TF | AUROC, pairs of two unseen TFs | AUROC, pairs with one unseen TF, caliper-matched |
|---|---|---|---|
| MeVD-GRN | 0.753 | *0.529* | 0.487 |
| Logistic regression on the two link counts (4 parameters) | **0.779** | 0.500 | 0.509 |
| beta-model (one parameter per TF) | **0.779** | 0.500 | 0.508 |
| Latent-space model, 2 dimensions | 0.772 | 0.500 | *0.527* |
| Latent-space model, 4 dimensions | 0.776 | 0.500 | **0.528** |
| Latent-space model, 8 dimensions | 0.777 | 0.500 | 0.522 |
| Degree-corrected latent-space model, 2 dimensions | 0.776 | 0.500 | 0.515 |
| Degree-corrected latent-space model, 8 dimensions | *0.778* | 0.500 | 0.514 |
| LDA on cell-data pair features (Gaussian) | 0.574 | 0.493 | 0.494 |
| QDA on cell-data pair features | 0.583 | 0.528 | 0.481 |
| Gaussian naive Bayes on cell-data pair features | 0.558 | 0.493 | 0.481 |
| Gaussian mixture classifier on cell-data pair features | 0.562 | **0.545** | 0.493 |

## 9. Label-free, non-GNN models on BEAR-GRN, ChIP (item 5)

- No labels at all, one model per dataset scores every TF, as the released BEAR-GRN methods do.

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | covariate-adjusted | degree-corrected |
|---|---|---|---|---|---|---|---|
| Gaussian graphical model (ridge neighbourhood regression) | 3 | 0.511 | pending | 0.525 | *0.512* | *0.513* | 0.508 |
| RegDiffusion (diffusion model with a learned gene x gene adjacency) | 3 | 0.504 | pending | **0.537** | 0.503 | 0.504 | 0.503 |
| DeepSEM (sparse structural-equation VAE) | – | pending | pending | pending | pending | pending | pending |
| Pearson | 2 | *0.520* | pending | 0.521 | 0.507 | 0.507 | 0.513 |
| GRNBoost2 | 2 | **0.563** | pending | 0.506 | 0.504 | 0.510 | **0.552** |
| SCENIC+ (released) | 3 | 0.503 | pending | 0.504 | 0.504 | 0.503 | 0.504 |
| LINGER (released) | 3 | 0.507 | pending | *0.532* | **0.514** | **0.517** | *0.544* |

## 10. Conclusions and next steps

- To be written from the final numbers.
