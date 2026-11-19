# Degree-controlled evaluation, ATAC value, simple models and non-GNN models

Advisor meeting, follow-up on the feedback of 19 November 2026.

## 1. Summary

- **Citations (item 1).** Degree-matched test negatives are published (Aiyappa et al., ICML 2025; Cappelletti et al.,
  Bioinformatics Advances 2024; BioGraphBench-PPI, Frontiers in Bioinformatics 2026). The caliper test is the
  covariate-adjusted ROC of biostatistics (Janes & Pepe, Biometrika 2009) with degree as the covariate. No
  link-prediction paper uses a caliper on degree; the statistical method is established.
- **BEAR-GRN under degree control (item 2).** On ChIP, the best caliper-matched AUROC of any method is
  0.561 (logistic regression, cell data only (20 parameters)); the best covariate-adjusted AUROC is 0.568 (logistic regression, all 20 features (21 parameters)). MeVD-GRN:
  0.524 caliper-matched,
  0.537 covariate-adjusted.
- **ATAC value (item 3), MeVD-GRN on BEAR-GRN, ChIP AUROC:** full 0.576,
  ATAC removed 0.577, ATAC shuffled
  0.584.
- **Simple models (item 4) and non-GNN models (item 5):** slides 7-10.

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
| LINGER | 9 | 0.520 | 0.519 | 0.515 | 0.507 | 0.511 | 0.507 |
| TRIPOD | 9 | 0.532 | 0.508 | 0.509 | 0.509 | 0.510 | 0.507 |
| CellOracle | 9 | 0.519 | 0.512 | 0.510 | 0.510 | 0.510 | 0.509 |
| SCENIC+ | 9 | 0.501 | 0.500 | 0.501 | 0.501 | 0.501 | 0.502 |
| FigR | 9 | 0.503 | 0.504 | 0.500 | 0.502 | 0.502 | 0.500 |
| GRaNIE | 9 | 0.510 | 0.506 | 0.506 | 0.504 | 0.505 | 0.503 |
| Pando-GLM | 9 | 0.502 | 0.501 | 0.503 | 0.501 | 0.502 | 0.502 |
| Pando-XGB | 9 | 0.502 | 0.501 | 0.502 | 0.501 | 0.501 | 0.502 |
| In-degree (training TFs' links) | 9 | **0.610** | **0.717** | *0.516* | **0.525** | *0.532* | **0.580** |
| Pearson | 9 | 0.505 | 0.505 | 0.509 | 0.503 | 0.503 | 0.505 |
| GRNBoost2 | 9 | 0.508 | 0.501 | 0.504 | 0.504 | 0.505 | 0.513 |
| All measured pairs | 9 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 |
| MeVD-GRN | 9 | *0.576* | *0.662* | **0.523** | *0.524* | **0.537** | *0.553* |
| MeVD-GRN, L2 (no cell-type labels) | 9 | 0.508 | 0.516 | 0.512 | 0.497 | 0.497 | 0.492 |

## 4. BEAR-GRN under degree control: knockout ground truth (item 2)

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | covariate-adjusted | degree-corrected |
|---|---|---|---|---|---|---|---|
| LINGER | 5 | 0.514 | 0.459 | 0.512 | 0.489 | 0.483 | *0.522* |
| TRIPOD | 5 | **0.598** | 0.515 | 0.527 | 0.520 | 0.529 | **0.545** |
| CellOracle | 5 | 0.476 | 0.503 | 0.506 | 0.504 | 0.505 | 0.496 |
| SCENIC+ | 5 | 0.501 | 0.502 | 0.502 | 0.502 | 0.502 | 0.500 |
| FigR | 5 | 0.532 | 0.522 | 0.516 | 0.518 | 0.518 | 0.512 |
| GRaNIE | 5 | 0.490 | 0.490 | 0.498 | 0.498 | 0.498 | 0.502 |
| Pando-GLM | 5 | 0.511 | 0.512 | 0.512 | 0.510 | 0.510 | 0.505 |
| Pando-XGB | 5 | 0.510 | 0.511 | 0.510 | 0.509 | 0.510 | 0.505 |
| In-degree (training TFs' links) | 5 | 0.446 | 0.425 | 0.498 | 0.496 | 0.501 | 0.479 |
| Pearson | 5 | *0.540* | **0.559** | **0.558** | *0.545* | *0.547* | 0.511 |
| GRNBoost2 | 5 | 0.501 | 0.512 | 0.511 | 0.511 | 0.510 | 0.501 |
| All measured pairs | 5 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 |
| MeVD-GRN | 5 | 0.531 | 0.524 | *0.544* | **0.546** | **0.580** | 0.510 |
| MeVD-GRN, L2 (no cell-type labels) | 5 | 0.537 | *0.535* | 0.525 | 0.525 | 0.528 | 0.518 |

## 5. How much ATAC adds on BEAR-GRN (item 3)

- MeVD-GRN uses ATAC in four places: the ATAC features, promoter openness, the peak-derived motif features and the
  accessibility filter of the TF-candidate graph. "Removed" cuts all four; "shuffled" gives every gene another gene's
  ATAC values (same distribution, broken gene mapping). Seeds 42-44.
- LINGER is re-run on Macrophage_S1 with LINGER's own bulk reference: full data, ATAC peaks shuffled (genomic location
  broken, cell structure kept) and ATAC fully randomised (no ATAC information left).

ChIP-seq ground truth:

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | covariate-adjusted | degree-corrected |
|---|---|---|---|---|---|---|---|
| MeVD-GRN, full data | 9 | 0.576 | *0.662* | **0.523** | **0.524** | **0.537** | **0.553** |
| MeVD-GRN, ATAC removed (all four routes) | 7 | *0.577* | *0.662* | *0.522* | *0.522* | **0.537** | *0.549* |
| MeVD-GRN, ATAC shuffled across genes | 7 | **0.584** | **0.663** | *0.522* | 0.521 | *0.531* | **0.553** |
| LINGER, released network (BEAR-GRN) | 9 | 0.520 | 0.519 | 0.515 | 0.507 | 0.511 | 0.507 |
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
| In-degree | 9 | **0.610** | **0.717** | 0.516 | 0.525 | 0.532 | *0.580* |
| MeVD-GRN | 9 | 0.576 | 0.662 | 0.523 | 0.524 | 0.537 | 0.553 |
| logistic regression, all 20 features (21 parameters) | 9 | 0.576 | **0.717** | 0.526 | 0.534 | **0.568** | 0.576 |
| logistic regression, cell data only (20 parameters) | 9 | 0.523 | 0.584 | **0.579** | **0.561** | 0.562 | 0.549 |
| LDA, all 20 features (250 parameters) | 9 | 0.559 | 0.683 | 0.545 | 0.542 | 0.555 | 0.564 |
| LDA, cell data only (228 parameters) | 9 | 0.523 | 0.585 | *0.574* | 0.557 | 0.558 | 0.548 |
| QDA, all 20 features (460 parameters) | 9 | 0.557 | 0.631 | 0.560 | 0.548 | 0.557 | 0.568 |
| QDA, cell data only (418 parameters) | 9 | 0.528 | 0.569 | 0.568 | 0.551 | 0.555 | 0.557 |
| Gaussian naive Bayes, all 20 features (80 parameters) | 9 | 0.552 | 0.614 | 0.562 | 0.541 | 0.550 | 0.558 |
| Gaussian naive Bayes, cell data only (76 parameters) | 9 | 0.524 | 0.562 | 0.567 | 0.546 | 0.548 | 0.551 |
| Gaussian mixture classifier, all 20 features (1848 parameters) | 9 | 0.470 | 0.482 | 0.509 | 0.511 | 0.517 | 0.496 |
| Gaussian mixture classifier, cell data only (1680 parameters) | 9 | 0.498 | 0.509 | 0.504 | 0.506 | 0.509 | 0.500 |
| Gaussian-kernel logistic regression, all 20 features (301 parameters) | 9 | 0.561 | 0.685 | 0.532 | 0.537 | 0.547 | 0.575 |
| Gaussian-kernel logistic regression, cell data only (301 parameters) | 9 | 0.518 | 0.580 | 0.572 | *0.560* | 0.559 | 0.553 |
| MLP (32 units), all 20 features (705 parameters) | 9 | 0.575 | 0.700 | 0.528 | 0.534 | 0.553 | 0.574 |
| MLP (32 units), cell data only (673 parameters) | 9 | 0.527 | 0.582 | 0.569 | 0.556 | 0.555 | 0.542 |
| gradient boosting, all 20 features (200 trees) | 9 | *0.581* | *0.701* | 0.528 | 0.531 | 0.559 | **0.582** |
| gradient boosting, cell data only (200 trees) | 9 | 0.543 | 0.637 | 0.568 | 0.551 | *0.563* | 0.559 |

## 8. Simple and Gaussian models on BEAR-GRN, knockout ground truth (item 4)

- Training labels as for MeVD-GRN: binding (ChIP-derived) and knockout labels pooled. Binding labels dominate, and
  the models then rank knockout targets below chance within a TF. The second block trains the same models on the
  knockout labels only (6 datasets with knockout labels; same folds).

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | covariate-adjusted | degree-corrected |
|---|---|---|---|---|---|---|---|
| In-degree | 5 | 0.446 | 0.425 | 0.498 | 0.496 | 0.501 | 0.479 |
| Pearson | 5 | **0.540** | **0.559** | **0.558** | *0.545* | *0.547* | *0.511* |
| MeVD-GRN | 5 | 0.531 | 0.524 | *0.544* | **0.546** | **0.580** | 0.510 |
| logistic regression, all 20 features (21 parameters) | 5 | 0.419 | 0.406 | 0.489 | 0.481 | 0.428 | 0.473 |
| logistic regression, cell data only (20 parameters) | 5 | 0.373 | 0.328 | 0.428 | 0.405 | 0.389 | 0.465 |
| LDA, all 20 features (250 parameters) | 5 | 0.386 | 0.366 | 0.448 | 0.438 | 0.406 | 0.466 |
| LDA, cell data only (228 parameters) | 5 | 0.377 | 0.334 | 0.438 | 0.411 | 0.396 | 0.467 |
| QDA, all 20 features (460 parameters) | 5 | 0.383 | 0.373 | 0.467 | 0.454 | 0.443 | 0.467 |
| QDA, cell data only (418 parameters) | 5 | 0.377 | 0.368 | 0.462 | 0.442 | 0.431 | 0.466 |
| Gaussian naive Bayes, all 20 features (80 parameters) | 5 | 0.387 | 0.377 | 0.486 | 0.464 | 0.460 | 0.465 |
| Gaussian naive Bayes, cell data only (76 parameters) | 5 | 0.381 | 0.385 | 0.483 | 0.458 | 0.450 | 0.463 |
| Gaussian mixture classifier, all 20 features (1848 parameters) | 5 | 0.525 | 0.527 | 0.513 | 0.508 | 0.503 | *0.511* |
| Gaussian mixture classifier, cell data only (1680 parameters) | 5 | *0.537* | *0.537* | 0.519 | 0.514 | 0.521 | **0.514** |
| Gaussian-kernel logistic regression, all 20 features (301 parameters) | 5 | 0.416 | 0.420 | 0.507 | 0.492 | 0.497 | 0.475 |
| Gaussian-kernel logistic regression, cell data only (301 parameters) | 5 | 0.391 | 0.375 | 0.473 | 0.443 | 0.436 | 0.474 |
| MLP (32 units), all 20 features (705 parameters) | 5 | 0.433 | 0.437 | 0.518 | 0.504 | 0.518 | 0.480 |
| MLP (32 units), cell data only (673 parameters) | 5 | 0.387 | 0.368 | 0.465 | 0.440 | 0.430 | 0.475 |
| gradient boosting, all 20 features (200 trees) | 5 | 0.473 | 0.449 | 0.518 | 0.509 | 0.520 | 0.497 |
| gradient boosting, cell data only (200 trees) | 5 | 0.433 | 0.398 | 0.499 | 0.468 | 0.456 | 0.489 |

Trained on knockout labels only:

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | covariate-adjusted | degree-corrected |
|---|---|---|---|---|---|---|---|
| MeVD-GRN | 5 | 0.531 | 0.524 | 0.544 | 0.546 | 0.580 | 0.510 |
| In-degree rule, knockout labels only (no fit) | 5 | 0.682 | 0.773 | 0.687 | *0.693* | 0.706 | 0.500 |
| logistic regression, cell data only, trained on knockout labels only | 5 | 0.615 | 0.676 | 0.627 | 0.615 | 0.630 | 0.539 |
| logistic regression, all 20 features, trained on knockout labels only | 5 | **0.693** | *0.776* | **0.692** | **0.694** | *0.710* | 0.510 |
| LDA, cell data only, trained on knockout labels only | 5 | 0.620 | 0.684 | 0.630 | 0.618 | 0.635 | **0.542** |
| LDA, all 20 features, trained on knockout labels only | 5 | *0.684* | **0.777** | *0.689* | 0.692 | **0.711** | 0.511 |
| Gaussian naive Bayes, cell data only, trained on knockout labels only | 5 | 0.567 | 0.592 | 0.590 | 0.570 | 0.585 | 0.530 |
| Gaussian naive Bayes, all 20 features, trained on knockout labels only | 5 | 0.643 | 0.707 | 0.661 | 0.656 | 0.672 | 0.516 |
| gradient boosting, cell data only, trained on knockout labels only | 5 | 0.663 | 0.736 | 0.670 | 0.657 | 0.673 | *0.541* |
| gradient boosting, all 20 features, trained on knockout labels only | 5 | 0.679 | 0.761 | 0.681 | 0.679 | 0.697 | 0.511 |

- The matched tests control the in-degree counted over all label files. Trained on knockout labels only, the models
  pick up "genes that many knockouts affect", a hub of the knockout network that this covariate does not control;
  the in-degree rule on knockout labels alone measures how much of the gain that is.

Same knockout ground truth, matched tests re-drawn with the covariate = target in-degree over the knockout labels
of the training TFs (5 datasets):

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | covariate-adjusted | degree-corrected |
|---|---|---|---|---|---|---|---|
| In-degree rule, knockout labels only (= covariate) | 5 | 0.682 | 0.773 | 0.500 | 0.500 | 0.514 | 0.500 |
| In-degree, all labels | 5 | 0.446 | 0.425 | **0.618** | 0.477 | 0.476 | 0.479 |
| Pearson | 5 | 0.540 | 0.559 | 0.544 | 0.516 | 0.517 | 0.511 |
| MeVD-GRN | 5 | 0.531 | 0.524 | *0.608* | 0.501 | 0.503 | 0.510 |
| Logistic regression, cell data only, all labels | 5 | 0.373 | 0.328 | 0.526 | 0.469 | 0.463 | 0.465 |
| LDA, cell data only, all labels | 5 | 0.377 | 0.334 | 0.534 | 0.471 | 0.465 | 0.467 |
| logistic regression, cell data only, trained on knockout labels only | 5 | 0.615 | 0.676 | 0.587 | *0.538* | 0.548 | 0.539 |
| logistic regression, all 20 features, trained on knockout labels only | 5 | **0.693** | *0.776* | 0.520 | 0.508 | 0.548 | 0.510 |
| LDA, cell data only, trained on knockout labels only | 5 | 0.620 | 0.684 | 0.580 | **0.539** | *0.550* | **0.542** |
| LDA, all 20 features, trained on knockout labels only | 5 | *0.684* | **0.777** | 0.518 | 0.509 | **0.552** | 0.511 |
| Gaussian naive Bayes, cell data only, trained on knockout labels only | 5 | 0.567 | 0.592 | 0.593 | 0.527 | 0.534 | 0.530 |
| Gaussian naive Bayes, all 20 features, trained on knockout labels only | 5 | 0.643 | 0.707 | 0.551 | 0.517 | 0.536 | 0.516 |
| gradient boosting, cell data only, trained on knockout labels only | 5 | 0.663 | 0.736 | 0.580 | 0.530 | 0.549 | *0.541* |
| gradient boosting, all 20 features, trained on knockout labels only | 5 | 0.679 | 0.761 | 0.537 | 0.511 | 0.534 | 0.511 |

## 9. Parameter ladder on the scMultiomeGRN benchmark (item 4)

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

## 10. Label-free, non-GNN models on BEAR-GRN (item 5)

- No labels at all, one model per dataset scores every TF, as the released BEAR-GRN methods do.

ChIP-seq ground truth:

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | covariate-adjusted | degree-corrected |
|---|---|---|---|---|---|---|---|
| Gaussian graphical model (ridge neighbourhood regression) | 9 | 0.502 | 0.501 | 0.506 | 0.505 | *0.505* | 0.502 |
| RegDiffusion (diffusion model with a learned gene x gene adjacency) | 9 | 0.505 | 0.510 | *0.520* | **0.511** | **0.511** | *0.507* |
| DeepSEM (sparse structural-equation VAE) | 9 | *0.514* | **0.526** | **0.531** | **0.511** | **0.511** | **0.513** |
| Pearson | 9 | 0.505 | 0.505 | 0.509 | 0.503 | 0.503 | 0.505 |
| GRNBoost2 | 9 | 0.508 | 0.501 | 0.504 | 0.504 | *0.505* | **0.513** |
| SCENIC+ (released) | 9 | 0.501 | 0.500 | 0.501 | 0.501 | 0.501 | 0.502 |
| LINGER (released) | 9 | **0.520** | *0.519* | 0.515 | *0.507* | **0.511** | *0.507* |

Knockout ground truth:

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | covariate-adjusted | degree-corrected |
|---|---|---|---|---|---|---|---|
| Gaussian graphical model (ridge neighbourhood regression) | 5 | *0.523* | *0.534* | 0.528 | *0.527* | *0.526* | 0.508 |
| RegDiffusion (diffusion model with a learned gene x gene adjacency) | 5 | 0.489 | 0.487 | 0.502 | 0.493 | 0.493 | 0.497 |
| DeepSEM (sparse structural-equation VAE) | 5 | 0.506 | 0.505 | *0.540* | 0.507 | 0.508 | 0.503 |
| Pearson | 5 | **0.540** | **0.559** | **0.558** | **0.545** | **0.547** | *0.511* |
| GRNBoost2 | 5 | 0.501 | 0.512 | 0.511 | 0.511 | 0.510 | 0.501 |
| SCENIC+ (released) | 5 | 0.501 | 0.502 | 0.502 | 0.502 | 0.502 | 0.500 |
| LINGER (released) | 5 | 0.514 | 0.459 | 0.512 | 0.489 | 0.483 | **0.522** |

## 11. Conclusions and next steps

- To be written from the final numbers.
