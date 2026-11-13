# Baselines, ablations and hub-controlled evaluation on the scMultiomeGRN benchmark

Advisor meeting, work from 18 October to 14 November 2026.

Every scMultiomeGRN-benchmark number is an average over the 12 cell types (9 lung, 3 blood) and 10 folds, scored
with scMultiomeGRN's metric code on held-out pairs only. Brackets are 95% confidence intervals from resampling cell
types; p-values are paired Wilcoxon tests over cell type x fold. In each table **bold** marks the best value in a
column and *italics* the second best.

## 1. Summary of results

- About 50 simple and statistical models, the two GNN baselines of the scMultiomeGRN paper (DeepTFni, GENELink),
  six ablations of scMultiomeGRN and eleven ablations of MeVD-GRN were evaluated on scMultiomeGRN's datasets.
- No label-free statistical score recovers the ground-truth network: all score 0.50 to 0.54 AUROC on the official
  test.
- Counting each TF's training links is enough to beat scMultiomeGRN on the official test. A logistic regression on
  the two counts scores 0.931 AUPR against 0.919 for
  scMultiomeGRN and 0.935 for MeVD-GRN.
- On tests that remove the advantage of hub TFs, every method falls towards chance. On the caliper-matched test,
  MeVD-GRN scores 0.560 AUROC and scMultiomeGRN 0.524.
- On a new split that holds out whole TFs, the link-count rule is the best method when one TF of the pair is unseen
  (0.779 AUROC against 0.753 for MeVD-GRN), and
  every model is near chance once hubness is controlled (0.487 for MeVD-GRN).
- Conclusion: on this benchmark, the scores of all methods, MeVD-GRN included, are explained by hub TFs.
- BEAR-GRN (ChIP-seq and knockout ground truths, interim results): MeVD-GRN beats LINGER on most datasets, but a
  baseline that ranks target genes by their number of regulators is stronger still on ChIP (section 11).
- Correction: an earlier draft reported 0.878 AUPR for MeVD-GRN on the degree-matched test. That value came from an
  evaluation bug, fixed in this period; the corrected value is 0.779 (section 9).

## 2. Model families and evaluation protocol

- **Label-free statistical scores.** No labels: correlations, mutual information, GRNBoost2, GENIE3, ridge and lasso
  regression (the BEELINE and geneRNIB baselines).
- **Supervised models on pair features.** Logistic regression, linear regression and random forest on scMultiomeGRN's
  inputs (MAESTRO accessibility, GRNBoost2 profiles, expression histograms) and the Geneformer embedding.
- **Link-prediction heuristics on the training graph.** Degree, common neighbours, Adamic-Adar, Katz, random walks,
  low-rank completion and look-alike TFs, computed from the fold's training links only.
- **GNN baselines.** DeepTFni and GENELink, the baselines of the scMultiomeGRN paper, run from the copies shipped in
  the official repository on the same folds.
- **scMultiomeGRN ablations.** One part of the authors' code is switched off at a time; with every switch at its
  default the patched code reproduces the published code exactly.
- **MeVD-GRN ablations.** One component removed at a time from the full model, and each component added alone to a
  floor model that has none of them (3 seeds each).
- **Official test.** scMultiomeGRN's 10 folds split the links: each test link is paired with a random non-link.
- **Degree-matched test.** Each test link is paired with one of the three test non-links closest to it in the two
  TFs' training-link counts ("hidden sampling biases" preprint, December 2025).
- **Caliper-matched test.** Stricter: a test link is kept only if a test non-link has both TFs' log link counts
  within 0.1 of the link's, so the link-count rule scores 0.5 by construction.
- **Held-out-TF split (new).** The TFs, not the links, are split into 10 groups (Park and Marcotte 2012). C2 = one TF
  of the pair unseen in training, C3 = both unseen; a C2 caliper-matched test pairs each link of a hidden TF with a
  non-link of the same TF whose partner has the same link count.
- **Metrics.** Hub-controlled columns use AUROC, because scMultiomeGRN's trapezoidal AUPR gives a method with all
  scores tied (1 + share of positives) / 2, about 0.75, instead of 0.5.
- Nothing from the test fold is used for fitting; every model is checked for zero test-to-train leakage.

## 3. Label-free statistical scores

- These methods ask whether two TFs behave alike across cells, without any label.
- All of them sit between 0.50 and 0.54 AUROC on every test, which is chance.
- A link in the ground-truth network means one TF's motif lies in the other TF's open promoter; that says little
  about co-expression, so these scores cannot recover it.

| Method | AUPR | AUROC | AUROC, degree-matched | AUROC, caliper-matched |
|---|---|---|---|---|
| MeVD-GRN | **0.935** | **0.924** | **0.708** | **0.560** |
| scMultiomeGRN (official code) | *0.919* | *0.913* | *0.651* | 0.524 |
| Pearson correlation of expression | 0.519 | 0.537 | 0.513 | 0.521 |
| Spearman correlation of expression | 0.518 | 0.533 | 0.510 | 0.516 |
| Partial correlation of expression | 0.508 | 0.506 | 0.503 | 0.504 |
| Mutual information of expression | 0.518 | 0.537 | 0.519 | 0.520 |
| CLR on mutual information | 0.497 | 0.505 | 0.498 | 0.514 |
| GRNBoost2 between the two TFs | 0.508 | 0.511 | 0.505 | *0.526* |
| GENIE3 between the two TFs | 0.504 | 0.520 | 0.508 | 0.516 |
| Ridge regression | 0.505 | 0.500 | 0.493 | 0.494 |
| Lasso regression | 0.505 | 0.499 | 0.495 | 0.494 |
| Correlation of MAESTRO accessibility | 0.528 | 0.520 | 0.519 | 0.511 |
| Similarity of Geneformer embeddings | 0.529 | 0.529 | 0.525 | 0.516 |

## 4. Supervised models on pair features

- Learned from the training links, scMultiomeGRN's inputs become informative on the official test: GRNBoost2 alone
  reaches 0.729 AUPR, the two inputs side by side 0.788,
  a random forest 0.836.
- Giving the inputs side by side works better than multiplying them.
- Under the caliper-matched test all of them fall to between 0.49 and 0.52: what they learn is mostly which TFs are
  hubs.

| Method | AUPR | AUROC | AUROC, degree-matched | AUROC, caliper-matched |
|---|---|---|---|---|
| MeVD-GRN | **0.935** | **0.924** | **0.708** | **0.560** |
| scMultiomeGRN (official code) | *0.919* | *0.913* | *0.651* | *0.524* |
| MAESTRO accessibility only | 0.699 | 0.674 | 0.545 | 0.487 |
| GRNBoost2 profiles only | 0.729 | 0.702 | 0.560 | 0.510 |
| Geneformer embedding only | 0.695 | 0.695 | 0.526 | 0.503 |
| Joint expression histograms only | 0.634 | 0.641 | 0.538 | 0.521 |
| Promoter openness only | 0.586 | 0.592 | 0.521 | 0.506 |
| MAESTRO and GRNBoost2 side by side | 0.788 | 0.765 | 0.594 | 0.514 |
| MAESTRO multiplied by GRNBoost2 | 0.721 | 0.695 | 0.545 | 0.487 |
| Product of the two TFs' features only | 0.523 | 0.509 | 0.490 | 0.487 |
| Linear regression, side by side | 0.788 | 0.763 | 0.596 | 0.514 |
| Random forest, side by side | 0.836 | 0.814 | 0.606 | 0.502 |

## 5. Link-prediction heuristics on the training graph

- The product of the two link counts scores 0.924 AUPR and a logistic regression on
  the two counts 0.931; both beat scMultiomeGRN (0.919).
- Look-alike TFs (are TFs with similar training links linked to the partner?) score
  0.929 and keep 0.688 AUROC on the degree-matched
  test, close to MeVD-GRN (0.708): the family structure of the training links explains
  most of what link counts do not.
- On the caliper-matched test every heuristic is between 0.48 and 0.57.

| Method | AUPR | AUROC | AUROC, degree-matched | AUROC, caliper-matched |
|---|---|---|---|---|
| MeVD-GRN | **0.935** | **0.924** | **0.708** | *0.560* |
| scMultiomeGRN (official code) | 0.919 | 0.913 | 0.651 | 0.524 |
| Product of the two link counts | 0.924 | 0.912 | 0.653 | 0.506 |
| Sum of the two link counts | 0.928 | 0.917 | 0.677 | 0.508 |
| Common neighbours | 0.781 | 0.757 | 0.586 | 0.488 |
| Adamic-Adar index | 0.806 | 0.782 | 0.585 | 0.485 |
| Resource allocation index | 0.832 | 0.808 | 0.578 | 0.483 |
| Katz index | 0.887 | 0.866 | 0.642 | 0.509 |
| Random walk with restart | 0.924 | 0.913 | 0.670 | 0.508 |
| Low-rank SVD of the training links | 0.717 | 0.675 | 0.562 | **0.567** |
| Look-alike TFs (matrix completion) | 0.929 | 0.914 | *0.688* | 0.533 |
| Logistic regression on the link counts | 0.931 | *0.919* | 0.671 | 0.507 |
| Logistic regression on counts and neighbour scores | *0.932* | *0.919* | 0.682 | 0.514 |
| Link counts + MAESTRO + GRNBoost2 | 0.926 | 0.915 | 0.667 | 0.520 |
| Counts and neighbour scores + Geneformer | 0.931 | *0.919* | 0.681 | 0.534 |

## 6. GNN baselines of the scMultiomeGRN paper

- DeepTFni scores 0.882 AUPR and GENELink 0.902 on the official test,
  both below scMultiomeGRN, consistent with the scMultiomeGRN paper.
- Both are close to chance on the caliper-matched test.

| Method | AUPR | AUROC | AUROC, degree-matched | AUROC, caliper-matched |
|---|---|---|---|---|
| MeVD-GRN | **0.935** | **0.924** | **0.708** | **0.560** |
| scMultiomeGRN (official code) | *0.919* | *0.913* | *0.651* | 0.524 |
| DeepTFni | 0.882 | 0.871 | 0.647 | 0.517 |
| GENELink | 0.902 | 0.901 | 0.627 | *0.537* |

## 7. scMultiomeGRN ablations

- The known-link graph is what matters: without it AUPR drops from 0.919 to
  0.826; a STRING network in its place gives
  0.781.
- The pair histograms, either single modality and an added Geneformer view change nothing measurable.
- On the hub-controlled tests the ablations are as close to chance as the full model.

| Method | AUPR | AUROC | AUROC, degree-matched | AUROC, caliper-matched |
|---|---|---|---|---|
| MeVD-GRN | **0.935** | **0.924** | **0.708** | **0.560** |
| scMultiomeGRN (official code) | 0.919 | 0.913 | *0.651* | 0.524 |
| Without the TF-TF graph | 0.826 | 0.815 | 0.600 | 0.501 |
| STRING network instead of the known links | 0.781 | 0.770 | 0.580 | 0.496 |
| Without the pair histograms | *0.920* | *0.914* | *0.651* | 0.523 |
| ATAC features only | 0.916 | 0.912 | 0.649 | 0.522 |
| RNA features only | 0.916 | 0.910 | 0.650 | 0.511 |
| With Geneformer added | 0.917 | 0.911 | 0.649 | *0.527* |

## 8. MeVD-GRN ablations

- Official test: Geneformer is the component that matters. The floor model plus Geneformer alone reaches
  0.931 AUPR (full model 0.935); removing Geneformer gives
  0.922, and the smaller base size 0.925.
- Removing the ATAC tower, the graph network, either graph type or the training-link relation changes AUPR by at most
  0.001.
- Caliper-matched test, official split: every single removal stays between 0.547 and
  0.562 AUROC (full model 0.560), while the floor model is at
  0.501: the small signal beyond link counts comes from the
  components together, and no single one is needed for it.
- Held-out-TF split, C2 caliper-matched: the full model and every variant are between 0.482
  and 0.508 AUROC, i.e. chance.

| Method | AUPR | AUROC | AUROC, degree-matched | AUROC, caliper-matched |
|---|---|---|---|---|
| MeVD-GRN | **0.935** | **0.924** | **0.708** | *0.560* |
| Without Geneformer | 0.922 | 0.910 | 0.686 | 0.559 |
| Without the ATAC tower | **0.935** | *0.923* | *0.707* | 0.557 |
| Without the graph network | *0.934* | *0.923* | 0.705 | **0.562** |
| Co-expression graph only | *0.934* | 0.922 | 0.704 | 0.552 |
| TF-candidate graph only | *0.934* | *0.923* | 0.701 | 0.553 |
| Base size (smaller model) | 0.925 | 0.912 | 0.691 | *0.560* |
| Without the training-link relation | *0.934* | *0.923* | 0.702 | 0.547 |
| Floor model (none of the four parts) | 0.633 | 0.635 | 0.533 | 0.501 |
| Floor model + ATAC tower | 0.833 | 0.822 | 0.619 | 0.522 |
| Floor model + Geneformer | 0.931 | 0.921 | 0.697 | 0.552 |
| Floor model + graph network | 0.875 | 0.862 | 0.642 | 0.528 |

Held-out-TF split:

| Method | AUROC, C2 | AUROC, C3 | AUROC, C2 caliper-matched |
|---|---|---|---|
| MeVD-GRN | 0.753 | 0.529 | 0.487 |
| Without Geneformer | 0.716 | *0.542* | 0.492 |
| Without the ATAC tower | 0.750 | 0.528 | 0.482 |
| Without the graph network | **0.761** | *0.542* | 0.493 |
| Co-expression graph only | 0.754 | 0.517 | 0.494 |
| TF-candidate graph only | *0.756* | 0.531 | 0.492 |
| Base size (smaller model) | 0.718 | 0.526 | 0.496 |
| Without the training-link relation | 0.750 | 0.534 | 0.493 |
| Floor model (none of the four parts) | 0.576 | *0.542* | 0.499 |
| Floor model + ATAC tower | 0.643 | **0.549** | **0.508** |
| Floor model + Geneformer | **0.761** | 0.497 | 0.491 |
| Floor model + graph network | 0.640 | 0.529 | *0.500* |

## 9. Hub-controlled tests on the official split, and a correction

- Degree-matched test: MeVD-GRN 0.708 AUROC [0.682, 0.738],
  scMultiomeGRN 0.651, best simple model 0.688, from Look-alike TFs (matrix completion).
  MeVD-GRN is best, but only 0.020 above the best simple model.
- The link-count rule still scores 0.653 on the degree-matched test,
  because non-links between two hubs are rare in these dense graphs. The caliper-matched test removes that.
- Caliper-matched test: MeVD-GRN 0.560 [0.541, 0.577],
  scMultiomeGRN 0.524, best other method 0.567, from Low-rank SVD of the training links.
- MeVD-GRN's interval lies above 0.5, so a small signal beyond link counts remains on the official split, of the
  same size as low-rank SVD's (0.567 [0.548, 0.585])
  and GENELink's (0.537 [0.519, 0.560]). These sets are
  small (a few pairs per fold in the densest cell types).
- **Correction.** The degree-matched value of 0.878 AUPR reported for MeVD-GRN in an earlier draft was wrong. The
  training script reduced each fold to one row per TF pair before scoring, and the matched tests counted a TF's links
  from the first column only, so MeVD-GRN was scored on different pairs from every other method. The fix counts links
  on unique pairs; every method now uses identical test pairs on every fold (checked on all 55,000 fold results), and
  the corrected value is 0.779 AUPR.

## 10. Held-out-TF split

- With one TF unseen (C2), the link-count rule is the best method: it still knows the partner TF's link count.
  MeVD-GRN scores 0.753 [0.741, 0.765], scMultiomeGRN
  0.744, the link-count rule 0.779.
- With both TFs unseen (C3), every supervised model is at chance; the best score is 0.575,
  from Joint expression histograms only. C3 sets are small (4 to 102 links per fold).
- On the C2 caliper-matched test no method exceeds 0.52: MeVD-GRN 0.487
  [0.481, 0.494], best method 0.520 [0.510, 0.530],
  from Mutual information of expression.
- Neighbourhood heuristics and look-alike TFs score exactly 0.5 here, because a hidden TF has no training links.

| Method | AUROC, C2 | AUROC, C3 | AUROC, C2 caliper-matched |
|---|---|---|---|
| MeVD-GRN | 0.753 | 0.529 | 0.487 |
| scMultiomeGRN (official code) | 0.744 | 0.526 | 0.483 |
| DeepTFni | 0.576 | 0.483 | 0.501 |
| GENELink | 0.661 | 0.518 | 0.492 |
| Product of the two link counts | **0.779** | 0.500 | 0.509 |
| Logistic regression on the link counts | **0.779** | 0.500 | 0.509 |
| Look-alike TFs (matrix completion) | 0.500 | 0.500 | 0.500 |
| Geneformer look-alikes linked to the partner | 0.731 | 0.500 | 0.507 |
| Mutual information of expression | 0.536 | *0.573* | **0.520** |
| GRNBoost2 between the two TFs | 0.509 | 0.540 | *0.517* |
| Joint expression histograms only | 0.593 | **0.575** | 0.502 |
| Random forest, MAESTRO and GRNBoost2 | 0.695 | 0.535 | *0.517* |
| Link counts + Geneformer | *0.774* | 0.539 | 0.504 |

## 11. BEAR-GRN: ChIP-seq and knockout ground truths

- BEAR-GRN (Nat Methods 2025) scores TF-to-gene networks against ChIP-seq, knockout (KO), union and intersection
  ground truths on 9 datasets (human K562, macrophage, iPS; mouse embryo and ESC). None of them comes from a motif rule.
- MeVD-GRN follows the pre-registered protocol (section 12 of `docs/benchmarks/bear_grn.md`): TF-disjoint 5-fold
  cross-fitting, so every scored TF is unseen in training; the configuration was fixed on K562 before any other
  dataset was run. Released outputs of the BEAR methods are re-scored with the same code.
- **Interim:** 9 of 9 datasets, 2-3 seeds per dataset (5 seeds and the strict L2 regime
  are still running).
- Over the four ground truths (52 dataset x metric cells), MeVD-GRN beats LINGER in 33, every released method in 21, and the target
  in-degree baseline in 23. The pre-registered win (above LINGER, in-degree and coverage) holds in 16.
- On ChIP, MeVD-GRN beats LINGER in AUPRC on 9 of 9 datasets, but the in-degree
  baseline is above MeVD-GRN on 7 of 9 (AUROC) and also above LINGER on 9 of 9 (AUPRC).
- So BEAR-GRN shows the same pattern as scMultiomeGRN: ranking target genes by how many TFs regulate them beats
  the published methods, and MeVD-GRN does not yet clear that baseline on ChIP.

ChIP-seq ground truth (in-degree = mean of 5 seeds; best of the three per metric in bold):

| Dataset | AUROC MeVD-GRN | AUROC LINGER | AUROC in-degree | AUPRC MeVD-GRN | AUPRC LINGER | AUPRC in-degree | AUPRC random |
|---|---|---|---|---|---|---|---|
| K562 | 0.603 | 0.536 | **0.652** | 0.473 | 0.430 | **0.484** | 0.331 |
| Macrophage_S1 | **0.625** | 0.596 | 0.609 | **0.373** | 0.321 | 0.365 | 0.272 |
| Macrophage_S2 | 0.588 | 0.548 | **0.610** | 0.367 | 0.324 | **0.374** | 0.271 |
| Naive_mESC | 0.541 | 0.476 | **0.598** | 0.285 | 0.160 | **0.308** | 0.153 |
| iPS | **0.617** | 0.496 | 0.604 | **0.156** | 0.108 | 0.150 | 0.091 |
| mESC_E7.5_rep1 | 0.549 | **0.652** | 0.604 | 0.296 | 0.241 | **0.323** | 0.155 |
| mESC_E7.5_rep2 | 0.563 | **0.644** | 0.596 | 0.294 | 0.225 | **0.306** | 0.153 |
| mESC_E8.5_rep1 | 0.544 | **0.658** | 0.607 | 0.283 | 0.227 | **0.312** | 0.154 |
| mESC_E8.5_rep2 | 0.539 | **0.656** | 0.609 | 0.277 | 0.224 | **0.313** | 0.155 |

Hub control (pre-registered, seed 42): AUROC when MeVD-GRN's scores are ranked only within groups of target
genes with the same in-degree (0 alone, then 20 quantile bins), so the in-degree ranking scores 0.5.
- ChIP: a small signal beyond in-degree remains on the human datasets (0.530 to 0.573); on the mouse datasets it is 0.475 to 0.542, above 0.52 only on mESC_E7.5_rep2 (1 of 5).
- KO: 0.516 to 0.657 on 5 datasets, the clearest signal beyond in-degree so far.

| Dataset | ChIP AUROC | ChIP AUROC within in-degree groups | KO AUROC | KO AUROC within in-degree groups |
|---|---|---|---|---|
| K562 | 0.610 | 0.530 | 0.559 | 0.516 |
| Macrophage_S1 | 0.630 | 0.573 | – | – |
| Macrophage_S2 | 0.575 | 0.530 | – | – |
| Naive_mESC | 0.539 | 0.498 | – | – |
| iPS | 0.551 | 0.541 | – | – |
| mESC_E7.5_rep1 | 0.561 | 0.509 | 0.487 | 0.524 |
| mESC_E7.5_rep2 | 0.586 | 0.542 | 0.547 | 0.611 |
| mESC_E8.5_rep1 | 0.532 | 0.477 | 0.522 | 0.565 |
| mESC_E8.5_rep2 | 0.536 | 0.475 | 0.597 | 0.657 |

## 12. Conclusions and next steps

- The scMultiomeGRN benchmark measures how well a method knows which TFs are hubs. Its ground-truth network comes
  from a motif rule on the same ATAC data the models read, and its random non-links are easy to separate by link count.
- Simple link-count models match or beat the published GNNs; label-free statistics are at chance.
- Once hubness is controlled (caliper matching, held-out TFs), no method, MeVD-GRN included, keeps a meaningful signal.
- BEAR-GRN shows the same effect on target genes: an in-degree ranking beats every published method on ChIP AUPRC,
  and MeVD-GRN beats LINGER but not yet that ranking. Within groups of genes with the same in-degree, MeVD-GRN
  keeps a small signal on human ChIP and on knockout ground truths, and little on mouse ChIP.
- These results are a finding about how GRN benchmarks are scored and belong in the paper as an evaluation
  contribution, alongside MeVD-GRN.
- Next: complete BEAR-GRN (5 seeds, the strict L2 regime, ablations) with a hub control against the in-degree
  baseline (within-in-degree-stratum scores); then PBMC10k with LINGER's own candidate genes.
