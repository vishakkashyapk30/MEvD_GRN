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
- **Floor model.** MeVD-GRN with the four components removed: no Geneformer embedding, no ATAC tower, no graph
  network and no training-link relation, at the smaller base size (128 wide). Each TF is described only by three
  numbers from its RNA (mean expression, variance, share of cells that detect it), passed through a small MLP; a pair
  is scored by the symmetric decoder plus a co-expression term (correlation of the two TFs' expression signatures).
  It is trained on the fold's training links like the full model.
- **Official test** [1]. scMultiomeGRN's 10 folds split the links: each test link is paired with a random non-link.
- **Degree-matched test** (own variant of degree-aware negatives [2, 3]). Each test link is paired with one of the
  three test non-links closest to it in the two TFs' training-link counts.
- **Caliper-matched test** (own construction: caliper matching [4, 5] applied to link counts). Stricter: a test link
  is kept only if a test non-link has both TFs' log link counts within 0.1 of the link's, so the link-count rule
  scores 0.5 by construction.
- **Held-out-TF split (new)** [6]. The TFs, not the links, are split into 10 groups. Test pairs fall in
  two classes: **one TF unseen** in training (Park and Marcotte's "C2") and **both TFs unseen** ("C3"). The
  one-TF-unseen caliper-matched test (own construction) pairs each link of a hidden TF with a
  non-link of the same TF whose partner has the same link count.
- **Metrics.** AUROC; AUPR with scMultiomeGRN's trapezoidal code [1]; average precision [7]. Hub-controlled columns
  use AUROC, because trapezoidal AUPR gives a method with all scores tied (1 + share of positives) / 2, about 0.75,
  instead of 0.5 (linear interpolation in precision-recall space overestimates the area [7]).
- **Statistics** (standard methods). 95% CIs by bootstrap over cell types; paired Wilcoxon signed-rank tests over
  cell type x fold. Sources of every test: slide 13.
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

- Floor model = MeVD-GRN without Geneformer, the ATAC tower, the graph network and the training-link relation: each
  TF is reduced to three RNA numbers (mean, variance, detection rate) plus a co-expression term (definition on slide 2).
- Official test: Geneformer is the component that matters. The floor model plus Geneformer alone reaches
  0.931 AUPR (full model 0.935); removing Geneformer gives
  0.922, and the smaller base size 0.925.
- Removing the ATAC tower, the graph network, either graph type or the training-link relation changes AUPR by at most
  0.001.
- Caliper-matched test, official split: every single removal stays between 0.547 and
  0.562 AUROC (full model 0.560), while the floor model is at
  0.501: the small signal beyond link counts comes from the
  components together, and no single one is needed for it.
- Held-out-TF split, one TF unseen, caliper-matched: the full model and every variant are between 0.482
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

| Method | AUROC, one TF unseen | AUROC, both TFs unseen | AUROC, one TF unseen, caliper-matched |
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

- With one TF unseen, the link-count rule is the best method: it still knows the partner TF's link count.
  MeVD-GRN scores 0.753 [0.741, 0.765], scMultiomeGRN
  0.744, the link-count rule 0.779.
- With both TFs unseen, every supervised model is at chance; the best score is 0.575,
  from Joint expression histograms only. These sets are small (4 to 102 links per fold).
- On the one-TF-unseen caliper-matched test no method exceeds 0.52: MeVD-GRN 0.487
  [0.481, 0.494], best method 0.520 [0.510, 0.530],
  from Mutual information of expression.
- Neighbourhood heuristics and look-alike TFs score exactly 0.5 here, because a hidden TF has no training links.

| Method | AUROC, one TF unseen | AUROC, both TFs unseen | AUROC, one TF unseen, caliper-matched |
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

- BEAR-GRN (Nat Commun 2026) [8] scores TF-to-gene networks against ChIP-seq, knockout (KO), union and intersection
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

**Baselines on BEAR-GRN** (9 released methods + 4 simple baselines, all on the 9 datasets):
- **Released methods** (CellOracle, DIRECT-NET, FigR, GRaNIE, LINGER,
  Pando-GLM, Pando-XGB, SCENIC+, TRIPOD): the networks published with BEAR-GRN, re-scored with its code.
- **In-degree (supervised hub prior).** Score of an edge TF -> gene = how many *training* TFs (same TF-disjoint
  folds as MeVD-GRN) have that gene as a target in the ground truth. Every TF gets the same ranking of genes; no cell
  data is used. It asks only "is this gene a target of many TFs?" (Stock et al. [2] call it the in-degree sorter).
- **Pearson.** |correlation| of the TF's and the gene's expression across cells. No labels.
- **GRNBoost2.** GRNBoost2 importance of the TF for the gene, from expression only. No labels.
- **All measured pairs.** Every TF-gene pair present in the data gets the same score: no ranking at all. BEAR-GRN
  scores pairs a method leaves out as 0, so this measures how much AUPRC comes from just emitting a dense network.
- **Random.** Expected AUPRC of a random ranking (AUROC 0.5).
- On ChIP AUPRC, Pearson beats LINGER on 9 of 9 datasets and "all measured pairs" on 9 of 9.

ChIP-seq ground truth, AUPRC (MeVD-GRN = mean over seeds; in-degree = mean of 5 fold seeds; bold = best and italics = second best in each dataset):

| Model | K562 | Macrophage_S1 | Macrophage_S2 | Naive_mESC | iPS | mESC_E7.5_rep1 | mESC_E7.5_rep2 | mESC_E8.5_rep1 | mESC_E8.5_rep2 |
|---|---|---|---|---|---|---|---|---|---|
| MeVD-GRN | *0.473* | **0.373** | *0.367* | *0.285* | **0.156** | *0.296* | *0.294* | *0.283* | *0.277* |
| LINGER | 0.430 | 0.321 | 0.324 | 0.160 | 0.108 | 0.241 | 0.225 | 0.227 | 0.224 |
| CellOracle | 0.390 | 0.293 | 0.300 | 0.163 | 0.110 | 0.191 | 0.182 | 0.185 | 0.188 |
| DIRECT-NET | – | – | – | – | – | – | – | – | – |
| FigR | 0.376 | 0.291 | 0.296 | 0.171 | 0.094 | 0.178 | 0.172 | 0.173 | 0.171 |
| GRaNIE | 0.368 | 0.269 | 0.268 | 0.182 | 0.099 | 0.188 | 0.179 | 0.184 | 0.184 |
| Pando-GLM | 0.348 | 0.289 | 0.275 | 0.164 | 0.091 | 0.159 | 0.155 | 0.157 | 0.160 |
| Pando-XGB | 0.346 | 0.289 | 0.280 | 0.165 | 0.091 | 0.159 | 0.155 | 0.158 | 0.161 |
| SCENIC+ | 0.332 | 0.285 | 0.278 | 0.153 | 0.091 | 0.156 | 0.155 | 0.155 | 0.156 |
| TRIPOD | 0.405 | 0.332 | 0.336 | 0.192 | 0.119 | 0.193 | 0.186 | 0.181 | 0.183 |
| In-degree | **0.484** | *0.365* | **0.374** | **0.308** | *0.150* | **0.323** | **0.306** | **0.312** | **0.313** |
| Pearson | 0.434 | 0.343 | 0.351 | 0.266 | 0.139 | 0.273 | 0.267 | 0.258 | 0.259 |
| GRNBoost2 | 0.380 | 0.362 | 0.361 | 0.203 | 0.131 | 0.215 | 0.210 | 0.209 | 0.205 |
| All measured pairs | 0.431 | 0.331 | 0.342 | 0.269 | 0.137 | 0.274 | 0.268 | 0.265 | 0.263 |
| *Random* | 0.331 | 0.272 | 0.271 | 0.153 | 0.091 | 0.155 | 0.153 | 0.154 | 0.155 |

DIRECT-NET's released network has no edge scores, so it has no AUROC or AUPRC (BEAR-GRN reports only its
top-10k precision).

Same, AUROC ("all measured pairs" ties every pair, so it is 0.5 by construction):

| Model | K562 | Macrophage_S1 | Macrophage_S2 | Naive_mESC | iPS | mESC_E7.5_rep1 | mESC_E7.5_rep2 | mESC_E8.5_rep1 | mESC_E8.5_rep2 |
|---|---|---|---|---|---|---|---|---|---|
| MeVD-GRN | *0.603* | **0.625** | 0.588 | 0.541 | **0.617** | 0.549 | 0.563 | 0.544 | 0.539 |
| LINGER | 0.536 | 0.596 | 0.548 | 0.476 | 0.496 | **0.652** | **0.644** | **0.658** | **0.656** |
| CellOracle | 0.539 | 0.588 | *0.605* | 0.508 | 0.478 | 0.529 | 0.503 | 0.481 | 0.511 |
| DIRECT-NET | – | – | – | – | – | – | – | – | – |
| FigR | 0.504 | 0.529 | 0.529 | 0.511 | 0.468 | 0.509 | 0.500 | 0.496 | 0.490 |
| GRaNIE | 0.512 | 0.571 | 0.462 | 0.539 | 0.568 | 0.570 | 0.594 | 0.588 | *0.621* |
| Pando-GLM | 0.576 | 0.575 | 0.500 | 0.454 | 0.506 | 0.455 | 0.394 | 0.382 | 0.401 |
| Pando-XGB | 0.518 | 0.557 | 0.535 | 0.534 | 0.437 | 0.508 | 0.515 | 0.516 | 0.532 |
| SCENIC+ | 0.465 | 0.599 | 0.513 | *0.556* | 0.478 | 0.409 | 0.461 | 0.386 | 0.427 |
| TRIPOD | 0.511 | 0.553 | 0.539 | 0.516 | 0.563 | 0.492 | 0.515 | 0.494 | 0.525 |
| In-degree | **0.652** | *0.609* | **0.610** | **0.598** | *0.604* | *0.604* | *0.596* | *0.607* | 0.609 |
| Pearson | 0.506 | 0.523 | 0.517 | 0.494 | 0.531 | 0.498 | 0.496 | 0.487 | 0.492 |
| GRNBoost2 | 0.496 | 0.567 | 0.562 | 0.484 | 0.465 | 0.511 | 0.506 | 0.504 | 0.497 |
| All measured pairs | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 |
| *Random* | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 |

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

## 13. Sources of the evaluation tests

| Test or metric | Source | Published, or own construction |
|---|---|---|
| Official test: random non-links, 10 folds, AUROC and trapezoidal AUPR | scMultiomeGRN [1] | Published (Nucleic Acids Research) |
| Degree-matched test | Idea: degree-preserving negatives [3], degree-aware negatives [2]; nearest-3 matching on training-link counts is own | Idea published (Bioinformatics); closest GRN version is a preprint; this variant is own |
| Caliper-matched test (official split, and held-out-TF with one TF unseen) | Caliper matching from causal inference [4, 5], applied here to TF link counts | Own construction; no published link-prediction paper using it was found |
| Held-out-TF split, one-TF-unseen and both-TFs-unseen classes (C2 and C3 in the paper) | Park and Marcotte [6] | Published (Nature Methods) |
| Average precision; overestimate of trapezoidal AUPR | Davis and Goadrich [7] | Published (ICML) |
| BEAR-GRN: AUROC and AUPRC on ChIP, KO, union and intersection ground truths, random-baseline AUPRC | BEAR-GRN [8], scored with its own code (ported, checked against the R code) | Published (Nature Communications) |
| BEAR-GRN in-degree baseline | In-degree sorter of [2] | Preprint |
| BEAR-GRN hub control: AUROC within in-degree groups | Stratified evaluation, pre-registered in `docs/benchmarks/bear_grn.md` | Own construction |
| Degree as a shortcut in pair-input prediction | [3]; Chatterjee et al. [9]; hard negatives for GNN link prediction, Li et al. [10] | Published (Bioinformatics, Nat Commun, NeurIPS) |

1. Xu J. et al. Deep learning-based cell-specific gene regulatory networks inferred from single-cell multiome data.
   *Nucleic Acids Research* 53(5), gkaf138 (2025). https://doi.org/10.1093/nar/gkaf138
2. Stock M. et al. Hidden sampling biases inflate performance in gene regulatory network inference. *bioRxiv*
   (December 2025, preprint, not peer reviewed). https://doi.org/10.64898/2025.12.19.695616
3. Yu J. et al. Simple sequence-based kernels do not predict protein-protein interactions. *Bioinformatics* 26(20),
   2610-2614 (2010). https://doi.org/10.1093/bioinformatics/btq483
4. Rosenbaum P. R., Rubin D. B. Constructing a control group using multivariate matched sampling methods that
   incorporate the propensity score. *The American Statistician* 39(1), 33-38 (1985).
5. Austin P. C. Optimal caliper widths for propensity-score matching when estimating differences in means and
   differences in proportions in observational studies. *Pharmaceutical Statistics* 10(2), 150-161 (2011).
   https://doi.org/10.1002/pst.433
6. Park Y., Marcotte E. M. Flaws in evaluation schemes for pair-input computational predictions. *Nature Methods* 9,
   1134-1136 (2012). https://doi.org/10.1038/nmeth.2259
7. Davis J., Goadrich M. The relationship between precision-recall and ROC curves. *ICML* (2006).
8. Karamveer et al. (Uzun lab). BEAR-GRN. *Nature Communications* (2026). https://doi.org/10.1038/s41467-026-77838-w
9. Chatterjee A. et al. Improving the generalizability of protein-ligand binding predictions with AI-Bind.
   *Nature Communications* 14, 1989 (2023). https://doi.org/10.1038/s41467-023-37572-z
10. Li J. et al. Evaluating graph neural networks for link prediction: current pitfalls and new benchmarking.
    *NeurIPS* (2023). https://arxiv.org/abs/2306.10453
