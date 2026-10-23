# Simple baselines and scMultiomeGRN ablations on the scMultiomeGRN benchmark

Advisor meeting, week of 18 to 24 October 2026.

Every number is an average over the 12 cell types (9 lung, 3 blood) and the 10 official folds, scored on the held-out test links with scMultiomeGRN's own metric code. In each table, **bold** marks the best value in a column and *italics* the second best.

## 1. Summary: label-free and supervised baselines, scMultiomeGRN ablations, degree-controlled evaluation

- 38 simple and statistical models were evaluated on scMultiomeGRN's datasets and splits, and six ablations of scMultiomeGRN itself were run with the authors' code.
- None of the classic network-inference scores from the cell data finds the answer key. Correlation, mutual information, GRNBoost2, GENIE3 and regression all score between 0.50 and 0.53 AUPR, which is chance on this test.
- Counting each TF's known training links is enough to beat scMultiomeGRN. A logistic regression on the two counts scores 0.931 against 0.918, and MeVD-GRN scores 0.935.
- scMultiomeGRN depends on the known links it passes messages over. Without them it drops to 0.826, and with a STRING network in their place to 0.781. Adding Geneformer to it changes nothing.
- Every model was also scored on a harder test where each true link is paired with a non-link whose TFs have similar numbers of links.
- On that test, scMultiomeGRN and every simple model that counts links fall to between 0.70 and 0.76 AUPR, close to the counting rule at 0.724. MeVD-GRN keeps 0.878.

## 2. Baseline families and evaluation protocol

- **Statistical scores from the cell data.** These use no labels at all. They are the methods from the BEELINE and geneRNIB benchmarks: correlations, mutual information, GRNBoost2, GENIE3, and ridge and lasso regression.
- **Simple supervised models on features.** These are logistic regressions on scMultiomeGRN's own inputs, MAESTRO accessibility and GRNBoost2 profiles, plus the Geneformer embedding. The inputs are given one at a time, side by side, and multiplied together.
- **Counting rules on the known links.** These are the classic link-prediction heuristics: degree, common neighbours, Adamic-Adar, Katz and random walks. All of them use only the fold's training links.
- **scMultiomeGRN ablations.** One part of the authors' code is changed at a time, and with no change the patched code gives exactly their numbers.
- **A degree-matched test.** This follows a December 2025 preprint, which showed that hub TFs inflate scores in this kind of benchmark. Each test link is paired with the test non-link closest to it in how many training links its two TFs have, so knowing the hubs no longer helps much.
- Nothing from the test fold is used for fitting, and the regularisation of every logistic regression is chosen on the training pairs only.

## 3. Label-free co-expression, mutual-information and regression scores (BEELINE family)

- None of these methods sees a label. They ask whether two TFs behave alike across cells.
- All of them land near 0.5, which is what a random guess gets on a test with as many links as non-links.
- The reason is the answer key itself. A link means that one TF's motif sits in the other TF's open promoter, and that says little about whether the two TFs are co-expressed.
- GRNBoost2 and GENIE3, the standard tools for gene networks, do no better than plain correlation here.
- This matches the scMultiomeGRN paper, where GENIE3 also came last among its baselines.

| Method | AUPR | AUROC | AUPR, degree-matched | AUROC, degree-matched |
|---|---|---|---|---|
| MeVD-GRN | **0.935** | **0.924** | **0.878** | **0.857** |
| scMultiomeGRN (unchanged) | *0.918* | *0.912* | *0.705* | *0.649* |
| Pearson correlation of expression | 0.519 | 0.537 | 0.534 | 0.512 |
| Spearman correlation of expression | 0.518 | 0.533 | 0.534 | 0.510 |
| Partial correlation of expression | 0.508 | 0.506 | 0.530 | 0.504 |
| Mutual information of expression | 0.518 | 0.537 | 0.539 | 0.518 |
| CLR on mutual information | 0.497 | 0.505 | 0.528 | 0.497 |
| GRNBoost2 between the two TFs | 0.508 | 0.511 | 0.539 | 0.503 |
| GENIE3 between the two TFs | 0.504 | 0.520 | 0.533 | 0.506 |
| Ridge regression | 0.505 | 0.500 | 0.528 | 0.494 |
| Lasso regression | 0.505 | 0.499 | 0.527 | 0.495 |
| Correlation of MAESTRO accessibility | 0.528 | 0.520 | 0.561 | 0.519 |
| Similarity of Geneformer embeddings | 0.529 | 0.529 | 0.561 | 0.524 |

## 4. Supervised logistic regression, linear regression and random forest on symmetric pair features

- Once a logistic regression can learn from the training links, the same inputs become useful. MAESTRO alone reaches 0.699 AUPR and GRNBoost2 alone 0.729.
- Giving the two inputs side by side works better than multiplying them. Side by side reaches 0.788, the per-TF product 0.721, and the product of the two TFs' features only 0.523.
- Geneformer on its own reaches 0.695, below the cell data, even though it knows nothing about these cells.
- A plain linear regression on the 0/1 label gives the same as the logistic regression. A random forest is the best of this group at 0.836.
- All of these models stay well below scMultiomeGRN, so the cell data alone does not explain its score.

| Method | AUPR | AUROC | AUPR, degree-matched | AUROC, degree-matched |
|---|---|---|---|---|
| MeVD-GRN | **0.935** | **0.924** | **0.878** | **0.857** |
| scMultiomeGRN (unchanged) | *0.918* | *0.912* | *0.705* | *0.649* |
| MAESTRO accessibility only | 0.699 | 0.674 | 0.600 | 0.545 |
| GRNBoost2 profiles only | 0.729 | 0.702 | 0.612 | 0.560 |
| Geneformer embedding only | 0.695 | 0.695 | 0.564 | 0.526 |
| Joint expression histograms only | 0.634 | 0.641 | 0.579 | 0.538 |
| Promoter openness only | 0.586 | 0.592 | 0.560 | 0.521 |
| MAESTRO and GRNBoost2 side by side | 0.788 | 0.765 | 0.642 | 0.593 |
| MAESTRO multiplied by GRNBoost2 | 0.721 | 0.695 | 0.605 | 0.545 |
| Product of the two TFs' features only | 0.523 | 0.509 | 0.534 | 0.490 |
| Linear regression, side by side | 0.788 | 0.763 | 0.644 | 0.596 |
| Random forest, side by side | 0.836 | 0.814 | 0.650 | 0.606 |

## 5. Topological link-prediction heuristics on the training graph

- These rules use only the links in the fold's training part, never the cell data.
- The simplest one wins. The product of the two TFs' link counts scores 0.924, and a logistic regression on the two counts 0.931. Both beat scMultiomeGRN.
- Rules based on shared neighbours do worse. Common neighbours scores 0.781, and the Jaccard index ranks pairs almost backwards.
- Adding the cell data or Geneformer to the counts does not help. Every combination stays between 0.926 and 0.931, no better than the counts alone.
- This is the pattern the recent literature warns about: on this benchmark, knowing which TFs are hubs is worth more than any feature.

| Method | AUPR | AUROC | AUPR, degree-matched | AUROC, degree-matched |
|---|---|---|---|---|
| MeVD-GRN | **0.935** | **0.924** | **0.878** | **0.857** |
| scMultiomeGRN (unchanged) | 0.918 | 0.912 | 0.705 | 0.649 |
| Product of the two link counts | 0.924 | 0.912 | 0.724 | 0.653 |
| Sum of the two link counts | 0.928 | 0.917 | 0.750 | 0.677 |
| Common neighbours | 0.781 | 0.757 | 0.662 | 0.587 |
| Adamic-Adar index | 0.806 | 0.782 | 0.662 | 0.586 |
| Resource allocation index | 0.832 | 0.808 | 0.651 | 0.579 |
| Katz index | 0.887 | 0.866 | 0.720 | 0.643 |
| Random walk with restart | 0.924 | 0.913 | 0.737 | 0.671 |
| Jaccard index | 0.410 | 0.365 | 0.522 | 0.487 |
| Logistic regression on the link counts | 0.931 | *0.919* | 0.745 | 0.672 |
| Logistic regression on counts and neighbour scores | *0.932* | *0.919* | *0.759* | *0.683* |
| Link counts + MAESTRO + GRNBoost2 | 0.926 | 0.915 | 0.734 | 0.667 |
| Link counts + Geneformer | 0.927 | 0.917 | 0.730 | 0.665 |
| Counts and neighbour scores + Geneformer | 0.931 | *0.919* | 0.754 | 0.681 |

## 6. Ablations of scMultiomeGRN: message-passing graph, modalities and edge features

- The authors' code was patched so that each part can be switched off. With every switch at its default, the patched code gives exactly the published code's numbers.
- The known-link graph is what matters. Passing messages over no links at all drops AUPR from 0.918 to 0.826. The drop is small on lung but large on blood, where it falls to 0.585.
- A STRING protein network in place of the known links is worse still, at 0.781. STRING contains only 2 to 6% of the answer-key links.
- The pair histograms, the ATAC view and the RNA view can each be removed without a measurable change.
- Adding the Geneformer embedding as a third view also changes nothing, at 0.917. On blood almost all of scMultiomeGRN's score comes from the known links, and none of its extra inputs adds anything measurable.

| Method | AUPR | AUROC | AUPR, degree-matched | AUROC, degree-matched |
|---|---|---|---|---|
| MeVD-GRN | **0.935** | **0.924** | **0.878** | **0.857** |
| scMultiomeGRN (unchanged) | 0.918 | 0.912 | 0.705 | 0.649 |
| Without the TF-TF graph | 0.826 | 0.815 | 0.667 | 0.600 |
| STRING network instead of the known links | 0.781 | 0.770 | 0.624 | 0.580 |
| Without the pair histograms | *0.920* | *0.914* | *0.709* | *0.651* |
| ATAC features only | 0.916 | 0.912 | 0.695 | 0.649 |
| RNA features only | 0.916 | 0.910 | 0.705 | 0.650 |
| With Geneformer added | 0.917 | 0.911 | 0.705 | 0.649 |

## 7. Degree-matched negative sampling as a hub-controlled evaluation

- On the official test, a non-link is any random pair, so it usually involves TFs with few links. A model can do well just by knowing the hubs.
- On the degree-matched test, each test link is paired with one of the three test non-links closest to it in the two TFs' training-link counts. The pairs come from the official test set, so nothing new is revealed to any model.
- The counting rule falls from 0.924 to 0.724. In these dense graphs non-links between two hubs are rare, so it cannot reach exactly 0.5.
- scMultiomeGRN falls to 0.705 and the best simple model to 0.759. Both end up close to the counting rule.
- MeVD-GRN keeps 0.878. This is the clearest sign so far that it uses something beyond hubness.
- The test can only reuse the official test non-links, so hard non-links are scarce, and the next step is a split with more hard non-links, such as holding out whole TFs.

| Method | AUPR | AUROC | AUPR, degree-matched | AUROC, degree-matched |
|---|---|---|---|---|
| MeVD-GRN | **0.935** | **0.924** | **0.878** | **0.857** |
| scMultiomeGRN (unchanged) | 0.918 | 0.912 | 0.705 | 0.649 |
| Product of the two link counts | 0.924 | 0.912 | 0.724 | 0.653 |
| Logistic regression on the link counts | 0.931 | *0.919* | 0.745 | 0.672 |
| Logistic regression on counts and neighbour scores | *0.932* | *0.919* | *0.759* | *0.683* |
| Random walk with restart | 0.924 | 0.913 | 0.737 | 0.671 |
| Random forest on MAESTRO and GRNBoost2 | 0.836 | 0.814 | 0.650 | 0.606 |
| scMultiomeGRN without the TF-TF graph | 0.826 | 0.815 | 0.667 | 0.600 |

## 8. Conclusions and next steps: held-out-TF split

- On the official test, this benchmark mostly rewards knowing the hub TFs. A two-number logistic regression beats the published model, as recent papers found in other GRN benchmarks.
- Label-free network-inference methods cannot recover a motif-based answer key, so they are not a useful baseline for this task.
- scMultiomeGRN's features add almost nothing on top of its known-link graph, and Geneformer does not help it.
- MeVD-GRN is the only model that keeps most of its score once hubs are controlled for, at 0.878 against 0.705 for scMultiomeGRN.
- Next step: confirm that gap on a split that holds out whole TFs, where no model can count a TF's links.
- After that: identify which part of MeVD-GRN carries the gap, by scoring the existing ablations under the degree-matched test.
