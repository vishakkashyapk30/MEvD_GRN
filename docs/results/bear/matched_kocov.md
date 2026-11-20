# BEAR-GRN: degree-controlled evaluation of every method

Mean over seeds within a dataset, then over datasets. Pairs a method does not output score 0. Covariate = target in-degree among the training TFs of the fold holding the scored TF out (seed 42). In-space pairs = (ground-truth TF, measured ground-truth target). Matched sets are identical for every method.

## KO

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | caliper, pair win rate | covariate-adjusted | degree-corrected (Aiyappa) |
|---|---|---|---|---|---|---|---|---|
| TRIPOD | 5 | 0.598 | 0.515 | 0.573 | 0.507 | 0.510 | 0.511 | 0.545 |
| FigR | 5 | 0.532 | 0.522 | 0.515 | 0.505 | 0.506 | 0.506 | 0.512 |
| Pando-GLM | 5 | 0.511 | 0.512 | 0.514 | 0.503 | 0.503 | 0.504 | 0.505 |
| Pando-XGB | 5 | 0.510 | 0.511 | 0.513 | 0.503 | 0.503 | 0.503 | 0.505 |
| CellOracle | 5 | 0.476 | 0.503 | 0.511 | 0.502 | 0.503 | 0.502 | 0.496 |
| SCENIC+ | 5 | 0.501 | 0.502 | 0.502 | 0.501 | 0.501 | 0.500 | 0.500 |
| GRaNIE | 5 | 0.490 | 0.490 | 0.506 | 0.499 | 0.500 | 0.499 | 0.502 |
| LINGER | 5 | 0.514 | 0.459 | 0.544 | 0.498 | 0.495 | 0.494 | 0.522 |
| LDA, cell data only, trained on perturbation labels only | 5 | 0.620 | 0.684 | 0.580 | 0.539 | 0.550 | 0.550 | 0.542 |
| logistic regression, cell data only, trained on perturbation labels only | 5 | 0.615 | 0.676 | 0.587 | 0.538 | 0.547 | 0.548 | 0.539 |
| gradient boosting, cell data only, trained on perturbation labels only | 5 | 0.663 | 0.736 | 0.580 | 0.530 | 0.544 | 0.549 | 0.541 |
| Gaussian naive Bayes, cell data only, trained on perturbation labels only | 5 | 0.567 | 0.592 | 0.593 | 0.527 | 0.534 | 0.534 | 0.530 |
| Gaussian naive Bayes, all features, trained on perturbation labels only | 5 | 0.643 | 0.707 | 0.551 | 0.517 | 0.534 | 0.536 | 0.516 |
| Pearson | 5 | 0.540 | 0.559 | 0.544 | 0.516 | 0.518 | 0.517 | 0.511 |
| gradient boosting, all features, trained on perturbation labels only | 5 | 0.679 | 0.761 | 0.537 | 0.511 | 0.531 | 0.534 | 0.511 |
| Gaussian graphical model (ridge, label-free) | 5 | 0.523 | 0.534 | 0.536 | 0.509 | 0.509 | 0.511 | 0.508 |
| LDA, all features, trained on perturbation labels only | 5 | 0.684 | 0.777 | 0.518 | 0.509 | 0.545 | 0.552 | 0.511 |
| logistic regression, all features, trained on perturbation labels only | 5 | 0.693 | 0.776 | 0.520 | 0.508 | 0.543 | 0.548 | 0.510 |
| DeepSEM (label-free) | 5 | 0.506 | 0.505 | 0.550 | 0.503 | 0.503 | 0.503 | 0.503 |
| MeVD-GRN, seeds 42-44 only (the seeds of the ATAC and negative-sampling runs) | 5 | 0.537 | 0.525 | 0.607 | 0.501 | 0.505 | 0.504 | 0.513 |
| MeVD-GRN | 5 | 0.531 | 0.524 | 0.608 | 0.501 | 0.503 | 0.503 | 0.510 |
| In-degree rule, knockout labels only | 5 | 0.682 | 0.773 | 0.500 | 0.500 | 0.507 | 0.514 | 0.500 |
| MeVD-GRN, same-TF training negatives (simple models' rule) | 5 | 0.527 | 0.496 | 0.639 | 0.500 | 0.502 | 0.502 | 0.516 |
| RegDiffusion (label-free) | 5 | 0.489 | 0.487 | 0.484 | 0.496 | 0.496 | 0.496 | 0.497 |
| In-degree | 5 | 0.446 | 0.425 | 0.618 | 0.477 | 0.478 | 0.476 | 0.479 |
| LDA, cell data only | 5 | 0.377 | 0.334 | 0.534 | 0.471 | 0.467 | 0.465 | 0.467 |
| logistic regression, cell data only | 5 | 0.373 | 0.328 | 0.526 | 0.469 | 0.464 | 0.463 | 0.465 |

