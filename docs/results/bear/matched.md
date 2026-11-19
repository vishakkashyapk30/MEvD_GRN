# BEAR-GRN: degree-controlled evaluation of every method

Mean over seeds within a dataset, then over datasets. Pairs a method does not output score 0. Covariate = target in-degree among the training TFs of the fold holding the scored TF out (seed 42). In-space pairs = (ground-truth TF, measured ground-truth target). Matched sets are identical for every method.

## ChIP

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | caliper, pair win rate | covariate-adjusted | degree-corrected (Aiyappa) |
|---|---|---|---|---|---|---|---|---|
| CellOracle | 9 | 0.519 | 0.512 | 0.510 | 0.510 | 0.510 | 0.510 | 0.509 |
| TRIPOD | 9 | 0.532 | 0.508 | 0.509 | 0.509 | 0.510 | 0.510 | 0.507 |
| LINGER | 9 | 0.520 | 0.519 | 0.515 | 0.507 | 0.512 | 0.511 | 0.507 |
| GRaNIE | 9 | 0.510 | 0.506 | 0.506 | 0.504 | 0.506 | 0.505 | 0.503 |
| FigR | 9 | 0.503 | 0.504 | 0.500 | 0.502 | 0.502 | 0.502 | 0.500 |
| Pando-GLM | 9 | 0.502 | 0.501 | 0.503 | 0.501 | 0.501 | 0.502 | 0.502 |
| Pando-XGB | 9 | 0.502 | 0.501 | 0.502 | 0.501 | 0.501 | 0.501 | 0.502 |
| SCENIC+ | 9 | 0.501 | 0.500 | 0.501 | 0.501 | 0.501 | 0.501 | 0.502 |
| logistic regression, cell data only | 9 | 0.523 | 0.584 | 0.579 | 0.561 | 0.567 | 0.562 | 0.549 |
| Gaussian-kernel logistic regression, cell data only | 9 | 0.518 | 0.580 | 0.572 | 0.560 | 0.565 | 0.559 | 0.553 |
| LDA, cell data only | 9 | 0.523 | 0.585 | 0.574 | 0.557 | 0.564 | 0.558 | 0.548 |
| MLP (32 units), cell data only | 9 | 0.527 | 0.582 | 0.569 | 0.556 | 0.561 | 0.555 | 0.542 |
| gradient boosting, cell data only | 9 | 0.543 | 0.637 | 0.568 | 0.551 | 0.571 | 0.563 | 0.559 |
| QDA, cell data only | 9 | 0.528 | 0.569 | 0.568 | 0.551 | 0.559 | 0.555 | 0.557 |
| QDA, all features | 9 | 0.557 | 0.631 | 0.560 | 0.548 | 0.560 | 0.557 | 0.568 |
| Gaussian naive Bayes, cell data only | 9 | 0.524 | 0.562 | 0.567 | 0.546 | 0.552 | 0.548 | 0.551 |
| LDA, all features | 9 | 0.559 | 0.683 | 0.545 | 0.542 | 0.563 | 0.555 | 0.564 |
| Gaussian naive Bayes, all features | 9 | 0.552 | 0.614 | 0.562 | 0.541 | 0.548 | 0.550 | 0.558 |
| Gaussian-kernel logistic regression, all features | 9 | 0.561 | 0.685 | 0.532 | 0.537 | 0.567 | 0.547 | 0.575 |
| MLP (32 units), all features | 9 | 0.575 | 0.700 | 0.528 | 0.534 | 0.577 | 0.553 | 0.574 |
| logistic regression, all features | 9 | 0.576 | 0.717 | 0.526 | 0.534 | 0.600 | 0.568 | 0.576 |
| gradient boosting, all features | 9 | 0.581 | 0.701 | 0.528 | 0.531 | 0.583 | 0.559 | 0.582 |
| In-degree rule, all label files (script 25 check; = covariate) | 6 | 0.611 | 0.727 | 0.511 | 0.527 | 0.615 | 0.523 | 0.599 |
| MeVD-GRN, ATAC features zeroed (old ablation) | 3 | 0.584 | 0.685 | 0.523 | 0.525 | 0.553 | 0.543 | 0.556 |
| In-degree | 9 | 0.610 | 0.717 | 0.516 | 0.525 | 0.578 | 0.532 | 0.580 |
| MeVD-GRN | 9 | 0.576 | 0.662 | 0.523 | 0.524 | 0.547 | 0.537 | 0.553 |
| MeVD-GRN, ATAC removed (all routes) | 7 | 0.577 | 0.662 | 0.522 | 0.522 | 0.548 | 0.538 | 0.550 |
| MeVD-GRN, ATAC shuffled across genes | 7 | 0.584 | 0.663 | 0.522 | 0.522 | 0.541 | 0.532 | 0.554 |
| DeepSEM (label-free) | 9 | 0.514 | 0.526 | 0.531 | 0.511 | 0.511 | 0.511 | 0.513 |
| RegDiffusion (label-free) | 9 | 0.505 | 0.510 | 0.520 | 0.511 | 0.511 | 0.511 | 0.507 |
| Gaussian mixture classifier, all features | 9 | 0.470 | 0.482 | 0.509 | 0.511 | 0.519 | 0.517 | 0.496 |
| MeVD-GRN without Geneformer | 3 | 0.493 | 0.598 | 0.506 | 0.509 | 0.528 | 0.522 | 0.504 |
| Gaussian mixture classifier, cell data only | 9 | 0.498 | 0.509 | 0.504 | 0.506 | 0.508 | 0.509 | 0.500 |
| Gaussian graphical model (ridge, label-free) | 9 | 0.502 | 0.501 | 0.506 | 0.505 | 0.505 | 0.505 | 0.502 |
| GRNBoost2 | 9 | 0.508 | 0.501 | 0.504 | 0.504 | 0.505 | 0.505 | 0.513 |
| Pearson | 9 | 0.505 | 0.505 | 0.509 | 0.503 | 0.503 | 0.503 | 0.505 |
| All measured pairs | 9 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 |
| MeVD-GRN, L2 (no cell-type labels) | 9 | 0.508 | 0.516 | 0.512 | 0.497 | 0.497 | 0.497 | 0.492 |
| Gaussian naive Bayes, cell data only, trained on perturbation labels only | 6 | 0.495 | 0.469 | 0.489 | 0.485 | 0.483 | 0.486 | 0.481 |
| MeVD-GRN, sequential curriculum | 2 | 0.525 | 0.608 | 0.471 | 0.482 | 0.469 | 0.459 | 0.514 |
| logistic regression, cell data only, trained on perturbation labels only | 6 | 0.477 | 0.444 | 0.458 | 0.460 | 0.456 | 0.461 | 0.472 |
| LDA, cell data only, trained on perturbation labels only | 6 | 0.477 | 0.445 | 0.456 | 0.459 | 0.453 | 0.458 | 0.477 |
| gradient boosting, cell data only, trained on perturbation labels only | 6 | 0.490 | 0.463 | 0.437 | 0.443 | 0.437 | 0.442 | 0.484 |
| Gaussian naive Bayes, all features, trained on perturbation labels only | 6 | 0.484 | 0.459 | 0.436 | 0.442 | 0.439 | 0.444 | 0.475 |
| gradient boosting, all features, trained on perturbation labels only | 6 | 0.495 | 0.472 | 0.418 | 0.427 | 0.423 | 0.428 | 0.485 |
| LDA, all features, trained on perturbation labels only | 6 | 0.478 | 0.458 | 0.405 | 0.418 | 0.414 | 0.419 | 0.474 |
| In-degree rule, knockout labels only | 6 | 0.479 | 0.466 | 0.405 | 0.417 | 0.413 | 0.418 | 0.478 |
| logistic regression, all features, trained on perturbation labels only | 6 | 0.480 | 0.459 | 0.403 | 0.416 | 0.413 | 0.419 | 0.470 |

## KO

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | caliper, pair win rate | covariate-adjusted | degree-corrected (Aiyappa) |
|---|---|---|---|---|---|---|---|---|
| TRIPOD | 5 | 0.598 | 0.515 | 0.527 | 0.520 | 0.529 | 0.529 | 0.545 |
| FigR | 5 | 0.532 | 0.522 | 0.516 | 0.518 | 0.518 | 0.518 | 0.512 |
| Pando-GLM | 5 | 0.511 | 0.512 | 0.512 | 0.510 | 0.510 | 0.510 | 0.505 |
| Pando-XGB | 5 | 0.510 | 0.511 | 0.510 | 0.509 | 0.509 | 0.510 | 0.505 |
| CellOracle | 5 | 0.476 | 0.503 | 0.506 | 0.504 | 0.506 | 0.505 | 0.496 |
| SCENIC+ | 5 | 0.501 | 0.502 | 0.502 | 0.502 | 0.501 | 0.502 | 0.500 |
| GRaNIE | 5 | 0.490 | 0.490 | 0.498 | 0.498 | 0.498 | 0.498 | 0.502 |
| LINGER | 5 | 0.514 | 0.459 | 0.512 | 0.489 | 0.483 | 0.483 | 0.522 |
| logistic regression, all features, trained on perturbation labels only | 5 | 0.693 | 0.776 | 0.692 | 0.694 | 0.707 | 0.710 | 0.510 |
| In-degree rule, knockout labels only | 5 | 0.682 | 0.773 | 0.687 | 0.693 | 0.703 | 0.706 | 0.500 |
| LDA, all features, trained on perturbation labels only | 5 | 0.684 | 0.777 | 0.689 | 0.692 | 0.709 | 0.711 | 0.511 |
| gradient boosting, all features, trained on perturbation labels only | 5 | 0.679 | 0.761 | 0.681 | 0.679 | 0.695 | 0.697 | 0.511 |
| gradient boosting, cell data only, trained on perturbation labels only | 5 | 0.663 | 0.736 | 0.670 | 0.657 | 0.671 | 0.673 | 0.541 |
| Gaussian naive Bayes, all features, trained on perturbation labels only | 5 | 0.643 | 0.707 | 0.661 | 0.656 | 0.670 | 0.672 | 0.516 |
| LDA, cell data only, trained on perturbation labels only | 5 | 0.620 | 0.684 | 0.630 | 0.618 | 0.634 | 0.635 | 0.542 |
| logistic regression, cell data only, trained on perturbation labels only | 5 | 0.615 | 0.676 | 0.627 | 0.615 | 0.629 | 0.630 | 0.539 |
| Gaussian naive Bayes, cell data only, trained on perturbation labels only | 5 | 0.567 | 0.592 | 0.590 | 0.570 | 0.584 | 0.585 | 0.530 |
| MeVD-GRN, ATAC shuffled across genes | 4 | 0.535 | 0.537 | 0.548 | 0.551 | 0.581 | 0.585 | 0.518 |
| MeVD-GRN, ATAC features zeroed (old ablation) | 2 | 0.515 | 0.564 | 0.543 | 0.551 | 0.585 | 0.585 | 0.507 |
| MeVD-GRN, sequential curriculum | 2 | 0.606 | 0.671 | 0.542 | 0.550 | 0.619 | 0.619 | 0.509 |
| MeVD-GRN | 5 | 0.531 | 0.524 | 0.544 | 0.546 | 0.576 | 0.580 | 0.510 |
| MeVD-GRN, ATAC removed (all routes) | 4 | 0.538 | 0.535 | 0.544 | 0.546 | 0.585 | 0.587 | 0.514 |
| Pearson | 5 | 0.540 | 0.559 | 0.558 | 0.545 | 0.546 | 0.547 | 0.511 |
| Gaussian graphical model (ridge, label-free) | 5 | 0.523 | 0.534 | 0.528 | 0.527 | 0.528 | 0.526 | 0.508 |
| MeVD-GRN, L2 (no cell-type labels) | 5 | 0.537 | 0.535 | 0.525 | 0.525 | 0.529 | 0.528 | 0.518 |
| Gaussian mixture classifier, cell data only | 5 | 0.537 | 0.537 | 0.519 | 0.514 | 0.520 | 0.521 | 0.514 |
| GRNBoost2 | 5 | 0.501 | 0.512 | 0.511 | 0.511 | 0.509 | 0.510 | 0.501 |
| gradient boosting, all features | 5 | 0.473 | 0.449 | 0.518 | 0.509 | 0.513 | 0.520 | 0.497 |
| Gaussian mixture classifier, all features | 5 | 0.525 | 0.527 | 0.513 | 0.508 | 0.510 | 0.503 | 0.511 |
| MeVD-GRN without Geneformer | 2 | 0.516 | 0.529 | 0.505 | 0.507 | 0.524 | 0.524 | 0.515 |
| DeepSEM (label-free) | 5 | 0.506 | 0.505 | 0.540 | 0.507 | 0.507 | 0.508 | 0.503 |
| MLP (32 units), all features | 5 | 0.433 | 0.437 | 0.518 | 0.504 | 0.506 | 0.518 | 0.480 |
| All measured pairs | 5 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 |
| In-degree rule, all label files (script 25 check; = covariate) | 5 | 0.447 | 0.426 | 0.500 | 0.497 | 0.481 | 0.506 | 0.480 |
| In-degree | 5 | 0.446 | 0.425 | 0.498 | 0.496 | 0.486 | 0.501 | 0.479 |
| RegDiffusion (label-free) | 5 | 0.489 | 0.487 | 0.502 | 0.493 | 0.492 | 0.493 | 0.497 |
| Gaussian-kernel logistic regression, all features | 5 | 0.416 | 0.420 | 0.507 | 0.492 | 0.489 | 0.497 | 0.475 |
| logistic regression, all features | 5 | 0.419 | 0.406 | 0.489 | 0.481 | 0.424 | 0.428 | 0.473 |
| gradient boosting, cell data only | 5 | 0.433 | 0.398 | 0.499 | 0.468 | 0.454 | 0.456 | 0.489 |
| Gaussian naive Bayes, all features | 5 | 0.387 | 0.377 | 0.486 | 0.464 | 0.445 | 0.460 | 0.465 |
| Gaussian naive Bayes, cell data only | 5 | 0.381 | 0.385 | 0.483 | 0.458 | 0.450 | 0.450 | 0.463 |
| QDA, all features | 5 | 0.383 | 0.373 | 0.467 | 0.454 | 0.434 | 0.443 | 0.467 |
| Gaussian-kernel logistic regression, cell data only | 5 | 0.391 | 0.375 | 0.473 | 0.443 | 0.435 | 0.436 | 0.474 |
| QDA, cell data only | 5 | 0.377 | 0.368 | 0.462 | 0.442 | 0.433 | 0.431 | 0.466 |
| MLP (32 units), cell data only | 5 | 0.387 | 0.368 | 0.465 | 0.440 | 0.430 | 0.430 | 0.475 |
| LDA, all features | 5 | 0.386 | 0.366 | 0.448 | 0.438 | 0.399 | 0.406 | 0.466 |
| LDA, cell data only | 5 | 0.377 | 0.334 | 0.438 | 0.411 | 0.397 | 0.396 | 0.467 |
| logistic regression, cell data only | 5 | 0.373 | 0.328 | 0.428 | 0.405 | 0.389 | 0.389 | 0.465 |

## Union

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | caliper, pair win rate | covariate-adjusted | degree-corrected (Aiyappa) |
|---|---|---|---|---|---|---|---|---|
| CellOracle | 6 | 0.521 | 0.510 | 0.509 | 0.509 | 0.509 | 0.509 | 0.510 |
| GRaNIE | 6 | 0.517 | 0.506 | 0.507 | 0.505 | 0.507 | 0.506 | 0.505 |
| TRIPOD | 6 | 0.523 | 0.508 | 0.504 | 0.505 | 0.506 | 0.506 | 0.485 |
| LINGER | 6 | 0.531 | 0.510 | 0.506 | 0.503 | 0.508 | 0.507 | 0.493 |
| FigR | 6 | 0.510 | 0.504 | 0.502 | 0.501 | 0.501 | 0.501 | 0.497 |
| Pando-GLM | 6 | 0.504 | 0.502 | 0.502 | 0.501 | 0.501 | 0.501 | 0.502 |
| Pando-XGB | 6 | 0.504 | 0.502 | 0.501 | 0.501 | 0.501 | 0.501 | 0.502 |
| SCENIC+ | 6 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 |
| Gaussian-kernel logistic regression, cell data only | 6 | 0.560 | 0.567 | 0.587 | 0.567 | 0.571 | 0.563 | 0.577 |
| MLP (32 units), cell data only | 6 | 0.559 | 0.568 | 0.586 | 0.565 | 0.571 | 0.563 | 0.575 |
| logistic regression, cell data only | 6 | 0.550 | 0.559 | 0.587 | 0.562 | 0.566 | 0.558 | 0.563 |
| LDA, cell data only | 6 | 0.551 | 0.559 | 0.587 | 0.560 | 0.565 | 0.557 | 0.562 |
| gradient boosting, cell data only | 6 | 0.564 | 0.623 | 0.583 | 0.559 | 0.580 | 0.569 | 0.569 |
| QDA, cell data only | 6 | 0.541 | 0.552 | 0.574 | 0.552 | 0.558 | 0.551 | 0.561 |
| LDA, all features | 6 | 0.586 | 0.647 | 0.561 | 0.552 | 0.569 | 0.556 | 0.590 |
| QDA, all features | 6 | 0.570 | 0.609 | 0.566 | 0.549 | 0.558 | 0.554 | 0.577 |
| Gaussian-kernel logistic regression, all features | 6 | 0.604 | 0.674 | 0.545 | 0.548 | 0.590 | 0.561 | 0.612 |
| MLP (32 units), all features | 6 | 0.607 | 0.690 | 0.539 | 0.544 | 0.600 | 0.563 | 0.614 |
| Gaussian naive Bayes, cell data only | 6 | 0.531 | 0.542 | 0.568 | 0.542 | 0.546 | 0.541 | 0.546 |
| logistic regression, all features | 6 | 0.605 | 0.694 | 0.533 | 0.541 | 0.622 | 0.574 | 0.609 |
| Gaussian naive Bayes, all features | 6 | 0.560 | 0.593 | 0.563 | 0.536 | 0.540 | 0.544 | 0.559 |
| gradient boosting, all features | 6 | 0.593 | 0.686 | 0.530 | 0.535 | 0.601 | 0.567 | 0.591 |
| In-degree | 6 | 0.603 | 0.698 | 0.517 | 0.531 | 0.603 | 0.535 | 0.592 |
| MeVD-GRN | 6 | 0.570 | 0.675 | 0.524 | 0.530 | 0.556 | 0.541 | 0.558 |
| MeVD-GRN, ATAC removed (all routes) | 4 | 0.566 | 0.678 | 0.524 | 0.530 | 0.562 | 0.545 | 0.558 |
| MeVD-GRN, ATAC features zeroed (old ablation) | 2 | 0.592 | 0.688 | 0.519 | 0.527 | 0.558 | 0.544 | 0.570 |
| MeVD-GRN, ATAC shuffled across genes | 4 | 0.568 | 0.678 | 0.521 | 0.527 | 0.549 | 0.534 | 0.555 |
| In-degree rule, all label files (script 25 check; = covariate) | 6 | 0.603 | 0.698 | 0.511 | 0.526 | 0.612 | 0.523 | 0.592 |
| RegDiffusion (label-free) | 6 | 0.504 | 0.510 | 0.511 | 0.514 | 0.515 | 0.513 | 0.508 |
| MeVD-GRN without Geneformer | 2 | 0.559 | 0.612 | 0.508 | 0.510 | 0.535 | 0.527 | 0.542 |
| Gaussian mixture classifier, all features | 6 | 0.489 | 0.463 | 0.516 | 0.510 | 0.516 | 0.512 | 0.504 |
| DeepSEM (label-free) | 6 | 0.508 | 0.515 | 0.528 | 0.509 | 0.509 | 0.507 | 0.509 |
| GRNBoost2 | 6 | 0.495 | 0.505 | 0.504 | 0.505 | 0.505 | 0.505 | 0.503 |
| Gaussian graphical model (ridge, label-free) | 6 | 0.497 | 0.504 | 0.497 | 0.502 | 0.502 | 0.501 | 0.500 |
| Pearson | 6 | 0.495 | 0.510 | 0.506 | 0.502 | 0.502 | 0.502 | 0.501 |
| Gaussian mixture classifier, cell data only | 6 | 0.511 | 0.505 | 0.506 | 0.500 | 0.502 | 0.502 | 0.494 |
| All measured pairs | 6 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 |
| MeVD-GRN, L2 (no cell-type labels) | 6 | 0.521 | 0.514 | 0.491 | 0.498 | 0.497 | 0.496 | 0.497 |
| Gaussian naive Bayes, cell data only, trained on perturbation labels only | 6 | 0.513 | 0.482 | 0.491 | 0.487 | 0.485 | 0.488 | 0.489 |
| MeVD-GRN, sequential curriculum | 2 | 0.475 | 0.640 | 0.476 | 0.487 | 0.483 | 0.473 | 0.491 |
| logistic regression, cell data only, trained on perturbation labels only | 6 | 0.492 | 0.475 | 0.464 | 0.466 | 0.463 | 0.468 | 0.478 |
| LDA, cell data only, trained on perturbation labels only | 6 | 0.494 | 0.477 | 0.462 | 0.465 | 0.461 | 0.465 | 0.484 |
| gradient boosting, cell data only, trained on perturbation labels only | 6 | 0.503 | 0.506 | 0.448 | 0.453 | 0.449 | 0.454 | 0.487 |
| Gaussian naive Bayes, all features, trained on perturbation labels only | 6 | 0.503 | 0.497 | 0.444 | 0.451 | 0.451 | 0.455 | 0.481 |
| gradient boosting, all features, trained on perturbation labels only | 6 | 0.504 | 0.517 | 0.430 | 0.441 | 0.438 | 0.443 | 0.485 |
| LDA, all features, trained on perturbation labels only | 6 | 0.492 | 0.507 | 0.418 | 0.432 | 0.430 | 0.435 | 0.477 |
| In-degree rule, knockout labels only | 6 | 0.490 | 0.513 | 0.419 | 0.431 | 0.429 | 0.433 | 0.478 |
| logistic regression, all features, trained on perturbation labels only | 6 | 0.492 | 0.507 | 0.416 | 0.430 | 0.430 | 0.434 | 0.471 |

## Intersection

| Method | Datasets | AUROC | AUROC per TF | degree-matched | caliper-matched | caliper, pair win rate | covariate-adjusted | degree-corrected (Aiyappa) |
|---|---|---|---|---|---|---|---|---|
| TRIPOD | 6 | 0.628 | 0.519 | 0.529 | 0.518 | 0.521 | 0.524 | 0.556 |
| CellOracle | 6 | 0.468 | 0.507 | 0.512 | 0.512 | 0.512 | 0.512 | 0.500 |
| GRaNIE | 6 | 0.523 | 0.496 | 0.511 | 0.508 | 0.509 | 0.510 | 0.515 |
| FigR | 6 | 0.512 | 0.517 | 0.512 | 0.507 | 0.507 | 0.507 | 0.506 |
| LINGER | 6 | 0.548 | 0.489 | 0.527 | 0.506 | 0.507 | 0.506 | 0.536 |
| Pando-XGB | 6 | 0.510 | 0.506 | 0.513 | 0.504 | 0.504 | 0.505 | 0.504 |
| Pando-GLM | 6 | 0.509 | 0.506 | 0.513 | 0.504 | 0.504 | 0.505 | 0.503 |
| SCENIC+ | 6 | 0.500 | 0.501 | 0.501 | 0.501 | 0.501 | 0.501 | 0.500 |
| Gaussian naive Bayes, all features, trained on perturbation labels only | 6 | 0.548 | 0.621 | 0.577 | 0.573 | 0.580 | 0.579 | 0.490 |
| logistic regression, all features, trained on perturbation labels only | 6 | 0.574 | 0.669 | 0.572 | 0.571 | 0.574 | 0.574 | 0.471 |
| gradient boosting, cell data only, trained on perturbation labels only | 6 | 0.588 | 0.648 | 0.585 | 0.570 | 0.576 | 0.576 | 0.515 |
| LDA, all features, trained on perturbation labels only | 6 | 0.570 | 0.670 | 0.570 | 0.570 | 0.575 | 0.574 | 0.474 |
| In-degree rule, knockout labels only | 6 | 0.560 | 0.665 | 0.564 | 0.567 | 0.571 | 0.568 | 0.466 |
| LDA, cell data only, trained on perturbation labels only | 6 | 0.545 | 0.613 | 0.571 | 0.558 | 0.569 | 0.568 | 0.513 |
| logistic regression, cell data only, trained on perturbation labels only | 6 | 0.540 | 0.607 | 0.570 | 0.556 | 0.567 | 0.565 | 0.510 |
| gradient boosting, all features, trained on perturbation labels only | 6 | 0.561 | 0.657 | 0.556 | 0.555 | 0.564 | 0.560 | 0.469 |
| Gaussian naive Bayes, cell data only, trained on perturbation labels only | 6 | 0.529 | 0.562 | 0.564 | 0.549 | 0.559 | 0.559 | 0.519 |
| Pearson | 6 | 0.526 | 0.529 | 0.544 | 0.522 | 0.524 | 0.525 | 0.507 |
| MeVD-GRN, L2 (no cell-type labels) | 6 | 0.546 | 0.530 | 0.542 | 0.520 | 0.522 | 0.524 | 0.514 |
| MeVD-GRN without Geneformer | 2 | 0.499 | 0.523 | 0.515 | 0.518 | 0.531 | 0.535 | 0.476 |
| MeVD-GRN, ATAC removed (all routes) | 4 | 0.519 | 0.572 | 0.518 | 0.518 | 0.542 | 0.544 | 0.522 |
| MeVD-GRN | 6 | 0.512 | 0.556 | 0.521 | 0.518 | 0.526 | 0.527 | 0.506 |
| MLP (32 units), all features | 6 | 0.492 | 0.518 | 0.526 | 0.515 | 0.528 | 0.535 | 0.503 |
| Gaussian graphical model (ridge, label-free) | 6 | 0.511 | 0.520 | 0.518 | 0.513 | 0.514 | 0.514 | 0.506 |
| Gaussian mixture classifier, cell data only | 6 | 0.561 | 0.513 | 0.518 | 0.513 | 0.519 | 0.516 | 0.518 |
| Gaussian mixture classifier, all features | 6 | 0.500 | 0.503 | 0.513 | 0.512 | 0.514 | 0.509 | 0.503 |
| MeVD-GRN, ATAC shuffled across genes | 4 | 0.495 | 0.559 | 0.505 | 0.512 | 0.520 | 0.518 | 0.511 |
| gradient boosting, all features | 6 | 0.544 | 0.522 | 0.516 | 0.510 | 0.523 | 0.526 | 0.522 |
| DeepSEM (label-free) | 6 | 0.506 | 0.515 | 0.542 | 0.508 | 0.508 | 0.510 | 0.506 |
| Gaussian-kernel logistic regression, all features | 6 | 0.461 | 0.496 | 0.513 | 0.505 | 0.514 | 0.517 | 0.495 |
| RegDiffusion (label-free) | 6 | 0.509 | 0.502 | 0.522 | 0.505 | 0.505 | 0.504 | 0.508 |
| GRNBoost2 | 6 | 0.506 | 0.511 | 0.504 | 0.503 | 0.502 | 0.505 | 0.501 |
| All measured pairs | 6 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 | 0.500 |
| In-degree rule, all label files (script 25 check; = covariate) | 6 | 0.496 | 0.507 | 0.500 | 0.499 | 0.497 | 0.507 | 0.501 |
| In-degree | 6 | 0.496 | 0.505 | 0.499 | 0.499 | 0.500 | 0.505 | 0.500 |
| gradient boosting, cell data only | 6 | 0.491 | 0.479 | 0.519 | 0.495 | 0.488 | 0.495 | 0.515 |
| logistic regression, all features | 6 | 0.465 | 0.491 | 0.498 | 0.495 | 0.480 | 0.485 | 0.495 |
| MeVD-GRN, ATAC features zeroed (old ablation) | 2 | 0.397 | 0.547 | 0.487 | 0.493 | 0.492 | 0.490 | 0.434 |
| MeVD-GRN, sequential curriculum | 2 | 0.570 | 0.677 | 0.488 | 0.491 | 0.514 | 0.519 | 0.506 |
| Gaussian naive Bayes, cell data only | 6 | 0.393 | 0.452 | 0.503 | 0.487 | 0.486 | 0.488 | 0.478 |
| Gaussian naive Bayes, all features | 6 | 0.409 | 0.445 | 0.503 | 0.486 | 0.480 | 0.491 | 0.483 |
| Gaussian-kernel logistic regression, cell data only | 6 | 0.441 | 0.452 | 0.506 | 0.486 | 0.484 | 0.486 | 0.496 |
| MLP (32 units), cell data only | 6 | 0.453 | 0.450 | 0.499 | 0.482 | 0.476 | 0.481 | 0.502 |
| QDA, all features | 6 | 0.417 | 0.444 | 0.488 | 0.479 | 0.476 | 0.483 | 0.485 |
| QDA, cell data only | 6 | 0.413 | 0.444 | 0.490 | 0.478 | 0.477 | 0.479 | 0.485 |
| LDA, all features | 6 | 0.428 | 0.448 | 0.471 | 0.467 | 0.454 | 0.458 | 0.484 |
| LDA, cell data only | 6 | 0.427 | 0.424 | 0.485 | 0.462 | 0.458 | 0.460 | 0.492 |
| logistic regression, cell data only | 6 | 0.425 | 0.421 | 0.480 | 0.459 | 0.454 | 0.458 | 0.491 |

