# MeVD-GRN under the BEAR-GRN protocol (SC-MO-GRN-DB baseline)

Status: **pipeline built and verified; Ada chain 587 -> 588 -> 615 queued; first K562 result in s10.2** (living document; written first, updated as the build
progresses so nothing is lost if a session is interrupted). All decisions
below were made without user input (the user asked not to be consulted before
18:00 on 2026-10-01). Each one is marked **[decision]** with its reason.

Goal: BEAR-GRN (Karamveer, Moeller, Valensi, Manful, Uzun, *Nat Commun*,
18 Sep 2026, doi 10.1038/s41467-026-77838-w) benchmarks single-cell multiome
GRN methods on SC-MO-GRN-DB's own datasets and ground truths. It is from the
SC-MO-GRN-DB lab and is the new published baseline for our SC-MO-GRN-DB
numbers. We score MeVD-GRN with BEAR-GRN's own metric code and slot the
numbers into its tables, using a supervision protocol that is fair to its
unsupervised methods.

File lanes (this experiment only): `scripts/20-24_bear_*`, `scripts/_bear_pget.py`,
`configs/bear/`, `slurm/bear_*.sh`, `src/benchmarks/bear_*.py`, this file.
Local scratch: `~/.cache/bear/`. Ada data root: `/share1/$USER/mevd_grn/bear`.

---

## 1. Sources

| What | Where | Notes |
|---|---|---|
| Paper | nature.com/articles/s41467-026-77838-w | Accelerated Article Preview: **no main text online yet**; the abstract + 12 supplementary files are. |
| Supplementary files | `~/.cache/bear/paper/MOESM1-12_ESM.*` | SI (figs 1-15), Reporting Summary, Supp Data 1-7 (xlsx), Source Data, **peer-review file (MOESM11)** which holds most method details |
| Tidy numbers (extracted from Supp Data 2-5, Source Data) | `~/.cache/bear/paper/bear_grn_auroc_auprc_tidy.csv`, `bear_grn_top10k_prf_chip_tidy.csv`, `bear_grn_modality_shuffle_tidy.csv` | every number in section 5 comes from these |
| Code (R package `BEARGRN` v1.0.0) | Zenodo 10.5281/zenodo.22015920 = GitHub UzunLab/BEAR-GRN @033a2cd; extracted at `~/.cache/bear/paper/code/UzunLab-BEAR-GRN-033a2cd/` | MIT licence. Scoring: `R/reproduce_ROC_PR_full.R` (paper numbers), `R/reproduce_early_metrics.R`, `R/reproduce_stability.R`; `benchmark_new_method*.R` = same code for one new method |
| Data | Zenodo 10.5281/zenodo.20704929 (21.7 GB total) | `INPUT.DATA.zip` 707 MB (the exact filtered RNA+ATAC matrices every method was run on), `INFERRED.GRNS.zip` 1.44 GB (**released outputs of all 9 methods x 9 datasets**), `GROUND.TRUTHS*.zip` (ChIP, KO, union, intersection, core, cell-type-exclusive), `INPUT.DATA.STABILITY.zip` 1.67 GB, `STABILITY_GRNS.zip` 12.5 GB, `GENOME_AND_ANNOTATION.zip` 5.3 GB |

Zenodo throttles one connection to ~10-50 KB/s (measured 2026-10-01 both
locally and on ada-gw1), so `scripts/20_bear_download.py` fetches only the
needed zip members with parallel HTTP range requests into a sparse local copy
of the zip (Deflate64 members of INFERRED.GRNS.zip are unpacked with Info-ZIP
`unzip`, or `zipfile-deflate64` if installed).

## 2. Datasets (Supp Data 1)

| BEAR name | Cell type | Paired? | Cells raw -> after QC | Genes after QC | Accession | SC-MO-GRN-DB DS (inferred, unverified) | ChIP GT |
|---|---|---|---|---|---|---|---|
| K562 | K562 (human) | yes | 571 -> 411 | 10,012 | GSE178707 | ~DS025 | RN117 |
| Macrophage_S1 | macrophage buffer 1 (human) | yes | 892 -> 474 | 10,050 | PRJNA1102756 / Zenodo 11001642 | DS026 | RN204 |
| Macrophage_S2 | macrophage buffer 2 (human) | yes | 855 -> 528 | 9,171 | same | DS026 | RN204 |
| iPS | iPSC WT_UNTREATED_D13 (human) | yes | 6,577 -> 5,710 | 4,808 | GSE278751 | not in DB list ("RN000") | RN000 |
| mESC_E7.5_rep1 | mouse embryo E7.5 r1 | yes | 7,416 -> 5,064 | 9,374 | GSE205117 | ~DS014 | RN111 |
| mESC_E7.5_rep2 | mouse embryo E7.5 r2 | yes | 10,772 -> 1,447 | 7,712 | GSE205117 | ~DS014 | RN111 |
| mESC_E8.5_rep1 | mouse embryo E8.5 r1 | yes | 8,826 -> 6,557 | 9,296 | GSE205117 | ~DS014 | RN111 |
| mESC_E8.5_rep2 | mouse embryo E8.5 r2 | yes | 6,918 -> 6,252 | 9,882 | GSE205117 | ~DS014 | RN111 |
| Naive_mESC | naive mESC | **no** (RNA 3,959 / ATAC 5,000 cells, CCA-paired) | -> 3,568 | 8,049 | GSE198730 | ~DS011 | RN111 |

`INPUT.DATA/<dataset>/` holds the post-QC matrices as CSV, features x cells
(K562: 10,012 genes x 411 cells RNA; 177,439 peaks x 411 cells ATAC, peak
names `chr1:115460-115953`). Preprocessing: `inst/scripts/DATA.PREPROCESSING/
Step1.QC.R` -> `Step2.filter_cells_genes.R` -> `Step3.select_common_cells.R`
(-> `Step4.subsample_cells.R` for stability). **[decision]** We do not redo
their QC: MeVD-GRN reads `INPUT.DATA` directly, so every method sees the same
cells and genes.

Our existing local data (`data/raw/single_cell/`) is SC-MO-GRN-DB's DS019 (K562
scRNA, unpaired with DS025 ATAC h5), DS026 (Macrophage, both buffers merged)
and DS027 (MCF7, not in BEAR). None of it is BEAR's filtered input, so for
this benchmark all inputs come from Zenodo. ESC DS010/DS012 are not in BEAR
(BEAR's mouse data are embryo DS014-like and naive mESC).

## 3. Ground truths ("four complementary definitions")

From the GT zips (counts after upper-casing + dedup, K562 measured locally):

| Definition | K562 file (edges / TFs / targets) | Mouse | Macrophage / iPS |
|---|---|---|---|
| **ChIP** (`GROUND.TRUTHS/filtered_<RN>_<ds>.tsv`) | RN117: 1,377,816 / 150 / 27,747 | RN111 (243 TFs) | RN204 (23 TFs) / RN000 (75 TFs) |
| **KO** (`GROUND.TRUTHS.KO/`, = `filtered_RN118_K562.tsv`) | RN118: 229,092 / 96 / 14,927 | RN112 (55 TFs) | none |
| **Union** ChIP or KO (`GROUND.TRUTHS.UNION/`) | 1,587,015 / 226 / 28,765 | 272 TFs | none |
| **Intersection** ChIP and KO, computed by BEAR (not RN119) (`GROUND.TRUTHS.INETRSECTION/`, sic) | 19,893 / 20 / 8,305 | 26 TFs | none |
| core (edges in >50% of human cell types) / cell-type-exclusive (Supp Fig 12) | 59,573 / 24 / 23,698 ; 1,648,555 / 226 / 28,822 | none | yes |

GT files are "filtered based on the genes expressed in each sample" (Supp
Data 1 footnote), but the target sets are far larger than the post-QC gene
lists (27,747 vs 10,012 for K562), so most GT targets cannot be scored by any
method and carry score 0 in the AUPRC universe.

These map almost one-to-one onto our curriculum tiers for K562: ChIP =
localization (RN117), KO = perturbation (RN118), intersection ~ dual evidence
(BEAR's own ChIP and KO, not RN119).

## 4. Scoring protocol (exact, from `R/reproduce_ROC_PR_full.R`, `R/reproduce_early_metrics.R`, `R/benchmark_new_method_stability.R`)

1. **Input edge list:** (TF, target, score) per method; `score <- abs(score)`,
   names upper-cased, duplicates collapsed with `max(score)`.
2. **Tested space:** tested TFs = unique GT sources, tested targets = unique GT
   targets. Edges outside are dropped (`use_filtered_approach = TRUE`).
3. **AUROC:** on the method's **own predicted edges** in the tested space only.
   Classes balanced 1:1 by down-sampling the larger one (`set.seed(42 + i)`,
   `sample_n`), then `pROC::roc(levels = c(0,1), direction = "<")`.
   Consequence: sparse methods are scored only on the edges they chose to emit.
4. **AUPRC:** on the **complete universe** tested_TFs x tested_targets;
   unpredicted pairs get score 0; all positives + a random sample of
   min(10 x positives, all) negatives; `PRROC::pr.curve(...)$auc.integral`.
   For dense GTs the 10x cap is not binding (K562 ChIP: 1.38 M positives,
   2.78 M negatives), so AUPRC_random equals the GT density (K562 ChIP 0.331).
5. **Random baseline:** the same with `runif` scores over the universe.
6. **Early metrics (ChIP GT only in the paper):** sort by |score| desc, keep
   edges with GT TF and GT target, `head(10000)`; precision = TP/(TP+FP),
   recall = TP/|GT|, F1.
7. **Stability:** each method re-run on cell subsamples (`INPUT.DATA.STABILITY`);
   per network keep the top `round(10%)` edges by score (a random 10% as
   control); pairwise Jaccard over networks; median JI per dataset x method.
8. **Modality perturbation (Source Data, Supp Figs 13-14):** AUROC/AUPRC on
   `Original`, `Cell_Shuffle`, `Gene_Shuffle` (RNA) and `ATAC_Gene_Shuffle`
   inputs; finding: methods are predominantly RNA-driven.
9. **Mixed multiome (Supp Fig 15):** RNA from one cell type, ATAC from another.

Known code hazards (from the helper's read of the code):
- `benchmark_new_method*.R` map `Macrophage_S2` to the Buffer1 GT;
  `reproduce_ROC_PR_full.R` (the paper) uses Buffer2. **[decision]** follow the paper (Buffer2).
- GT-file matching uses regex on file names; point `ground_truth_dir` at a
  directory holding one GT type only.
- DIRECT-NET's released network has no score column, so it has no AUROC/AUPRC
  in the paper (early metrics only, as an unsorted head(10000)).

**Our scorer.** `src/benchmarks/bear_metrics.py` is a line-by-line Python port.
`prroc_auc_integral` and the pROC AUC were checked against R (PRROC 1.4,
pROC 1.18) on a tied-score test vector: identical to 7 digits
(0.2702473 / 0.7199483). The only unavoidable difference is the RNG of the
1:1 / 1:10 subsamples (R's `set.seed` + `sample_n` cannot be reproduced in
numpy); `n_rep` averages several draws. (Running the official R functions
directly was attempted locally, but tidyverse binaries are not available for
this laptop's Ubuntu 25.04 and Ada has no R; the port is instead verified on
the released GRNs, section 8.)

## 5. Numbers to beat (BEAR-GRN Supp Data 2-5; 3 decimals)

AUROC / AUPRC per dataset (rows) x method (cols). DIRECT-NET has no score
column and is not in these tables. AUPRC_random = GT density of the universe.

**ChIP GT, AUROC**

| dataset | CellOracle | FigR | GRaNIE | LINGER | Pando-GLM | Pando-XGB | SCENIC+ | TRIPOD |
|---|---|---|---|---|---|---|---|---|
| K562 | 0.539 | 0.504 | 0.512 | 0.536 | **0.575** | 0.518 | 0.466 | 0.512 |
| Macrophage_S1 | 0.584 | 0.529 | **0.650** | 0.596 | 0.574 | 0.558 | 0.573 | 0.553 |
| Macrophage_S2 | **0.606** | 0.533 | 0.473 | 0.547 | 0.507 | 0.536 | 0.514 | 0.539 |
| iPS | 0.481 | 0.468 | **0.568** | 0.501 | 0.503 | 0.442 | 0.440 | 0.564 |
| E7.5_rep1 | 0.529 | 0.507 | 0.570 | **0.652** | 0.458 | 0.501 | 0.415 | 0.492 |
| E7.5_rep2 | 0.504 | 0.500 | 0.596 | **0.644** | 0.392 | 0.514 | 0.459 | 0.515 |
| E8.5_rep1 | 0.481 | 0.496 | 0.589 | **0.658** | 0.382 | 0.520 | 0.391 | 0.494 |
| E8.5_rep2 | 0.512 | 0.490 | 0.622 | **0.656** | 0.401 | 0.534 | 0.427 | 0.525 |
| Naive_mESC | 0.508 | 0.511 | 0.539 | 0.473 | 0.455 | 0.534 | **0.556** | 0.515 |

**ChIP GT, AUPRC** (random = 0.331 K562, 0.273 / 0.271 Mac S1 / S2, 0.091 iPS, 0.153-0.156 mouse)

| dataset | CellOracle | FigR | GRaNIE | LINGER | Pando-GLM | Pando-XGB | SCENIC+ | TRIPOD |
|---|---|---|---|---|---|---|---|---|
| K562 | 0.390 | 0.376 | 0.368 | **0.430** | 0.348 | 0.346 | 0.332 | 0.405 |
| Macrophage_S1 | 0.293 | 0.291 | 0.269 | 0.321 | 0.289 | 0.289 | 0.285 | **0.332** |
| Macrophage_S2 | 0.301 | 0.296 | 0.268 | 0.324 | 0.276 | 0.281 | 0.279 | **0.336** |
| iPS | 0.110 | 0.094 | 0.099 | 0.108 | 0.091 | 0.091 | 0.091 | **0.119** |
| E7.5_rep1 | 0.191 | 0.178 | 0.188 | **0.241** | 0.159 | 0.159 | 0.156 | 0.193 |
| E7.5_rep2 | 0.182 | 0.172 | 0.179 | **0.225** | 0.155 | 0.155 | 0.155 | 0.186 |
| E8.5_rep1 | 0.185 | 0.173 | 0.184 | **0.227** | 0.157 | 0.158 | 0.155 | 0.181 |
| E8.5_rep2 | 0.188 | 0.171 | 0.184 | **0.224** | 0.160 | 0.161 | 0.156 | 0.183 |
| Naive_mESC | 0.163 | 0.171 | 0.182 | 0.160 | 0.164 | 0.165 | 0.153 | **0.192** |

**KO GT** (K562 random 0.159; mouse ~0.090; naive 0.110)

| dataset | metric | CellOracle | FigR | GRaNIE | LINGER | Pando-GLM | Pando-XGB | SCENIC+ | TRIPOD |
|---|---|---|---|---|---|---|---|---|---|
| K562 | AUROC | 0.492 | 0.450 | 0.507 | 0.402 | 0.472 | 0.496 | **0.573** | 0.494 |
| K562 | AUPRC | 0.147 | 0.149 | 0.152 | 0.139 | 0.156 | 0.156 | **0.159** | 0.149 |
| E7.5_rep1 | AUROC / AUPRC | 0.541 / 0.090 | 0.550 / 0.102 | 0.302 / 0.094 | 0.463 / 0.095 | 0.526 / 0.095 | 0.509 / 0.095 | 0.535 / 0.091 | 0.571 / 0.113 |
| E7.5_rep2 | AUROC / AUPRC | 0.503 / 0.089 | 0.531 / 0.099 | 0.408 / 0.088 | 0.455 / 0.097 | 0.573 / 0.093 | 0.490 / 0.093 | 0.503 / 0.091 | 0.520 / 0.108 |
| E8.5_rep1 | AUROC / AUPRC | 0.589 / 0.088 | 0.530 / 0.097 | 0.404 / 0.089 | 0.483 / 0.095 | 0.527 / 0.095 | 0.579 / 0.095 | 0.619 / 0.092 | 0.540 / 0.106 |
| E8.5_rep2 | AUROC / AUPRC | 0.602 / 0.090 | 0.544 / 0.100 | 0.311 / 0.089 | 0.475 / 0.097 | 0.559 / 0.098 | 0.545 / 0.098 | 0.564 / 0.092 | 0.579 / 0.118 |
| Naive_mESC | AUROC / AUPRC | 0.559 / 0.110 | 0.471 / 0.111 | 0.492 / 0.113 | 0.596 / 0.112 | 0.463 / 0.113 | 0.506 / 0.114 | 0.292 / 0.110 | 0.523 / 0.121 |

Every K562 KO AUPRC is at or below the 0.159 random baseline.

**Union GT** (K562 random 0.244; mouse 0.128-0.134)

| dataset | metric | CellOracle | FigR | GRaNIE | LINGER | Pando-GLM | Pando-XGB | SCENIC+ | TRIPOD |
|---|---|---|---|---|---|---|---|---|---|
| K562 | AUROC | 0.521 | 0.497 | 0.526 | 0.525 | **0.570** | 0.522 | 0.473 | 0.510 |
| K562 | AUPRC | 0.293 | 0.288 | 0.281 | **0.344** | 0.261 | 0.259 | 0.245 | 0.322 |
| E7.5_rep1 | AUROC / AUPRC | 0.517 / 0.162 | 0.516 / 0.152 | 0.555 / 0.159 | 0.633 / 0.205 | 0.468 / 0.135 | 0.517 / 0.134 | 0.433 / 0.132 | 0.507 / 0.168 |
| E7.5_rep2 | AUROC / AUPRC | 0.492 / 0.153 | 0.505 / 0.146 | 0.585 / 0.151 | 0.629 / 0.192 | 0.418 / 0.130 | 0.526 / 0.130 | 0.453 / 0.130 | 0.516 / 0.160 |
| E8.5_rep1 | AUROC / AUPRC | 0.470 / 0.159 | 0.501 / 0.149 | 0.594 / 0.158 | 0.644 / 0.195 | 0.412 / 0.135 | 0.529 / 0.135 | 0.425 / 0.133 | 0.506 / 0.157 |
| E8.5_rep2 | AUROC / AUPRC | 0.505 / 0.161 | 0.497 / 0.148 | 0.610 / 0.157 | 0.637 / 0.192 | 0.429 / 0.137 | 0.537 / 0.138 | 0.471 / 0.133 | 0.534 / 0.160 |
| Naive_mESC | AUROC / AUPRC | 0.521 / 0.144 | 0.496 / 0.151 | 0.552 / 0.162 | 0.475 / 0.142 | 0.459 / 0.145 | 0.531 / 0.146 | 0.478 / 0.134 | 0.523 / 0.171 |

**Intersection GT** (K562 random 0.119; mouse ~0.090; naive 0.094)

| dataset | metric | CellOracle | FigR | GRaNIE | LINGER | Pando-GLM | Pando-XGB | SCENIC+ | TRIPOD |
|---|---|---|---|---|---|---|---|---|---|
| K562 | AUROC | 0.443 | 0.442 | 0.472 | 0.359 | **0.529** | 0.445 | 0.444 | 0.498 |
| K562 | AUPRC | 0.106 | 0.111 | **0.120** | 0.099 | 0.114 | 0.114 | 0.119 | 0.115 |
| E7.5_rep1 | AUROC / AUPRC | 0.518 / 0.091 | 0.556 / 0.113 | 0.348 / 0.109 | 0.515 / 0.112 | 0.526 / 0.094 | 0.469 / 0.094 | 0.540 / 0.092 | 0.553 / 0.140 |
| E7.5_rep2 | AUROC / AUPRC | 0.476 / 0.087 | 0.530 / 0.102 | 0.476 / 0.087 | 0.506 / 0.108 | 0.537 / 0.092 | 0.433 / 0.092 | 0.542 / 0.092 | 0.506 / 0.122 |
| E8.5_rep1 | AUROC / AUPRC | 0.428 / 0.082 | 0.497 / 0.091 | 0.555 / 0.093 | 0.585 / 0.104 | 0.473 / 0.093 | 0.585 / 0.094 | 0.637 / 0.091 | 0.523 / 0.103 |
| E8.5_rep2 | AUROC / AUPRC | 0.570 / 0.090 | 0.527 / 0.097 | 0.484 / 0.093 | 0.564 / 0.110 | 0.540 / 0.100 | 0.542 / 0.100 | 0.529 / 0.092 | 0.591 / 0.139 |
| Naive_mESC | AUROC / AUPRC | 0.613 / 0.097 | 0.483 / 0.099 | 0.503 / 0.104 | 0.747 / 0.106 | 0.428 / 0.105 | 0.511 / 0.108 | 0.500 / 0.096 | 0.523 / 0.124 |

**Top-10k precision / recall / F1, ChIP GT** (`bear_grn_top10k_prf_chip_tidy.csv`, includes DIRECT-NET)

| dataset | CellOracle | DIRECT-NET | FigR | GRaNIE | LINGER | Pando-GLM | Pando-XGB | SCENIC+ | TRIPOD |
|---|---|---|---|---|---|---|---|---|---|
| K562 precision | 0.675 | 0.575 | 0.562 | 0.553 | 0.587 | **0.682** | 0.570 | 0.504 | 0.575 |
| K562 F1 | **0.0097** | 0.0083 | 0.0081 | 0.0080 | 0.0085 | 0.0029 | 0.0082 | 0.0022 | 0.0083 |
| Mac-S1 precision | 0.448 | 0.647 | 0.403 | 0.034 | 0.477 | 0.552 | 0.524 | **0.882** | 0.550 |
| Mac-S1 F1 | 0.0506 | **0.0732** | 0.0456 | 0.0006 | 0.0539 | 0.0284 | 0.0426 | 0.0161 | 0.0622 |
| Mac-S2 precision | 0.511 | 0.640 | 0.444 | 0.067 | 0.460 | **0.688** | 0.554 | 0.487 | 0.562 |
| Mac-S2 F1 | 0.0579 | **0.0725** | 0.0503 | 0.0017 | 0.0521 | 0.0058 | 0.0193 | 0.0185 | 0.0637 |

(Full 9-dataset P/R/F1 tables are in the CSV; F1 at 10k edges is bounded by
recall, which is <= 10,000/|GT|.)

Reading of the bar: on ChIP, LINGER has the best AUPRC on K562 and all 4
mouse embryo sets, TRIPOD on macrophage / iPS / naive. No method reaches
AUROC 0.66 anywhere; K562 KO and intersection are at or below random for
almost everyone. **"Beating BEAR-GRN" = per (dataset, GT): AUPRC above the
best listed method and above AUPRC_random, AUROC above 0.5 and the best listed,
mean over 5 seeds with the seed sd reported.**

## 6. Fair-comparison design for a supervised model (pre-registered before any MeVD-GRN result)

BEAR-GRN excluded supervised methods on purpose (peer-review file, response to
Reviewer 2 Comment 1): evaluating against the same label class used for
training "would inflate performance estimates", and "fairly evaluating these
supervised approaches would require a dedicated design (e.g., leave-one-TF-out
cross-validation using multi-cell-type ChIP-seq compendia)". Our design is
exactly that.

**6.1 Edge lists, not splits.** MeVD-GRN must emit one ranked edge list per
dataset, like every unsupervised method, so BEAR's code can score it
unchanged. A supervised model can do that only by **cross-fitting**:

- **TF-disjoint K-fold cross-fitting (K = 5)** over the union of TFs of all of
  the dataset's GTs (ChIP ∪ KO), degree-stratified and seeded
  (`bear_protocol.tf_disjoint_folds`). The same partition is used for every GT
  of a dataset (union/intersection share TFs with ChIP/KO).
- For fold k, **every label whose source TF is in fold k is removed from every
  label source** (all tiers, all networks), so a held-out TF contributes no
  positive and no negative anywhere in training. An assertion
  (`assert_tf_disjoint`) checks this on every run.
- The fold-k model scores all (fold-k TF, expressed gene) pairs. The 5 blocks
  are concatenated into one edge list over all TFs and scored once by BEAR's
  code. Every scored edge was predicted by a model that never saw its TF.
- **[decision]** K = 5 (not leave-one-TF-out): 150 K562 ChIP TFs would mean
  150 trainings per seed; 5 folds keeps 80% of TFs for training and is
  cheap. Macrophage has only 23 ChIP TFs, so 5 folds of 4-5 TFs.
- **[decision]** Raw sigmoid scores are concatenated across folds without
  re-calibration (all fold models share the same loss, negative ratio and
  schedule). A rank-within-fold variant is reported as a robustness check.

**6.2 Label regimes.**
- **L1 (headline): in-cell-type, TF-disjoint.** Training labels = the
  SC-MO-GRN-DB networks of the same cell type, restricted to the training
  TFs: K562 = RN117 (localization) + RN118 (perturbation), Macrophage =
  RN204 (localization only, single-tier), mouse = RN111 + RN112. These are
  the BEAR GT files themselves with the held-out TFs removed (BEAR's GTs are
  the RN networks filtered to the sample). Protocol = the project's
  recommended `all_at_once` over the available tiers.
- **L2 (robustness, "no cell-type labels"): compendium labels.** Training
  labels = SC-MO-GRN-DB's non-specific ChIP network (RN005 human / RN006 mouse)
  minus every TF that appears in any GT of the dataset. One model scores all
  GT TFs (no cross-fitting needed). This matches BEAR's stated use-case
  (cell-type-specific ChIP unavailable). Built: `configs/bear/mevd_fm_h384_L2.yaml`,
  `scripts/21 --regime L2` -> `processed/<ds>__L2` (TF nodes = RNA genes that
  are a compendium or GT source). K562: RN005 has 386,277 edges / 2,159 TFs;
  in the RNA universe 57,376 edges / 874 TFs; after removing the 204 GT TFs
  the training set is 18,072 positives over 680 TFs (+113 validation TFs).
  Note: RN005.tsv separates Target and Relationship with two spaces, not a
  tab, on many lines; `bear_data.read_gt_pairs` keeps the first whitespace token.

**6.3 Negatives.** **Degree-matched** (Stock et al. 2025): a training negative
(tf, g) is drawn with tf ~ out-degree and g ~ in-degree (+1 smoothing) of the
training positives, excluding every known positive of any label source.
Uniform TF x gene negatives are an ablation. The trainer's leak fix stays on
(`restrict_negative_pool` is called; the pool is built per fold from training
TFs only, so no held-out pair can be in it).

**6.4 Validation / model selection.** A TF-disjoint inner split: 15% of the
training TFs (degree-stratified) are validation TFs, with degree-matched
negatives at 1:5. Early stopping and checkpoint selection use only these.
The held-out fold is touched once, at scoring. **[decision]** Their labels are
otherwise not trained on (they are carved out of the training TFs, not added).

**6.5 Label-free inputs (no second leak).** `02_preprocess.py` builds the
TF-candidate graph with `exclude_positives=True` over **all** splits, so a
pair's absence from a TF's candidate list correlates with it being a held-out
positive (generalization doc s10). `21_bear_prepare.py` builds the
co-expression kNN graph and the TF-candidate graph with
`exclude_positives=False` and **no evidence argument at all**; the gene
universe is the post-QC RNA gene list; the TF node list is label-free with
respect to edges (genes that are a source in any SC-MO-GRN-DB network of the
species, ∩ RNA genes). The motif relation is off. A unit check
(`bear_prepare` summary `label_free: true`) asserts the graphs are identical
whether or not labels are loaded.

**6.6 What MeVD-GRN outputs.** Scores for every (GT-or-training TF in the RNA
gene list) x (RNA gene) pair. This is "dense" output: BEAR's AUROC is then the
full-space AUROC (on 1:1 subsamples), not a conditional one over a chosen
sparse edge set. **[decision]** Headline = dense. A sparse variant (top-N
edges, N = median in-space edge count of the 8 scored methods on that dataset)
is also scored, because BEAR's AUROC depends on output size.

**6.7 Mandatory trivial baselines, same cross-fitting and scoring.**
- **Target in-degree ranker:** score(tf, g) = in-degree of g among the
  training TFs' labels. It is the "supervised hub prior" that TF-disjoint
  supervision can still exploit; if MeVD-GRN does not beat it, the gain is
  not from the data.
- **|Pearson r|** of TF and target on the INPUT.DATA RNA (unsupervised).
- **GRNBoost2** (arboreto, dask-free per-target path of
  `src/benchmarks/scmgrn_features.py`) on the INPUT.DATA RNA, TFs as regulators.
- The 9 BEAR methods are **not re-run**: their released outputs
  (INFERRED.GRNS) are re-scored with our scorer (that is also the port's
  validation, section 8).

**6.8 Seeds and variants.** 5 seeds (42-46); the seed sets the fold
partition, the inner validation split, negatives and initialisation.
Variants: `fm_h384` (Geneformer + h384/l2, recommended config; headline),
`base` (h128/l2, no FM), `fm_h384_rna_only` (ATAC zeroed: modality ablation
matching BEAR's finding that methods are RNA-driven), `fm_h384_uniformneg`
(uniform negatives).

**6.9 Stability and modality shuffles (phase 2).** Stability: rerun the
headline variant on `INPUT.DATA.STABILITY` subsamples and compute median JI of
the top 10% edges with the ported code. Modality shuffles: apply
Cell_Shuffle / Gene_Shuffle / ATAC_Gene_Shuffle to INPUT.DATA and rescore.

## 7. Build log
- 2026-10-01 morning: doc started; Ada SSH (ada.iiit.ac.in) timed out.
- 2026-10-01: protocol + numbers extracted (helper agent; files in
  `~/.cache/bear/paper/`). Metric port written and checked against R on a
  synthetic tied-score vector. K562 INPUT.DATA + all K562 GT variants fetched
  locally (234 MB). Ada is back at `ada-gw1` (see `docs/reference/ada.md`
  top); Zenodo is ~10-50 KB/s per connection from both sites, hence the
  parallel sparse-zip downloader. Phase-1 fetch (K562, Mac-S1, Mac-S2: input,
  GT, inferred) started on the gateway into `/share1/$USER/mevd_grn/bear`.

## 8. Port validation against the released GRNs (2026-10-01)

The released `INFERRED.GRNS` files were re-scored with `scripts/24_bear_score.py`
(single draw) and compared with Supp Data 2-5 (`results/port_validation.csv`).

| dataset | GT | method | AUROC ours / paper | AUPRC ours / paper | evaluable edges ours / paper |
|---|---|---|---|---|---|
| K562 | ChIP | CellOracle | 0.5389 / 0.5388 | 0.3895 / 0.3895 | 215,375 / 215,375 |
| K562 | ChIP | LINGER | 0.5360 / 0.5358 | 0.4297 / 0.4297 | 707,264 / 707,264 |
| K562 | KO | CellOracle | 0.4868 / 0.4919 | 0.1467 / 0.1467 | 63,481 / 63,481 |
| K562 | KO | LINGER | 0.4027 / 0.4019 | 0.1390 / 0.1390 | 171,511 / 171,511 |
| K562 | Union | CellOracle | 0.5214 / 0.5212 | 0.2927 / 0.2927 | 256,948 / 256,948 |
| K562 | Union | LINGER | 0.5247 / 0.5248 | 0.3443 / 0.3443 | 794,485 / 794,485 |
| K562 | Intersection | CellOracle | 0.4567 / 0.4433 | 0.1058 / 0.1058 | 15,728 / 15,728 |
| K562 | Intersection | LINGER | 0.3574 / 0.3595 | 0.0989 / 0.0989 | 50,932 / 50,932 |
| Macrophage_S1 | ChIP | CellOracle | 0.5887 / 0.5844 | 0.2925 / 0.2925 | 28,611 / 28,611 |
| Macrophage_S1 | ChIP | LINGER | 0.5971 / 0.5957 | 0.3214 / 0.3214 | 101,232 / 101,232 |

- **AUPRC and every edge / positive count match exactly** (4 decimals) in all
  10 cases. For these GTs the 1:10 negative cap is not binding, so AUPRC is
  deterministic (sd over 30 draws = 0.0000).
- **AUROC** matches to 3 decimals in 6/10 cases; the rest differ by
  0.001-0.013. AUROC depends on the random 1:1 down-sampling (R `set.seed(42+i)`
  + `sample_n`, not reproducible in numpy). Over 30 draws the paper value lies
  within the port's draw distribution: z = (paper - mean)/sd = 1.77 (CellOracle
  KO, sd 0.0043), -0.23 (CellOracle Int., sd 0.0095), 0.83 (CellOracle ChIP),
  0.75 / 1.35 / -0.68 (LINGER KO / Int. / ChIP), and Mac-S1 CellOracle
  0.5882 +- 0.0016 vs 0.5844 (z = -2.4), LINGER 0.5961 +- 0.0010 vs 0.5957.
  **[decision]** The port is accepted as verified. Our methods are reported
  with `n_rep = 20` (mean of 20 draws) to remove this sampling noise; released
  methods keep the paper's single-draw numbers.
- The random-baseline AUPRC equals the paper's to 3 decimals for every GT
  (K562: ChIP 0.331, KO 0.160, Union 0.244, Intersection 0.119), confirming
  the universe construction.

## 9. Runbook

Data (gateway ada-gw1, resumable; ~2.4 GB for phase 1, done 2026-10-01):
```bash
cd ~/mevd_grn
python3 scripts/20_bear_download.py --root /share1/$USER/mevd_grn/bear \
    --datasets K562 Macrophage_S1 Macrophage_S2 --what gt input inferred --conns 32
# phase 2 (running): iPS + 4 mouse embryo + naive mESC; INPUT.DATA.STABILITY for phase-1 sets
# GTFs: /share1/$USER/mevd_grn/bear/gtf/gencode.v38.annotation.gtf.gz (hg38), gencode.vM25 (mm10)
```
Jobs (gateway, repo root):
```bash
bash slurm/bear_submit_all.sh        # prep (CPU) -> work[0-3]%2 (GPU, 2080 Ti) -> score (CPU)
# if the score job is rejected by the 20-submitted-jobs QoS cap, submit it when slots free up:
sbatch --dependency=afterany:<work job id> slurm/bear_score.sh
```
- `bear_prep.sh`: per dataset, `21_bear_prepare.py --fm` (label-free processed
  dir + Geneformer), `23_bear_baselines.py` (indegree x 5 seeds, pearson,
  coverage, grnboost2), `24_bear_score.py --released all --ours all` (re-scores
  all 9 BEAR methods; writes `results/port_validation.csv`).
- `bear_work.sh`: units = variant x seed x dataset, headline `fm_h384` first;
  `22_bear_train_eval.py --skip_existing`; each unit pushes
  `results/<ds>/mevd_<variant>/seed<s>/{grn.tsv.gz,meta.json}` to /share1.
- `bear_score.sh`: `24_bear_score.py --ours all --n_rep 20 --sparse_match`, then
  `--compile_only` -> `results/summary_tables.md`, `summary_agg.csv`.
- Local equivalent (laptop smoke): same scripts with `--root ~/.cache/bear/data`
  (`configs/bear/smoke.yaml`, `--max_epochs`, `--folds 0`).

Phase 2 (iPS, 4 mouse embryo sets, naive mESC; data on /share1 since
2026-10-01 16:10). **Status: cancelled, to be resubmitted after phase 1**
(see s11, 2026-10-01 ~18:30). `slurm/bear_queue_phase2.sh` (nohup on ada-gw1, log
`~/bear_tmp/queue_phase2.log`) waits for free QoS submit slots and submits
prep -> work array 0-1%1 -> score with `VARIANTS="fm_h384 fm_h384_L2"`,
`SEEDS="42 43 44"`. **[decision]** Phase 2 uses 3 seeds and 2 variants (~36
GPU-h) until phase-1 results justify the full grid; the mouse sets are as
large as K562 in TFs x genes and have 5-12x more cells.
Genomes (from BEAR's GENOME_AND_ANNOTATION.zip listing): human hg38 (GENCODE
v38, EnsDb v86), mouse **mm10** (EnsDb v79 mm10) -> we use GENCODE v38 / vM25.
Geneformer is a human model: for mouse, genes are matched by upper-cased
symbol (orthologs with identical symbols), so the FM channel is partial there.

Cost (measured locally on the RTX 4060): K562 fm_h384 ~20 s/epoch per fold
(3.26 M params, ~0.9 M training pairs), i.e. <= ~100 min per unit at the
60-epoch cap; base h128 ~5 s/epoch. Macrophage units are ~5x smaller (23 TFs).

## 10. Results

### 10.1 Trivial baselines change how BEAR's AUPRC must be read (K562, port)

indegree = mean (sd) over seeds 42-46 (TF-disjoint folds re-drawn per seed),
each scored with 20 AUROC draws; coverage and |Pearson r| are deterministic
(1 run, 20 draws). BEAR column = best method in Supp Data 2-5.

| K562 | coverage (all measured pairs = 1) | indegree (TF-disjoint, 5 seeds) | \|Pearson r\| | best BEAR method (paper) | random |
|---|---|---|---|---|---|
| ChIP AUROC | 0.500 | **0.652 (0.000)** | 0.506 | 0.575 (Pando-GLM) | - |
| ChIP AUPRC | 0.431 | **0.484 (0.000)** | 0.434 | 0.430 (LINGER) | 0.331 |
| KO AUROC | 0.500 | 0.569 (0.001) | 0.491 | 0.573 (SCENIC+) | - |
| KO AUPRC | 0.122 | 0.134 (0.000) | 0.121 | 0.159 (SCENIC+) | 0.160 |
| Union AUROC | 0.500 | **0.626 (0.000)** | 0.497 | 0.570 (Pando-GLM) | - |
| Union AUPRC | 0.303 | 0.344 (0.000) | 0.302 | 0.344 (LINGER) | 0.244 |
| Intersection AUROC | 0.500 | 0.550 (0.001) | 0.467 | 0.529 (Pando-GLM) | - |
| Intersection AUPRC | 0.094 | 0.100 (0.000) | 0.092 | 0.120 (GRaNIE) | 0.119 |

Findings (full computations, not smoke numbers):
1. **Coverage artifact.** Emitting a constant score for every measured
   (TF, gene) pair already gives K562 ChIP AUPRC 0.431, equal to the best
   method in the paper (LINGER 0.430). BEAR's universe contains all GT
   targets (27,747), but only 10,012 genes are measured; GT density among
   measured pairs is 46% vs 33% overall, so "being measured" is itself
   predictive. Any dense method inherits this floor. **[decision]** MeVD-GRN
   is only credited for AUPRC above the coverage baseline, and the sparse-
   matched variant (s6.6) is reported alongside.
2. **Hub prior.** The TF-disjoint in-degree ranker (no expression, no ATAC)
   beats every BEAR method on K562 ChIP (AUROC 0.652, AUPRC 0.484) and Union
   AUROC. This is exactly the InfoSEM / Stock et al. concern: supervision
   transfers target hubness across TFs. MeVD-GRN must beat **indegree**, not
   just the BEAR methods, for a data-driven claim.
3. On KO and Intersection nothing (BEAR methods or our baselines) is above
   random AUPRC on K562.
4. Top-10k precision of the coverage baseline (0.909) is meaningless: with all
   scores tied, `head(10000)` takes the file order. Tied outputs must not be
   compared on early metrics.

### 10.2 MeVD-GRN, first complete run (K562, headline fm_h384, seed 42 only; laptop RTX 4060)

Full budget, all 5 TF-disjoint folds, 4,201 s. Per-fold best val AUPR
(inner validation TFs) 0.461 / 0.589 / 0.599 / 0.515 / 0.504, stopped at
epochs 35 / 45 / 40 / 60 / 25. Scored with the port, 20 AUROC draws.
**n = 1 seed: a direction, not a result.** The 5-seed Ada runs (588) are the
real numbers.

| K562 | MeVD-GRN dense | MeVD-GRN rank-within-fold | MeVD-GRN sparse-matched* | indegree (5 seeds) | coverage | best BEAR (paper) | random |
|---|---|---|---|---|---|---|---|
| ChIP AUROC | 0.602 | 0.606 | 0.543 | **0.652** | 0.500 | 0.575 | - |
| ChIP AUPRC | 0.468 | 0.470 | 0.426 | **0.484** | 0.431 | 0.430 | 0.331 |
| KO AUROC | 0.539 | 0.541 | 0.432 | 0.569 | 0.500 | **0.573** | - |
| KO AUPRC | 0.123 | 0.124 | 0.146 | 0.134 | 0.122 | **0.159** | 0.160 |
| Union AUROC | 0.634 | **0.637** | 0.538 | 0.626 | 0.500 | 0.570 | - |
| Union AUPRC | 0.347 | **0.349** | 0.314 | 0.344 | 0.303 | 0.344 | 0.244 |
| Intersection AUROC | 0.355 | 0.346 | 0.434 | **0.550** | 0.500 | 0.529 | - |
| Intersection AUPRC | 0.089 | 0.089 | 0.095 | 0.100 | 0.094 | **0.120** | 0.119 |
| Core AUROC / AUPRC | 0.666 / 0.169 | | | 0.576 / 0.148 (seed 42) | 0.500 / 0.134 | (Supp Fig 12, not tabulated) | 0.105 |

\* top-N in-space edges, N = median in-space edge count of the released
methods scored locally so far (only LINGER and CellOracle; the Ada score job
recomputes it over all 8).

Reading (preliminary, 1 seed):
- On ChIP and Union, MeVD-GRN (dense) is above every BEAR method on both
  AUROC and AUPRC and above the coverage floor. **But on ChIP it is below the
  TF-disjoint in-degree ranker** (0.602 vs 0.652 AUROC, 0.468 vs 0.484 AUPRC).
  So its advantage over BEAR's methods is mostly the transferable target-hub
  prior that any TF-disjoint supervised model gets; the expression/ATAC part
  adds nothing measurable on top of it here. It is only marginally above
  indegree on Union (0.634 vs 0.626 AUROC, 0.347 vs 0.344 AUPRC).
- On KO and Intersection, MeVD-GRN is at or below random AUPRC, like every
  BEAR method; its Intersection AUROC (0.355) is among the worst. Training
  mixes ChIP-dominated labels (447k of ~450k positives come from RN117), so
  the model learns "binding hubs", which is anti-correlated with BEAR's
  ChIP-and-KO intersection.
- Re-calibrating scores within folds (rank-within-fold) changes nothing
  material, so the fold concatenation is not driving the numbers.
- Pre-registered next checks (already in the Ada grid): base model, L2
  compendium labels (no cell-type labels at all), RNA-only, uniform negatives.
  A planned but not yet built control is an indegree-residual model
  (MeVD-GRN score with the indegree rank regressed out), to test whether
  any signal beyond hubness exists.


### 10.3 Hub-controlled check: is there signal beyond target in-degree? (2026-10-08, K562, CPU)

**Question.** §10.2 showed MeVD-GRN below the TF-disjoint in-degree ranker on ChIP. Does it carry any
signal beyond that prior?

**Inputs (all local, no retraining).**
- MeVD-GRN `results/K562/mevd_fm_h384/seed42/grn.tsv.gz`, the §10.2 run. It is the only MeVD-GRN
  K562 run on disk; the 5-seed Ada runs are not reachable. **n = 1 seed.**
- `baseline_indegree/seed42`. Its fold partition is identical to the MeVD-GRN run's (`meta.json` folds
  equal), so each edge's in-degree is the target's in-degree among that fold's *training* TFs.
- `baseline_pearson/seed0`, and uniform random scores, as controls.

**Method.** Scripts are in `~/.cache/local_runs/scripts/`: `bear_hub_control.py`, `bear_per_tf.py`,
`bear_tf_level.py`. Outputs are in `~/.cache/local_runs/bear_hub/`.
1. **Residual scores, scored with BEAR's own metric** (`bear_metrics.score_roc_pr`, the scripts/24
   call, 5 draws, rng 43).
   - `resid_lin`: per fold, logit(score) minus its OLS fit on a cubic in log1p(in-degree).
   - `within_stratum`: per fold, the score's percentile within its in-degree stratum (0, then 20
     quantile bins of the non-zero in-degrees).
   - Both are shifted to be > 0, because BEAR takes |score| and unscored pairs are 0.
   - The same transforms are applied to |Pearson r| and to random scores. The random rows are the
     dense-output floor, i.e. the "coverage" floor of §10.1.
2. **Per-TF target ranking.** For each GT TF with ≥10 positives and ≥10 negatives among its in-space
   pairs, the AUROC of its own row: raw, and as a percentile within 20 in-degree bins of that row.
   This removes TF-level score differences, which the pooled BEAR metric keeps.

**Pooled BEAR metric (ROC on the method's in-space edges, PR on the full universe):**

| K562 GT | MeVD raw AUROC / AUPRC | MeVD in-degree residual (lin) | MeVD within in-degree stratum | in-degree ranker | \|Pearson r\| within stratum | random scores (dense floor) |
|---|---|---|---|---|---|---|
| ChIP | 0.602 / 0.468 | 0.512 / 0.440 | 0.519 / 0.435 | **0.652 / 0.484** | 0.502 / 0.432 | 0.500 / 0.431 |
| KO | 0.540 / 0.123 | 0.477 / 0.120 | 0.487 / 0.121 | **0.569 / 0.134** | 0.490 / 0.121 | 0.499 / 0.122 |
| Union | **0.634 / 0.347** | 0.565 / 0.322 | 0.572 / 0.320 | 0.626 / 0.343 | 0.495 / 0.300 | 0.500 / 0.303 |
| Intersection | 0.356 / 0.089 | 0.264 / 0.088 | 0.275 / 0.088 | **0.553 / 0.100** | 0.469 / 0.092 | 0.505 / 0.094 |
| Core | **0.665 / 0.169** | 0.642 / 0.159 | 0.644 / 0.154 | 0.577 / 0.148 | 0.473 / 0.128 | 0.498 / 0.134 |
| CellTypeExclusive | **0.634 / 0.348** | 0.565 / 0.324 | 0.572 / 0.322 | 0.626 / 0.345 | 0.495 / 0.303 | 0.500 / 0.305 |

(The raw MeVD-GRN row reproduces §10.2 exactly, e.g. ChIP 0.6017 / 0.4677. Uniform-universe random
AUPRC is 0.331 / 0.160 / 0.244 / 0.120 / 0.105 / 0.253.)

**Per-TF target ranking (mean AUROC over TFs):**

| K562 GT | TFs | MeVD raw | in-degree ranker | \|Pearson r\| | MeVD within in-degree bins | share of TFs > 0.5 (within bins) |
|---|---|---|---|---|---|---|
| ChIP | 129 / 150 | 0.720 | **0.765** | 0.510 | 0.523 | 67% |
| KO | 93 / 96 | 0.611 | **0.613** | 0.491 | 0.519 | 55% |
| Union | 203 / 226 | 0.707 | **0.725** | 0.505 | 0.543 | 67% |
| Intersection | 17 / 20 | 0.544 | **0.644** | 0.492 | 0.388 | 18% |
| Core | 17 / 24 | 0.633 | **0.652** | 0.499 | 0.530 | 65% |

**Within in-degree strata, pooled over all in-space pairs** (n-weighted mean AUROC over 20-21 strata):
- ChIP: MeVD 0.510, Pearson 0.502. MeVD's signal is concentrated in the highest-in-degree strata
  (0.535 and 0.548 in the top two; 0.497-0.515 elsewhere).
- Union: MeVD 0.571 (0.53-0.60 in every stratum), Pearson 0.494.
- Core: MeVD 0.645, Pearson 0.472.
- KO: MeVD 0.500.
- Intersection: MeVD 0.285.

**Reading (n = 1 seed):**
1. **For ranking each TF's targets, MeVD-GRN is a slightly noisier copy of the in-degree ranker.**
   - Per TF it is below in-degree on every GT.
   - Within in-degree bins it keeps a small, consistent positive residual on ChIP, KO, Union and Core
     (0.52-0.54 AUROC; about two thirds of TFs above 0.5). That is real but small.
   - On ChIP, the pooled AUPRC beyond in-degree is 0.435-0.440, against a dense floor of 0.431.
2. **The larger pooled residuals on Union / Core / CellTypeExclusive (AUROC 0.57-0.64)** come mostly
   from how the score scale differs *between* TFs: per TF the residual is only 0.53-0.54.
   - A TF's mean score does not track its GT density (Spearman -0.20 to +0.31, none significant), so
     this is not a simple TF-hub prior. The mechanism is unresolved.
   - It should not be described as target-level regulatory signal.
3. **Intersection is anti-correlated** beyond in-degree (pooled 0.26-0.28; per TF 0.39, only 3 of 17 TFs
   above 0.5). It is the clearest failure. The training labels are 99% ChIP (§10.2), so the model ranks
   "binding but not KO-responsive" targets above the ChIP∩KO ones.
4. |Pearson r| carries no signal beyond in-degree on any GT (0.47-0.50).
5. **Claim that survives:** on K562, MeVD-GRN's margin over BEAR's methods is the transferable
   target-hub prior, plus a small within-hub residual (per-TF AUROC +0.02 to +0.04 over chance within
   in-degree bins). It does not beat a TF-disjoint in-degree ranker at ranking any TF's targets. Repeat
   this on the 5-seed Ada outputs, and on Macrophage, before quoting it.

## 12. Pre-registration: "fairly beat LINGER" (written 2026-10-09 ~01:00 IST, before any new run)

New mandate (user via coordinator, 2026-10-08/09): BEAR-GRN is THE benchmark for MeVD-GRN; the
goal is to properly and fairly beat LINGER, BEAR's overall #1. Ada is locked, so everything runs on
the laptop (RTX 4060 8 GB shared with the PBMC/K562 queue, 22 GB RAM, 16 cores, ~6.7 GB free disk).
Everything below is fixed **now**, before any of the runs it governs. Changes after this point are
logged in s12.10 with a reason and are flagged in the results.

### 12.1 Coverage (item 1)
All 9 BEAR datasets: K562, Macrophage_S1, Macrophage_S2, iPS, mESC_E7.5_rep1/rep2,
mESC_E8.5_rep1/rep2, Naive_mESC. Every GT BEAR defines: ChIP for all; KO, Union, Intersection for
K562 and the 5 mouse sets. Core and CellTypeExclusive (human only, Supp Fig 12, no tables) are
scored as secondary. If a dataset cannot be processed locally (disk/RAM), it is reported as
missing with the reason, never silently dropped.

### 12.2 Development vs held-out datasets
- **K562 is the development dataset.** It has already been looked at (s10). The model choice in
  s12.6 is made on K562 only, and only with TF-disjoint inner-validation metrics (never BEAR scores).
- **The other 8 datasets are held out.** No MeVD-GRN run is started on them until the configuration
  is frozen (s12.6). K562 numbers are always reported, but flagged "dev".

### 12.3 Scoring (item 2)
The verified port (s8) of BEAR's code: AUROC (own in-space edges, 1:1), AUPRC (full GT-TF x GT-target
universe, unscored = 0, 1:10), AUPRC_random, top-10k precision/recall/F1 (on every GT, not only ChIP).
**[decision]** Every method, including the released ones, is scored with `n_rep = 20` draws of the
1:1 / 1:10 subsamples (mean), so AUROC differences are not RNG artefacts. The paper's single-draw
values are shown alongside for the released methods. Stability (Jaccard of the top 10% across
BEAR's 5 cell subsamples, ported code) and a modality-perturbation check are run where feasible
(s12.8); the stability inputs are `.rds` and need `pyreadr` (not installed); if that fails they are
reported as not done.

### 12.4 Comparison (item 3)
Against the **released** outputs of all 9 BEAR methods (Zenodo `INFERRED.GRNS`), scored by the
same code. LINGER is the named target; "best other method" = best of the other 8 per dataset x GT x metric.

### 12.5 No label overlap with evaluated TFs (item 4)
- **L1 (main):** TF-disjoint 5-fold cross-fitting (s6.1). Every scored TF is unseen in training;
  its labels are removed from every tier. Inner validation = 15% of the training TFs.
- **L2 (strict):** training labels = SC-MO-GRN-DB non-specific ChIP compendium only (RN005 human /
  RN006 mouse), with every TF of every GT of the dataset removed; one model scores all GT TFs.
  L2 is the closest analogue to LINGER (external data only, no cell-type labels).

### 12.6 Model (items 7, 8) and how the headline is chosen
Candidates, all with label-free graphs (`prior_exclude_positives: false`), the negative-leak fix
(`exclude_eval_negatives: true`), all_at_once over the dataset's label tiers, Geneformer + h384/l2:
- **M0** = `fm_h384` (s6, degree-matched negatives) - the configuration of s10.2.
- **M1** = M0 with **natural (uniform) negatives**. Reason: BEAR scores the natural TF x gene
  distribution, where target hubness is predictive; degree-matched negatives train that signal
  away (s10.3 shows MeVD-GRN below the in-degree ranker). The in-degree baseline (s12.7) is the
  guard against crediting hubness to the model.
- **M2** = M1 + **hub term**: the decoder adds `w * log1p(indeg_train(target))`, where indeg_train
  is the target's in-degree among the *current fold's training TFs* (the same quantity as the
  in-degree baseline; leak-free by construction).
- **M3** = M2 + **TF-motif pair features**: for each (TF, gene), the best JASPAR 2024 PWM log-odds
  score of the TF (relative to the PWM maximum) in the gene's top-5 RP-weighted accessible peaks,
  the number of those peaks with a hit >= 0.8 x max, and a has-motif indicator. Every BEAR method
  uses TF motifs; ChIP/KO GTs are not motif-derived, so this is not circular.
- M2/M3 are implemented by subclassing MEvDGRN in `src/benchmarks/bear_model.py` (no shared-file edits).

**Selection rule (fixed now).** On K562, seed 42, all 5 folds, compare M0-M3 by the mean over folds of
the **TF-disjoint inner-validation AUPR on the natural distribution** (all validation-TF x gene pairs;
labels = union of the label tiers). Choose the most complex model whose inner-val AUPR exceeds the
next simpler one by >= 0.005; otherwise keep the simpler. The chosen model is the **headline** and
is frozen for all datasets, seeds and regimes. Early stopping in every variant uses the same
natural-distribution inner-val AUPR. (This replaces s6.4's 1:5 degree-matched validation negatives,
which made inner-val AUPR incomparable across negative schemes.)

### 12.7 Mandatory baselines (item 5), on every dataset x GT
in-degree (5 seeds, the headline's fold partitions), coverage, |Pearson r|, GRNBoost2. **A win over
LINGER counts only if MeVD-GRN also beats in-degree and coverage** on the same metric. Hub control
(s10.3 method): MeVD-GRN's in-degree residual (`resid_lin`) and within-in-degree-stratum scores,
pooled BEAR metric, plus per-TF AUROC within in-degree bins.

### 12.8 Output network size (item 6), ablations (item 8), modality check
- **Primary rule (fixed): dense** - score every (TF node in the RNA gene list) x (RNA gene) pair,
  self-pairs excluded, as LINGER does (LINGER emits a dense TF x gene matrix).
- **Sensitivity:** the top-N in-space edges with N = LINGER's in-space edge count for that dataset x
  GT, and with N = the median of the 8 released methods.
- **Ablations** (headline vs one change, 3 seeds 42-44, on K562, Macrophage_S1, mESC_E7.5_rep1):
  without ATAC (`zero_atac`: ATAC features + openness zeroed), without Geneformer (`use_fm: false`,
  same h384/l2), sequential curriculum (ChIP stage then KO stage; K562 and mouse only) vs all_at_once.
- **Modality check analogous to BEAR's shuffles:** ATAC_Gene_Shuffle = permute the ATAC feature
  and openness rows across genes; Gene_Shuffle = permute the RNA feature/signature/FM rows; both at
  inference on the trained headline models. Cell_Shuffle is a no-op for MeVD-GRN by construction
  (each modality is summarised per gene before pairing), which is stated, not run.

### 12.9 Seeds, success criterion, honesty
- 5 seeds (42-46) for the headline (L1) and L2 on all datasets; seed sets the fold partition, inner
  split, negatives and initialisation. Mean (sd) reported.
- **Fair win over LINGER on a dataset x GT x metric** = MeVD-GRN 5-seed mean > LINGER, > in-degree
  mean, > coverage (and AUPRC > AUPRC_random). Wins, ties and losses are counted over all cells
  and all are reported, including losses.
- **LINGER's advantage, stated:** LINGER is pre-trained on external ENCODE bulk data (and uses motif
  priors). MeVD-GRN L1 instead uses the cell type's own SC-MO-GRN-DB labels for *other* TFs; L2 uses
  only a non-cell-type-specific ChIP compendium. Neither sees any label of an evaluated TF. The two
  information sources differ, and the comparison is reported as such.
- Compute order (GPU shared, one BEAR GPU job at a time, started only when free memory >= measured
  peak + 1 GB): (1) K562 M0-M3 selection, seed 42; (2) headline seeds 42-44 on all 9 datasets;
  (3) seeds 45-46; (4) L2; (5) ablations. Whatever is not finished is reported as not finished.

### 12.9b Implementation of s12 (2026-10-09, before any s12 GPU run)
- `src/benchmarks/bear_model.py`: `BearMEvDGRN` (subclass; hub term + motif pair term in `decode`).
- `src/benchmarks/bear_motif.py`: UCSC hg38 / mm10 streamed once into a 2-bit store (no 3 GB FASTA
  on disk); JASPAR 2024 CORE vertebrates (redundant set; `A::B` matrices count for both partners);
  vectorised torch conv1d scan, both strands, window-validity masked; peaks centre-cropped to 2 kb.
  Unit-tested on a planted motif (forward and reverse-complement hits = 1.0, short sequence = -1).
- `scripts/22`: natural-distribution inner validation; per-fold hub vector (fit TFs for training and
  validation, fit + validation TFs when scoring the held-out fold); sequential-curriculum ablation
  (one stage per label tier, config order); `--modality_check` (s12.8); records peak GPU memory and
  the zero-leak guard status per fold. CPU smokes passed for M2 (hub), sequential, modality check.
- `scripts/20 --no_extract` + `bear_data.read_csv_sparse('zip:...')`: INPUT.DATA members stay
  compressed and are parsed in row chunks (a 2.5 GB mouse ATAC CSV never lands on disk).
- `scripts/24`: `--headline` head-to-head table with the s12.9 WIN rule; sparse sensitivity at
  N = median of released and N = LINGER's in-space edge count; modality-shuffled variants scored
  as `<variant>@atacshuf` / `@rnashuf`.
- `src/benchmarks/bear_hub.py`: the s10.3 hub control, generalised (pooled residual / within-
  stratum BEAR scores + per-TF AUROC within in-degree bins).
- Drivers: `slurm/bear_local_dataset.sh <DS> [all|prep|baselines|released|motif]` (disk-frugal,
  aborts below 3 GB free), `slurm/bear_gpu_queue.sh` (one BEAR GPU job at a time through the
  shared gate).

### 12.10 Deviations log
- 2026-10-09 01:20: the s10.2 MeVD-GRN K562 run (old 1:5 validation, old processed dir) was moved
  to `~/.cache/bear/archive/K562_s10.2/` so that M0 is re-run under the s12.6 validation rule;
  its numbers stay in s10.2 / s10.3 and are not used for selection.
- GPU jobs go through the coordinator's shared gate (`~/.cache/local_runs/scripts/gpu_gate.sh`,
  flock-serialised) via `slurm/bear_gpu_queue.sh` (queue `~/.cache/bear/gpu_queue.txt`).

## 13. Results under the s12 protocol (laptop, 2026-10-09 onward)

Cells are AUROC / AUPRC / top-10k precision, scored with the port, 20 draws (s12.3). Released
methods are their Zenodo outputs. DIRECT-NET has no score column, so only its unsorted head(10000)
early metric exists (as in the paper). in-degree = mean over seeds 42-46 (sd in brackets).
Top-10k precision of `baseline_coverage` is not meaningful (all scores tied: head(10000) is file order).

### 13.1 K562 (development dataset): released methods and baselines

| method | ChIP | KO | Union | Intersection | Core | CellTypeExclusive |
|---|---|---|---|---|---|---|
| *random AUPRC* | 0.331 | 0.160 | 0.244 | 0.120 | 0.105 | 0.253 |
| LINGER | 0.536 / 0.430 / 0.587 | 0.401 / 0.139 / 0.067 | 0.525 / 0.344 / 0.530 | 0.355 / 0.099 / 0.038 | 0.507 / 0.134 / 0.153 | 0.525 / 0.348 / 0.529 |
| CellOracle | 0.539 / 0.390 / 0.675 | 0.482 / 0.147 / 0.050 | 0.522 / 0.293 / 0.521 | 0.444 / 0.106 / 0.046 | 0.491 / 0.123 / 0.158 | 0.522 / 0.299 / 0.521 |
| SCENIC+ | 0.465 / 0.332 / 0.503 | 0.580 / 0.159 / 0.048 | 0.477 / 0.245 / 0.436 | 0.578 / 0.119 / 0.014 | 0.518 / 0.104 / 0.035 | 0.476 / 0.254 / 0.434 |
| Pando-GLM | 0.576 / 0.348 / 0.682 | 0.477 / 0.156 / 0.051 | 0.570 / 0.261 / 0.620 | 0.508 / 0.114 / 0.047 | 0.519 / 0.105 / 0.110 | 0.570 / 0.269 / 0.620 |
| Pando-XGB | 0.518 / 0.346 / 0.570 | 0.490 / 0.156 / 0.052 | 0.522 / 0.259 / 0.527 | 0.468 / 0.114 / 0.045 | 0.434 / 0.105 / 0.110 | 0.522 / 0.267 / 0.527 |
| FigR | 0.504 / 0.376 / 0.562 | 0.457 / 0.149 / 0.050 | 0.497 / 0.288 / 0.473 | 0.437 / 0.111 / 0.073 | 0.474 / 0.113 / 0.127 | 0.497 / 0.295 / 0.473 |
| TRIPOD | 0.511 / 0.405 / 0.575 | 0.501 / 0.149 / 0.081 | 0.510 / 0.322 / 0.547 | 0.499 / 0.115 / 0.103 | 0.513 / 0.138 / 0.202 | 0.510 / 0.327 / 0.547 |
| GRaNIE | 0.512 / 0.368 / 0.553 | 0.509 / 0.152 / 0.076 | 0.526 / 0.281 / 0.538 | 0.466 / 0.120 / 0.126 | 0.500 / 0.112 / 0.153 | 0.526 / 0.288 / 0.538 |
| DIRECT-NET | - / - / 0.575 | - / - / 0.103 | - / - / 0.570 | - / - / 0.142 | - / - / 0.245 | - / - / 0.569 |
| baseline_indegree | 0.652 / 0.484 / 0.696 (sd 0.000/0.000) | 0.569 / 0.134 / 0.275 (sd 0.001/0.000) | 0.626 / 0.343 / 0.582 (sd 0.000/0.000) | 0.550 / 0.100 / 0.099 (sd 0.001/0.000) | 0.577 / 0.148 / 0.205 (sd 0.002/0.000) | 0.626 / 0.345 / 0.582 (sd 0.000/0.000) |
| baseline_coverage | 0.500 / 0.431 / 0.909 | 0.500 / 0.122 / 0.207 | 0.500 / 0.303 / 0.254 | 0.500 / 0.094 / 0.034 | 0.500 / 0.134 / 0.049 | 0.500 / 0.305 / 0.254 |
| baseline_pearson | 0.506 / 0.434 / 0.548 | 0.491 / 0.121 / 0.081 | 0.497 / 0.301 / 0.323 | 0.467 / 0.092 / 0.053 | 0.475 / 0.128 / 0.102 | 0.497 / 0.304 / 0.323 |
| baseline_grnboost2 | 0.496 / 0.380 / 0.481 | 0.502 / 0.135 / 0.092 | 0.494 / 0.269 / 0.307 | 0.494 / 0.101 / 0.060 | 0.464 / 0.112 / 0.110 | 0.494 / 0.275 / 0.307 |

Port re-check on all 8 scorable released methods x 4 GTs (`results/port_validation.csv`): AUPRC
and edge counts match Supp Data 2-5 exactly (GRaNIE edge counts differ by 10 of ~182k, from "NA"
gene-name parsing; its AUPRC still matches to 4 decimals); AUROC differences are subsampling noise,
largest where there are very few positives (SCENIC+ Intersection: 3 positives).

Reading: on K562 ChIP the coverage floor (0.431) already equals LINGER (0.430) and the TF-disjoint
in-degree ranker (0.652 / 0.484) beats every released method on ChIP and Union AUROC/AUPRC. On KO
and Intersection every method, released or trivial, is at or below random AUPRC.

## 11a. Ada status, 2026-10-09 (v2, GPU work moved to Ada, then blocked again)
What was done on Ada on 2026-10-09 ~15:30-16:05 IST, while ssh worked, and nothing else:
- **No jobs were submitted.** No job IDs exist for the v2 plan.
- `/share1/$USER` is at its **file-count quota** (3189 files used, soft 3000 / hard 3200, no grace), so
  `mkdir` there fails with "Disk quota exceeded". Space is fine (55 of 100 GB). The new work therefore
  uses **`~/bear_v2` (home, visible on compute nodes)**, not /share1. The old `/share1/$USER/mevd_grn/bear`
  (894 files) and `pbmc` (1264 files) were not touched; tarring or deleting files there is the way to
  free inodes if /share1 is needed again.
- Created `~/bear_v2/{processed,results,claims,done,failed,logs}` and started
  `rsync -a ~/.cache/bear/data/processed -> ada-gw1:bear_v2/` (1.1 GB, 18 dirs). The laptop rebooted
  before it was confirmed to finish (43 MB were on Ada at the last check), so **the transfer is
  incomplete and must be re-run** (`rsync -a --partial ~/.cache/bear/data/processed vishakkashyap.k@ada-gw1:bear_v2/`);
  it also has to be repeated after the motif features of the other datasets are built.
- Synced the lane files to `~/mevd_grn` (scripts 20-24, `configs/bear/`, `src/benchmarks/bear_*.py`,
  this runbook, `slurm/bear_ada_*.sh`); `scripts/22` is byte-identical to the laptop copy (md5 9a2d3f6f...).
- Written but never run on Ada: `slurm/bear_ada_common.sh` (verified 5x-retried copies, no gateway hop),
  `bear_ada_dev.sh` (array 0-3 = M0-M3 on K562 seed 42), `bear_ada_steal.sh` (work-stealing worker over
  `~/bear_v2/units.txt`, atomic `mkdir` claims), `bear_ada_reset.sh`, `bear_ada_score.sh`.
- Home usage was 41 of 50 GB on 2026-10-09 (envs 5.8 GB, miniforge 12 GB, another project 7.1 GB).
- Ada is blocked again from ~16:00 ("Your password has expired", no TTY). Do not retry in a loop.
  Plan B is the laptop queue below; the pre-registration (s12) is unchanged either way.

## 11b. Local GPU queue (the active compute path while Ada is blocked)
Relaunch (always under an inhibitor; the laptop suspends itself otherwise), from the repo root:
`setsid nohup systemd-inhibit --what=sleep:idle:handle-lid-switch --who="MeVD-GRN BEAR" --why="BEAR GPU queue" --mode=block bash slurm/bear_gpu_queue.sh ~/.cache/bear/gpu_queue.txt > ~/.cache/bear/gpu_queue.out 2>&1 < /dev/null &`
The queue file lists one job per line (`<need_MiB> <command>`); finished lines are in `gpu_queue.txt.done`
(a line is added even when the job failed, so check `rc=` in `gpu_queue.out` and remove a failed line from
`.done` to retry). `scripts/22 --skip_existing` makes completed runs a no-op.

## 11. Ada status
- 2026-10-01 ~15:40 IST: first submission 572/573, cancelled while pending to
  add regime L2.
- 2026-10-01 ~16:05 IST: **resubmitted** from ada-gw1 `~/mevd_grn`
  (`bash slurm/bear_submit_all.sh`): `bear-prep` **587** (CPU, 10 cores),
  `bear-work` **588** (array 0-3, max 2 running, 1x 2080 Ti + 10 cores each,
  `afterok:587`; units = {fm_h384, base, fm_h384_L2, fm_h384_rna_only,
  fm_h384_uniformneg} x seeds 42-46 x {K562, Macrophage_S1, Macrophage_S2}).
  `bear-score` was rejected by `QOSMaxSubmitJobPerUserLimit` (the user's
  20-submitted-jobs cap is full with the PBMC / leak-fix arrays). A watcher on
  ada-gw1 (`~/bear_tmp/submit_score_when_free.sh`, log `~/bear_tmp/submit_score.log`)
  submits `sbatch --dependency=afterany:588 slurm/bear_score.sh` as soon as a
  slot frees. Manual fallback: the same command.
  Jobs pend on `QOSMaxCpuPerUserLimit` (40 CPUs per user, used by the PBMC
  jobs, which have GPU priority).
- 2026-10-01 16:11: the score watcher submitted `bear-score` **615**
  (`afterany:588`) and exited.
- 2026-10-01 16:39-18:19: the phase-2 watcher submitted chain 2:
  prep **660** -> work **774** (array 0-1%1) -> score **816**.
- 2026-10-01 ~18:30: **cancelled chain 2 (660, 774_[0-1], 816)** at the
  coordinator's request. Chain 2 was not a duplicate (phase 2 = iPS + mouse,
  disjoint from chain 1's K562 / Macrophage), but together the two chains held
  10 of the user's 20 QoS submit slots, which made the PBMC LINGER baseline
  submission fail with `QOSMaxSubmitJobPerUserLimit`. The rule from now on:
  **BEAR holds at most 6 submitted job records** (array tasks count one by one).
  Chain 1 (587 prep, 588_0-3 work, 615 score = exactly 6 records) is kept.
  Both gateway watchers have exited; nothing will resubmit automatically.
  Phase 2 is to be resubmitted only after chain 1 has left the queue, keeping
  the same 6-record limit, e.g. (on ada-gw1, repo root):
  `export DATASETS="iPS mESC_E7.5_rep1 mESC_E7.5_rep2 mESC_E8.5_rep1 mESC_E8.5_rep2 Naive_mESC" VARIANTS="fm_h384 fm_h384_L2" SEEDS="42 43 44"; WORK_CONC=2 bash slurm/bear_submit_all.sh`
  (1 prep + 4 work + 1 score = 6 records; set the prep time with
  `sbatch --time=1-00:00:00` if the mouse ATAC parsing needs more than 8 h).
- 2026-10-01 evening: the coordinator put `bear-prep` 587 on `scontrol hold`
  (with PBMC 474_3 / 475) so the PBMC LINGER baseline (job 819) gets the next
  free CPU slot; `~/linger_first.sh` on ada-gw1 (log `~/logs/linger_first.log`)
  releases them automatically when 819 starts. **Do not release or resubmit
  587 by hand.**
- Local laptop: the first K562 fm_h384 seed-42 run died when the laptop
  rebooted (~18:23); rerun detached, finished 21:40 (s10.2). **[decision]**
  Local results are NOT pushed to /share1: Ada recomputes everything so all
  numbers come from one processed dir and one environment.

## 14. Dev selection finished (2026-10-10 03:00 IST)

All four candidates trained on K562, seed 42, 5 TF-held-out folds, zero-leak guard passing on every fold.
M0 and M1 ran on the laptop RTX 4060 (about 70 min each); M2 and M3 ran on Ada 2080 Ti nodes (jobs 14388_2/3,
about 35 min each; peak GPU about 4.8 GB). `python -m src.benchmarks.bear_select` applied the pre-registered rule
(section 12.6):

| Candidate | Mean inner-validation AUPR | vs predecessor |
|---|---|---|
| M0 fm_h384 | 0.5255 | n/a |
| M1 uniformneg | 0.5447 | +0.0192 pass |
| M2 hub | 0.5387 | -0.0060 fail |
| M3 hubmotif | 0.5532 | +0.0145 pass |

**HEADLINE = M3.** BEAR-scoring (20 draws, our verified port), K562, seed 42:

| GT | M3 AUROC / AUPRC | M2 AUROC / AUPRC | random AUPRC |
|---|---|---|---|
| ChIP | 0.6099 / 0.4761 | 0.6053 / 0.4740 | 0.3311 |
| Union | 0.6573 / 0.3609 | 0.6481 / 0.3570 | 0.2441 |
| KO | 0.5590 / 0.1257 | 0.5611 / 0.1265 | 0.1600 |
| Intersection | 0.2829 / 0.0879 | 0.2765 / 0.0878 | 0.1205 |

Reference (same GT): LINGER ChIP 0.536 / 0.430, Union 0.525 / 0.344; target-count baseline ChIP 0.652 / 0.484, Union 0.626 / 0.343.
M3 has hub terms derived from training-TF ChIP labels, so it belongs to the supervised reference track.
Incidents on the way: first Ada submissions died on a full node `/scratch` (gnode082) and on a stage-in bug
(rsync did not create the parent dir; fixed in f745c2f); the local queue sat blocked for 4 h on a 6800 MiB gate
threshold the desktop's GPU use made unreachable (now 5800).
