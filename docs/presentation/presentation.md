# MeVD-GRN: progress, 27 Sept to 9 Oct 2026

Advisor meeting, 10 Oct 2026, 11:00. Slides are separated by `---`.

1. **Part 1:** MeVD-GRN vs scMultiomeGRN (zero-leakage versions)
2. **Part 2:** 10x PBMC multiome results, and why we ran them
3. **Part 3:** What the BEAR-GRN paper found, and what we found auditing it

Figures are in `figs/` (regenerate with `python make_plots.py`; every number is read from the result files).
Everything shown is **zero test-to-train leakage** unless a row says "with leaks".

---

# PART 1. MeVD-GRN vs scMultiomeGRN

---

## 1.1 The problem, and the model in one picture

**Which transcription factor (TF) regulates which gene?** Predict TF → gene links from paired single-cell RNA + ATAC.

![MeVD-GRN architecture](../figures/architecture_simple_2026-10-09.png)

- Two graph encoders (RNA tower, ATAC tower) over label-free gene graphs, plus a pretrained Geneformer prior.
- A gated score for each TF-gene pair: chromatin openness acts as a gate.
- Trained on ChIP-seq, then knockout evidence; the overlap of the two is held out as a zero-shot test.

---

## 1.2 First, a correction: two leaks, both closed

![Effect of the leaks](figs/fig2_leak_effect.png)

- **Leak 1:** the trainer drew negatives from a pool that also held the validation and test negatives.
- **Leak 2:** input graphs were built with held-out positives removed, so a pair's absence from the graph hinted at its label.
- Training now **refuses to run** unless both are closed, and every result file records it.
- Switching the leaks back on reproduces the old paper numbers exactly, so the drops come from the fixes alone.
- **Every number in the rest of this deck is zero-leak.** The old headline model lost 0.11 AUPR on perturbation.

---

## 1.3 The two models side by side

![Comparison](../figures/comparison_scmultiomegrn_vs_mevd.png)

- Shared: separate RNA and ATAC pathways (GraFRank idea), class-weighted loss, regulatory-potential ATAC features.
- Biggest differences: **what is predicted** (TF-TF vs TF → any gene), **what the graph is** (known edges vs label-free), **where modalities meet** (encoder vs decoder), **what the labels are** (motif scans of the same ATAC vs experiments).

---

## 1.4 Where each idea comes from

| Component | Source | Status |
|---|---|---|
| Separate RNA / ATAC pathways | GraFRank (WWW 2021), used by scMultiomeGRN | **borrowed** |
| Regulatory-potential ATAC score | MAESTRO / BETA | **borrowed** (scMultiomeGRN uses it too) |
| GraphSAGE over several relations | GraphSAGE, R-GCN | adapted |
| Geneformer gene prior | Geneformer (Nature 2023); scRegNet did it first | adapted, not new |
| ATAC-gated bilinear decoder | PECA / CellOracle biology | adapted; no accuracy gain, so offered for interpretability |
| Label-free input graphs | contrast with scMultiomeGRN, GENELink, scRegNet, which pass messages over known edges | design choice |
| **Evidence-tier protocol:** train on ChIP and knockout, hold out their overlap | no prior method found that trains on separate evidence types | **our main claimable idea** |
| Geneformer helps in-domain but hurts transfer to a new cell type | scRegNet reported only the in-domain gain | **finding**; numbers predate the fixes, to be re-run |

---

## 1.5 K562 results: same data, same splits, zero leakage

![K562 results](figs/fig1_k562_zero_leak.png)

AUPR / AUROC, MeVD-GRN = mean of 4 seeds (std at most 0.012):

| Model | Localization | Perturbation | Dual evidence (zero-shot) |
|---|---|---|---|
| **MeVD-GRN, evidence curriculum** | 0.730 / 0.751 | **0.814 / 0.960** | **0.954 / 0.990** |
| MeVD-GRN, single stage | **0.953 / 0.947** | 0.582 / 0.883 | 0.902 / 0.979 |
| scMultiomeGRN (our re-run, 1 seed) | 0.906 / 0.914 | 0.655 / 0.919 | 0.847 / 0.971 |
| GRNBoost2 | 0.535 / 0.499 | 0.206 / 0.550 | 0.173 / 0.525 |

- The curriculum wins perturbation (+0.16 AUPR) and the zero-shot tier (+0.11), and loses localization (it forgets the first tier).
- The single-stage model wins localization but loses perturbation.

---

## 1.6 Where scMultiomeGRN is stronger, and the caveats

**What it does that we don't**
- **Per-cell features and pairwise edge features** (a 16 x 16 joint expression histogram); we reduce each gene to a few numbers.
- **Uses the known network as input**; we withhold it on purpose.
- **Outputs a final network** by a 10-run vote, with robustness and biology case studies (SPI1 in monocytes, etc.); we output rankings only.
- **Trained to convergence** (up to 2,000 epochs) with released code; we train 45 epochs.

**Caveats on our comparison**
- The scMultiomeGRN row is **our re-implementation, one seed**, not their released code. It was not run on their own benchmark (fetal lung), whose labels are motif scans of the same ATAC it reads.
- Dual evidence is **not a clean test set**: some design choices used it.
- Both models are leak-free here. The old paper claim of winning 5 of 6 metrics falls to 3 of 6.

---

# PART 2. 10x PBMC multiome

---

## 2.1 Why we ran on PBMC

- **K562 alone is not enough.** All K562 numbers come from one resource (SC-MO-GRN-DB). A reviewer asks: does it hold on data and ground truth we did not build?
- **A ~130-paper survey** (journals, ML venues, preprints) found one benchmark used most: the 10x PBMC multiome, scored against Cistrome ChIP-seq as in LINGER (17-22 papers). scTFBridge (Nat Commun 2025), KEGNI (Genome Biol 2025) and regX reuse LINGER's exact setup.
- **It is truly paired** (RNA and ATAC from the same nuclei) and **independent** of SC-MO-GRN-DB.
- **It compares us with LINGER**, the strongest published method, on LINGER's own ground truth.
- **It forced a fair supervised protocol**: our model uses labels, LINGER does not, so test TFs must be kept out of training.
- scMultiomeGRN's own benchmark was not usable for this: circular labels and data probably not truly paired.

---

## 2.2 PBMC10k: the task and the bar

| | |
|---|---|
| **Data** | 10x PBMC granulocyte-sorted 10k Multiome; LINGER's 9,543 labelled cells in 4 types: classical monocytes 1,848, naive CD4 T 1,373, naive B 282, myeloid DC 232 |
| **Ground truth** | 20 Cistrome ChIP-seq datasets = 10 TFs (MYC, RUNX1, IRF4, STAT1, IRF1, ETS1, FOXP3, CTCF, REST, SPI1) x 4 cell types |
| **Metric** | positives = a TF's top-1,000 genes by ChIP regulatory potential; score its gene ranking with AUROC and AUPR ratio; mean over the 19 evaluable datasets |
| **The bar** | LINGER (Nat Biotechnol 2024): **AUROC 0.7143, AUPR ratio 2.2526**. SCENIC+ 0.548 / 1.29, GENIE3 0.539 / 1.17 |

---

## 2.3 PBMC10k: our fair protocol (fixed before any result existed)

- **Training labels:** CollecTRI (literature-curated), DoRothEA A-B as a second source. **All 10 evaluation TFs removed as regulators**: a test TF is never seen in training.
- **Negatives:** degree-matched, so a model cannot win by learning popular targets.
- **Inputs:** gene graphs built without labels; motif relation off; early stopping on validation only; 5 seeds; one model per cell type; zero-leak guard on every run.
- **Variants (220 runs on Ada):** full model; no ATAC; no Geneformer; neither; uniform negatives; DoRothEA labels; random-edge and held-out-target splits.
- **Baselines under the same code:** target-count ranker, gene-ID regression, Pearson. LINGER's published numbers are used (our re-run failed: it needs ~120 GB of memory and a missing package).
- **Win rule written in advance:** beat 0.714 and the simple baselines. *(Too weak: see 2.5.)*

---

## 2.4 PBMC10k: results

![PBMC10k ablation](figs/fig4_pbmc10k.png)

| 5 seeds, TFs held out | AUROC | AUPR ratio |
|---|---|---|
| **MeVD-GRN** | **0.734 ± 0.005** | 2.135 |
| LINGER (published) | 0.714 | **2.253** |
| MeVD-GRN, no ATAC | 0.720 ± 0.030 | 2.178 |
| MeVD-GRN, no Geneformer | 0.670 ± 0.008 | 1.722 |
| MeVD-GRN, neither | 0.544 ± 0.011 | 1.337 |
| MeVD-GRN, uniform negatives | 0.597 ± 0.014 | 1.698 |
| MeVD-GRN, DoRothEA labels | 0.658 ± 0.018 | 1.933 |
| Target count / gene ID / Pearson | 0.579 / 0.578 / 0.584 | 1.4-1.6 |

- ATAC adds **+0.126** AUROC without Geneformer, only **+0.015** with it (within noise).
- Degree-matched negatives are essential (0.734 vs 0.597); CollecTRI beats DoRothEA.

![By cell type](figs/fig5_pbmc_per_celltype.png)

- The average win comes from **monocytes** (10 of 19 datasets): 0.760 vs 0.711. On naive CD4 T and naive B, MeVD-GRN is below LINGER.

---

## 2.5 PBMC10k: is the signal TF-specific? (no)

![TF specificity](figs/fig6_pbmc_tf_specificity.png)

- A ranker that gives **every TF the same gene ranking**, using ATAC accessibility alone, scores **0.814**, above MeVD-GRN and LINGER. The metric rewards open chromatin.
- **A dataset's own TF row scores no better than a different TF's row** (0.690 vs 0.685; seed-42 local models). Within accessibility strata it falls to 0.599.
- So **"0.734 vs 0.714" is not evidence of TF-specific regulation**, for either method. LINGER's own rows are unchecked: that needs the LINGER re-run.
- Held-out target genes: the original split scored 0.409 (below random) because of its negative design; the fixed split scores 0.745.
- **Open:** a local rerun of the headline gives 0.690, not 0.734 (cause not found); scoring used all expressed genes, not LINGER's own candidates.

---

# PART 3. BEAR-GRN

---

## 3.1 What BEAR-GRN is

**"Systematic assessment of single-cell multi-omics-based gene regulatory network inference methods"**, Karamveer, Moeller, Valensi, Manful, Uzun. *Nature Communications*, 18 Sept 2026 (accelerated article preview: main text not online yet; we read the abstract, all 12 supplementary files, the peer-review file, and the released code and data).

- Same lab that built **SC-MO-GRN-DB**, our training resource.
- Benchmarks **9 methods** on **9 paired multiome datasets** against **four kinds of ground truth**.
- Everything is released: data (Zenodo 20704929), code (GitHub UzunLab/BEAR-GRN, Zenodo 22015920) and an R package.
- **No paper cites it yet** (Crossref, Semantic Scholar, Europe PMC: 0, three weeks after release).

---

## 3.2 Datasets and ground truths

| Dataset | Cells after filtering | Source |
|---|---|---|
| K562 | 411 | GSE178707 |
| Macrophage (2 replicates) | 474 / 528 | PRJNA1102756 |
| iPSC | 5,710 | GSE278751 |
| Mouse embryo E7.5 (2 reps) | 5,064 / 1,447 | GSE205117 |
| Mouse embryo E8.5 (2 reps) | 6,557 / 6,252 | GSE205117 |
| Naive mESC (unpaired, paired computationally) | 3,568 | GSE198730 |

**Ground truths:** ChIP-seq (all datasets); knockout/knockdown, the **union** and the **intersection** of the two (K562 and mouse only); plus a core network shared across human cell types.

- ChIP and knockout agree on only **1.3-1.4%** of edges, so the choice of ground truth changes the story.
- Data-quality issues found in the released ground truths: GFP and epitope tags listed as TFs, a stray header row.

---

## 3.3 How methods are scored

- **Methods:** LINGER, CellOracle, TRIPOD, GRaNIE, FigR, SCENIC+, Pando (two variants), DIRECT-NET (has no score, so no AUROC/AUPRC). scGLUE and STREAM were run but excluded.
- **AUROC:** only on the edges a method itself predicted, among ground-truth TFs and targets, balanced 1:1.
- **AUPRC:** over *all* ground-truth TF x target pairs, with **unpredicted pairs scored 0**, so a method is rewarded for covering many measured pairs.
- **Also:** top-10k precision / recall / F1, stability across subsamples (Jaccard), runtime and memory, and shuffling RNA or ATAC to test whether the second modality matters.
- Our Python port reproduces BEAR's published AUPRC **exactly** on the methods we re-scored.

---

## 3.4 What BEAR-GRN found

- **Accuracy is low across the board.** On ChIP no method exceeds 0.66 AUROC. On the K562 knockout and intersection ground truths **every method is at or below random AUPR.**
- **LINGER ranks first overall** (overall score 5.29), then DIRECT-NET (5.17) and CellOracle (4.93). SCENIC+ ranks last.
- **LINGER is also the most stable**, but the heaviest: about 8.5 hours and 89 GB of memory at 5,000 cells.
- **The second modality barely matters.** Shuffling RNA or ATAC changes AUPR by about 0.02-0.03 for most methods. Only LINGER drops clearly (K562 0.430 to 0.342).
- **Ground-truth choice changes rankings substantially**, and methods are described as predominantly RNA-driven.

---

## 3.5 Our audit: what a trivial baseline scores on BEAR-GRN

![BEAR-GRN ChIP AUPR](figs/fig3_bear_chip_auprc.png)

- **Predicting every measured pair** (no model, no labels) gets within 0.001 of LINGER on K562 (0.431 vs 0.430).
- **A ranker that scores genes by how many TFs regulate them** beats every published method on all 9 datasets. It uses ChIP labels from other TFs, so it is a supervised reference, not an unsupervised competitor.
- **First MeVD-GRN result (K562, model M0, one seed, 5 TF-held-out folds):**

| K562 | MeVD-GRN M0 | LINGER | Target-count baseline |
|---|---|---|---|
| ChIP AUROC / AUPRC | 0.602 / 0.468 | 0.536 / 0.430 | 0.652 / 0.484 |
| Union AUROC / AUPRC | 0.635 / 0.347 | 0.525 / 0.344 | 0.626 / 0.343 |

---

## 3.6 Scope: BEAR-GRN is for methods that use no ChIP labels

From the authors' reply to reviewers (peer-review file):

> "BEAR-GRN targets settings where cell-type-specific ChIP-seq labels are unavailable."
> "Fairly evaluating these supervised approaches would require a dedicated design (e.g., leave-one-TF-out cross-validation using multi-cell-type ChIP-seq compendia), which we leave for future work."

- **Our current MeVD-GRN is supervised, so it is outside BEAR's stated scope**, and so is the target-count baseline.
- Our leave-TF-out protocol is the design they name, but our main regime trains on the same cell type's ChIP for other TFs. The in-scope regime (train on non-specific compendium only) is planned, not run.
- Safe wording: *"under the leave-TF-out design BEAR's authors propose"*, not *"we beat LINGER on BEAR"*.
- Not verified: the final Methods text, since the paper has no body online yet.

---

# Where this leaves us

---

## 4.1 What we can and cannot claim

**Supported**
1. Two leaks found and closed, enforced in code.
2. With zero leakage, the evidence curriculum beats our re-run of scMultiomeGRN on perturbation and the zero-shot tier on K562.
3. Benchmark audit: trivial per-gene rankers match or beat published methods on both PBMC10k (ATAC accessibility alone: 0.814) and BEAR-GRN (predict every pair; target count).

**Not supported**
- A PBMC10k win that reflects TF-specific signal (a different TF's row scores the same).
- MeVD-GRN beating LINGER on BEAR-GRN in the in-scope regime.
- Anything on scMultiomeGRN's own benchmark.

---

## 4.2 Proposed direction: self-supervised, TF-specific by construction

- **Objective:** predict each gene's expression **in each cell** from candidate-TF expression, gated by that cell's chromatin accessibility around the gene. The sparse TF → target weights are the network. A gene-level popularity prior cannot explain cell-to-cell variation, so this forces TF-specific signal.
- **Label-free evidence tiers** replace ChIP/knockout: motif in an accessible peak, then peak-gene co-accessibility, then expression dependence. This keeps the curriculum idea.
- **Pretrain on larger multiome sets** (PBMC, bone marrow) and transfer; K562 has only 411 cells, and LINGER uses atlas-scale pretraining.
- **Keep:** graph encoders, ATAC-gated decoder, Geneformer (declared as external pretraining), the zero-leak rigor and the TF-swap / baseline controls.
- **Evaluate:** BEAR-GRN on all 9 datasets and 4 ground truths, plus the trivial baselines and the TF-swap test, plus a second label-free benchmark (geneRNIB).

---

## 4.3 Status and timeline

| Job | Status | ETA |
|---|---|---|
| BEAR-GRN candidates M0, M1 (laptop) | done, ~70 min each | done |
| M2 (hub term), M3 (hub + motif features) on **Ada**, one GPU each | running, 2 of 5 folds at 22:20 | ~23:00 |
| Model selection (K562 validation only) | after M2 and M3 | ~23:10 |
| K562 seed 46 and PBMC held-out-target seeds (laptop) | running | overnight |

| Plan to mid-January | Dates |
|---|---|
| Label-free baselines and TF-swap audit on all 9 BEAR datasets (CPU) | Oct 10-24 |
| Build the self-supervised model; develop on K562 and one macrophage set only | Oct 24-Nov 21 |
| **Decision gate:** TF-specific signal beyond priors, or fall back to an audit-centred PLOS ONE paper | Nov 7 |
| Freeze; 5 seeds x 9 datasets on Ada | Nov 21-Dec 12 |
| Ablations, second benchmark, writing (ICML 2027 mid-January; PLOS ONE by February) | Dec 12-Jan 15 |

---

## 4.4 Decisions for the advisor

1. **Scope:** report BEAR-GRN results only in the leave-TF-out framing (supervised reference), and pursue a self-supervised model for the in-scope claim. Agree?
2. **Venue:** ICML 2027 needs ML novelty (benchmark audit + a TF-specific self-supervised objective). PLOS ONE judges soundness and fits what we have already. Which first?
3. **scMultiomeGRN's own benchmark:** the pipeline is built; run it for a head-to-head, or skip? (Circular labels, probably unpaired data.)
4. **Headline model if we keep the supervised track:** curriculum (wins perturbation and zero-shot) or single stage (wins localization).

---

## Appendix A. Reproducibility

- Branch `leak-fix-and-benchmarks`; key commits: `c26c03d` (leak 1), `d854fb6` (zero-leak guard), `8ed78ea` (BEAR pipeline), `0a46cf3` (selection tools).
- Runbooks: `docs/experiments/` (`leakfix_rerun.md`, `labelfree_graph_rerun.md`, `pbmc10k_linger_benchmark.md`, `bear_grn_benchmark.md`); comparison study: `docs/mevd_vs_scmultiomegrn.md`.
- Plots: `docs/presentation/make_plots.py`. Diagrams: `docs/figures/make_architecture_simple.py`, `make_comparison_diagram.py`.
- Figures are drawn locally: the BioRender connector is waiting on organisation approval.
