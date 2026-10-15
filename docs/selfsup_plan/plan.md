# Label-free GRN inference for BEAR-GRN: decision plan (10 Oct 2026)

Status: complete draft, 10 Oct 2026 ~04:30 IST. The executive summary is final; s1-s10 and Appendix A are filled. Items still open are marked "week-1 task".
Evidence labels used throughout: **[V]** = verified, with the primary source given; **[C]** = computed by us
(the exact command is in Appendix A); **[I]** = inference or judgement; **[G]** = guess, with low confidence.

---

## 0. Executive summary (one page)

**What we found before planning anything.** BEAR-GRN's headline numbers can be matched with *no
TF-specific information at all*. Any plan that does not start from this will produce a claim we
cannot defend.
- **AUPRC mostly measures coverage.** Emitting every measured TF x gene pair with one constant score
  ("coverage") gives ChIP AUPRC at or above the best published method on 8 of 9 datasets (the exception
  is Macrophage S1: 0.331 vs TRIPOD's 0.332) [C, `results/summary_tables.md`]. Pairs a method does not
  emit score 0, and measured genes are about twice as dense in positives as the full ground-truth
  universe (K562: 0.498 vs 0.331) [C].
- **Gene-level priors already beat LINGER on K562.** Ranking targets by their mean expression, the same
  ranking for every TF, scores K562 ChIP AUROC/AUPRC 0.554/0.452. LINGER scores 0.536/0.430 [C].
- **On K562, the published methods show almost no TF-specific signal.** Permuting TF identities in the
  released K562 networks leaves BEAR's metrics essentially unchanged (LINGER ChIP AUROC 0.536 → 0.537,
  AUPRC 0.430 → 0.432) [C].
- **LINGER's lead on the mouse sets comes from how it spreads score across TFs, not from picking the
  right targets.** On E7.5 rep1, LINGER's ChIP AUROC of 0.653 falls to 0.509 when TF identities are
  permuted, and to 0.541 when its scores are rank-normalised within each TF [C].
  - We also ran a "TF-swap" test: does TF *t*'s own score row predict *t*'s targets better than other
    TFs' rows do? LINGER's gap is only +0.005. CellOracle's is +0.033, GRaNIE's +0.024, and a label-free
    motif x accessibility score's +0.004 [C].
  - As a diagnostic only, ranking edges by each TF's ground-truth out-degree reaches ChIP AUROC 0.87-0.92
    [C]. BEAR's AUROC therefore mostly rewards knowing *which TFs have many ChIP peaks*.
  - Simple label-free TF-level proxies (TF expression, motif hit rate) do not reproduce that signal
    (≤0.52 on mouse) [C]. LINGER's TF-level signal is real, but it is not TF-target specificity.
  - **The comparison at equal coverage.** If LINGER's E7.5 network is filled out to every measured pair
    (missing pairs ranked last), its ChIP AUROC is **0.531**, against 0.532 for motif x accessibility.
    Its AUPRC is 0.313, against 0.296 for motif x accessibility and 0.274 for the constant-score floor
    [C]. On K562, coverage-completed LINGER gives 0.560/0.452, and target mean expression alone gives
    0.554/0.452.
  - BEAR's AUROC rewards emitting a selective subset. Its AUPRC rewards emitting everything. A single
    output cannot win both cleanly.

**Recommendation.**
1. **Replace MeVD-GRN for the in-scope track, and keep it as the supervised reference track.** Its GNN
   smoothing and gene-level inputs are exactly the hub and coverage channels described above.
2. **Build a small, label-free, cell-level model: accessibility-gated sparse regression (AGR).**
   - The core is a per-dataset metacell regression: target expression ~ Σ_t w_gt · TF_t expression ·
     gate_cgt.
   - The gate gate_cgt is each cell's accessibility at peaks near gene *g* that carry TF *t*'s motif.
   - Edge scores come from conditional, held-out predictive gain plus stability selection.
   - The score is factorised into a per-TF "activity" term and a within-TF target term, and each is
     reported and ablated separately. This keeps any gain from between-TF allocation visible rather than
     hidden.
   - It is TF-specific by construction, cheap (GPU-hours, not weeks), and in scope.
   - What is new is not the regression idea (Pando, LINGER and DIRECT-NET are close relatives). It is
     the **evaluation contract**: coverage-matched comparisons, label-free marginal baselines, the TF-swap
     test, TF-normalised AUROC, and BEAR's own shuffles, stability and runtime checks.
3. **Make the paper "what BEAR-GRN rewards, plus a method that wins only where it should".** Do not make
   the paper "we beat LINGER". That second claim is available to a constant score.
4. **Do not** pretrain a large foundation model, read GRNs from attention, build diffusion/flow models,
   or tune anything against BEAR labels. These are explained in s3.6.

**Odds (our calibrated judgement, [I]).**

| Outcome by 30 Nov | Odds |
|---|---|
| Beats every published method on BEAR ChIP AUPRC on ≥7/9 datasets (BEAR's own scoring) | ~85%, and mostly a coverage effect |
| Beats every published method on ChIP AUROC on ≥5/9 (LINGER's mouse 0.64-0.66 is mostly TF-level signal) | ~25% |
| Beats all label-free baselines (coverage, target expression, motif x accessibility) on ChIP AUROC **and** AUPRC on ≥6/9 | ~35% |
| TF-swap gap > 0 (bootstrap CI) **and** larger than motif x accessibility on ≥6/9 ChIP sets | ~30% |
| Above-random AUPRC on K562 knockout or intersection | ~15% |

| Venue outcome | Odds |
|---|---|
| ICML 2027 acceptance, if submitted | ~5% |
| PLOS ONE eventual acceptance of an audit-plus-method paper | ~65% (decision likely after Feb 2027) |

**Gates.**
- **G1, 23 Oct:** the model must recover planted edges on simulated multiome data.
- **G2, 31 Oct, K562 development set:** the model must beat the label-free baselines and show a
  positive swap gap. If it fails, we switch to the PLOS ONE audit paper, with AGR included as an honest
  baseline.
- **G3, 23 Nov:** ICML go or no-go.

The compute need is at most about 330 GPU-hours (about 2 GPU-weeks, an upper-end guess; s7) of the
3-4 GPU-weeks available. About 1-2 GPU-weeks stay in reserve.

**Decisions we need from you.**
1. **Accept the success criterion.** Success means beating the label-free baselines and passing the
   TF-swap test on frozen datasets, with BEAR's own metrics reported alongside. It does not mean
   "beat LINGER".
2. **Choose the venue strategy.**
   - ICML's dual-submission rule (2026 text) forbids "substantially similar" work submitted in parallel
     to a journal.
   - So ICML in January and PLOS ONE in February must be *different papers*: method vs audit. The
     alternative is to choose one venue.
   - We recommend PLOS ONE as the primary venue, with ICML only if G3 passes.
3. **Set the external-knowledge policy and author contact.**
   - We propose that only JASPAR/CIS-BP motifs, a published TF census, and the genome are allowed.
     ENCODE bulk pretraining would be optional and declared. No curated TF-target database is used.
   - We propose emailing the BEAR authors this week with the questions in s1.6.

---

## 1. Verified BEAR-GRN requirements

Sources:
- `~/.cache/bear/paper/` (SI = `MOESM1_ESM.txt`; peer-review file = `MOESM11_ESM.txt`; reporting summary
  = `MOESM2_ESM.txt`; Supp Data dump = `all_xlsx_dump.txt`).
- Code snapshot `code/UzunLab-BEAR-GRN-033a2cd/` (cited as `R/<file>:<line>`).
- Zenodo listings in `paper/zenodo/`.

The main text is still not online (accelerated preview; `article_fulltext.txt` holds only the abstract).
Anything that lives only in the main text is therefore unverified.

### 1.1 What a new method must supply
- **One edge list per dataset: TF, target, and a numeric score**, in any column names, which are passed
  as `tf_column`, `target_column` and `score_column` [V README l.336-343]. Columns must not be quoted
  [V README l.174].
  - A method with no score gets early metrics only. DIRECT-NET is excluded from AUROC/AUPRC for that
    reason: it "lacks confidence scores or rankable edge metrics" [V MOESM11 l.273-276].
- **Inputs are BEAR's post-QC matrices** (`INPUT.DATA`: gene x cell RNA counts and peak x cell ATAC).
  The README says these exist "to train your new GRN inference method" and "ensure fair comparison"
  [V README l.304-318]. We must not redo the QC.
- **Stability** requires running the method on `INPUT.DATA.STABILITY`.
  - Released stability networks are 10 per method per dataset, named `70_percent_subsampled_<k>`
    [V `zenodo/listing_STABILITY_GRNS.zip.txt`].
  - They cover K562, Macrophage S1 and S2, and the four embryo sets, not iPS or naive mESC [V same
    listing].
  - Scoring: top 10% of edges by score vs a random 10%, pairwise Jaccard, median
    [V `R/reproduce_stability.R` l.13-14, 368].
- **Runtime and memory** are measured on 1k-5k-cell subsamples of E8.5 rep1 and rep2
  [V SI Figs 18-19 captions]. Subsampling uses `set.seed(123)` and `sizes <- c(1000, ..., 5000)`
  [V `inst/scripts/DATA.PREPROCESSING/Step4.subsample_cells.R`]. Hardware is not stated.
  - Read from the figure (approximate): at 5k cells LINGER takes about 30k s (about 8.3 h) elapsed with
    about 85 GB peak RSS. SCENIC+ takes about 25k s with about 91 GB, so SCENIC+, not LINGER, has the
    highest memory at 5k [V by eye, SI Fig 18a-b].
- **Modality shuffles** (Supp Data 7) cover four conditions (Original, Cell_Shuffle, Gene_Shuffle,
  ATAC_Gene_Shuffle) on 8 datasets (no naive mESC) [V tidy CSV]. A mixed-multiome analysis (RNA from one
  cell type, ATAC from another) is in Supp Figs 15-16 [V].

### 1.2 Scoring, exactly (from `R/reproduce_ROC_PR_full.R`)

| Step | What the code does | Lines |
|---|---|---|
| Standardise | `abs(score)`, upper-case names, `max` over duplicate (TF, target) | 278-286 |
| Tested space | GT TFs x GT targets; edges outside are dropped | 203-204, 292-294 |
| AUROC | on the method's **own** edges in the tested space; larger class down-sampled 1:1 with `set.seed(42 + i)` (i = the method's row in `method_info`, so a new method's draw depends on where it is inserted); `pROC::roc(direction = "<")` | 243, 308-331 |
| AUPRC | complete universe GT-TF x GT-target; **unpredicted = 0**; all positives + `min(10 x pos, all)` negatives; `PRROC::pr.curve(...)$auc.integral` | 215-223, 348-369 |
| Random | `runif` over the same universe, same sampling | 386-397 |
| Early metrics | top 10,000 edges after filtering to GT TFs/targets; precision, recall, F1; ChIP only in the paper (Supp Data 6) | `R/reproduce_early_metrics.R` |

Our Python port (`src/benchmarks/bear_metrics.py`) matches AUPRC exactly. AUROC differs only through the
RNG (`docs/experiments/bear_grn_benchmark.md` s4, s13.1). We re-verified this today: the fast
re-implementation in `analysis/labelfree_baselines.py` reproduces the coverage baseline's K562 ChIP
AUPRC of 0.431 [C].

**New port bug [C]:** `GROUND.TRUTHS.KO/filtered_RN111_Naive_mESC.tsv` (and `GROUND.TRUTHS/filtered_RN112_Naive_mESC.tsv`)
have 62,992 target names with a trailing space ("AATF "). R's `readr` trims whitespace by default, but
our `bm.load_gt` does not. As a result, naive-mESC KO is currently unscorable in the port (it shows as
`nan` in `results/summary_tables.md`). Fix: add `.str.strip()` in `bm.load_gt` and in
`bear_data.read_gt_pairs`. This is a one-line change; I did not make it, because I may not edit repo
files. No other GT file has trailing whitespace (checked with `grep -c -P "[ \t]$"` over all GT files).

### 1.3 Our beliefs vs the sources

| Belief | Verdict | Evidence |
|---|---|---|
| 9 datasets; K562 411 cells, Mac 474/528, iPSC 5,710, embryo 5,064/1,447/6,557/6,252, naive mESC 3,568 (CCA-paired) | **Correct** | Supp Data 1 [V] |
| Genes after QC | 10,050 / 9,171 (Mac), 10,012 (K562), **4,808 (iPSC)**, 7,712-9,882 (mouse), 8,049 (naive) | Supp Data 1 [V] |
| 9 methods; DIRECT-NET has no scores; scGLUE and STREAM excluded | **Correct** | MOESM11 l.273-276; reporting summary "Two methods (STREAM and scGLUE) were excluded from accuracy evaluation due to technical incompatibilities" [V] |
| "LINGER ranks first overall" | **Partly.** The abstract says "LINGER and DIRECT-NET show the highest overall accuracy and stability". The Fig 6 ranking is in the unpublished main text | abstract [V] |
| Ground-truth types per dataset | **Correct** (ChIP everywhere; KO/Union/Intersection on K562 and all 5 mouse sets; Core and cell-type-exclusive on the 4 human sets) | GT dirs [V] |
| ChIP and KO edges overlap about 1.3% | **Correct**: 1.3% K562, 1.4% mouse embryo; target overlap about 49% | SI Fig 1 caption [V] |
| "GT filtered to expressed genes" | **Stated, but not to the post-QC gene list**: K562 has 27,747 GT targets vs 10,012 measured genes. Most universe pairs can never be scored by anyone | Supp Data 1, summary.json [V] |
| Core GT | One file shared by all 4 human datasets (identical size, 743,318 bytes) | `ls` [V] |
| Max ChIP AUROC about 0.66 | **Correct** (0.658, LINGER E8.5 rep1). On other ground truths: 0.747 (LINGER, naive Intersection, 976 edges), 0.650 (GRaNIE Mac S1 ChIP, **100 edges**) | Supp Data 2-5 [V] |
| K562 KO and Intersection at or below random AUPRC for all methods | **Correct** | Supp Data 3, 5 [V] |
| Shuffling changes AUPRC by about 0.02-0.03 for most methods | **Wrong in magnitude.** For 6 of 8 scored methods the median AUPRC change is at most 0.005 under all three shuffles. GRaNIE is the exception under the ATAC-gene shuffle (-0.029). LINGER drops 0.04-0.09 AUPRC and **0.15-0.20 AUROC** on mouse. But LINGER's evaluable edges also collapse under shuffles (E7.5 rep1: 496,595 → 68,400 under cell shuffle), so part of that drop is a coverage effect | tidy CSV [C, Appendix A.4]; Supp Data 7 [V] |
| LINGER most stable | **Supported for cross-replicate edge stability** (LINGER and TRIPOD, Jaccard 0.65-0.87). The subsampling-stability heatmap is in the main text and unverified | MOESM11 l.437-440 [V] |
| LINGER about 8.5 h / 89 GB at 5k cells | **Roughly**: about 8.3 h, about 85 GB. SCENIC+ has higher memory | SI Fig 18 [V by eye] |
| Trivial baselines | Coverage beats or ties every method's ChIP AUPRC on 8/9 datasets; target mean expression beats it on 9/9 | Appendix A [C] |

### 1.4 Scope and eligibility: exact text and what it permits

Verbatim from the authors' reply to Reviewer 2, Comment 1 (`MOESM11_ESM.txt` l.330-346):

> "Circularity: BEAR-GRN evaluates methods against ChIP-seq-derived ground truth networks. Merlin+P and MTLRank require ChIP-seq-derived TF-gene labels as part of their training input, evaluating them against the same data class used to train them would inflate performance estimates and invalidate comparisons with the unsupervised methods benchmarked here."
>
> "Different tasks and data requirements: All methods in BEAR-GRN infer GRNs de novo from paired scRNA-seq and scATAC-seq data. scGREAT and DeepIMAGER rely solely on scRNA-seq, lacking the chromatin accessibility signals central to BEAR-GRN's multi-omics framework. [...]"
>
> "Scope: BEAR-GRN targets settings where cell-type-specific ChIP-seq labels are unavailable. Methods that require such labels for training are inapplicable in precisely the scenarios this benchmark is designed to address."
>
> "Fairly evaluating these supervised approaches would require a dedicated design (e.g., leave-one-TF-out cross-validation using multi-cell-type ChIP-seq compendia), which we leave for future work."

Related precedent from the same file: LINGER "incorporate[s] large-scale bulk regulatory data during
neural network training", and "this pre-training on bulk data effectively gives LINGER a regulatory
prior" (l.237-243). LINGER uses "HOMER's curated motif database" (l.231). CellOracle uses "prior network
structure" (l.242). The README lists the GT sources as "ChIP-seq experiments, BEELINE, Literature-curated
databases (ChIP-Atlas, ESCAPE, ENCODE, CistromeDB), Experimental knockout validation studies (KnockDB,
ESCAPE)" (README l.709-713).

| Input or practice | In scope? | Basis |
|---|---|---|
| Self-supervised or unsupervised fitting to the dataset's own RNA + ATAC | **Yes** [I] | "infer GRNs de novo from paired scRNA-seq and scATAC-seq data". The word "self-supervised" never appears, but every benchmarked method fits a model to the data itself (regression, GBDT, NN) |
| Motif databases (JASPAR, HOMER, CIS-BP, HOCOMOCO, cisTarget) | **Yes** [V by precedent] | Used by every scored method (SI Table 1). Caveat [I]: many PWMs were themselves derived from ChIP-seq peaks. Declare this, and run a sensitivity check with a non-ChIP motif subset (CIS-BP PBM/SELEX) |
| External **unlabeled** pretraining (bulk ENCODE RNA/ATAC, atlases, Geneformer) | **Yes, if declared** [I, LINGER precedent] | l.237-243. It must not include TF-target labels |
| Curated TF-target databases (DoRothEA, CollecTRI, TRRUST, RegNetwork) | **No** [I + our zero-leakage rule] | Not addressed by the authors. These resources aggregate ChIP-seq and literature evidence that overlaps BEAR's GT sources (ChIP-Atlas, ENCODE, CistromeDB). CellOracle's "prior network structure" is motif-derived in BEAR's runs (SI Table 1), so it is not a precedent |
| Non-cell-type-specific ChIP compendia (ReMap, SC-MO-GRN-DB RN005/RN006), even with BEAR TFs removed | **No for the in-scope track** [I] | The authors describe exactly this as the supervised "dedicated design ... for future work". This is our MeVD-GRN supervised reference track |
| Choosing the TF list from GT sources | **Avoid** [I] | `bear_data.py` l.128-143 defines TF nodes as GT sources. For the label-free track, take the TF list from a census (Lambert et al. 2018 human TFs; AnimalTFDB mouse). This does not change scores, because scoring filters to GT TFs, but it removes any GT dependence from the model |
| Any use of BEAR GT for tuning, selection, early stopping or thresholds | **No** | Our rule. Also implied by "de novo" |

### 1.5 Consequences for design [I]
1. **AUPRC is dominated by coverage.** Emit **all** measured TF-census x gene pairs (dense output).
   Otherwise we lose AUPRC to coverage and win or lose for reasons unrelated to biology. Always compare
   against the coverage floor.
2. **AUROC is computed on our own edges.** The dense output fixes the evaluated set, which makes our
   AUROC comparable across our variants. Published methods are each scored on different edge sets, so we
   also report their **coverage-completed** versions (s4.2).
3. **AUROC is dominated by between-TF allocation.** We must report raw, TF-normalised, and TF-permuted
   AUROC side by side.

### 1.6 Questions to send to the BEAR-GRN authors (this week)
1. Is a self-supervised method that uses only `INPUT.DATA`, motif PWMs and the genome in scope? Is
   declared **unlabeled** external pretraining (as in LINGER) acceptable? Do you consider curated
   TF-target databases (DoRothEA/CollecTRI) out of scope?
2. The GT targets (e.g. K562 27,747) far exceed the post-QC genes (10,012). Which expression filter
   produced the "filtered for expressed genes" GT? Do you intend unmeasured targets to stay in the AUPRC
   universe?
3. Should a new method output a dense or a thresholded network? Given that AUROC is computed on each
   method's own edges and AUPRC scores unpredicted pairs as 0, is there a recommended output size?
4. `benchmark_new_method*.R` map Macrophage_S2 to the Buffer1 GT, while `reproduce_ROC_PR_full.R` uses
   Buffer2. Which is canonical?
5. Will you accept numbers from a faithful Python port, or must new methods use the R functions? (The
   AUROC draw depends on `set.seed(42 + i)` and on the method's position.)
6. What are the stability subsample parameters (70%, 10 subsamples, seeds)? What are the 7 files in
   `STABILITY_GRNS/K562/new_method_stability/`? Are iPS and naive mESC excluded from stability by
   design?
7. On what hardware were runtime and memory measured? Was a GPU allowed or used (LINGER can use one)?
8. Are the shuffled inputs (cell, gene, ATAC-gene) released? If not, which seeds should we use to
   regenerate them?
9. What is the provenance of the iPS "RN000" ground truth? It does not appear in SC-MO-GRN-DB's list.
10. The naive-mESC KO GT has trailing whitespace in target names. Was it trimmed (readr default) in your
    scoring? We assume yes.
11. Do you plan a public leaderboard, and would you list a self-supervised method together with a
    trivial-baseline panel?
12. The ChIP GTs are very dense (K562: 1.38 M edges over 150 TFs, about 9,200 targets per TF). Which
    peak-to-gene assignment and thresholds produced them? This determines what "TF-specific" can mean.

## 2. Literature position (2023-2026)

**Method.** One helper search on 10 Oct 2026 covered PubMed metadata, web search, arXiv and OpenReview,
about 30 queries. DOIs below were confirmed against PubMed unless marked. "[mem]" means the one-line
summary comes from prior knowledge of the paper and was not re-read today. This is a targeted survey, not
a systematic review. Expand it in December for the related-work section.

### 2.1 Closest prior art to AGR (label-free, multiome)

| Method | Reference | What it does | Relation to AGR |
|---|---|---|---|
| **Pando** | Fleck et al., *Nature* 2023, [10.1038/s41586-022-05279-8](https://doi.org/10.1038/s41586-022-05279-8) | Per-gene regression of target expression on TF expression x accessibility of motif-bearing peaks [mem] | **Closest. AGR's regression form is essentially Pando's.** We must credit it, and the ablation "AGR → Pando-like" is mandatory |
| **LINGER** | Yuan & Duren, *Nat Biotechnol* 2024/25, [10.1038/s41587-024-02182-7](https://doi.org/10.1038/s41587-024-02182-7) | Per-gene NN from TF expression and RE accessibility, bulk-ENCODE lifelong-learning prior, motif manifold regularisation, Shapley scores [mem; also MOESM11 l.212-243] | Benchmark leader; Candidate D mirrors its bulk prior |
| **TRIPOD** | Jiang et al., *Cell Systems* 2022, [10.1016/j.cels.2022.08.004](https://doi.org/10.1016/j.cels.2022.08.004) | Nonparametric TF-peak-gene trio tests on matched metacells | Conditional-dependence idea, metacells |
| **CellOracle** | Kamimoto et al., *Nature* 2023, [10.1038/s41586-022-05688-9](https://doi.org/10.1038/s41586-022-05688-9) | Motif base GRN, then per-cluster regularised regression, plus in-silico knockout | Motif restriction and in-silico perturbation, but cluster-level, not cell-gated |
| **DIRECT-NET** | Zhang et al., *Sci Adv* 2022, [10.1126/sciadv.abl7393](https://doi.org/10.1126/sciadv.abl7393) | XGBoost on kNN pseudo-cells, CRE-gene links | Pseudo-cells; no scores in BEAR |
| **SCENIC+** | Bravo González-Blas et al., *Nat Methods* 2023, [10.1038/s41592-023-01938-4](https://doi.org/10.1038/s41592-023-01938-4) | cisTopic, cisTarget, then GRNBoost2 | Baseline |
| **GRaNIE / GRaNPA** | Kamal et al., *Mol Syst Biol* 2023, [10.15252/msb.202311627](https://doi.org/10.15252/msb.202311627) | Correlation-based eGRN. GRaNPA evaluates GRNs by predictive power for differential expression | Predictive gain used for *evaluation*. AGR uses it for *inference* |
| **FigR** | Kartha et al., *Cell Genomics* 2022, [10.1016/j.xgen.2022.100166](https://doi.org/10.1016/j.xgen.2022.100166) | DORCs plus motif enrichment x TF-RNA correlation | Baseline |
| **Dictys** | Wang et al., *Nat Methods* 2023, [10.1038/s41592-023-01971-3](https://doi.org/10.1038/s41592-023-01971-3) | ATAC footprinting plus a stochastic expression model; dynamic GRNs [mem] | Not in BEAR |
| **scGLUE** | Cao & Gao, *Nat Biotechnol* 2022, [10.1038/s41587-022-01284-4](https://doi.org/10.1038/s41587-022-01284-4) | Graph-linked VAE; regulatory scores from embeddings | Excluded by BEAR ("technical incompatibilities") |
| **scTFBridge** | Wang et al., *Nat Commun* 2025, [10.1038/s41467-025-64227-y](https://doi.org/10.1038/s41467-025-64227-y) | Disentangled hybrid VAE, motif-constrained ATAC decoder, explainability-based edges | Label-free with a motif prior, so **in scope**. It would be a strong extra comparator if it runs on BEAR inputs (week 2-3 option) |
| **scMultiomeGRN** | Xu et al., *Nucleic Acids Res* 2025 (DOI not confirmed) | Graph model with cross-modal attention; cell-specific GRNs | Supervision status not confirmed. Check before citing it as label-free |
| **EpiAwareNet** | Xu et al., arXiv 2606.00685 (Jun 2026), venue not confirmed | Transformer with gene-peak cross-attention; bulk GRN prior as noisy positives | **Weakly supervised by a GRN prior, so out of scope** for the label-free track |
| **GET** | Fu et al., *Nature* 2025, [10.1038/s41586-024-08391-z](https://doi.org/10.1038/s41586-024-08391-z) | Accessibility + sequence foundation model predicting expression; regulation via attribution [mem] | Possible external pretraining source. Heavy |
| seq2PRINT / PRINT | Hu et al., *Nature* 2025, [10.1038/s41586-024-08443-4](https://doi.org/10.1038/s41586-024-08443-4) | Multiscale ATAC footprints | Could sharpen motif-in-peak evidence later. Not a GRN method |

### 2.2 Foundation models and GRN evaluation
- Ahlmann-Eltze, Huber & Anders, *Nat Methods* 2025,
  [10.1038/s41592-025-02772-6](https://doi.org/10.1038/s41592-025-02772-6): deep perturbation-effect
  predictors "do not yet outperform simple linear baselines".
- Kedzierska et al., *Genome Biol* 2025,
  [10.1186/s13059-025-03574-x](https://doi.org/10.1186/s13059-025-03574-x): zero-shot limitations of
  single-cell foundation models.
- A 2026 arXiv preprint (Kendiukhov, [2602.17532](https://arxiv.org/abs/2602.17532)) reports that
  scGPT/Geneformer attention edges add no value over gene-level baselines for perturbation prediction.
  This was reported by the helper and **not read by us**.
- geneRNIB (bioRxiv [10.1101/2025.02.25.640181](https://www.biorxiv.org/content/10.1101/2025.02.25.640181v1)):
  simple models often beat complex pipelines, and expression-only methods beat multimodal ones.
- **Implication:** this supports rejecting Candidate B (attention-derived networks) for this cycle.

### 2.3 Benchmarks and trivial-baseline precedent
- **BEAR-GRN** (Nat Commun 2026). No citing papers were found as of 10 Oct 2026 (helper search,
  consistent with the brief).
- **GRETA** (Badia-i-Mompel et al., bioRxiv 2024,
  [10.1101/2024.12.20.629764](https://doi.org/10.1101/2024.12.20.629764)): a multimodal GRN benchmark
  reporting modest performance and heterogeneous networks. **This is the most important neighbour for the
  audit paper**; read it in full in week 1.
- BEELINE (Pratapa et al., *Nat Methods* 2020, [10.1038/s41592-019-0690-6](https://doi.org/10.1038/s41592-019-0690-6)).
- PEREGGRN (Kernfeld et al., *Genome Biol* 2025, [10.1186/s13059-025-03840-y](https://doi.org/10.1186/s13059-025-03840-y)).
- CausalBench (Chevalley et al., arXiv [2210.17283](https://arxiv.org/abs/2210.17283); journal version
  not confirmed).
- **Degree and hub precedent:** Gillis & Pavlidis, *PLoS ONE* 2011,
  [10.1371/journal.pone.0017258](https://doi.org/10.1371/journal.pone.0017258): node degree alone gives
  strong gene-function prediction, and they propose computational controls. Also *PLoS Comput Biol* 2012,
  [10.1371/journal.pcbi.1002444](https://doi.org/10.1371/journal.pcbi.1002444). It is fitting for a PLOS ONE
  audit to cite this lineage.
- InfoSEM (ICML 2025, [PMLR v267 cui25e](https://proceedings.mlr.press/v267/cui25e.html)) reportedly
  notes that degree-exploiting baselines do well on supervised GRN benchmarks. This is from a search
  summary and must be verified in the paper.

### 2.4 ML-venue precedent
- ICLR 2024 (oral): Dai et al., "Gene Regulatory Network Inference in the Presence of Dropouts: a Causal
  View" (scRNA only) [helper-confirmed OpenReview].
- ICML 2025: InfoSEM (scRNA only, label-free text-embedding priors) [PMLR].
- **No multiome (RNA + ATAC) GRN paper was confirmed in an ICML, NeurIPS or ICLR main track** (helper
  search; "not found" is not proof of absence).
- Implication: ICML is not hostile to GRN inference, but accepted papers brought a clear ML contribution,
  such as causal identifiability under dropouts, or informative priors in a generative model.

### 2.5 What would be novel [I]
- **For PLOS ONE (soundness-first):** the first decomposition of a multiome GRN benchmark into coverage,
  target-marginal, TF-marginal and TF-specific components. We also show that published methods' BEAR
  scores are matched by label-free gene-level baselines at equal coverage, and introduce coverage
  completion, TF-normalised AUROC and the TF-swap test as standard controls.
  - No paper doing this for multiome GRN benchmarks was found, though the hub-control *idea* goes back to
    Gillis & Pavlidis 2011. This is novel enough for PLOS ONE, and PLOS does not require novelty anyway.
- **For ICML:** novelty would have to come from AGR. Specifically:
  - (i) cell-level, motif-specific gating as a source of TF identifiability, with a formal statement and
    synthetic evidence;
  - (ii) conditional held-out gain with stability selection as an edge score;
  - (iii) the controls framework.
  - Given Pando and LINGER, (i)-(ii) are incremental unless the results are clearly better on the
    TF-specific channel. That is why ICML is conditional on G3.

## 3. Candidate approaches and recommendation

### 3.0 What the diagnostics say a method must do [C → I]
Three facts drive the choice of method. All rows are on ChIP; the evidence is in Appendix A.

| Signal channel | How much BEAR rewards it | Can a label-free method get it? |
|---|---|---|
| Coverage (emit every measured pair) | AUPRC floor 0.27 on mouse, against a random baseline of 0.155; at or above every published method on 8/9 | Yes, trivially |
| Target marginal (hub genes) | ChIP AUROC 0.52-0.55 label-free; 0.59-0.65 with labels (LOO in-degree oracle) | Partly (expression and accessibility) |
| TF marginal (which TFs are "big") | AUROC up to **0.87-0.93** for the GT oracle. LINGER's mouse AUROC is mostly this (0.653 → 0.541 TF-normalised) | Simple label-free proxies fail (≤0.52). LINGER gets it from somewhere |
| TF x target interaction (true TF-specific edges) | Small: swap gaps of 0.001-0.033 for every method incl. motif x accessibility (0.001-0.008) | This is the scientifically interesting channel. Its effect on headline metrics will be small |

So the plan has to be honest about two things. First, a genuinely TF-specific method may move BEAR's
headline numbers only modestly. Second, the TF-marginal channel can move them a lot. Any method we build
must therefore (a) be dense, (b) separate its per-TF activity from its within-TF target ranking, and
(c) be judged on both.

### 3.1 Candidate A: accessibility-gated sparse regression with a perturbation scoring head (AGR). **Recommended.**

**Objective (self-supervised, per dataset).**
- **Metacells.** Build them with kNN in a joint RNA PCA + ATAC LSI space, with k about 20-50 and capped
  overlap. These are the LINGER and DIRECT-NET defaults (LINGER k=20; DIRECT-NET k=50, overlap 0.5;
  MOESM11 l.221-224).
- **Model.** Predict held-out metacell expression of gene *g*:

  ŷ_cg = b_g + Σ_{t ∈ C_g} w_gt · x_ct · a_cgt

  - x_ct is TF *t*'s expression in metacell *c*.
  - a_cgt = Σ_p RP(p,g) · motif_t(p) · acc_cp is the **cell-level accessibility of peaks near *g* that
    contain TF *t*'s motif** (RP = the exponential TSS-distance weights already in
    `preprocess_scatac`).
- **Candidate set C_g.** Defined by label-free evidence tiers:
  - (1) a motif in an accessible peak within ±100 kb;
  - (2) plus peak-promoter co-accessibility across metacells;
  - (3) expression-only TFs, which enter with a stronger penalty.
- **Loss.** Masked-metacell MSE (5-fold over metacells), plus an L1 / group penalty on w, plus a
  low-rank-plus-sparse factorisation w_gt = u_tᵀv_g + s_gt. The factorisation shares TF activity across
  genes and gives an explicit TF-activity term A_t.
- **Edge score.**
  - The core score is **conditional held-out predictive gain** (in-silico knockout): the increase in
    held-out loss when TF *t*'s term is removed from gene *g*'s model.
  - It is multiplied by its **stability-selection frequency** over about 30 metacell bootstraps.
  - It is reported in factorised form: score_gt = A_t x r_gt, where r_gt is the within-TF percentile.

**Why it should produce TF-specific signal [I].**
- The gate a_cgt differs between TFs only through *their own* motifs. Two co-expressed TFs with
  different motifs get different regressors.
- The conditional gain partials out shared hub variation.
- b_g absorbs per-gene level, so target popularity does not inflate w.
- The cell-level (metacell-level) pairing is exactly what BEAR's cell shuffle destroys, so the shuffle
  controls become a direct test.
- *This is a hypothesis.* A pure hub model can still leak in through high-variance genes. The swap test
  (s4.4) is the guard.

**Reuse from our codebase.**
- `bear_data.py`: INPUT.DATA readers and the RP-weighted peak-gene incidence.
- `bear_motif.py`: 2-bit genome, JASPAR loader, torch conv1d scanner. It needs extending from "top-5
  peaks per gene" to all peaks in the window, and to CIS-BP/HOCOMOCO motifs for mouse. Today only 49-57
  of about 230 mouse TF nodes have a JASPAR match [C].
- `bear_metrics.py` with the strip fix; `24_bear_score.py`; Slurm wrappers; leak guard.
- Not reused: GraphSAGE towers, the bilinear decoder, supervised samplers, and the Geneformer features
  (allowed only as an ablation).

**Compute [I/G].**
- One dataset is about 10k genes x about 30-150 candidate TFs per gene x 200-700 metacells. A batched
  GPU solve takes minutes, and 30 bootstraps take under 1 GPU-hour on a 2080 Ti.
- A full motif scan of about 180k peaks x about 800 PWMs takes about 1-3 GPU-hours per genome, run once.
- Total for 9 datasets x 5 seeds plus controls: about 3-5 GPU-days (s7).

**Main failure modes.**
1. TF collinearity: co-expressed paralogs share motifs, giving unstable w. Mitigations: group penalties
   over motif clusters, and stability selection.
2. Small datasets: K562 has 411 cells and the macrophage sets about 500, so only about 20-25 metacells at
   k=20. Fall back to kNN-smoothed single cells with heavier regularisation.
3. Motif coverage gaps for mouse TFs.
4. The very dense ChIP GTs (about 9,200 targets per TF in K562) bound how much TF-specific ranking can
   show.
5. The method may beat the baselines only through A_t, the TF-activity channel. We would report that
   honestly; it would not be a TF-specificity win.

**Honest odds [I].**

| Outcome | Odds |
|---|---|
| Working pipeline on all 9 datasets by 6 Nov | 80% |
| Beats all label-free baselines on dev K562 ChIP AUROC by ≥0.02 | 45% |
| Swap gap CI > 0 and above motif x accessibility on ≥6/9 frozen datasets | 30% |
| Matches LINGER's raw mouse AUROC (≥0.64) | 15% |
| Matches LINGER's **TF-normalised** mouse AUROC (0.541) | 50% |

**Novelty [I].** The regression form is *not* new. Pando fits gene expression ~ TF x motif-peak
interactions; LINGER, DIRECT-NET and CellOracle are regression-based (MOESM11 l.212-220; SI Table 1).
The defensible novelties are:
- the cell-level motif-specific gate combined with conditional, stability-selected predictive gain,
  with an explicit TF-activity / within-TF factorisation;
- the first evaluation of a multiome GRN method against coverage, marginal and TF-swap controls on
  BEAR.

For PLOS ONE that is enough. For ICML it is thin (s8).

### 3.2 Candidate B: masked multimodal pretraining with attention- or perturbation-derived networks. **Not now.**
- **Idea.** A transformer over gene (and peak) tokens is pretrained on public multiome atlases with
  masked expression and accessibility reconstruction. The GRN is read from TF→gene attention, or from
  in-silico TF knockdown.
- **In scope?** Yes: unlabeled pretraining, which must be declared.
- **Why not.**
  - At least 10k gene tokens plus about 10^5 peaks does not fit 11 GB GPUs without heavy approximation.
  - Two species are involved.
  - Atlas wrangling alone would take weeks.
  - Attention weights are not causal regulatory weights. Evidence that foundation-model attention or
    in-silico perturbation recovers GRNs well is weak (s2).
  - Most importantly, nothing in this design targets the TF-specific channel more directly than A does.
- **Odds** of beating label-free baselines *with* a positive swap gap by 30 Nov: about 10%. This is
  future work, after A shows signal.

### 3.3 Candidate C: in-silico perturbation scoring. **Adopt as A's scoring head, not as a separate model.**
- Removing a TF's term (or setting x_ct to 0) and measuring the held-out loss change *is* the edge score
  in A.
- The same head can be applied to any predictor, which is useful as an ablation: "marginal correlation
  vs conditional gain".

### 3.4 Candidate D: LINGER-style external bulk prior with a lighter model. **Optional add-on after G2.**
- **Idea.**
  - Fit A's TF-activity term and/or w on bulk paired RNA + ATAC (ENCODE, the resource LINGER uses for
    its "regulatory prior"; MOESM11 l.237-243).
  - Then regularise the per-dataset fit toward it, with L2 to the pretrained weights or EWC.
- **Why consider it.** Our diagnostics suggest LINGER's mouse AUROC is mostly TF-level allocation. A bulk
  prior is the most plausible label-free source of that.
- **In scope.** Yes, declared, matching LINGER.
- **Cost and risks.**
  - We do not know LINGER's exact bulk resource format. LINGER's repository is said to ship
    preprocessed bulk data [G, unverified]. If wrangling takes more than 3 days, drop it.
  - Risk: the bulk prior improves raw AUROC only via A_t. That is still worth knowing and reporting, but
    it is not TF-specificity.
- **Odds** that D lifts raw mouse AUROC by ≥0.03 over A: about 35% [G].

### 3.5 Candidate E: repurpose MeVD-GRN self-supervised (masked features or link prediction on the label-free graphs). **No.**
- Its objectives would learn symmetric co-expression and gene "popularity". GraphSAGE over kNN graphs
  propagates hubness, and the model has no cell-level input.
- On PBMC, an accessibility-only ranker already beat it (0.814 vs 0.734).
- Keep MeVD-GRN as the **supervised reference track** (M3 headline), clearly labelled as out of BEAR's
  scope.

### 3.6 Also rejected for this cycle
- **Diffusion or flow models; causal discovery (NOTEARS/DCDI-style).** Not identifiable from snapshot
  observational data at 10k nodes. Heavy, and no route to BEAR's TF-specific channel within 7 weeks.
- **Disentangled VAEs (scTFBridge-style) and prior-guided Transformers (EpiAwareNet-style).** Their
  regulatory readouts are latent-factor or attention based. Whether they are supervised by curated
  TF-target priors must be checked case by case (s2); if they use curated priors they are out of scope.
- **Contrastive objectives.** Contrasting matched vs mismatched RNA-ATAC cells learns cell pairing, not
  TF→target edges. It could serve as a pretraining auxiliary for A's metacell embedding, but it is not
  worth the time now.

### 3.7 Recommendation

**Replace MeVD-GRN's model for the in-scope track with Candidate A. Reuse MeVD-GRN's data, motif,
scoring and cluster infrastructure. Keep MeVD-GRN as the supervised reference.**

Reasons:
- A is the only candidate that targets the TF x target channel by construction.
- It runs in GPU-hours, so the budget goes to controls (seeds, stability, shuffles) rather than to one
  big model.
- It is in scope with no grey areas.
- If it shows no TF-specific signal, we still have a complete, honest PLOS ONE paper: the audit plus a
  well-controlled modern baseline.

### 3.8 What we should not do
1. Do not headline "beats all published methods": AUPRC wins come from coverage.
2. Do not emit a thresholded network whose size or threshold was chosen by looking at BEAR scores.
3. Do not define TFs from GT sources. Do not use curated TF-target databases or ChIP compendia anywhere
   in the in-scope model.
4. Do not tune on K562 beyond the pre-registered peek budget (s5). Do not look at the 8 frozen datasets
   until code freeze.
5. Do not let a TF-activity term silently drive the headline. Always report TF-normalised AUROC and the
   swap gap next to raw AUROC.
6. Do not use Geneformer or other gene-level embeddings in the scoring path except as a labelled
   ablation. They are a hub channel.
7. Do not spend more than about 20% of GPU time on pretraining (Candidate D) before G2 passes.
8. Do not chase K562 KO or Intersection AUPRC: every method, trivial or published, is at or below random
   there. Report it, but do not optimise for it.
9. Do not re-run the nine published methods. Their released networks, shuffle tables and stability
   networks are enough.

## 4. Evaluation and controls

Every claim in the paper must map to one row of the table in s4.9. All scoring uses the port, with the
strip fix and `n_rep = 20` draws for AUROC. For the camera-ready and leaderboard numbers we also run
BEAR's own R functions on a machine where tidyverse installs (Ada has no R today, so a conda R env is a
week-1 task).

### 4.1 BEAR's own metrics, run exactly as published
- AUROC, AUPRC and AUPRC_random on every ground truth for each dataset: ChIP, KO, Union, Intersection,
  Core, cell-type-exclusive.
- Top-10k precision, recall and F1 (ChIP).
- Stability: 10 x 70% subsamples on the 7 datasets with released stability networks; median Jaccard of
  the top 10% vs a random 10%.
- Modality shuffles: cell, gene, ATAC-gene.
- Runtime and peak memory at 1k-5k cells, E8.5 rep1 and rep2.

We submit two pre-registered outputs and score each as its own BEAR entry:
- **AGR-dense (primary).** The full TF-census x measured-gene list.
- **AGR-selected (secondary).** Only edges with stability-selection frequency ≥ 0.6. The threshold is
  fixed now and never tuned.

The reason for two entries is structural. BEAR computes AUROC on a method's own edges, which rewards
selective output (LINGER E7.5: 0.653 on its own 24% of pairs vs 0.531 coverage-completed). It computes
AUPRC over the universe, which rewards dense output. We report both entries rather than choosing after
seeing scores. We also report a LINGER-size truncation of the dense list (`--sparse_match` already
exists in `scripts/24`).

### 4.2 Coverage-matched comparisons (mandatory)
1. **Coverage floor.** Report the constant-score dense list next to every AUPRC.
2. **Coverage-completed published methods.** Append every measured TF x gene pair a method did not emit,
   with a score below its minimum, then re-score. This gives each published method the same coverage as
   ours, so AUPRC differences reflect ranking. (Computed for K562 and E7.5: Appendix A.5.)
3. **Size-matched truncation.** Score our list truncated to each published method's evaluable edge
   count, which matters for AUROC on own edges.

### 4.3 Baseline panel (fixed definitions; all 9 datasets; same dense output; script `analysis/labelfree_baselines.py`)

| Baseline | What it uses | Already computed |
|---|---|---|
| coverage | nothing (constant) | yes |
| target mean expression; target detection rate | RNA, gene-level | yes |
| target ATAC activity; target proximal RP; target peak count | ATAC, gene-level | yes |
| TF expression x target ATAC | rank-1, label-free | yes |
| motif x accessibility | JASPAR + peaks, TF-specific | yes (JASPAR only; extend to CIS-BP) |
| \|Pearson r\| (TF, target); GRNBoost2 | RNA only | yes, from the earlier `scripts/23` runs |
| TF-level probes (TF expression, motif hit rate) | label-free TF marginal | yes, 4 datasets |
| *Supervised reference (out of scope):* TF-disjoint in-degree; MeVD-GRN M3 | GT labels of other TFs | in-degree yes; M3 from the supervised track |
| *Diagnostics, never methods:* oracle TF out-degree, LOO in-degree, degree product | GT | yes |

### 4.4 TF-specificity battery (the core contribution)
For every method (ours, the released ones, and the baselines), on every dataset and ground truth:
1. **TF-identity permutation.** Permute TF labels among the method's TFs, 100 times. Report ΔAUROC and
   ΔAUPRC with a permutation p-value. This removes TF-marginal and TF-specific signal together.
2. **TF-normalised AUROC/AUPRC.** Convert scores to within-TF percentiles before scoring. This removes the
   TF marginal and keeps target hubness and TF-specificity.
3. **TF-swap gap.**
   - For each TF *t*, compute AUROC(own row, labels of *t*) − mean over *u* ≠ *t* of AUROC(row *u*,
     labels of *t*), across targets.
   - Report the mean over TFs, a 95% bootstrap CI over TFs, and the fraction of TFs whose own row beats
     the 95th percentile of other rows.
   - A TF-agnostic score gives exactly 0. Target hubs and TF marginals cannot create a gap.
   - This is the test from our PBMC lesson.
4. **Target-marginal projection.** Replace each edge's score by its target's mean score over TFs. If this
   matches the original, the method's information is gene-level.
5. **Within-target AUROC** (rank TFs for each gene). This is reported, but it is contaminated by TF
   marginals (oracle 0.91-0.93), so it is not a pass/fail criterion.

**Pre-registered TF-specificity criterion:** on a dataset x GT, the swap gap's 95% CI excludes 0 **and**
the point estimate exceeds that of motif x accessibility. The current reference values on ChIP are
+0.001 to +0.008 for motif x accessibility, and for the released methods +0.001 to +0.012 on K562 and
+0.005 to +0.033 on E7.5 [C].

### 4.5 Robustness, as in BEAR
- **Shuffles.** We regenerate cell, gene and ATAC-gene shuffles with fixed seeds, unless BEAR releases
  theirs (question 8).
  - Our output is dense, so coverage is constant under shuffles and any drop is a ranking effect. That is
    cleaner than LINGER, whose evaluable edges fall from 496,595 to 68,400.
  - Expected behaviour if AGR works: a cell shuffle should hurt more than for the published median
    (because the gate is cell-level), and a gene shuffle should push the swap gap to 0.
- **Mixed multiome** (RNA from K562, ATAC from macrophage, and the reverse): optional, as in Supp
  Figs 15-16.
- **Stability.** INPUT.DATA.STABILITY uses `.rds` files and needs `pyreadr`, or an R env to convert them.
- **Runtime and memory.** Report the same 1k-5k subsamples on (a) 40 CPU cores, no GPU, and (b) one
  2080 Ti. State the hardware; BEAR's is unknown.

### 4.6 Ablations (reported, never used for selection)
1. Remove the gate (x_ct only), i.e. RNA-only regression.
2. Remove the motif prior (all TFs are candidates).
3. Single cells instead of metacells.
4. Marginal correlation instead of conditional gain.
5. Remove A_t (within-TF ranks only).
6. Remove stability selection.
7. JASPAR only vs JASPAR + CIS-BP.
8. Add the bulk prior (Candidate D).
9. Add Geneformer gene embeddings (the hub-channel check).

### 4.7 Seeds and statistics
- 5 seeds (metacell construction, folds, bootstraps) on all 9 datasets. Report mean ± sd.
- Cross-dataset claims use paired one-sided Wilcoxon tests over datasets vs each baseline. Report the
  count of wins out of 9, with the K562 development set reported separately and excluded from the
  confirmatory count.
- Within a dataset, bootstrap CIs resample **TFs**, not edges, because edges of a TF are dependent.
- Holm correction across the pre-registered primary comparisons (s5.4).

### 4.8 Synthetic sanity check (allowed for tuning)
- Simulated paired RNA + ATAC with a planted TF→peak→gene network: our own lightweight simulator
  (metacell-level linear-Gaussian with motif-gated accessibility), or scMultiSim if it installs cleanly
  (its existence, version and licence are unverified [G]).
- Purpose: confirm that AGR recovers planted edges and that the swap test detects them. It is not
  evidence about biology.

### 4.9 Claims and the evidence each needs

| Claim we might make | Evidence required |
|---|---|
| "Beats published methods on BEAR" | Raw BEAR metrics **plus** coverage-completed published methods **plus** the coverage floor, on the 8 frozen datasets |
| "Signal beyond gene-level priors" | Beats every label-free baseline on AUROC **and** AUPRC on ≥6/8 frozen datasets |
| "TF-specific regulation" | Swap-gap criterion (s4.4) on ≥6/8, and target-marginal projection < original |
| "Uses chromatin accessibility at the cell level" | Gate ablation, plus a cell shuffle that hurts AUROC and the swap gap more than the gene-level baselines' (≈0) |
| "Stable / efficient" | BEAR stability JI and runtime vs released numbers |
| "BEAR's metrics reward coverage and marginals" | Appendix A tables, reproduced with BEAR's R code |

## 5. Experimental protocol and pre-registration

### 5.1 Development and test split
- **Development (GT-scored, limited peeks): K562 only.** Its ground truth (RN117) is not shared with any
  other dataset. We have also already inspected it heavily, so it cannot serve as a clean test set.
- **Engineering-only (inputs used, GT never scored before freeze): mESC E8.5 rep1.** This is the
  largest mouse set, used to size metacells and check runtime and memory. Using its inputs without
  labels leaks nothing.
- **Frozen test set (8 datasets):** Macrophage S1 and S2, iPS, E7.5 rep1 and rep2, E8.5 rep1 and rep2,
  naive mESC. Our method's outputs on these are scored **once**, after code freeze.
- **Why not a mouse development set:** all five mouse sets share the RN111 ChIP GT, and both macrophage
  sets share RN204. Tuning on any one of them would leak into its siblings.
- **Disclosure [I]:** our design is informed by benchmark-level diagnostics computed today on all 9
  datasets (label-free baselines and released methods; Appendix A). No output of our model has been
  scored. The paper must say this.

### 5.2 Hyperparameter policy
- Defaults come from the literature or BEAR's methods: metacell k=20 (LINGER), ±100 kb window, RP decay
  10 kb (already ours), and top-k motif peaks.
- They are tuned only on (a) the synthetic data, and (b) label-free criteria computed on the inputs of
  K562 and E8.5 rep1.
- At most 12 configurations reach the K562 GT peek. Each peek is logged with a timestamp, config hash
  and result in `docs/selfsup_plan/peek_log.md`.
- The **peek budget is 5 scorings of K562** before freeze. Peeks are for catching bugs, not for choosing
  between configurations unless the label-free rule ties.

### 5.3 Label-free model-selection rule (decided now)
Rank configurations by the sum of their ranks on three criteria, computed without any GT:
1. Held-out metacell predictive R² (masked-metacell CV).
2. Edge stability: median Jaccard of the top 10% across 10 x 70% subsamples (BEAR's own stability
   definition, computable label-free).
3. **Self-consistency with TF identity.** The model's held-out R² with true TF expression minus its R²
   when TF expression vectors are permuted among TFs. A model that does not care which TF is which scores
   0.

Ties are broken by fewer parameters. The selected configuration is frozen with its git hash, config
file and a SHA-256 of the output edge list.

### 5.4 Pre-registration checklist
Commit this to the repo before the first frozen-set scoring.
- [ ] Code commit hash, config file, and environment lock file.
- [ ] TF census list and motif database versions (JASPAR 2024 release; CIS-BP version) with file hashes.
- [ ] Any external data (ENCODE bulk, if Candidate D is used): accessions, preprocessing and hashes.
- [ ] Statement that no ChIP, KO or curated TF-target resource touches training, selection or
      thresholds.
- [ ] Output definition: AGR-dense (TF census x measured genes; score = A_t x r_gt; factorised parts
      also reported) and AGR-selected (stability frequency ≥ 0.6). Both thresholds are fixed now.
- [ ] Datasets: dev = K562; frozen = the 8 listed in s5.1; the K562 peek log is attached.
- [ ] Primary endpoints:
  - (P1) ChIP AUROC vs the best label-free baseline;
  - (P2) ChIP AUPRC vs the coverage-completed best published method;
  - (P3) ChIP swap-gap criterion.
  - Success = P1 and P3 on ≥6/8 frozen datasets.
  - BEAR raw metrics vs published methods are *secondary*.
- [ ] Secondary endpoints: Union, KO, Intersection, Core and cell-type-exclusive metrics; top-10k;
      stability; shuffles; runtime; ablations.
- [ ] Statistics: 5 seeds; Wilcoxon over datasets; TF-bootstrap CIs; Holm over P1-P3.
- [ ] Deviation log. Any change after freeze is reported as a deviation, with both versions scored.
- [ ] What we will report if it fails: everything, including the null.

## 6. Timeline and decision gates (10 Oct to 30 Nov experiments; writing to the deadlines)

| Week | Dates | Work | Gate / output |
|---|---|---|---|
| 1 | 10-16 Oct | Port strip fix. Baseline panel on all GTs (done for 9 datasets). Coverage-completed and swap tests for released methods on all 9 (download the remaining released GRNs, about 1.4 GB in total, to Ada `/share1`). Email BEAR authors. TF census lists. Full-window motif scan (JASPAR + CIS-BP), human and mouse. Metacell builder. Synthetic simulator. Write s5.4 v0 | Baseline/audit tables complete for all 9 |
| 2 | 17-23 Oct | AGR v0 (gate, sparse regression, conditional-gain head) on synthetic data, K562 inputs and E8.5 rep1 inputs. Label-free selection criteria implemented | **G1 (23 Oct): AGR on synthetic data recovers planted edges (AUPRC ≥2x motif x accessibility) and the swap test detects them.** If it fails, fix the model for at most 1 week, then go to fallback F1 |
| 3 | 24-31 Oct | AGR v1: stability selection, factorised A_t x r_gt, collinearity handling. Label-free selection over ≤12 configs. K562 peeks 1-3 | **G2 (31 Oct): on K562 dev, ChIP AUROC ≥ best label-free baseline + 0.02 *and* swap gap > motif x accessibility (CI > 0).** If it fails, F1: PLOS ONE audit paper, with AGR as a controlled baseline. Writing starts 1 Nov |
| 4 | 1-8 Nov | If G2 passed, Candidate D (bulk prior), timeboxed to 4 days; drop it if the data are not usable. Runtime and memory profiling. Freeze config | **Code freeze, 8 Nov.** Pre-registration committed |
| 5 | 9-15 Nov | Frozen runs: 9 datasets x 5 seeds, scored once. TF-specificity battery | Headline tables |
| 6 | 16-22 Nov | Stability (7 datasets x 10 subsamples x 1 seed); shuffles (9 x 3 x 3 seeds); ablations (9 x 9 x 1 seed); mixed multiome if time | Robustness tables |
| 7 | 23-30 Nov | **G3 (23 Nov): ICML go/no-go** (s8.1 criteria). Figures, reruns, buffer, R-port cross-check of headline numbers | **All results frozen 30 Nov** |
| Dec | 1-31 Dec | Write. ICML draft (if go) by 20 Dec, internal review 21-31 Dec. PLOS ONE draft in parallel; it is a different paper if ICML is a go (s8.3) | |
| Jan | 1 Jan to the ICML deadline | ICML polish, appendix, reproducibility package. The 2027 deadline is **not yet announced**; ICML 2026's was 28 Jan AoE, and a tracker estimates about 22 Jan 2027 [V/G]. Plan for 15 Jan | ICML submission |
| Feb | to 28 Feb | PLOS ONE submission: data and code on Zenodo/GitHub, DAS, cover letter | PLOS ONE submission |

**Fallbacks.**
- **F1 (G1 or G2 fails):** PLOS ONE paper, *"What does BEAR-GRN reward? Coverage, marginal and
  TF-specificity controls for multiome GRN inference"*. It contains the full audit on all 9 datasets x 6
  GTs, coverage-completed rescoring of all released methods, the swap battery, and AGR plus motif x
  accessibility as controlled label-free baselines. No ICML submission.
- **F2 (G2 passes, G3 fails):** PLOS ONE method-plus-audit paper. ICML is skipped.
- **F3 (G3 passes):** ICML method paper, and a separate PLOS ONE audit paper with no AGR claims, subject
  to the dual-submission check in s8.3.

## 7. Compute budget

These are guesses [G] calibrated against our MeVD-GRN runs. AGR's per-run cost is the main
uncertainty, and we will re-estimate it after week 2. "Run" means one dataset x one seed, including 30
stability bootstraps. Ada limits: 4 x 2080 Ti, 40 CPUs, 125 GB RAM, 20 queued jobs, 4-day walltime.
The login node is never used for compute.

| Item | Runs | GPU-h each | GPU-h | CPU-h | Where |
|---|---|---|---|---|---|
| Released-GRN download (about 1.4 GB) and audit on 9 datasets (swap, coverage-completed, probes) | 9 | 0 | 0 | 15 | Ada CPU job (compute nodes have internet) |
| Full-window motif scans, human + mouse, JASPAR + CIS-BP | 2 | 4 | 8 | 4 | Ada GPU |
| Metacells and gate tensors (precompute) | 9 | 0 | 0 | 20 | Ada CPU |
| Synthetic development | ~20 | 0.5 | 10 | - | laptop 4060 / Ada |
| Label-free selection, ≤12 configs x (K562 + E8.5 r1 inputs) | 24 | 1 | 24 | - | Ada GPU |
| Candidate D bulk prior (only if G2 passes; 4-day timebox) | ~6 | 4 | 24 | 30 | Ada GPU/CPU |
| Frozen headline runs, 9 datasets x 5 seeds x 2 outputs (one fit) | 45 | 1 | 45 | - | Ada GPU |
| Stability: 7 datasets x 10 subsamples | 70 | 0.7 | 49 | - | Ada GPU |
| Modality shuffles: 9 x 3 types x 3 seeds | 81 | 1 | 81 | - | Ada GPU |
| Ablations: 9 x 9 variants x 1 seed | 81 | 1 | 81 | - | Ada GPU |
| Runtime and memory profile, 1k-5k cells x 2 sets x (CPU, GPU) | 20 | 0.5 | 10 | 20 | Ada (exclusive node) |
| Port vs official R cross-check of headline tables | - | 0 | 0 | 10 | laptop or Ada conda R |
| **Total** | | | **≈332** | **≈100** | |

**Split of the 3-4 GPU-weeks (504-672 GPU-h):**
- about 50% to AGR and its controls, as above;
- about 10% to Candidate D;
- about 40% reserve, for reruns after bugs, a possible second seed set for ablations, and the supervised
  M3 track if the humans want more seeds there.

With 4 GPUs, the frozen-phase work (headline + stability + shuffles + ablations, about 256 GPU-h) takes
about 3 days of wall-clock time. Weeks 5-6 therefore have slack.

## 8. Paper plans (ICML 2027, PLOS ONE)

### 8.0 Venue facts

**ICML 2027.**
- The call for papers is **not published**: `icml.cc/Conferences/2027/CallForPapers` and
  `icml.cc/Conferences/2027` both returned 404 on 10 Oct 2026 [V]. A third-party tracker
  (opencurious.com) estimates abstracts on 16 Jan and papers on 22 Jan 2027, in South America, and labels
  this an estimate [G].
- Everything below is the **ICML 2026** call (icml.cc/Conferences/2026/CallForPapers) [V], used as a
  proxy:
  - Abstract deadline 23 Jan 2026 AoE; full paper 28 Jan 2026 AoE; no separate supplementary deadline.
    The 2026 conference was in Seoul, 6-11 July.
  - **8 pages** of main text, "with unlimited pages for references, the impact statement, and
    appendices"; one extra page for camera-ready.
  - **Impact statement required.**
  - Claims "must be clearly stated and supported by reproducible experiments and/or sound theoretical
    analysis". There is no mandatory code release.
  - At least one author must serve as a reciprocal reviewer.
  - LLM assistance is allowed; LLMs are not authors; prompt injection means desk rejection.
  - Application-driven ML is in scope ("... driven by the needs of end-users in applications such as
    healthcare, physical sciences, biosciences ...").
  - **Dual submission:** "Authors may not submit papers that are identical, or substantially similar, to
    versions that have been previously published, accepted for publication, or submitted in parallel to
    other conferences or journals." Workshops without proceedings are exempt.
  - Notification dates are not on the call page [not verified].

**PLOS ONE** [V, journals.plos.org/plosone/s/criteria-for-publication, /s/submission-guidelines,
/s/data-availability, /s/materials-software-and-code-sharing].
- Seven criteria, with no novelty criterion. These include original research, "Results reported have
  not been published elsewhere", "high technical standard", and "Conclusions ... supported by the data".
- Methods and tools papers: "must present a proven advantage over existing alternatives, where
  applicable". Software "must be open source, deposited in an appropriate archive".
- Data: authors must "make all data necessary to replicate their study's findings publicly available
  without restriction". "Available on request" is not acceptable for non-sensitive data. A Data
  Availability Statement is mandatory.
- Code: "We expect all researchers ... in which author-generated code underpins the findings ... to make
  all author-generated code available without restrictions upon publication". It may be required as a
  condition of publication.
- Preprints: "does not impact consideration".
- No length limit; abstract ≤300 words; LaTeX submissions as PDF.
- **Not verified:** the publication fee and India pricing (plos.org/publish/fees returned 403). Check this
  before February.

### 8.1 ICML paper (only if G3 passes on 23 Nov)

**G3 criteria:**
- P1 and P3 (s5.4) hold on ≥6/8 frozen datasets;
- the gate ablation shows the cell-level gate matters;
- synthetic recovery holds.

Otherwise we do not submit. A rejection of a thin ICML paper costs December.

**Working title:** *Popular or specific? Accessibility-gated self-supervised regression and
specificity controls for multiome gene-regulatory inference.*

**Outline (8 pages).**
1. Introduction: why multiome GRN benchmarks reward popularity; our contributions.
2. Setup: BEAR-GRN and its metrics.
3. Pathologies and controls, with small formal statements:
   - AUROC on own edges is not invariant to output selection;
   - the swap gap is 0 for any TF-agnostic scorer;
   - the decomposition into coverage, target marginal, TF marginal and interaction.
4. AGR: the model, motif-specific gating as a source of TF identifiability, and conditional-gain
   scoring.
5. Experiments: synthetic recovery; BEAR (9 datasets x 6 GTs); the TF-specificity battery; shuffles and
   stability; ablations; runtime.
6. Related work.
7. Limitations.

**Claims we could honestly make (if the data allow).**
- (a) BEAR's headline metrics are dominated by coverage and node marginals.
- (b) AGR is the first label-free method on BEAR with a significant TF-swap gap on ≥6/8 frozen datasets,
  above motif x accessibility.
- (c) At equal coverage, AGR beats all published methods and all label-free baselines on ChIP
  AUROC/AUPRC, if that is what the data show.

**Figures.**
1. Decomposition bars (coverage, target, TF, interaction) per dataset.
2. The swap-gap forest plot across methods and datasets.
3. AGR schematic.
4. Raw vs coverage-completed vs TF-normalised AUROC heatmap.
5. Ablations and shuffles.
6. Stability and runtime.

**Likely reviewer attacks and how we pre-empt them.**

| Attack | Pre-emption |
|---|---|
| "Incremental over Pando/LINGER" | Head-to-head ablation that turns AGR into a Pando-like model (no gate, marginal scores) |
| "One benchmark" | BEAR is 9 datasets, 2 species, 6 GTs; plus synthetic data. *Our memory note says PBMC was dropped. ICML reviewers may still ask for a second benchmark, and we accept that risk* |
| "Effects are tiny" | Report effect sizes with TF-bootstrap CIs; argue that the size of the TF-specific channel is itself the finding |
| "Ground truth is noisy" | Agree. Show results hold on the Core and Union GTs |
| "Why ICML?" | Frame as application-driven ML: evaluation methodology plus a model, matching the call's wording |

**Odds [I]:** about 5% acceptance if submitted. Conditional on G3 passing strongly, about 12%.

### 8.2 PLOS ONE paper (primary; submit by 28 Feb 2027)

**Working title:** *What does a multiome GRN benchmark reward? Coverage, marginal and TF-specificity
controls on BEAR-GRN* (F1/F3 version). In the F2 version, the subtitle adds "... and a label-free
accessibility-gated baseline".

**Outline.**
- Introduction.
- Methods: datasets and GTs; scoring port with R cross-check; baseline panel; coverage completion; the
  TF-specificity battery; AGR (F2 only); statistics.
- Results:
  1. Coverage explains AUPRC.
  2. Gene-level priors match the published methods.
  3. Published methods' TF-specificity, with LINGER's TF-level allocation.
  4. Equal-coverage re-ranking of methods.
  5. AGR (F2 only).
  6. Stability, shuffles and runtime.
- Discussion: concrete recommendations for benchmark designers, i.e. report coverage-completed AUPRC,
  TF-normalised AUROC and the swap gap; limitations.
- Data and code availability: GitHub plus a Zenodo DOI for code, all scores, and the derived tables.
  BEAR inputs are cited, not redistributed.

**Main claims:** (a) and the equal-coverage re-ranking, which are near-certain given Appendix A; (b) and
(c) only if supported.

**Reviewer attacks.**
- "Adversarial to BEAR's authors." Pre-empt: contact them first, share the tables, and invite comment.
- "Port fidelity." Pre-empt: R cross-check.
- "LINGER's TF-level signal is legitimate." Pre-empt: we agree, and we report raw scores too.
- "Utility." Pre-empt: the advantage over existing alternatives is the control panel, which PLOS's
  methods-paper criterion accepts "where applicable".

**Odds [I]:** about 65% eventual acceptance. A first decision is likely after 2-4 months [G].

### 8.3 Dual-submission plan [I, verify with the ICML PCs once the 2027 call is out]
- If F3 (ICML go): the ICML paper is the **method** (AGR plus the controls it needs). The PLOS ONE paper
  is the **audit of existing methods**, without AGR. The two papers share the control definitions, so
  they must cite each other's preprints. The audit preprint goes on bioRxiv before the ICML deadline, so
  it is prior public work rather than a parallel submission of the same content.
- Ask the ICML 2027 PCs whether a journal submission of the audit during ICML review counts as
  "substantially similar". If the answer is unclear, delay the PLOS submission until ICML decisions,
  which would miss the "before end of February" target.
- If F1 or F2: PLOS ONE only.

## 9. Risk register

Likelihoods are our judgement [I].

| # | Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|---|
| R1 | AGR shows no TF-specific signal beyond motif x accessibility (G2 fails) | 55-60% | High for ICML, low for PLOS | F1 audit paper is already mostly computed; AGR becomes a controlled baseline |
| R2 | AGR's gains come only through the TF-activity term A_t | 40% | Medium | Factorised reporting; TF-normalised AUROC; report it as a TF-level finding, not specificity |
| R3 | Small datasets (K562 411 cells, macrophage about 500) give too few metacells | 50% | Medium | Overlapping kNN metacells or smoothed single cells; report per-dataset cell counts next to results |
| R4 | Mouse motif coverage too low (JASPAR matches only about 50/230 TF nodes) | 60% unless fixed | Medium | Add CIS-BP / HOCOMOCO v12 in week 1; expression-only tier with a stronger penalty |
| R5 | We overfit to K562 through peeks | 25% | High (credibility) | Peek budget of 5, a logged peek file, label-free selection rule, frozen 8 datasets |
| R6 | Port vs R discrepancies, e.g. the trailing-space bug found today | 30% that another appears | Medium | Cross-check headline tables with BEAR's R functions in a conda R env before submission |
| R7 | BEAR authors object to the audit framing, or change GT files or scoring | 20% | Medium | Contact them first; pin Zenodo versions and hashes; frame constructively |
| R8 | Ada unavailable or blocked (it happened on 9 Oct) | 30% for ≥3 days | Medium | AGR is light enough for the laptop 4060 for everything except the stability, shuffle and ablation sweeps; front-load Ada jobs |
| R9 | Laptop disk (6.7 GB free) | High | Low | Keep released GRNs and outputs on Ada `/share1`; laptop holds only processed features and scores |
| R10 | ICML 2027 dates or rules differ from 2026 | 30% | Low-medium | Re-check when the call appears (expected about Nov-Dec) |
| R11 | Dual-submission conflict between ICML and PLOS | 30% if both are pursued | High | Separate papers (s8.3); ask the PCs; default to PLOS-only |
| R12 | Coverage-completed rescoring is criticised as unfair to sparse methods | 40% | Low | Report raw too; frame it as an additional view, not a replacement |
| R13 | Candidate D data wrangling overruns | 50% | Low | 4-day timebox; optional |
| R14 | Writing load: two papers in Dec-Feb with two people | 50% | Medium | Decide at G3; the PLOS audit writing starts 1 Nov, independent of AGR |

## 10. The first 7 days (10-16 Oct), in order

Person A is the model owner, Person B the evaluation owner, and the advisor (Adv) makes decisions and
contacts the BEAR authors.

1. **Fri 10 Oct (all).** Present this plan and take the 3 decisions in s0. Adv sends the s1.6 questions
   to the BEAR authors (yuzun@pennstatehealth.psu.edu per SI).
2. **Fri 10 Oct (B).**
   - Apply the one-line `.str.strip()` fix in `bm.load_gt` and `bear_data.read_gt_pairs`.
   - Re-score naive-mESC KO for all released methods and baselines.
   - Commit `docs/selfsup_plan/analysis/` (scripts plus CSVs).
3. **Sat-Sun 11-12 Oct (B).**
   - Ada CPU job (compute node, internet): download all released GRNs for the 8 non-K562 datasets to
     `/share1`.
   - Run `labelfree_baselines.py released` and `coverage_completed.py` for all 9 datasets. Add the
     100-permutation version and the TF-bootstrap CIs.
   - This completes the audit tables (Appendix A, extended).
4. **Sat 11 Oct (A).** TF census files: human from Lambert et al. 2018; mouse from AnimalTFDB. Record
   versions. Check that every GT TF is in the census, and list those that are not.
5. **Sun-Mon 12-13 Oct (A).**
   - Extend `bear_motif.py` from the top-5 peaks per gene to all peaks within ±100 kb.
   - Add CIS-BP (and/or HOCOMOCO v12) PWMs. Check their licences.
   - One Ada GPU job per genome. Output: a sparse peak x TF motif matrix.
6. **Mon 13 Oct (A).**
   - Metacell builder: joint PCA (RNA) + LSI (ATAC) kNN, k ∈ {20, 50}, overlap cap.
   - Gate tensor a_cgt as a sparse product of RP x motif x accessibility.
   - Test on K562 and E8.5 rep1 inputs, no GT.
7. **Tue 14 Oct (B).**
   - Synthetic simulator v0: a planted TF→peak→gene network, motif-gated accessibility, metacell noise.
   - Unit tests: swap gap = 0 for TF-agnostic scores; gap > 0 on planted edges; coverage-completion
     invariants.
8. **Wed 15 Oct (A).**
   - AGR v0: batched PyTorch proximal-gradient sparse regression; masked-metacell CV; conditional-gain
     head.
   - Run on synthetic data, then on K562 *inputs* (no scoring).
9. **Wed 15 Oct (B).** Label-free selection criteria (held-out R², stability JI, TF-permutation
   self-consistency) and the `peek_log.md` template.
10. **Thu 16 Oct (all).**
    - Write pre-registration v0 (s5.4) and commit it.
    - Set up a conda R env (r-base, tidyverse, pROC, PRROC) on Ada or the laptop and reproduce one
      published table with BEAR's own R function.
    - Week-1 review: is G1 on track for 23 Oct?

## Appendix A. Numbers we computed (commands and outputs)

**How these were produced.**
- All analyses ran on the laptop CPU on 10 Oct 2026, each run under 1 min and under 2.4 GB RSS.
- Scripts: `docs/selfsup_plan/analysis/`. Raw CSVs sit next to them, and `labelfree_summary.md` holds
  every GT.
- Scoring follows BEAR's rules, using the vectorised re-implementation in `labelfree_baselines.py`
  (`Scorer`). AUROC is the mean of 3 draws.
- **Check:** the coverage baseline reproduces the K562 ChIP AUPRC of 0.431 from `results/summary_tables.md`.
- Inputs are the label-free processed features in `~/.cache/bear/data/processed/<ds>/`, built earlier by
  `scripts/21`.
- Released networks: K562 from `~/.cache/bear/data/INFERRED.GRNS/K562`. E7.5 rep1 (LINGER, CellOracle,
  TRIPOD, GRaNIE) was downloaded to the scratchpad with `scripts/20_bear_download.py`.

**Commands** (run from `/home/vishak/research/mevd_grn`):

```
for ds in K562 Macrophage_S1 Macrophage_S2 iPS mESC_E7.5_rep1 mESC_E7.5_rep2 mESC_E8.5_rep1 mESC_E8.5_rep2 Naive_mESC; do
  python docs/selfsup_plan/analysis/labelfree_baselines.py baselines --dataset $ds --n_rep 3; done   # ~15-50 s each, 2.4 GB peak (K562, /usr/bin/time -v)
python docs/selfsup_plan/analysis/labelfree_baselines.py released --dataset K562 --n_rep 3
python scripts/20_bear_download.py --root $SCRATCH/bear_dl --datasets mESC_E7.5_rep1 --what inferred --methods LINGER TRIPOD GRaNIE CellOracle
python docs/selfsup_plan/analysis/labelfree_baselines.py released --dataset mESC_E7.5_rep1 --inf_root $SCRATCH/bear_dl --n_rep 3
python docs/selfsup_plan/analysis/coverage_completed.py K562
python docs/selfsup_plan/analysis/coverage_completed.py mESC_E7.5_rep1 $SCRATCH/bear_dl
python docs/selfsup_plan/analysis/tf_marginal_probe.py K562 mESC_E7.5_rep1 mESC_E8.5_rep1 Macrophage_S1
(cd docs/selfsup_plan/analysis && python summarize.py)     # -> labelfree_summary.md
```

**Definitions.**
- All label-free scores are 1 + percentile, so every measured pair ranks above unmeasured ones, and they
  cover TF nodes x all measured genes.
- `motif x atac` uses the best relative JASPAR 2024 PWM score in the gene's top-5 RP peaks x the target's
  ATAC activity percentile. TFs without a motif fall back to 0.5 x the target's ATAC activity.
- `oracle_*` rows are computed **from the GT**, as diagnostics only:
  - *TF out-degree* = the TF's number of GT targets minus the edge itself;
  - *LOO in-degree* = the target's number of GT TFs minus the edge itself.

### A.1 ChIP: label-free baselines vs published methods (published numbers from Supp Data 2)

| dataset | random | coverage AUPRC | tgt mean expr AUROC / AUPRC | motif x atac AUROC / AUPRC | motif x atac swap gap | LINGER AUROC / AUPRC | best published AUROC | best published AUPRC | oracle TF out-degree AUROC* | oracle LOO in-degree AUROC* |
|---|---|---|---|---|---|---|---|---|---|---|
| K562 | 0.331 | 0.431 | 0.554 / 0.452 | 0.558 / 0.454 | +0.008 | 0.536 / 0.430 | 0.575 (Pando-GLM) | 0.430 (LINGER) | 0.870 | 0.654 |
| Macrophage_S1 | 0.272 | 0.331 | 0.552 / 0.350 | 0.510 / 0.333 | +0.007 | 0.596 / 0.321 | 0.650 (GRaNIE, 100 edges) | 0.332 (TRIPOD) | 0.876 | 0.562 |
| Macrophage_S2 | 0.270 | 0.342 | 0.553 / 0.360 | 0.522 / 0.348 | +0.008 | 0.547 / 0.324 | 0.606 (CellOracle) | 0.336 (TRIPOD) | 0.866 | 0.565 |
| iPS | 0.091 | 0.137 | 0.526 / 0.140 | 0.478 / 0.134 | +0.002 | 0.501 / 0.108 | 0.568 (GRaNIE) | 0.119 (TRIPOD) | 0.864 | 0.572 |
| E7.5 rep1 | 0.155 | 0.274 | 0.527 / 0.286 | 0.532 / 0.296 | +0.004 | 0.652 / 0.241 | 0.652 (LINGER) | 0.241 (LINGER) | 0.922 | 0.596 |
| E7.5 rep2 | 0.153 | 0.268 | 0.519 / 0.276 | 0.528 / 0.285 | +0.003 | 0.644 / 0.225 | 0.644 (LINGER) | 0.225 (LINGER) | 0.929 | 0.585 |
| E8.5 rep1 | 0.155 | 0.265 | 0.525 / 0.275 | 0.531 / 0.283 | +0.004 | 0.658 / 0.227 | 0.658 (LINGER) | 0.227 (LINGER) | 0.924 | 0.597 |
| E8.5 rep2 | 0.155 | 0.263 | 0.529 / 0.275 | 0.530 / 0.282 | +0.004 | 0.656 / 0.224 | 0.656 (LINGER) | 0.224 (LINGER) | 0.923 | 0.600 |
| Naive mESC | 0.153 | 0.269 | 0.521 / 0.277 | 0.514 / 0.281 | +0.001 | 0.473 / 0.160 | 0.556 (SCENIC+) | 0.192 (TRIPOD) | 0.926 | 0.587 |

\* GT-derived diagnostics, not methods.

Counts:
- coverage AUPRC ≥ best published on **8/9** (Macrophage S1: 0.331 vs 0.332);
- target mean expression AUPRC > best published on **9/9**;
- no single label-free baseline beats the best published ChIP AUROC on any dataset. The closest is K562,
  with 0.558 vs 0.575.

Motif coverage: JASPAR matched only 104/204 (K562), 19/28 and 18/27 (macrophage), 12/30 (iPS) and 48-57
of 226-242 (mouse) TF nodes. This is a known weakness of the current scan.

Other GTs (KO, Union, Intersection, Core, cell-type-exclusive) are in `analysis/labelfree_summary.md`.
On KO and Intersection, label-free baselines are at or below random, like every published method.

### A.2 TF-specificity of released networks (ChIP)

Definitions:
- **TF-permuted:** TF labels shuffled among the method's TFs (3 draws).
- **TF-normalised:** within-TF percentile ranks.
- **Target projection:** each edge gets its target's mean score.
- **Swap gap:** mean over TFs of AUROC(own row) − mean AUROC(other rows) on the TF's labels, across the
  method's tested targets, with unpredicted pairs scored 0.

| dataset | method | AUROC raw | TF-permuted | TF-normalised | target projection | AUPRC raw | TF-permuted | TF-normalised | swap gap | frac TFs own > 95th pct |
|---|---|---|---|---|---|---|---|---|---|---|
| K562 | CellOracle | 0.539 | 0.531 | 0.523 | 0.536 | 0.390 | 0.376 | 0.388 | +0.010 | 0.32 |
| K562 | SCENIC+ | 0.467 | 0.460 | 0.459 | 0.477 | 0.332 | 0.332 | 0.332 | +0.001 | 0.16 |
| K562 | Pando-GLM | 0.576 | 0.547 | 0.562 | 0.568 | 0.348 | 0.347 | 0.348 | +0.001 | 0.08 |
| K562 | LINGER | 0.536 | 0.537 | 0.553 | 0.549 | 0.430 | 0.432 | 0.436 | +0.007 | 0.38 |
| K562 | FigR | 0.503 | 0.506 | 0.504 | 0.511 | 0.376 | 0.378 | 0.376 | +0.002 | 0.10 |
| K562 | TRIPOD | 0.512 | 0.509 | 0.516 | 0.524 | 0.405 | 0.409 | 0.406 | +0.012 | 0.23 |
| K562 | GRaNIE | 0.512 | 0.514 | 0.521 | 0.529 | 0.368 | 0.376 | 0.369 | +0.008 | 0.25 |
| K562 | Pando-XGB | 0.518 | 0.501 | 0.515 | 0.530 | 0.346 | 0.344 | 0.346 | +0.001 | 0.09 |
| E7.5 rep1 | CellOracle | 0.530 | 0.462 | 0.503 | 0.503 | 0.191 | 0.174 | 0.190 | +0.033 | 0.34 |
| E7.5 rep1 | LINGER | 0.653 | 0.509 | 0.541 | 0.538 | 0.241 | 0.215 | 0.220 | +0.005 | 0.25 |
| E7.5 rep1 | TRIPOD | 0.491 | 0.510 | 0.495 | 0.493 | 0.193 | 0.197 | 0.194 | +0.005 | 0.11 |
| E7.5 rep1 | GRaNIE | 0.570 | 0.557 | 0.535 | 0.533 | 0.188 | 0.179 | 0.185 | +0.024 | 0.33 |

Reading:
- On K562, nothing the methods emit depends on TF identity in a way BEAR's metrics can see.
- On E7.5, LINGER's lead is between-TF allocation: its swap gap is no larger than motif x accessibility's
  (+0.004).
- CellOracle and GRaNIE have the largest genuine swap gaps, which is consistent with motif-restricted
  regression.
- Caveat: there is no CI yet. The 100-permutation version with TF-bootstrap CIs is a week-1 task.

### A.3 Coverage-completed released networks (ChIP)
Every measured TF-node x gene pair a method did not emit is appended with score 0.5 x its minimum
positive score.

| dataset | method | evaluable edges (raw) | AUROC raw → completed | AUPRC raw → completed |
|---|---|---|---|---|
| K562 | CellOracle | 215,375 | 0.539 → 0.542 | 0.390 → 0.454 |
| K562 | FigR | 261,330 | 0.503 → 0.516 | 0.376 → 0.439 |
| K562 | GRaNIE | 182,297 | 0.512 → 0.515 | 0.368 → 0.439 |
| K562 | LINGER | 707,264 | 0.536 → 0.560 | 0.430 → 0.452 |
| K562 | Pando-GLM | 64,367 | 0.576 → 0.505 | 0.348 → 0.436 |
| K562 | Pando-XGB | 59,356 | 0.518 → 0.504 | 0.346 → 0.435 |
| K562 | SCENIC+ | 3,023 | 0.467 → 0.500 | 0.332 → 0.431 |
| K562 | TRIPOD | 416,415 | 0.512 → 0.547 | 0.405 → 0.449 |
| E7.5 rep1 | CellOracle | 127,612 | 0.530 → 0.520 | 0.191 → 0.295 |
| E7.5 rep1 | GRaNIE | 89,103 | 0.570 → 0.517 | 0.188 → 0.296 |
| E7.5 rep1 | LINGER | 496,595 | 0.653 → 0.531 | 0.241 → 0.313 |
| E7.5 rep1 | TRIPOD | 282,998 | 0.491 → 0.513 | 0.193 → 0.282 |

The completed universe has 1,188,090 evaluable pairs on K562 and 2,043,450 on E7.5.

Reading:
- At equal coverage, the bar for a dense label-free method on K562 ChIP is about 0.56 AUROC and 0.454
  AUPRC.
- On E7.5 it is about 0.53 AUROC and 0.313 AUPRC. The 0.313 belongs to LINGER, which is the one place a
  published method still clearly leads the label-free baselines (0.296).

### A.4 BEAR modality shuffles (from `bear_grn_modality_shuffle_tidy.csv`, Supp Data 7)
Median change vs Original across the 8 datasets.

| Method | ΔAUPRC cell / gene / ATAC-gene | ΔAUROC cell / gene / ATAC-gene |
|---|---|---|
| CellOracle | 0.000 / -0.002 / -0.002 | -0.000 / -0.038 / -0.039 |
| FigR | -0.003 / -0.005 / +0.001 | +0.002 / -0.009 / +0.000 |
| GRaNIE | +0.002 / -0.003 / -0.029 | -0.031 / -0.070 / -0.073 |
| LINGER | -0.064 / -0.064 / -0.064 | -0.092 / -0.149 / -0.151 |
| Pando-GLM | 0.000 / -0.002 / -0.002 | -0.002 / +0.073 / +0.038 |
| Pando-XGB | -0.000 / -0.002 / -0.003 | +0.001 / -0.039 / -0.030 |
| SCENIC+ | -0.001 / -0.000 / -0.001 | +0.005 / +0.030 / +0.025 |
| TRIPOD | 0.000 / -0.001 / -0.002 | +0.007 / -0.005 / -0.004 |

LINGER on mouse loses 0.15-0.20 AUROC under every shuffle. Its evaluable edges also shrink (E7.5 rep1:
496,595 → 68,400, cell shuffle).

Command: pivot the tidy CSV by (`dataset`, `Method`) x `condition`, subtract `Original`, and take the
median per method. This is a one-off pandas snippet; the numbers above are its output.

### A.5 Label-free TF-level probes (ChIP AUROC; `tf_marginal_probe.csv`)

| score | K562 | Macrophage S1 | E7.5 rep1 | E8.5 rep1 |
|---|---|---|---|---|
| TF expression only | 0.457 | 0.592 | 0.486 | 0.478 |
| TF detection only | 0.458 | 0.603 | 0.487 | 0.481 |
| TF motif hit rate only | 0.535 | 0.524 | 0.508 | 0.519 |
| TF motif hit rate x target expression | 0.552 | 0.556 | 0.522 | 0.531 |
| TF expression x target expression | 0.508 | 0.601 | 0.512 | 0.503 |

None comes close to the oracle TF-marginal (0.87-0.93) or to LINGER's raw mouse AUROC. Macrophage has
only 28 TF nodes, so those values are noisy.
