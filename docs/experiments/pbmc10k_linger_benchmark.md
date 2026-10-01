# MeVD-GRN on PBMC10k Multiome, scored against LINGER-style Cistrome ChIP

Status: **design pre-registered (2026-10-01, before any result existed); code being built.**
This is a living document. It is updated as the build progresses, so nothing is lost if a session is interrupted.

Why this benchmark: see [../reference/multiome_grn_benchmark_consensus.md](../reference/multiome_grn_benchmark_consensus.md).
10x `pbmc_granulocyte_sorted_10k` Multiome with LINGER's blood Cistrome ChIP-seq evaluation is the most-used
paired-multiome GRN benchmark. LINGER (Nat Biotechnol 2024), KEGNI (Genome Biol 2025), scTFBridge
(Nat Commun 2025) and regX reuse it head to head.

---

## 1. Sources

Local copies of the small reference files are in `~/.cache/pbmc/ref/` and `~/.cache/pbmc/raw/` (not in git).

| What | Where |
|---|---|
| 10x data | `https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5` (192,125,528 bytes; 11,898 barcodes; 36,601 genes + 143,887 peaks). LINGER's tutorial uses the same matrix, as the `.tar.gz` MatrixMarket version. |
| Cell labels | LINGER's `PBMC_label.txt`, Google Drive id `17PXkQJr8fk0h90dCkTi3RGPmFNtDqHO_` (362,958 bytes). 9,543 barcodes, 14 cell types. Counts: classical monocytes 1,848, naive CD4 T 1,373, naive B 282, myeloid DC 232. These are exactly KEGNI's counts. |
| ChIP ground truth | 20 Cistrome DB "gene score" files, `<CistromeID>_gene_score_5fold.txt` (BETA-style RP, decay 10 kb, hg38 refGene; 112k transcript rows each). Downloaded by the helper to `~/.cache/pbmc/ref/groundtruth/<ID>_<TF>_gene_score_5fold.txt`; source URLs in the download script. |
| Dataset list | LINGER Supplementary Table 1 (`LINGER_SuppTables_MOESM2.xlsx`, sheet T1) |
| LINGER per-dataset numbers | LINGER Supplementary Table 7 (trans-regulation, cell-type-specific): LINGER / GENIE3 / PCC / PIDC / SCENIC+ AUC and AUPR ratio for each ChIP dataset |
| LINGER evaluation code | `LingerGRN/Benchmk.py::bm_trans` (package `LingerGRN==1.110`, repo github.com/Durenlab/LINGER, docs/Benchmark.md). KEGNI's `benchmark_cistrome.ipynb` copies it verbatim. |
| LINGER general GRN (`data_bulk.tar.gz`) | Google Drive id `1jwRgRHPJrKABOk7wImKONTtUupV7yJ9b`, **20,812,483,490 bytes**. Needed only to re-run LINGER, so it goes to /share1 on Ada, never the laptop. |
| TSS annotation | GENCODE v38 `gencode.v38.annotation.gtf.gz` (46,556,621 bytes) |
| Training labels | CollecTRI via OmniPath `interactions?datasets=collectri&genesymbols=yes&organisms=9606` (64,515 rows); DoRothEA A/B via OmniPath `datasets=dorothea&dorothea_levels=A,B` (15,266 rows) |

**LINGER released no PBMC output networks.** Neither the repo tree nor the docs contain them. The
LINGER numbers in section 2 are therefore quoted from Table S7. The LINGER baseline on our identical cells
comes from re-running `LingerGRN==1.110` on Ada (section 8).

## 2. Numbers to beat

**LINGER Supplementary Table 7** (cell-type-specific trans-regulation, PBMC10k). 19 of the 20 datasets are
evaluable; FOXP3/44098 is N/A in the table. The means are computed by us from the table's rows:

| Method | mean AUC (19) | mean AUPR ratio (19) |
|---|---|---|
| **LINGER** | **0.7143** | **2.2526** |
| SCENIC+ | 0.5481 | 1.2905 |
| PCC | 0.5408 | 1.2287 |
| GENIE3 | 0.5387 | 1.1686 |
| PIDC | 0.5292 | 1.1717 |

- The 0.7143 matches the "0.714" in the consensus doc (KEGNI's number), and the AUPR ratio 2.2526 matches
  its "about 2.25". That was a figure-derived value; it is now exact.
- LINGER per cell type (mean AUC): classical monocytes 0.7111 (n=10), naive CD4 T 0.7071 (n=5, FOXP3
  excluded), naive B 0.7220 (n=3), myeloid DC 0.7583 (n=1).
- Per-dataset LINGER AUC values are in Table S7. Examples: STAT1/41287 0.7573, RUNX1/81223 0.5988 (the
  lowest), RUNX1/8481 0.7583 (the highest).
- KEGNI: mean AUROC 0.699 (KEGNI paper). Its notebook prints RUNX1/44097 in naive CD4 T: AUC 0.7138,
  AUPR ratio 2.5004.
- scTFBridge (published, STAT1 in CD14 monocytes): AUC 0.693 / AUPR ratio 2.221. It quotes LINGER at
  0.719 / 2.105, which differs from Table S7's four STAT1 datasets. It used its own 17-TF subset and gene
  intersection.

Secondary protocol (EpiAwareNet, KDD 2026): a different PBMC file, DoRothEA labels and a random 80/20 edge
split. AUROC 0.8218 (±5 kb). Not implemented (section 9).

## 3. Protocol pinned down from code (LINGER `Benchmk.bm_trans`, verbatim)

- **Datasets:** 20 ChIP-seq sets (Table S1). The cell type is the *PBMC cell type whose cell-type-specific
  GRN is scored*.

| CistromeID | TF | Cell type |
|---|---|---|
| 5967 | MYC | naive B cells |
| 8481 | RUNX1 | myeloid DC |
| 40215 | IRF4 | naive B cells |
| 41287, 41288, 41289, 41290 | STAT1 | classical monocytes |
| 41301, 41302, 41303 | IRF1 | classical monocytes |
| 44092, 44093, 44094 | ETS1 | naive CD4 T cells |
| 44097 | RUNX1 | naive CD4 T cells |
| 44098 | FOXP3 | naive CD4 T cells (N/A in LINGER Table S7) |
| 45178 | CTCF | naive B cells |
| 45444 | CTCF | classical monocytes |
| 47435 | REST | naive CD4 T cells |
| 81223 | RUNX1 | classical monocytes |
| 85986 | SPI1 | classical monocytes |

  That is **10 distinct TFs** (MYC, RUNX1, IRF4, STAT1, IRF1, ETS1, FOXP3, CTCF, REST, SPI1) in **4 cell
  types** (classical monocytes, naive CD4 T, naive B, myeloid DC).
- **Target definition:**
  - Read the Cistrome file (`skiprows=5`, tab-separated) and take the **max `score` per `symbol`** over
    RefSeq transcripts.
  - Sort in descending order and take the **top N = 1000 symbols** as positives. The code comment says
    "top 500", but the code uses 1000.
  - The top 1000 is chosen from the *whole* Cistrome gene list, before any intersection with the method's
    genes.
- **Candidate space:** the TG rows of the method's own trans matrix (`TGset = data2.index`). Every TG row is
  scored, and a row is labelled 1 if its symbol is in the top-1000 list. LINGER's cell-type-specific TG set
  is the set of genes with ≥1 cis link in `cell_type_specific_cis_regulatory_<ct>.txt`, which depends on
  `data_bulk`. KEGNI uses all genes in its embedding. **The candidate space is therefore method-dependent in
  the published tables.**
- **Metrics:**
  - AUC = `sklearn.metrics.roc_auc_score(d1, Score)`
  - AUPR ratio = `average_precision_score(d1, Score) * len(d1) / sum(d1)`, i.e. AP over the positive
    fraction
- **Cells (LINGER `preprocess.get_adata` + tutorial):**
  - barcodes in `PBMC_label.txt`
  - RNA `pct_counts_mt < 5`
  - `filter_cells(min_genes=200)` on RNA and on ATAC
  - `filter_genes(min_cells=3)` on both
  - the intersection of the surviving barcodes
  - Genes are symbols (`features[1]`) with `var_names_make_unique`.
- **Scored GRN:** the *cell-type-specific* trans-regulatory matrix of the dataset's cell type.

## 4. Pre-registered design (fixed 2026-10-01, before any PBMC result)

MeVD-GRN is supervised. LINGER and the other baselines are not trained on TF-target labels. The design
therefore has to guarantee that MeVD-GRN never sees the evaluation TFs' labels, and that its gain is not
something a trivial ranker would get too.

### 4.1 Labels
- **Primary training labels: CollecTRI** (literature-curated, not motif- or accessibility-derived).
  - Multi-member complexes (AP-1 / NF-kB heteromers) are dropped.
  - A single-member `COMPLEX:<uniprot>` maps to its gene.
- **Secondary: DoRothEA levels A-B only.** About 14% of these edges carry TFBS (motif) support among other
  evidence, which is one reason the motif relation is off everywhere (4.4).
- **TF-disjoint from the evaluation.** Every TF in the LINGER Cistrome evaluation set is removed as a
  *regulator* from the training labels, in every split regime. Such a TF may still appear as a *target* of
  another TF, as in InfoSEM's definition. The removed edge count is logged per build.
- The evaluation ChIP set is never used for training, model selection or early stopping.

### 4.2 Splits of the training label set (70 / 15 / 15, split seed 0)
| Regime | Unit | Role |
|---|---|---|
| `tf` | regulator | **headline**. Val/test TFs never appear as regulators in training, just like the Cistrome TFs. |
| `target` | target gene | train pairs, including negatives, only touch train targets |
| `random` | edge | comparability with EpiAwareNet-style numbers only |

In every regime the **headline metric is the LINGER-protocol Cistrome evaluation**, which is always
TF-held-out (4.1). The internal val/test parts are used for model selection (val) and as a
supplementary report (test).

### 4.3 Negatives
- **Train:** a degree-matched pool.
  - Each train positive (t, g) spawns 4 x `neg_ratio` negatives (t, g'). They keep the same TF, so TF
    out-degree is preserved.
  - g' is drawn with probability proportional to its in-degree among the train positives, mixed with 20%
    uniform draws, so genes with no known regulator are still seen as negatives.
  - The trainer resamples `neg_ratio` (5) x n_pos negatives from this pool every epoch.
- **Val:** uniform same-TF negatives (20:1, fixed). This mirrors the headline per-TF ranking over all genes.
- **Test (internal):** both uniform (20:1) and degree-matched 1:1 (Stock et al. 2025).
- **General rules:**
  - Negatives never include a known positive of the label source or a self pair.
  - `MEvDTrainer.restrict_negative_pool` stays on (`exclude_eval_negatives: true`).
  - `train_stage` refuses to run without it.

### 4.4 Graphs and features (label-free)
- **Node set:** one gene universe shared by all cell types, so indices are comparable (section 3 records the
  exact rule). Per-cell-type features come from that cell type's cells:
  - RNA: mean / var / detection and the co-expression signature
  - ATAC: RP-weighted gene activity plus the 4-dim locus descriptor, from GENCODE v38 TSS
  - FM: Geneformer V2-104M token embeddings, which are cell-type-invariant
- **The co-expression kNN and the TF-candidate graph are built label-free** (`exclude_positives=False`).
  They are never built from label sets. This avoids the `02_preprocess.py` leak documented in
  `scmultiomegrn_generalization.md` section 10.
- **The motif relation is off** (`use_motif: false`) for every variant.
  - DoRothEA A/B is partly TFBS-supported.
  - The ChIP evaluation measures binding, which a motif scan of the same peaks would partly encode.
- **No per-gene ID embeddings exist in MeVD-GRN.** An unseen TF is represented only by its features, the FM
  embedding and its graph neighbourhood.

### 4.5 Training and model selection
- One stage (single tier), all-at-once. Hard negatives are off, since they are undefined with one tier.
- `neg_ratio` 5, lr 1e-3, AdamW, eval every 5 epochs.
- Early stopping on **validation AUPR only** (MeVD-GRN's native rule), patience 20 checks, max 300 epochs.
- The best-val checkpoint is used. Test and Cistrome are touched once, after training.
- **5 model seeds (42-46).** The split seed is fixed at 0, so every seed sees identical splits.
- One model per evaluated cell type, trained on that cell type's features. The labels are not
  cell-type-specific; only the features are.

### 4.6 Variants
| Variant | What changes | Purpose |
|---|---|---|
| `fm_h384` | FM + hidden 384 / 2 GNN layers | **headline** (the project's recommended config) |
| `fm_h384_rnaonly` | the same, with ATAC features, locus descriptor and accessibility gate of the candidate graph all zeroed or removed | **key ablation: does ATAC add anything?** |
| `base` | h128, no FM | smaller model |
| `base_rnaonly` | base, RNA only | ablation at small size |
| `fm_h384_uniformneg` | uniform instead of degree-matched training negatives | negative-sampling sensitivity |

The full grid: headline + ablation x 3 regimes x 5 seeds, and the other variants on the `tf` regime x 5 seeds.
All of it runs on CollecTRI. DoRothEA A-B runs only the headline + ablation, `tf` regime.

### 4.7 Baselines (identical evaluation set, identical cells)
- **Mandatory trivial baselines:**
  - **degree-only:** score(t, g) = in-degree of g in the train labels. It ignores the TF, which is all
    any model can know about an unseen TF from labels alone.
  - **gene-ID-only (InfoSEM-style):** logistic regression on one-hot(TF) + one-hot(target), trained on the
    same train positives and negatives. For an unseen TF it reduces to a learned target bias.
  - **Pearson:** |r| and signed r between TF and target expression within the cell type.
  - **GRNBoost2:** arboreto's SGBM routine, per cell type, with regulators = label TFs ∪ evaluation TFs.
- **Real baselines:**
  - **LINGER:** its released cell-type-specific trans-regulatory outputs, if available. Otherwise a re-run
    on Ada.
  - **SCENIC+, scTFBridge:** re-run on the same cells and genes if their released outputs or code allow.
    Otherwise their published numbers are quoted, clearly labelled as such.
  - All baselines are scored by the same evaluation code (section 3), over the same candidate genes.

### 4.8 Reporting
- For each variant: the mean ± std over 5 seeds of the per-dataset AUROC and AUPR ratio. Also the mean over
  datasets, which is the number that goes in the LINGER / KEGNI table.
- The same metrics for every baseline.
- Internal test results per regime:
  - per-TF mean AUROC over all genes
  - pooled AUROC/AUPR vs uniform negatives
  - pooled AUROC/AUPR vs degree-matched 1:1 negatives
- Results are also stratified by evaluation-TF target count.
- **Win criterion (pre-registered):** the headline variant's mean AUROC over the Cistrome datasets (5-seed
  mean) exceeds LINGER's 0.714 under the same code, and exceeds the degree-only and gene-ID baselines.
  The ATAC ablation is reported whatever its sign.

## 4a. Amendments (dated; LINGER's code wins for the Cistrome evaluation)

- **2026-10-01, A1 (cells):** MeVD-GRN's per-cell-type features use exactly LINGER's post-QC cells. These
  are the label file, mt < 5%, ≥200 genes in RNA and ≥200 peaks in ATAC, and the intersection. It replaces
  the "that cell type's cells" of s4.4, which did not specify QC.
- **2026-10-01, A1b (no mitochondrial filter):**
  - `LingerGRN==1.110`'s `get_adata` drops cells with RNA `pct_counts_mt >= 5`. On this dataset the
    median mt fraction is 9.7%, so that filter keeps only **514 of 9,543** labelled cells (4 naive B,
    4 mDC).
  - KEGNI's cell counts (1,848 / 1,373 / 282 / 232), which come from LINGER's processed data, equal the
    *unfiltered* label-file counts. So the published LINGER/KEGNI numbers did not use that filter; it was
    added to the package later.
  - Decision: **no mt filter.** QC is the label file plus ≥200 genes (RNA) and ≥200 peaks (ATAC).
  - The LINGER re-run must use the same cells, i.e. without the mt filter.
- **2026-10-01, A2 (gene universe):**
  - The node set is every gene that passes LINGER's `filter_genes(min_cells=3)` on the post-QC labelled
    cells. It is a superset of every TG set LINGER can output, plus the 10 evaluation TFs.
  - This replaces the "shared gene universe, rule in section 3" placeholder.
- **2026-10-01, A3 (candidate space):**
  - LINGER scores each method over its own TG rows. To compare methods fairly, **every method, including
    MeVD-GRN, all baselines and the LINGER re-run, is scored over the same gene set.**
  - **Primary set: `linger_tg`**, the LINGER re-run's cell-type-specific TG set. It is closest to the
    published LINGER numbers.
  - **Secondary set: `expressed`**, the gene universe of A2, which is available before the LINGER re-run
    exists.
  - Both are reported. Scores for all genes are saved (`eval_scores.npz`), so the choice does not require
    retraining.
- **2026-10-01, A4 (dataset set):**
  - The headline mean is over the **19 datasets LINGER reports** (FOXP3/44098 excluded, as in Table S7).
    This keeps it comparable with 0.7143 / 2.2526.
  - A 20-dataset mean is reported alongside when FOXP3 is evaluable.
- **2026-10-01, A5 (no change):** the evaluation TFs of s4.1 are exactly the 10 TFs of section 3.
  The label/split/negative design of s4.1-4.3 does not interact with LINGER's code, which uses no training
  labels, so it stands unchanged.
- **2026-10-01, A6 (scheduling only, no protocol change):** work-stealing unit claims, and GRNBoost2 as a
  single job. See the build log.
- **2026-10-01, A7 (scores storage):** eval-TF scores are saved as float32, not float16. float16 created
  ties that shifted AP in the 3rd decimal.

## 5. Files

| File | Purpose |
|---|---|
| `src/benchmarks/pbmc_labels.py` | label parsing, eval-TF exclusion, the 3 split regimes, degree-matched / uniform negatives, leak assertions |
| `src/benchmarks/pbmc_data.py` | single-pass 10x h5 loader; per-cell-type processed dirs (label-free co-expression kNN and TF-candidate graphs, plus an RNA-only candidate graph) |
| `src/benchmarks/pbmc_eval.py` | LINGER `bm_trans` scoring (verbatim semantics), candidate spaces, the 19-dataset summary, internal test metrics |
| `src/benchmarks/pbmc_baselines.py` | degree, gene-ID LR, Pearson (abs and signed), GRNBoost2 |
| `scripts/25_pbmc_download.sh` | `data` / `linger` / `fm` downloads on ada-gw1 into /share1 |
| `scripts/26_pbmc_build.py` | steps `cells,processed,fm,splits` |
| `scripts/27_pbmc_train_eval.py` | MeVD-GRN: one cell type x source x regime x seed x variant |
| `scripts/28_pbmc_baselines.py` | baselines, scored by the same code |
| `scripts/29_pbmc_compile.py` | re-scores every saved score file per candidate space, then summary tables |
| `configs/pbmc/benchmark.yaml` | dataset, QC, cell types, the 20 ChIP datasets, Table S7 targets |
| `configs/pbmc/mevd_base.yaml`, `mevd_fm_h384.yaml`, `mevd_fm_h384_uniformneg.yaml` | model configs (s4) |
| `src/benchmarks/pbmc_linger_run.py` | LINGER re-run driver, following LINGER's tutorial steps (`linger` env); exports `linger/tg_<ct>.txt` and eval-TF scores |
| `slurm/pbmc_common.sh`, `pbmc_build.sh`, `pbmc_work.sh` (array 0-3, work stealing), `pbmc_grnboost2.sh` (1 job), `pbmc_linger.sh`, `pbmc_submit_all.sh` | Ada jobs |

## 6. Build log
- 2026-10-01, morning: Ada SSH to the old host timed out. 10x h5, GENCODE v38 GTF, CollecTRI and DoRothEA
  A/B were downloaded to `~/.cache/pbmc/`. The design (s4) was written before any result.
- 2026-10-01, protocol extracted from LINGER code and tables (s1-s3). Amendments A1-A5 and A1b recorded.
- 2026-10-01, afternoon: Ada is back with a new layout (`ada-gw1`, /home2, partition `u22`). Downloads
  started on the gateway into `/share1/$USER/mevd_grn/pbmc`.
  - GitHub raw (ground truth) was very slow from Ada. The 20 Cistrome files and the two OmniPath label
    files were pushed from the laptop instead; md5 sums match.
  - The OmniPath label snapshot used is the laptop's 2026-10-01 download: CollecTRI md5
    `d03a4672cfb1532d4912516849fe9b97`, DoRothEA A/B md5 `0669f25ccb8406928f5c67804951bf6f`.
- 2026-10-01 15:01 (Ada clock): submitted build **473**, work **474** (array 0-3) and GRNBoost2 **475**.
  - **Decision (A6):** another agent's `leakfix` array (467, 12 tasks, %2) shares the per-user QoS (4 GPUs,
    40 CPUs, 20 submitted jobs). The work array therefore uses **work stealing**: an atomic `mkdir` claim
    under `/share1/.../claims/<jobid>/`. Whatever number of tasks the QoS lets run, they take units in
    priority order. GRNBoost2 is one job (all 4 cell types), to save submit slots.
  - Build 473 finished in about 8 min. Cells, universe and splits are identical to the local smoke build.
    Trivial baselines are done for all 4 cell types x {CollecTRI, DoRothEA A-B} x {tf, target, random}, plus
    Pearson.
  - Work tasks 474_0 and 474_1 started at 15:09 on gnode046 (2080 Ti); 474_2-3 and 475 were pending
    (`QOSMaxCpuPerUserLimit`).
  - Training: about 5.4 s/epoch. Val AUPR peaks early (0.398 at epoch 10 for classical monocytes) and
    drifts down while the training loss keeps falling, i.e. it overfits the training TFs. Early stopping at
    epoch 110 is the pre-registered rule (best-val checkpoint restored). About 9.3 min per unit.
- 2026-10-01: LINGER re-run prepared.
  - `scripts/25 linger-env` creates the conda env `linger` (LingerGRN==1.110, bedtools, pybedtools).
  - The driver is `src/benchmarks/pbmc_linger_run.py`; the job is `slurm/pbmc_linger.sh`.
  - The gateway watcher `~/pbmc_linger_watch.sh` (log `~/logs/pbmc_linger_watch.log`) submits it
    automatically once the 20.8 GB `data_bulk.tar.gz` is complete and the env imports cleanly.

## 7. Local smoke test (2026-10-01, RTX 4060 laptop, `~/.cache/pbmc/root`)

The smoke-test numbers below come from short runs and **are not results**.

**Build: `26 --steps cells,processed,fm,splits --cell_types naive_b`, 35 s**
- Cells: 9,543 labelled, all pass QC. Per cell type: classical monocytes 1,848, naive CD4 T 1,373,
  naive B 282, mDC 232. These equal KEGNI's counts.
- Universe: 25,477 genes. All 10 eval TFs are present.
- Processed dir (naive B):
  - 939 regulators in the candidate graph
  - co-expression kNN 973,000 edges; TF-candidate graph 469,500 edges (the same count for the RNA-only
    version)
  - TSS for 18,908 genes; ATAC non-zero for 70.8% of genes
  - Geneformer matched 15,804 / 25,477 genes
- Labels:
  - CollecTRI: 31,022 edges / 903 TFs / 5,091 targets after removing 2,130 eval-TF edges.
    TF regime: train 21,309 / val 4,161 / test 5,552 positives.
  - DoRothEA A-B: 8,438 edges / 282 TFs.
  - All leak assertions pass. Val/test positives overlap train in 0 pairs. Val/test TFs (tf regime) and
    targets (target regime) are disjoint from train. No eval TF is a regulator in any part.

**MeVD-GRN `fm_h384`, CollecTRI / tf, 10 epochs**
- 8.7 s/epoch. val AUPR rose from 0.357 to 0.437.
- Internal per-TF AUROC 0.897; degree-matched AUROC 0.592.
- Cistrome naive B (3 datasets): mean AUROC 0.711, AUPR ratio 1.73.

**`fm_h384_rnaonly`, target regime, 5 epochs:** runs end to end.

**Trivial baselines (naive B, CollecTRI / tf)**

| Baseline | Cistrome AUROC | AUPR ratio | Internal per-TF AUROC |
|---|---|---|---|
| degree | 0.562 | 1.22 | **0.924** |
| gene-ID | 0.560 | 1.35 | 0.920 |
| abs Pearson | 0.622 | 1.38 | n/a |
| signed Pearson | 0.498 | 1.12 | n/a |

The internal per-TF metric is dominated by target hubness: degree-only gets 0.924 on it. **So the
internal test is not evidence of anything without the degree-matched column.** The Cistrome evaluation is
where degree-only is weak.

**`29` compile:** works. Scores are saved as float32 so the AP re-score is exact.

## 8. Runbook for Ada (post-upgrade layout; see docs/reference/ada.md top section)

```bash
# laptop -> Ada (work is uncommitted; same relative paths, no git)
rsync -aR src/benchmarks/pbmc_*.py scripts/2[5-9]_pbmc_* configs/pbmc/ slurm/pbmc_*.sh \
      docs/experiments/pbmc10k_linger_benchmark.md vishakkashyap.k@ada-gw1:mevd_grn/
# gateway (ada-gw1), repo root ~/mevd_grn
bash scripts/25_pbmc_download.sh data          # ~0.3 GB into /share1/$USER/mevd_grn/pbmc
bash scripts/25_pbmc_download.sh fm            # Geneformer into ~/.cache/huggingface (home2 is visible on nodes)
GB2=1 bash slurm/pbmc_submit_all.sh            # build -> work[0-3] (+ GRNBoost2): 6 submit slots
bash scripts/25_pbmc_download.sh linger        # 20.8 GB LINGER data_bulk
bash scripts/25_pbmc_download.sh linger-env    # conda env `linger`
sbatch slurm/pbmc_linger.sh                    # LINGER re-run + export (or let ~/pbmc_linger_watch.sh do it)
# when results exist (gateway; reads /share1 directly)
python scripts/29_pbmc_compile.py --root /share1/$USER/mevd_grn/pbmc   # -> <root>/summary/summary.{md,tsv}
# resubmitting the work array resumes: finished units are skipped
SKIP_BUILD=1 bash slurm/pbmc_submit_all.sh
```

**Work grid** (`slurm/pbmc_work.sh`; 220 units, priority order). Each line is one variant x the 4 cell
types x 5 seeds:
1. fm_h384 / tf (headline)
2. fm_h384_rnaonly / tf
3. the same two variants on the target regime, then on the random regime
4. base and base_rnaonly / tf
5. fm_h384_uniformneg / tf
6. fm_h384 and fm_h384_rnaonly on DoRothEA A-B / tf

At about 9 min per unit on 2 concurrent GPUs, tier 1 (40 units) takes about 3 h, and the full grid about
17 h.

**First completed units (fm_h384, CollecTRI / tf, seed 42, candidate space `expressed`). These are early,
single-seed numbers, not the headline.**
- Classical monocytes: Cistrome mean AUROC 0.7629, AUPR ratio 2.569, over 10 datasets. LINGER Table S7 on
  the same 10 datasets: 0.7111.
- Naive CD4 T: 0.7087 / 1.705 over 5 datasets. LINGER: 0.7071.
- Internal test: degree-matched AUROC 0.684 / 0.693, against degree-only 0.576 and gene-ID 0.591 on the
  same split.
- Comparability caveat: the candidate space differs from LINGER's own TG set. The `linger_tg` re-score
  arrives with the LINGER re-run.


## 9. Not done yet, and known deviations

- **SCENIC+ and scTFBridge** have not been re-run. Until they are, only their published numbers are
  quoted: SCENIC+ 0.5481 / 1.2905 (Table S7), and scTFBridge STAT1 0.693 / 2.221 (its paper). Both are
  labelled as published, not re-run.
- **The EpiAwareNet protocol** (its own PBMC file, DoRothEA, 80/20 edge split) is not implemented. Our
  `random` regime is the closest analogue, but it is not the same file or label set, so it is not directly
  comparable.
- **The `linger_tg` candidate space** needs the LINGER re-run. Until then, every number is on the
  `expressed` space (25,477 genes), which differs from LINGER's own TG rows.
- **KnockTF and other supporting ground truths** (LINGER Table S8) are not built.
