# MeVD-GRN on scMultiomeGRN's benchmark (generalization experiment #1)

Status: **findings recorded, code being built** (living document; updated as
the build progresses so nothing is lost if a session is interrupted).

Goal: train and evaluate MeVD-GRN on the datasets and protocol used by
scMultiomeGRN (Xu *et al.*, *Nucleic Acids Research* 53(5):gkaf138, 2025;
PMID 40037709, PMC11879466), and compare against the numbers that paper reports.
SC-MO-GRN-DB results remain the bonus.

---

## 1. Sources

| What | Where |
|---|---|
| Paper | `docs/reference/scMultiomeGRN.pdf` (14 pp.) |
| Supplement (Table S1, Figs S1-S9) | Europe PMC bundle: `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11879466/supplementaryFiles` → `gkaf138_supplemental_file.pdf` (the OUP CDN link needs a signed URL) |
| Official code | Zenodo `10.5281/zenodo.14848389`, file `ScmultiomeGRN-main.zip` (8.4 MB), `https://zenodo.org/records/14848389/files/ScmultiomeGRN-main.zip?download=1`. No license is declared, so we **re-implement** their data steps with attribution rather than vendoring code. |
| Published processed ground truth / splits | **None.** The zip has code, motif/promoter resources, and one raw demo cell type (`demo/B_naive`, PBMC). The lung and hematopoiesis graphs, features, and splits are not published, so they have to be regenerated with their pipeline. |

Contents of the official zip that we use:
- `data_resource/HOCOMOCOv11.zip`: 769 `.meme` motifs covering 678 TF names (TF name = file name before `_`).
- `data_resource/gencode.v19.ProteinCoding_gene_promoter.txt` (lung, hg19) and `gencode.v38...` (PBMC/AD, hg38). Columns: `chr start_promoter end_promoter strand gene_name`, TSS ± 2 kb.
- `data_resource/download_genome/download.sh`: UCSC `hg19.fa.gz` / `hg38.fa.gz`.
- `demo/B_naive/{atac,scrna}/{matrix.mtx,barcodes.tsv,peaks.tsv|genes.tsv,var_features.tsv}`: 12,079 peaks × 382 cells (ATAC) and 15,029 genes × 817 cells (RNA). This is their input layout, and we use it for the local smoke test.

## 2. Datasets

| Dataset | Accession / URL | Species / build | Cell types used | Cells | Modalities | Role in paper |
|---|---|---|---|---|---|---|
| **Human fetal lung** (Domcke *et al.* 2020) | GEO **GSM4508936** (series GSE149683): `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4508nnn/GSM4508936/suppl/GSM4508936_lung_filtered.seurat.RDS.gz` (**3.01 GB** gz) | human, **hg19** | 9 (below) | 72,622 total (Table S1 sums to 72,662) | sci-ATAC-seq3. **Important:** the "scRNA" modality is the Seurat object's `RNA` assay, which GEO describes as a *gene-body(+2 kb upstream) accessibility count matrix used for annotation*, **not real RNA-seq**. | **Main benchmark** (Fig. 2B, Figs 3, S2-S9, Table S1) |
| Human adult hematopoiesis (Buenrostro *et al.* 2018) | scATAC: `github.com/pinellolab/scATAC-benchmarking/tree/master/Real_Data/Buenrostro_2018` (peaks `combined.sorted.merged.bed`, `metadata.tsv`; per-cell BAMs only via a Dropbox tarball, counted with `bedtools coverage`). scRNA: Buenrostro 2018 Data S2 `https://ars.els-cdn.com/content/image/1-s2.0-S009286741830446X-mmc4.zip` (138.6 MB) | human, hg19 | HSC (347 ATAC cells), CMP (502), GMP (402) | – | scATAC + scRNA (separate experiments) | Secondary benchmark (Supp. Fig. S1) |
| PBMC (10x) | scATAC `10k PBMCs ... v1.0.1`; scRNA `pbmc_10k_v3` (10x website) | human, hg38 | CD14 Mono, CD4 Memory/Naive, CD8 effector/Naive | 7,189 ATAC / 9,129 RNA | unpaired | Biology only (no metrics in paper). `demo/B_naive` comes from here. |
| Alzheimer's (Xiong *et al.* 2023) | Synapse syn52293417 | human, hg38 | microglia etc. | ~414k/437k nuclei | snRNA + snATAC | Biology only (no metrics) |

Lung cell types and Supplementary Table S1 (cell number, "Initial" = number of TF interactions in the motif ground truth, Accuracy):

| Cell type (dir name) | Cells | Initial TF interactions | Predicted | Accuracy |
|---|---:|---:|---:|---:|
| Stromal | 43,289 | 2,104 | 2,518 | 0.873 |
| Bronchiolar and alveolar epithelial | 18,988 | 2,712 | 3,504 | 0.849 |
| Vascular endothelial | 4,352 | 9,609 | 16,926 | 0.830 |
| Ciliated epithelial | 1,111 | 8,503 | 13,504 | 0.856 |
| Lymphatic endothelial | 662 | 1,889 | 2,734 | 0.847 |
| Lymphoid | 1,860 | 613 | 630 | 0.864 |
| Megakaryocytes | 568 | 930 | 1,312 | 0.831 |
| Myeloid | 1,288 | 542 | 674 | 0.867 |
| Neuroendocrine | 544 | 1,361 | 1,630 | 0.851 |

The official `src/experiment/train.py` lists the lung directories as
`Vascular_endothelial_cells, Ciliated_epithelial_cells, Bronchiolar_and_alveolar_epithelial_cells, Lymphatic_endothelial_cells, Lymphoid_cells, Megakaryocytes, Neuroendocrine_cells, Myeloid_cells` (+ Stromal cells from Fig. 2A),
built with `hg19` + `HOCOMOCOv11`.
The "Initial" counts let us check how faithfully the ground truth is regenerated.

## 3. Protocol (from the paper plus the official code; code wins where they differ)

**Ground truth ("initial adjacency matrix")**, per cell type (`src/utils.py`: `get_TFBS_from_ATAC`, `get_TFBS_from_promoter`, `build_Graph`, `build_cross_graph`):
1. Keep scATAC peaks detected in **>10 % of cells** (`rate > 0.1`), sorted BED.
2. Extract peak sequences (Bio::DB::Fasta `subseq(start=>end)`, 1-based inclusive; header `chr-start-end-`; **chrY skipped**).
3. FIMO (MEME 5.4.1) per HOCOMOCO v11 motif: `fimo --thresh 1e-4 --no-qvalue` (default max-stored-scores). Then keep hits with **p ≤ 1e-6**. TF name = `motif_id.split("_")[0]`, so redundant motifs of one TF are merged.
4. Promoters = GENCODE protein-coding TSS ± 2 kb for genes that have a motif.
5. Edge TF_a–TF_b if a TFBS of TF_b lies **fully inside** TF_a's promoter. It is **symmetrised**, so the task is an **undirected TF–TF** graph. Nodes = TFs whose promoter contains ≥1 TFBS. The diagonal is removed at load time.
6. Node set is intersected with TFs that have both node-feature types (MAESTRO RP and scRNA/GRNBoost2).

**Node / edge features (competitor's inputs):**
- ATAC: MAESTRO `scatac-genescore --model Enhanced --genedistance 1000`, species **GRCh38** (hard-coded default, even for hg19 lung). Feature vector = RP score of the TF in every cell (dim = #cells).
- RNA: GRNBoost2 (arboreto, seed 0) TF→gene importances over variable genes ∪ TFs (dim = #genes).
- Edge feature: 16×16 log joint histogram of the two TFs' log10(expr + 0.01) (paper Eq. 1-2).
- Model input graph = **training-positive adjacency** (+ self loops, GCN-normalised).

**Split** (`src/dataset.py: GraphDataModule.split_dataset`, defaults `dataset_seed=666, n_splits=10, train_split_rate=0.8, rest_all_train=True`):
- Unit = unordered TF pair (upper triangle). `KFold(10, shuffle=True, random_state=666)` over positive pairs.
- Per fold: test = that fold's positives (10 %) + an **equal number** of random non-edges.
- The remaining 90 % is split 80/20 into train (72 % of all positives) and val (18 %). Val negatives are drawn 1:1.
- Train negatives = **all** remaining non-edges. Quirk: the train/val negative pool includes the diagonal, which their loader relabels as 1.
- Each pair is evaluated in both directions. Scores are symmetric.

**Training:** Adam, lr 1e-5, max 2000 epochs, early stopping on **val loss** (patience 100), model seed 666, best checkpoint by val loss. The paper says the outputs of 10 runs are aggregated (≥6/10 votes) into the final network. Metrics are averaged over runs/folds.

**Metrics** (`src/utils.py: metric_fn`):
- AUROC = `auc(roc_curve)`
- **AUPR = trapezoidal `auc(recall, precision)`** (AP is logged separately)
- Accuracy uses threshold = **median** of the test-split scores (`pred >= thr`)

The code logs two evaluation sets:
- `test`: held-out fold, 1:1.
- `all`: every N×N entry, labels = full ground-truth adjacency. This includes the training edges and the diagonal (label 0). It uses the test-median threshold.

**Which set the paper's numbers come from:** the text says "on the test set", but in every ROC/PR figure (Fig. 3A, Figs S1-S9) the PR curves end at precision ≈ 0.13-0.37 at recall 1. That is the graph density, not the 0.5 of a 1:1 test set, and GENIE3's AUPR (0.13-0.24) equals that base rate. So the reported AUROC/AUPR are the **`all` (full-matrix, transductive) metrics**, mean ± std over folds. We will report **both**: `all` (directly comparable to their numbers) and `test` (the leakage-free held-out fold).

## 4. Target numbers (scMultiomeGRN, ATAC+scRNA fused model)

Lung, per cell type (Supp. Figs S2-S9 legends; Lymphoid from main Fig. 3A; Accuracy from Table S1):

| Cell type | AUROC | AUPR | Acc | Source |
|---|---|---|---|---|
| Bronchiolar & alveolar epithelial | 0.91 ± 0.01 | 0.76 ± 0.03 | 0.849 | Fig. S2 |
| Ciliated epithelial | 0.93 ± 0.00 | 0.75 ± 0.00 | 0.856 | Fig. S3 |
| Lymphatic endothelial | 0.93 ± 0.00 | 0.81 ± 0.00 | 0.847 | Fig. S4 |
| Vascular endothelial | 0.90 ± 0.00 | 0.67 ± 0.01 | 0.830 | Fig. S5 |
| Megakaryocytes | 0.92 ± 0.00 | 0.78 ± 0.01 | 0.831 | Fig. S6 |
| Myeloid | 0.95 ± 0.00 | 0.85 ± 0.00 | 0.867 | Fig. S7 |
| Neuroendocrine | 0.92 ± 0.01 | 0.82 ± 0.02 | 0.851 | Fig. S8 |
| Stromal | 0.93 ± 0.01 | 0.80 ± 0.03 | 0.873 | Fig. S9 |
| Lymphoid | 0.93 ± 0.01 | 0.86 ± 0.02 | 0.864 | Fig. 3A |
| **Mean of 9** | **0.924** | **0.789** | 0.852 | matches text "avg AUROC 0.924, AUPR 0.790" (Fig. 2B) |

Lung baselines (Fig. 2B, 9-cell-type averages, AUROC / AUPR):
GENIE3 0.540 / 0.242, DeepWalk 0.614 / 0.292, scMTNI 0.553 / 0.290, SCRIP 0.638 / 0.644, DeepTFni 0.871 / 0.721, GENELink 0.900 / 0.758.

Single-modality ablations of scMultiomeGRN (lymphoid, Fig. 3A): ATAC only 0.79 / 0.66, scRNA only 0.38 / 0.44, fused 0.93 / 0.86.

Hematopoiesis (Supp. Fig. S1, AUROC / AUPR):

| Cell type | scMultiomeGRN | DeepTFni | GENELink |
|---|---|---|---|
| CMP | 0.92 / 0.66 | 0.90 / 0.65 | 0.87 / 0.50 |
| GMP | 0.93 / 0.71 | 0.91 / 0.70 | 0.91 / 0.64 |
| HSC | 0.92 / 0.67 | 0.89 / 0.60 | 0.90 / 0.58 |

Official demo (PBMC B_naive, fold 1, 87 nodes / 946 edges): best-checkpoint `test` AUROC 0.811, AUPR 0.814, acc 0.745; `all` AUROC 0.836, AUPR 0.471.

## 5. Design decisions for MeVD-GRN (pre-registered before any results)

- **Single evidence tier.** The ground truth is one motif-in-accessible-promoter network, which is binding-type evidence, closest to SC-MO-GRN-DB's localization tier. There are no perturbation or dual tiers. We refuse to invent pseudo-tiers (for example by motif p-value strength), because that would be manufacturing evidence. The curriculum therefore degenerates to one stage trained on that tier (≡ `all_at_once` with one tier; hard negatives are undefined for tier 0).
- **Same splits, same metrics:** the 10 folds are regenerated exactly with their algorithm (seed 666), and their `metric_fn` is reproduced. Both directions of a pair go to the same fold. Scores are symmetrised, (σ(i→j) + σ(j→i))/2.
- **No leakage paths:**
  - (a) Negatives for training are drawn only from the fold's train negatives.
  - (b) The TF-candidate graph is built **label-free** (`exclude_positives=False`). The standard pipeline excludes *all* positives, including test ones, and on a ~300-node TF graph that exclusion would encode the test labels.
  - (c) MeVD-GRN's JASPAR motif relation is **disabled**, because the ground truth is itself a motif scan.
  - Model selection uses val only. Test is touched once.
- **Information parity:** scMultiomeGRN message-passes over the fold's *training* adjacency. Variant `obs` gives MeVD-GRN the same training adjacency as an extra relation (the 3rd relation slot). Variant `noobs` uses MeVD-GRN's label-free graphs only. Both are reported.
- **Budget:** matched to theirs, max 2000 epochs with early stopping at 100 epochs of no val improvement (MeVD-GRN's native val-AUPR criterion).
- **Headline config** = the project's recommended config (FM + h384/l2) with `obs`. The other variants (`noobs`, base h128) are reported alongside. Headline numbers use 5 model seeds (42-46) × 10 folds.
- **Reproduction anchor:** the official scMultiomeGRN code is re-run on our regenerated ground truth and splits. This checks that our regeneration matches their reported numbers, independently of MeVD-GRN.

## 6. Build log

- 2026-09-29: protocol and targets extracted (sections 1-4). Ada down ("major upgrade, back 1 Oct"; `Permission denied (hostbased)`).
- 2026-09-30: session resumed after reboot (scratch wiped). Code written (files below). Local smoke test passed (section 7). Ada still down at 03:12 IST (`Permission denied (hostbased)`), so **no jobs have been submitted**. Everything is ready for the launch in section 8.

Files (all new unless marked):
- `src/data/preprocessing.py` (**edited**, backward compatible): `load_mtx_dir()`; `load_features_by_cells` also accepts a MatrixMarket dir or `matrix.mtx`.
- `src/benchmarks/scmgrn_groundtruth.py`: ground-truth regeneration (peak filter, FASTA, FIMO, TFBS, promoter-containment graph, cross-graph).
- `src/benchmarks/scmgrn_features.py`: MAESTRO v1.2.1 Enhanced RP port; dask-free GRNBoost2 via `arboreto.core.infer_partial_network`.
- `src/benchmarks/scmgrn_protocol.py`: exact split replica, fold loader, official `metric_fn`, test/all evaluation.
- `configs/scmgrn/lung.yaml` (dataset, protocol, and the targets above), `pbmc_demo.yaml` (smoke), `mevd_base.yaml`, `mevd_fm_h384.yaml`.
- `scripts/14_download_scmultiomegrn_data.sh` (`env|data|fm`, login node), `15_scmgrn_split_seurat.R`, `16_scmgrn_build_benchmark.py`, `17_scmgrn_train_eval.py`, `18_scmgrn_official_repro.py`, `19_scmgrn_compile.py`.
- `slurm/scmgrn_common.sh`, `scmgrn_prep.sh`, `scmgrn_build.sh` (array 3), `scmgrn_work.sh` (array 4, %4), `scmgrn_submit_all.sh`.

## 7. Local smoke test (2026-09-30, RTX 4060 laptop)

**Data used:** the official demo cell type `demo/B_naive` (PBMC, hg38). Peaks were restricted to chr19 + chr22 (the laptop has no room for a whole genome): 1,222 peaks / 3.25 Mb. All 769 HOCOMOCO v11 motifs were scanned with real FIMO 5.4.1. Working files were kept in `~/.cache/mevd_scmgrn/` and deleted afterwards.

| Check | Result |
|---|---|
| Split replica vs official `GraphDataModule.split_dataset` (run from the downloaded repo, stubbed `pytorch_lightning`) | **byte-identical**: 30/30 files on B_naive, and also on a synthetic 60-node graph |
| Metric replica vs official `utils.metric_fn` (random scores with ties) | identical to 6 decimals on all of auroc, aupr, ap, acc, precision, recall, threshold |
| MAESTRO port vs MAESTRO v1.2.1 `ExtractGeneInfo` + `RP_AddExonRemovePromoter` (their code, run on the 12,079 demo peaks) | max abs diff **0.0** on 12 TFs. The port runs 3× faster. |
| GRNBoost2 | arboreto's own dask-based `grnboost2()` **fails** on current dask ("Must supply at least one delayed object"). Our per-target path calls the same `infer_partial_network` (same SGBM kwargs, early-stop window 25, seed 0). |
| `gt` step | 72,420 TFBS (p ≤ 1e-6); 461 promoter TFs; graph 16 → 13 nodes after the feature intersection; 16 undirected edges. 1.2 min. |
| `official`, `splits`, `verify`, `processed`, `fm` steps | all OK (GRNBoost2 13 × 2,011; Geneformer 13/13 matched) |
| `17` MeVD-GRN, base obs/noobs (2 folds) and fm_h384 obs (10 folds), short budgets | runs end to end; ~1 s per fold at this size |
| `18` official scMultiomeGRN code, unmodified, 2 folds × 30 epochs (lightning 2.5.5) | runs. Our recomputation from its saved `best_all_score_matrix` **equals its own logged** `test/best_test/*` and `test/best_all/*` metrics exactly. So MeVD-GRN and the reproduction are scored by the same code. |
| `19` compile | OK |
| R splitter `15` on a mock Seurat-like S4 RDS read without the class definitions | OK |
| `slurm/scmgrn_common.sh` `pull`/`push`/cleanup (local-mount mode) | OK. All shell scripts pass `bash -n`. |

The smoke numbers are meaningless as results: a 13-node graph gives 1-2 test pairs per fold.

## 8. Runbook for Ada (nothing has been run on Ada yet)

Once Ada is back, re-check the post-upgrade environment first:
`sinfo -o "%N %G %f"`, `module avail`, `sacctmgr show assoc user=$USER format=Account,QOS,DefaultQOS`.
The jobs assume the `mevd-grn` conda env at `$HOME/miniconda3`. Override with `CONDA_SH=...` if the path changed.
`slurm/scmgrn_submit_all.sh` auto-detects the 2080 Ti gres type or feature. If it finds neither, it excludes the 1080 Ti nodes (numbers 01-40).

```bash
# 0. (Ada login node) get the code: it is on branch leak-fix-and-benchmarks
cd <ada-repo-path> && git fetch origin && git checkout leak-fix-and-benchmarks && git pull
# 1. (Ada login node, repo root) one-time setup + downloads straight into /share1 (~6.5 GB)
bash scripts/14_download_scmultiomegrn_data.sh env    # conda env scmgrn-tools (MEME 5.4.1, R+Matrix); pip lightning, arboreto
bash scripts/14_download_scmultiomegrn_data.sh data   # $SCMGRN_ROOT=/share1/$USER/mevd_grn/scmgrn
bash scripts/14_download_scmultiomegrn_data.sh fm     # cache Geneformer for the offline compute nodes
# 2. launch everything (8 jobs = medium-QoS submit cap): prep -> build[0-2] -> work[0-3]%4
bash slurm/scmgrn_submit_all.sh
# 3. after the work jobs finish (login node, reads /share1 directly)
python scripts/19_scmgrn_compile.py --config configs/scmgrn/lung.yaml --root /share1/$USER/mevd_grn/scmgrn
```

What runs:
- **prep** (1 × 10 CPU, 30 GB): R splits the RDS into `sep_data/lung/<cell type>`. If the RDS is too big for 30 GB, resubmit with `--gres=gpu:2 -c 20` (still 1:10).
- **build** (3 tasks, 10 CPU each): every 3rd cell type. Runs FIMO (≈ 25 min per 30 Mb of peaks), MAESTRO RP, GRNBoost2 (the slowest part for Stromal, 43k cells), the splits, **the byte-identity check against the official split code**, and the processed dir.
- **work** (4 × 1 GPU, `%4`): 9 cell types × (official reproduction × 10 folds + 3 MeVD-GRN variants × 5 seeds × 10 folds) = 1,440 fold runs. Results are pushed to `/share1/.../results/lung/<cell type>/<method>/seed*/fold*.json` after every unit. Resubmitting skips finished folds.

Why independent single-GPU runs rather than DDP: each fold is a 50-500-node graph that trains in seconds to minutes, so DDP would only add communication overhead. A single 4-task array also stays inside the 8-submitted-jobs cap, because every array task counts toward it.

Sanity checks once jobs start:
- `logs/smg_build_*`: "cross graph: N nodes, E undirected edges (adj sum X; paper 'Initial' = Y)" should be of the same order as Table S1, and "split verification ... IDENTICAL".
- `logs/smg_work_*`: val AUPR rising over the first checks, and the official folds printing AUROC ≈ 0.9 on "all".

Results are **not** tracked by git as-is: `.gitignore` ignores `results/*`. To track the compiled summary, copy `summary.{json,md}` into a tracked location.

## 9. Deviations from scMultiomeGRN's protocol, and why

1. **The ground truth is regenerated, not downloaded.** Nothing processed is published. We re-implemented their steps with identical semantics, using the same FIMO binary/version/flags, the same HOCOMOCO v11 files and GENCODE promoters from their zip, and the same thresholds. Their perl `GetSequence.pl` and bedops steps are replaced by equivalent Python: containment implies bedops' ≥20% overlap. Table S1's "Initial" counts are the fidelity check, and the official-code reproduction (script 18) is the performance check.
2. **MAESTRO:** a port of v1.2.1's Enhanced RP model instead of the docker image (no docker on Ada). It was verified exact against their code on the demo. Species GRCh38 is kept, as their code hard-codes it.
3. **GRNBoost2:** a per-target `infer_partial_network` process pool instead of arboreto's dask client, which is broken on current dask. It is the same algorithm and seed. Only the TF × gene column order differs, and their own code has arbitrary set order there too.
4. **If the RDS lacks `var.features`,** GRNBoost2 runs on all genes, which is what their code does without a `var_features.tsv`. It is slower but faithful.
5. **Official code shim:** only `ReduceLROnPlateau(verbose=)` is dropped on torch ≥ 2.7 (a print flag). This is a no-op on Ada's torch 2.1.
6. **MeVD-GRN adaptations (pre-registered, section 5):**
   - single tier
   - symmetrised scores
   - self-pairs dropped from train/val (they are not TF-TF interactions, and their loader labels them 1)
   - model selection on val AUPR (their rule is val loss; both use val only)
   - label-free candidate graph
   - no motif relation
   - optional training-graph relation (`obs`)
7. **Metrics:** the paper's AUROC/AUPR are the transductive full-matrix `all` set, which contains the training edges (section 3). We report both `all` (comparable) and the held-out `test` fold (leakage-free). A win claim must hold on `all` to be comparable. It should also be shown on `test` to be meaningful.
8. **Seeds:** MeVD-GRN uses 5 seeds × 10 folds. The official reproduction uses its default seed 666 × 10 folds. The paper aggregates "ten runs", and its per-fold ± is the fold std.
9. **Hematopoiesis (Supp. Fig. S1) is not built yet (phase 2).** It needs Buenrostro's per-cell BAMs (Dropbox tarball, `bedtools coverage` over `combined.sorted.merged.bed`), and its scRNA source (mmc4.zip) is ambiguous. Once HSC/CMP/GMP exist in the same `<cell type>/{atac,scrna}` layout, the pipeline is dataset-agnostic: copy `lung.yaml` with the hg19 resources and the targets from section 4.

## 10. Findings about the existing SC-MO-GRN-DB pipeline

- **Test negatives are seen in training.** `MEvDTrainer._build_train_edges` samples random negatives from `data.negative_pool`, which is the *full* pool. `create_global_edge_splits` draws the val/test negatives from chunks of that same pool. For K562 localization (~818k train positives × `neg_ratio` 5 > the 1M pool cap), essentially every val/test negative is a label-0 training example every epoch. The fix is to sample from `splits["train"]["neg"]` (or the train chunk) only. This would change the existing reported numbers, so it needs a deliberate rerun. **Fixed 2026-09-30** on branch `leak-fix-and-benchmarks`. `MEvDTrainer.restrict_negative_pool` removes every tier's val/test negatives from the pool before training. Every training entry point (03, 06, 12, 17) calls it, `train_stage` refuses to run without it, and `scripts/00_sanity_check.py` checks both. Measured on the local K562 splits: 100% of every tier's val/test negatives were inside the sampling pool. Setting `training.exclude_eval_negatives: false` reproduces the pre-fix numbers for the before/after comparison. The numbers in results.md and the paper predate the fix.
- **The TF-candidate graph excludes all positives, including val/test.** `02_preprocess.py` calls `build_prior_graph(..., evidence=<all tiers, all splits>, exclude_positives=True)`, so a pair's *absence* from a TF's top-500 candidates is correlated with its being a held-out positive. The effect is weak at 22,943 genes, but it is label-dependent graph construction. The fix is to exclude only train positives, or build it label-free (as done here).
