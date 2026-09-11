# MEvD-GRN — Implementation Progress (Session 1)

> Handoff note for the next LLM/session. This captures everything done so far so you
> don't lose context. Date: 2026-07-01/02. Author: previous Claude session.

---

## 0. TL;DR

MEvD-GRN is a compact (~120K param) **dual-modality GNN** that predicts directed
**TF → target** regulatory edges from single-cell RNA + ATAC, trained via a
**3-stage evidence-aware curriculum** (localization → perturbation → dual-evidence)
over the SC-MO-GRN-DB reference networks. The full pipeline is **implemented and has
been run end-to-end locally on an RTX 4060** for two cell types (K562 human, ESC mouse),
including baselines, a 9-variant ablation study, and cross-cell-type transfer.
Status = working proof of concept. Headline: K562 dual-evidence AUPR ≈ 0.88, ESC ≈ 0.73,
vs unsupervised baselines ≈ 0.17.

Plans live in `MEvD_GRN_plan_v2.md` (canonical spec) and `MEvD_GRN_plan.md` (v1 w/ code).
Hardware ref: `ada.md` (IIIT Ada cluster — NOT used yet; everything ran locally).

---

## 1. Environment (IMPORTANT — read before running anything)

- **Use this Python:** `/home/vishak/miniforge3/bin/python` (miniforge base env).
  It has torch 2.9.1+cu128 and torch_geometric 2.8.0. Plain `python3` on PATH may
  resolve to linuxbrew Python which has **no torch_geometric** — always call the
  miniforge python explicitly, or `conda activate base`.
- GPU: NVIDIA RTX 4060 (8 GB). CUDA works. Model uses ~1.6 GB; **regdiffusion baseline
  OOMs** for >~14k genes on 8 GB.
- `scanpy`/`anndata` are NOT installed — preprocessing deliberately avoids them
  (custom loaders). `lightgbm` present; `arboreto`/`regdiffusion` installed but
  arboreto is broken with the modern dask (see §8). `h5py`, `sklearn`, `pandas`,
  `pyranges`(not needed) available.
- `batch_size: 8192` in configs (default was 2048; raised because the model re-encodes
  the whole graph every minibatch, so fewer/larger batches = big speedup at ~no VRAM cost).
- Disk was tight at times (~5–25 GB free on /home). Datasets were unzipped then zips deleted.

---

## 2. Data downloaded (SC-MO-GRN-DB)

Per-file URLs: networks `https://scmogrndb.psu.edu/RN_TSV/RN###.tsv`,
datasets `https://scmogrndb.psu.edu/DS_ZIP/DS###.zip`. Full single-cell ZIPs are
17.5 GB (human) / 7.2 GB (mouse) — DO NOT bulk download. Annotation tables:
`/Reference_Networks_Annotation.tsv`, `/Datasets_Annotation.tsv`.

**Only two cell types have all 3 evidence tiers** (required for the curriculum):

| Tier | K562 (human) | ESC (mouse) |
|---|---|---|
| localization | RN117 ChIP-seq | RN114 ChIP-chip (ChIPX) |
| perturbation | RN118 KO | RN115 LOGOF |
| dual-evidence | RN119 ChIPseq&KO | RN116 ChIPchip&LOGOF |

Networks live in `data/raw/networks/{K562,ESC}_{localization,perturbation,dual_evidence}.tsv`.
Network file format: 3-col TSV WITH header `Source\tTarget\tRelationship`, gene symbols UPPERCASE.

**Single-cell datasets used** (in `data/raw/single_cell/`):
- K562 RNA: **DS019** (`K562_scRNA/`) — `.out` file, **comma-separated** genes×cells CSV, 953 cells.
- K562 ATAC: **DS025** (`K562_multiome_h5/GSM5396329_..._raw_feature_bc_matrix.h5`) — a 10x
  CellRanger multiome h5 (Gene Expression + Peaks; we use Peaks, hg38, ~177k human peaks).
  (The 485 MB fragments file was skipped; only the .h5 was extracted.)
- ESC RNA+ATAC: **DS010** (`ESC_multiome/DS010_..._MOUSE_EB/`) — `scRNA_Expression.EB.tsv`
  (tab; row labels `coord|strand|ensembl|SYMBOL`) and `scATAC_PeakMatrix.EB.txt` (tab; `chr:start-end`).

**GTFs** (for ATAC peak→gene mapping) in `data/raw/gtf/`: `gencode.v38.gtf.gz` (hg38, K562),
`gencode.vM25.gtf.gz` (mm10, ESC — CONFIRMED correct build via 96.6% ATAC coverage),
`gencode.vM27.gtf.gz` (mm39, unused). Loader reads `.gz` directly.

### Data-format gotchas discovered (all handled in code)
- **Formats vary per dataset**: tab `.txt`, tab `.tsv`, comma `.csv/.out`, and 10x `.h5`.
  `load_features_by_cells` sniffs delimiter; `load_10x_h5` handles the CellRanger matrix.
- **Matrices are features×cells** (transposed vs scanpy) — loaders transpose.
- **Gene-symbol case & label formats differ**: networks UPPERCASE, RNA mixed-case, DS010 uses
  `chr..|strand|ENSG|Symbol`. `extract_symbol()` canonicalizes (takes last `|` field, uppercases).
- **10x peaks are barnyard** (`hg38.` / `mm10.` prefixes) — we filter to the target species prefix.

---

## 3. Repo structure (all implemented)

```
configs/           default.yaml (+ k562.yaml, esc.yaml inherit it; *_localtest/fasttest are smoke tests)
src/data/          preprocessing.py (loaders, QC, gene-activity, 10x h5, GTF TSS),
                   graph_builder.py (universe, prior graph, evidence edges, negatives, splits, nesting),
                   dataset.py (CellTypeData container)
src/models/        encoders.py, fusion.py (GatedFusion+ConcatFusion), gnn.py (GraphSAGE),
                   decoder.py (bilinear), mevd_grn.py (assembles; ablation flags use_atac/use_gnn/fusion)
src/training/      curriculum.py (stages), sampler.py (random+hard negatives),
                   losses.py (weighted BCE, focal), trainer.py (MEvDTrainer)
src/evaluation/    metrics.py (AUROC/AUPR/EarlyPrec/EPR), benchmarker.py (compile + scorer eval)
src/baselines/     common.py (expr loader + scorers), grnboost2_wrapper.py (LightGBM reimpl!),
                   regdiffusion_wrapper.py, gmfgrn_wrapper.py (run_gmf_gae self-contained GAE)
scripts/           00_sanity_check.py, 01_download_data.sh, 02_preprocess.py, 03_train.py,
                   04_evaluate.py (incl. transfer), 05_run_baselines.py, 06_ablation.py,
                   07_compile_results.py
slurm/             train.sh, ablation.sh, baselines.sh (for Ada; not used locally yet)
mevd_grn_architecture.html / .png   presentation diagram
results/           K562_results.json, ESC_results.json, *_transfer.json,
                   ablations/*_K562.json (9), baselines/*_{K562,ESC}.json, summary_table.csv, figures/*.png
data/processed/{K562,ESC}/  rna_features_aligned.npy, atac_features_aligned.npy, gene_index.json,
                   tf_indices.json, prior_edges.pt, negative_pool.pt, evidence_*.pt, summary.json
data/splits/       {cell}_{tier}_splits.pt
```

---

## 4. How to run (local)

```bash
PY=/home/vishak/miniforge3/bin/python
$PY scripts/00_sanity_check.py                                   # synthetic end-to-end check (no data/GPU)
$PY scripts/02_preprocess.py --config configs/k562.yaml          # + esc.yaml
$PY scripts/03_train.py      --config configs/k562.yaml --device cuda:0
$PY scripts/04_evaluate.py   --config configs/k562.yaml \
     --checkpoint results/checkpoints/K562/final_model.pt --transfer_config configs/esc.yaml
$PY scripts/06_ablation.py   --config configs/k562.yaml --ablation full_curriculum --device cuda:0
$PY scripts/05_run_baselines.py --config configs/k562.yaml --baselines gmf_gae,grnboost2 --device cuda:0
$PY scripts/07_compile_results.py                                # -> results/summary_table.csv + figures
```
**Do not run training and baselines concurrently** — both use the 8 GB GPU and will OOM.

---

## 5. Model & method (what's implemented)

- **Framing:** directed bipartite link prediction. Nodes = genes; per-gene features =
  [RNA mean, RNA var, ATAC mean, ATAC var]. ATAC = gene-activity score (sum peaks within
  ±100 kb of TSS, then mean/var across cells).
- **Architecture:** RNA encoder (2-layer MLP) + ATAC encoder (linear) → **gated fusion** →
  **2-layer GraphSAGE** over a prior candidate graph (bidirectional message passing) →
  **bilinear decoder** score(TF,tg)=σ(h_TF·W·h_tg). ~115–125K params.
- **Prior graph:** per TF, top-K=500 most-accessible target genes; EXCLUDES known positives
  (so val/test labels can't leak via message passing). It's a message-passing structure only —
  the decoder scores any pair.
- **Negatives:** random TF×gene pairs not positive in ANY tier (field standard). NOT the prior —
  restricting negatives to the accessibility-ranked prior inflates AUPR via an artifact (fixed).
- **Curriculum:** Stage1 localization (lr 1e-3, 30 ep) → Stage2 perturbation (3e-4, 15 ep) →
  Stage3 dual-evidence (1e-4, 10 ep, encoders frozen). Weighted BCE (pos_weight = actual neg/pos),
  hard negatives from adjacent tiers, AdamW + cosine LR, early stopping on val AUPR.
- **Eval:** per-tier held-out test (70/15/15 edge-level split) + random negatives.
  AUPR (primary), AUROC, Early Precision, EPR. `dual_evidence_mode: reconstruct` can rebuild
  dual = loc∩pert for a clean-nested variant.

---

## 6. Results obtained (preliminary, 2 cell types)

**MEvD-GRN (RNA+ATAC), per tier:**

| Cell | Tier | AUPR | AUROC | EP | EPR |
|---|---|---|---|---|---|
| K562 | localization | 0.710 | 0.639 | 0.635 | 17.3 |
| K562 | perturbation | 0.337 | 0.706 | 0.274 | 50.4 |
| K562 | **dual-evidence** | **0.878** | **0.976** | 0.793 | 949.8 |
| ESC | localization | 0.445 | 0.799 | 0.519 | 47.8 |
| ESC | perturbation | 0.358 | 0.666 | 0.389 | 39.6 |
| ESC | **dual-evidence** | **0.728** | **0.938** | 0.731 | 700.5 |

**Baselines (dual-evidence AUPR / AUROC):** all unsupervised, ~random on dual tier.
- K562: GRNBoost2 0.175/0.532, gmf_gae 0.167/0.530. RegDiffusion FAILED (data edge-case).
- ESC:  gmf_gae 0.177/0.523, GRNBoost2 0.171/0.529. RegDiffusion FAILED (OOM on 8 GB).
- (Baselines do better on the dense localization tier, e.g. K562 gmf_gae loc AUPR 0.541,
  because that tier is easy; they collapse on the sparse high-confidence tiers.)

**Ablations (K562, RNA+ATAC, dual-evidence):**

| Ablation | AUPR | AUROC | note |
|---|---|---|---|
| concat_fusion | 0.902 | 0.982 | ≈ full (gating not clearly better — honest caveat) |
| full_curriculum | 0.890 | 0.978 | the method |
| rna_only | 0.740 | 0.939 | **ATAC contributes +0.15 AUPR** |
| all_at_once | 0.709 | 0.880 | **curriculum > multi-task (+0.18)** |
| with_replay | 0.639 | 0.871 | replay doesn't help |
| loc_only | 0.638 | 0.809 | single tier worse |
| no_gnn | 0.634 | 0.929 | **GNN critical (big AUPR drop)** |
| dual_only | 0.623 | 0.851 | single tier worse |
| pert_only | 0.165 | 0.506 | fails without Stage-1 init |

(Main K562 dual 0.878 vs full_curriculum ablation 0.890 differ slightly — separate runs,
different negative sampling; both are "full curriculum + ATAC".)

**Cross-cell-type transfer (K562 model → ESC, no retraining):**
loc 0.282/0.701, pert 0.173/0.500, dual 0.317/0.770. Weak (cross-SPECIES, gene-identity-agnostic
features) but the ATAC model transfers better than the RNA-only one did (dual AUROC 0.59→0.77).

---

## 7. Key scientific decisions & findings

1. **Evidence nesting is verified, not assumed.** ESC: dual⊆loc and dual⊆pert = 100% (RN116 is
   the exact intersection), but pert⊆loc ≈ 11% (knockouts hit indirectly-regulated genes).
   K562: dual⊆loc ≈ 77% (RN119 is from a DIFFERENT study than RN117/RN118) — not strictly nested.
2. **Two bugs fixed early** (both would have faked good numbers): (a) negatives must be random
   TF×gene, not prior-restricted (else expression-artifact shortcut → fake AUPR≈0.99);
   (b) prior graph must exclude positives (else message-passing leaks test labels).
3. **ATAC helps** — K562 dual 0.767 (RNA-only) → 0.878 (with ATAC); ESC full 0.728 vs rna_only 0.708.
4. **Curriculum > multi-task > single-tier**, confirmed on real data. This is the paper's core claim.
5. **Baselines are unsupervised**; MEvD-GRN is supervised — the comparison shows the value of
   supervised evidence-curriculum learning, not a like-for-like architecture race. State this honestly.
6. Speed: fast because we train over ~20k GENE rows (pre-aggregated features), not millions of cells.

---

## 8. Benchmarking & baseline strategy (planned — not yet fully built)

### 8.1 Candidate baseline suite (4 families)
A broad suite was proposed (~14 methods); a lean, reviewer-friendly **8** is recommended
(most-cited/best-maintained per family, covering both supervision regimes):
- **Classical statistical (unsupervised):** GENIE3, GRNBoost2, PCC, PIDC — near-zero cost, universally expected.
- **scRNA GNNs:** GENELink, GNNLink, GNE (**supervised** link prediction) · DeepSEM, SCODE (**unsupervised** inference).
- **ATAC / motif (unsupervised):** DeepTFni, SCRIP.
- **Multiomic / deep generative:** SCENIC+ (unsup.), scMTNI (unsup./prior), LINGER (prior-informed).
- Lean 8 pick: GENIE3, GRNBoost2, PIDC, GENELink, DeepSEM, SCRIP, SCENIC+, LINGER.

### 8.2 The ONE rule that makes it apples-to-apples
Every method is scored by the **same harness on the identical held-out test set**: same positive
test edges, same sampled negatives, same gene universe (RNA ∩ network genes), same TF list, same
metrics (AUPR/AUROC/EP/EPR). This is already built: `evaluate_scorer_on_splits` (benchmarker.py)
+ the scorer-wrapper pattern (`DictScorer`/`DenseMatrixScorer`/`EmbeddingScorer` in
baselines/common.py). Each new baseline needs only a thin wrapper mapping its output → a score for
every candidate (TF,target) edge (0 for unranked pairs). GRNBoost2 + gmf_gae already run this way.

### 8.3 "Same training data" is necessary but NOT sufficient — stratify on TWO axes
**Axis 1 — supervision (the big one):**
- **Unsupervised** (GENIE3, GRNBoost2, PCC, PIDC, DeepSEM, SCODE, SCENIC+, scMTNI, DeepTFni, SCRIP):
  they never see ground-truth edges — they emit a full ranked network from data alone. They get **no
  train split**; we score their ranking on the test edges. This puts them at an inherent disadvantage
  vs supervised methods → **must be stated**; they are *reference points* ("what pure inference
  achieves"), not head-to-head competitors. (This is exactly how BEELINE frames them.)
- **Supervised link prediction** (GENELink, GNNLink, GNE — and MEvD-GRN): train on the **identical
  train edge split**, predict the **identical test split**. This is the *true* apples-to-apples group
  and the comparison reviewers trust most.

**Axis 2 — modality:** comparing full MEvD-GRN (RNA+ATAC) to an RNA-only tool conflates "better
method" with "more data". So report two brackets:
- **RNA-only bracket:** MEvD-GRN `rna_only` (ablation we already have) vs GENIE3/GRNBoost2/PCC/PIDC/
  DeepSEM/SCODE/GENELink/GNNLink/GNE.
- **Multiomic bracket:** MEvD-GRN full vs SCENIC+/LINGER/scMTNI/DeepTFni/SCRIP.

### 8.4 How each family is actually run
- **Unsupervised:** run on the **same input** (RNA ± ATAC per the method's design), restricted to the
  **same gene universe**; take its per-edge scores on the test split; same metrics. No training.
- **Supervised:** train on our train split, predict our test split. They have no evidence-tier notion,
  so train once (on the localization tier, or the union) and evaluate against all three tiers — same
  per-tier scoring everyone gets.

### 8.5 Where the curriculum claim is REALLY proven
External baselines **cannot** isolate the curriculum — none has one. The curriculum's value is proven
by our **ablations** (`full_curriculum` vs `all_at_once` vs `loc/pert/dual_only`), which hold
architecture/data/splits fixed and change only the training schedule. Keep two separate claims in the
paper: (1) *curriculum > multi-task/single-tier* → proven by ablations; (2) *architecture is
competitive with the field* → proven by external baselines.

### 8.6 Is there a BEELINE-style benchmark for our exact problem? — NO
- **BEELINE** standardizes datasets + ground-truth + metrics (AUPRC ratio, EPR) for **unsupervised
  scRNA** GRN inference — it has **no train/test split**. The split-based evaluation comes from the
  **supervised** papers (GENELink/GNNLink/GNE), which each roll **their own** splits on BEELINE data;
  there is **no canonical shared split** even for scRNA.
- For our setting (**supervised + multiomic RNA+ATAC + evidence-tiered**): **no equivalent benchmark
  exists.** Closest are review/benchmarking papers (Loers & Vermeirssen 2024, Brief Bioinform;
  Karamveer & Uzun 2024) — not reusable suites. **SC-MO-GRN-DB is the resource that fills this gap**
  (it's why it was built), but it's *data + tiered ground truth*, **not** a turnkey protocol
  (no prescribed splits/negatives/leaderboard). The evidence-tier framing is brand new.
- **Implication / opportunity:** we **define, freeze, and publish our own protocol** (edge-level
  split seeds, negative-sampling scheme, per-tier eval, gene-universe definition) → framed as "the
  first standardized supervised benchmark protocol on SC-MO-GRN-DB", this is a *contribution*, not a
  weakness. Borrow BEELINE's **metrics/conventions** for legibility; optionally add a BEELINE-style
  cross-check for MEvD-GRN `rna_only` (SC-MO-GRN-DB inherited BEELINE's ground-truth networks) to
  bridge to the established scRNA literature. All baselines run under **our** frozen protocol, since
  no external shared split exists to defer to.

### 8.7 Practical caveats when adding baselines
- **Output formats differ** (edge lists, adjacency, regulons, TF-activity) → one small adapter each →
  score-per-edge. Same 3 scorer types cover all.
- **Compute:** SCENIC+/LINGER/scMTNI are heavy/slow; RegDiffusion + LINGER **OOM on 8 GB** for 20k+
  genes → run these on **Ada**, not the laptop. Fix seeds; report mean ± std for the supervised group.
- **Same gene universe for all** so the AUPR denominator is identical (else fewer candidates = falsely
  better).
- Highest-priority to build first: **GENELink** (the fair supervised RNA-only competitor) and
  **LINGER/SCENIC+** (the fair multiomic competitors). Both slot into `evaluate_scorer_on_splits`.

---

## 9. Known issues / limitations

- **RegDiffusion**: OOMs on 8 GB for large gene sets; also a zero-variance edge-case on K562
  (denoise fix added but didn't fully cover). Works only on a bigger GPU / reduced gene set.
- **arboreto (original GRNBoost2)** is INCOMPATIBLE with modern dask (legacy DataFrame removed,
  "Nanny failed to start"). We reimplemented GRNBoost2 faithfully with **LightGBM** instead.
- **GMFGRN**: we run a self-contained graph-autoencoder stand-in (`run_gmf_gae`); the official
  repo hook (`run_gmfgrn`) is stubbed pending a clone.
- Cross-species transfer is weak (expected; no homolog mapping).
- Only 2 cell types validated so far.

---

## 10. Next steps / TODO

1. **More cell types** — within-species transfer (e.g. K562→GM12878/HepG2/MCF7/H1, all human,
   localization-only) using the existing transfer machinery; run partial curricula where only
   1–2 tiers exist (trainer already skips missing tiers); reconstruct dual=loc∩pert where possible.
2. **External validation** — evaluate against a held-out cell type or an independent resource
   (TRRUST, etc.) rather than only held-out edges of the same networks.
3. **Fair supervised baselines & full benchmark suite** — current baselines are unsupervised;
   see **§8** for the full strategy. Priority: build **GENELink** (fair supervised RNA-only
   competitor) + **LINGER/SCENIC+** (fair multiomic competitors), each as a scorer wrapper into
   `evaluate_scorer_on_splits`. Freeze & publish our benchmark protocol (§8.6).
4. **Run on Ada** for scale (slurm/ scripts ready; research account, 1 GPU, 1:10 GPU:CPU).
5. Figures/paper: `results/figures/*.png` + `summary_table.csv` are the raw material.
6. Consider seeds/repeats for error bars (single run currently).

---

## 11. Where the truth lives
- Config-driven; no hardcoded hyperparameters. All results reproducible via the scripts in §4.
- Memory files (cross-session): `~/.claude/projects/-home-vishak-research/memory/`
  (`mevd-grn-project.md`, `scmogrndb-download.md`).
