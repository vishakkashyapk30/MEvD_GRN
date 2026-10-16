# MeVD-GRN vs scMultiomeGRN on scMultiomeGRN's own benchmark

Advisor meeting, Saturday 17 Oct 2026, 11:00. Period covered: 11 to 17 Oct 2026. Slides are separated by `---`.
Branch `scmgrn-tftf`. Experiment write-up with every deviation: `docs/experiments/scmultiomegrn_generalization.md` (section 11).

> **Reading guide.** Every MeVD-GRN number here is zero-leakage: the held-out fold's links and non-links never enter training,
> and the input graphs are built without labels (checked in code on every fold, recorded in every result file).
> Sections marked **[UPDATE]** are appended as runs finish. A table shows "–" where a run has not finished.

---

## 1. Summary (as of the last update below)

1. **We rebuilt scMultiomeGRN's benchmark and reproduced its published numbers.** The regenerated TF-TF answer keys are within
   5% of the paper's link counts on all 9 lung cell types, our fold files are byte-identical to what their split code writes,
   and their own code, re-run on our data, scores what the paper reports (for example Neuroendocrine 0.920 / 0.820, Myeloid 0.943 / 0.849, HSC 0.927 / 0.69 on the paper's metric).
2. **We changed MeVD-GRN to predict undirected TF-TF links**, the task scMultiomeGRN solves (it used to predict directed TF-to-gene edges),
   and gave it scMultiomeGRN's own inputs.
3. **On the held-out test fold, MeVD-GRN is slightly better than scMultiomeGRN on 4 of the 5 lung cell types finished so far**
   (pooled mean +0.011 AUPR, +0.009 AUROC). The gain is small compared with the fold-to-fold spread, and individually significant for
   two cell types. See section 7.
4. **The ablations say the gain does not come from the scMultiomeGRN-specific inputs or from the training graph.** A variant that never sees
   the training links scores the same as the full model. This fits the main caveat of the benchmark: the answer key is generated from the
   same ATAC data the models read (section 9).
5. **Status:** 5 of 9 lung cell types finished; the other 4 and hematopoiesis (HSC, CMP, GMP) are running. See section 11 for the current table.

---

## 2. The datasets

**Human fetal lung** (Domcke et al., Science 2020; GEO GSM4508936): 72,662 cells, 9 cell types, sci-ATAC-seq3.
The object's "RNA" assay is gene-body accessibility (per the GEO description), not RNA, so both modalities come from ATAC.

| Cell type | Cells | TF nodes | TF-TF links (undirected) | Directed entries (ours / paper "Initial") | Link density |
|---|---|---|---|---|---|
| Myeloid | 1,288 | 46 | 266 | 532 / 542 | 0.26 |
| Lymphoid | 1,860 | 41 | 301 | 602 / 613 | 0.37 |
| Megakaryocytes | 568 | 64 | 459 | 918 / 930 | 0.23 |
| Neuroendocrine | 544 | 69 | 674 | 1,348 / 1,361 | 0.29 |
| Lymphatic endothelial | 662 | 94 | 936 | 1,872 / 1,889 | 0.21 |
| Stromal | 43,289 | 94 | 1,044 | 2,088 / 2,104 | 0.24 |
| Bronchiolar and alveolar epithelial | 18,988 | 115 | 1,346 | 2,692 / 2,712 | 0.21 |
| Ciliated epithelial | 1,111 | 228 | 4,254 | 8,508 / 8,503 | 0.16 |
| Vascular endothelial | 4,352 | 271 | 5,026 | 10,052 / 9,609 | 0.14 |

**Human adult hematopoiesis** (Buenrostro et al., Cell 2018; scATAC from GEO GSE96769, scRNA from the paper's Data S2; scATAC and scRNA are separate experiments, as in the scMultiomeGRN paper).

| Cell type | scATAC cells | scRNA cells | TF nodes | TF-TF links |
|---|---|---|---|---|
| HSC | 347 | 2,268 | 150 | 1,519 |
| CMP | 502 | 4,454 (cmp0+cmp1+cmp2, our reading) | 171 | 2,016 |
| GMP | 402 | 1,096 | 157 | 1,556 |

**The answer key is not experimental.** scMultiomeGRN (following DeepTFni) defines a link between two TFs when a binding-motif hit of one TF
(HOCOMOCO v11, FIMO p ≤ 1e-6) falls inside the promoter (TSS ± 2 kb) of the other, in peaks open in more than 10% of the cell type's cells.
We regenerated it with the same tools and thresholds. The task is therefore to recover a motif-in-promoter rule among TFs.

---

## 3. The two models

**scMultiomeGRN** (Xu et al., NAR 2025)
- Nodes are TFs. Inputs per TF: its accessibility score (MAESTRO regulatory potential, one value per cell) and its GRNBoost2 importance
  towards every gene; per TF pair, a 16 × 16 joint expression histogram.
- Two graph neural networks (one per modality) with attention over neighbours, a cross-modal attention layer that mixes them, and an
  inner-product decoder that rebuilds the adjacency. Trained on the known links with binary cross-entropy (up to 2,000 epochs, early stopping, lr 1e-5).
- Final network = majority vote over 10 runs. We use the official code unmodified.

**MeVD-GRN, as used here**
- A node-level encoder per source (RNA/accessibility summary features, ATAC gene activity, and a frozen Geneformer gene embedding), two relational
  GraphSAGE networks over label-free gene graphs (co-expression kNN, accessibility-based TF-candidate graph; optionally the training links),
  fused into one embedding per TF.
- **New this week:** an undirected pair decoder. The score of (u, v) is built from h_u*h_v, |h_u - h_v|, h_u + h_v, a diagonal bilinear term,
  an expression-correlation term and optional pair features, so it is identical whichever TF comes first (unit-tested).
- Inputs matched to scMultiomeGRN in the headline variant: its MAESTRO and GRNBoost2 features (compressed to 64 numbers each), its pair histogram,
  and the training links as message-passing edges. 5.07 M parameters (headline variant).
- Trains in about 3 seconds per fold; the official code takes about 1 to 5 minutes per fold depending on the cell type.

| Variant | What it is | Parameters |
|---|---|---|
| **H** (headline) | undirected + training links + scMultiomeGRN inputs | 5.07 M |
| A1 | undirected + training links, without the scMultiomeGRN input features | 4.88 M |
| A2 | undirected, label-free graphs only (never sees the training links) | 3.70 M |
| A3 | the old directed model, scores averaged over both directions | 4.44 M |

---

## 4. What "apples to apples" means here

| | scMultiomeGRN | MeVD-GRN |
|---|---|---|
| Data | Domcke lung (9 cell types), Buenrostro hematopoiesis | the same files |
| Answer key | motif-in-promoter graph | the same graph, regenerated with their method |
| Splits | 10 folds, seed 666: test = one fold of links + equal random non-links; val 18% of links | the same files (byte-identical to the official split code) |
| Task | undirected TF-TF links | undirected TF-TF links |
| Scoring code | their `metric_fn` (AUROC, AUPR, accuracy) | the same function |
| Inputs | MAESTRO, GRNBoost2, histograms, training graph | the same four (headline variant) |
| Repeats | official code, 10 folds each: 6 seeds (Myeloid, Lymphoid, Megakaryocytes, Neuroendocrine, Lymphatic), 2 seeds (the rest) | 5 seeds × 10 folds |

**Headline score = the held-out test fold.** The paper's published numbers use a different evaluation (the full N × N matrix, which includes
the links the model trained on). We report both and label the second.

---

## 5. How we kept it fair

- The held-out fold's links and non-links are never in the training data (asserted on every fold).
- Training non-links come only from the fold's training non-links; the best epoch is chosen on the validation split only.
- Input graphs are built without labels; the motif relation is off because the answer key is a motif scan.
- All MeVD-GRN result files record `zero_leak: true`; the compile step aborts otherwise.
- The official scMultiomeGRN code is re-run unmodified on our regenerated data, so any gap to the paper is visible.
- No hyper-parameter was tuned on test results: variants were fixed before any result existed.

---

## 6. Did we rebuild the benchmark correctly? (fidelity checks)

1. **Link counts:** within 5% of the paper's Table S1 on all 9 lung cell types (table in section 2). Largest gap: Vascular endothelial, +4.6%.
2. **Fold files:** identical to the output of the official split code, 30 of 30 files, for all 12 cell types.
3. **Official code, re-run on our data, on the paper's own metric (full matrix), AUROC / AUPR:**

| Cell type | Paper | Our re-run of the official code |
|---|---|---|
| Myeloid | 0.95 / 0.85 | 0.943 / 0.849 |
| Lymphoid | 0.93 / 0.86 | 0.929 / 0.861 |
| Megakaryocytes | 0.92 / 0.78 | 0.917 / 0.788 |
| Neuroendocrine | 0.92 / 0.82 | 0.920 / 0.820 |
| Lymphatic endothelial | 0.93 / 0.81 | 0.915 / 0.781 |
| HSC (hematopoiesis, folds 1 and 2) | 0.92 / 0.67 | 0.927 / 0.686 |

The benchmark therefore reproduces the paper; comparisons made on it are meaningful.

---

## 7. Headline result: held-out test fold (5 of 9 lung cell types so far)

AUROC / AUPR, mean over 10 folds then over seeds.

| Cell type | **MeVD-GRN H** | A1 | A2 | A3 | scMultiomeGRN (official code) |
|---|---|---|---|---|---|
| Lymphatic endothelial | **0.922 / 0.937** | 0.921 / 0.935 | 0.920 / 0.934 | 0.919 / 0.933 | 0.908 / 0.917 |
| Megakaryocytes | **0.905 / 0.923** | 0.907 / 0.925 | 0.902 / 0.923 | 0.901 / 0.920 | 0.898 / 0.915 |
| Myeloid | 0.927 / 0.933 | 0.928 / 0.936 | 0.928 / 0.937 | 0.930 / 0.939 | **0.930 / 0.943** |
| Neuroendocrine | **0.927 / 0.938** | 0.928 / 0.939 | 0.929 / 0.940 | 0.928 / 0.939 | 0.912 / 0.917 |
| Lymphoid | **0.926 / 0.935** | 0.922 / 0.932 | 0.927 / 0.937 | 0.921 / 0.926 | 0.917 / 0.918 |

![held-out test fold](figs_scmgrn/fig_test.png)

Fold-paired comparison, H minus official code (10 folds per cell type, Wilcoxon signed-rank):

| Cell type | AUPR difference | p | AUROC difference | p |
|---|---|---|---|---|
| Lymphatic endothelial | +0.020 | 0.002 | +0.014 | 0.010 |
| Neuroendocrine | +0.021 | 0.002 | +0.015 | 0.004 |
| Lymphoid | +0.017 | 0.105 | +0.009 | 0.375 |
| Megakaryocytes | +0.009 | 0.160 | +0.007 | 0.275 |
| Myeloid | -0.011 | 0.375 | -0.003 | 0.695 |
| **Pooled (50 fold pairs)** | **+0.011** | **0.0001** | **+0.009** | **0.0005** |

How to read this: the gains are small (about one point of AUPR), the fold-to-fold standard deviation is 2 to 5 points, only 2 of 5 cell types are
individually significant, and the pooled p-value treats folds within a cell type as independent although they share a graph, so it is optimistic.

---

## 8. Ablations: what produces the score?

On the held-out fold, the four MeVD-GRN variants are within about 0.005 of each other on every cell type finished.

- **A2** never sees the training links (only label-free graphs and features) and scores the same as H. The headline gain therefore does not come from
  message passing over known links.
- **A1 ≈ H:** the scMultiomeGRN-specific inputs (MAESTRO, GRNBoost2, pair histograms) add nothing measurable on top of our own features.
- **A3 ≈ H:** the old directed model scores the same, so the undirected decoder is needed for the task definition but is not what produces the gain.

![ablation](figs_scmgrn/fig_ablation.png)

Interpretation: a TF's accessibility and sequence-related features alone already predict whether it has a motif-promoter link,
because the answer key is a function of the same ATAC data (section 9).

---

## 9. Same benchmark, the paper's own protocol (full matrix; contains training links)

AUROC / AUPR. **This set includes the links the model was trained on, so it is not a leakage-free measure.** Shown to compare with the published numbers.

| Cell type | MeVD-GRN H | A3 | scMultiomeGRN (official code) | Paper (published) |
|---|---|---|---|---|
| Lymphatic endothelial | 0.939 / 0.850 | 0.942 / 0.862 | 0.915 / 0.781 | 0.93 / 0.81 |
| Megakaryocytes | 0.938 / 0.846 | 0.937 / 0.855 | 0.917 / 0.788 | 0.92 / 0.78 |
| Myeloid | 0.956 / 0.877 | 0.958 / 0.892 | 0.943 / 0.849 | 0.95 / 0.85 |
| Neuroendocrine | 0.951 / 0.893 | 0.952 / 0.900 | 0.920 / 0.820 | 0.92 / 0.82 |
| Lymphoid | 0.951 / 0.902 | 0.950 / 0.898 | 0.929 / 0.861 | 0.93 / 0.86 |

![paper protocol](figs_scmgrn/fig_paper.png)

The pooled difference on this set is +0.054 AUPR (50 pairs). It is larger than on the held-out fold because a higher-capacity model fits the training links better;
it should not be read as better generalisation.

---

## 10. Honest caveats

- **The answer key is built from the same ATAC data the models read** (motif hit in an open promoter). A model can score well by learning that
  rule, without learning regulation. The same holds for scMultiomeGRN; the ablations suggest it holds for us too. This benchmark cannot tell
  "learns regulation" from "learns the construction rule".
- **Lung "RNA" is gene-body accessibility**, not expression, so both modalities come from ATAC.
- **Small graphs:** the five finished cell types have 41 to 94 TFs and 266 to 936 links; each held-out fold has only 27 to 94 positive links.
- **Official code seeds:** 6 seeds for the first five cell types, 2 (666 and 42) for the rest, because the code is near-deterministic across seeds
  (Myeloid, 6 seeds: all within 0.001) and costs about 5 minutes per fold on the large cell types.
- **Deviation from the official preprocessing:** the Seurat object stores no variable-gene list, so GRNBoost2 uses genes detected in more than 10% of
  cells (plus all TFs), following the paper's own filter. Applies to both methods.
- **CMP scRNA** is our reading of three labels (cmp0, cmp1, cmp2); the paper does not say which it used.
- **Not tested here:** whether the result transfers to experimentally derived answer keys (that is what BEAR-GRN offers).

---

## 11. [UPDATE] Progress of the runs

| Dataset / cell type | Status at last update | MeVD-GRN (4 variants × 5 seeds × 10 folds) | Official code |
|---|---|---|---|
| Lung: Myeloid, Lymphoid, Megakaryocytes, Neuroendocrine, Lymphatic | **finished** | done | 6 seeds done |
| Lung: Bronchiolar | running | in progress | in progress (2 seeds) |
| Lung: Ciliated, Vascular | queued | – | – |
| Lung: Stromal (43k cells) | queued last; may hit the memory limit of the official code | – | – |
| Hematopoiesis: HSC, CMP, GMP | running | in progress | in progress (2 seeds) |

*(Further updates are appended below as results arrive.)*

---

## 12. Engineering notes (for the record)

Problems met and fixed while running on real data, all in `docs/experiments/scmultiomegrn_generalization.md` section 11.8:
- Ada's outbound internet is about 10 to 140 KB/s per connection; the 3 GB lung file was streamed from the laptop over the VPN (about 3.5 hours).
- The GEO lung file is gzipped twice; the Seurat object has no variable-gene list.
- Official code ran out of memory on the large cell types when workers were given 6 GB; fixed with 30 GB per worker.
- A bug in the worker queue (ssh swallowing its own input list) made the laptop worker stop early; fixed and re-tested.

---

## Appendix. Where things are

- Code: branch `scmgrn-tftf` (model: `src/models/mevd_grn.py`, `src/models/decoder.py`; features: `src/benchmarks/scmgrn_features.py`;
  scripts 14b/16/17/18/19; jobs `slurm/scmgrn_*.sh`; tests `tests/`).
- Results on Ada: `~/scmgrn/results/{lung,hema}/<cell type>/<method>/seed*/fold*.json`; compiled with `scripts/19_scmgrn_compile.py`.
- Figures: `docs/presentation/make_scmgrn_figs.py` (to `docs/presentation/figs_scmgrn/`).
