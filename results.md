# MEvD-GRN Results

Last updated: 2026-09-11. This file replaces the earlier results catalog,
which is preserved for the record at
`docs/archive/results_2026-09-07_pre_bugfix.md`. **Do not compare numbers
across the two files directly** — two real bugs (an EPR metric error and a
hard-negative train/test leakage bug, see Section 1) were fixed in between,
and several numbers moved substantially as a result. See `plan.md` for the
full narrative and rationale behind every change; this file is the numbers
reference.

Cell type: K562 (main), plus Macrophage and MCF7 (multi-cell-type / transfer,
Section 5). Primary metric: AUPR. Also reported: AUROC, early precision
(EP), EPR (early precision ratio vs. random).

---

## 0. Two bugs fixed before any of the numbers below

1. **EPR metric bug**: early precision was ranked over a small curated eval
   set but normalized against a genome-wide TF×gene candidate count,
   inflating EPR into the hundreds (previously up to ~450). Fixed to
   normalize against the actual eval set — EPR is now a sane 1-5x.
2. **Hard-negative leakage bug** (the dominant cause of the previously bad
   results): during Stage 2 training, "hard negative" candidates were drawn
   from a lower tier's *full* evidence (train+val+test), not just its train
   split — so the model was trained to score down exactly the edges later
   used to evaluate that tier. Fixed to restrict hard negatives to the lower
   tier's train split only. This alone took localization AUROC after Stage 2
   from **0.33 (below random)** to **0.51 (mild, ordinary forgetting)**, and
   zero-shot dual-evidence AUPR from 0.558 to 0.870.

See `src/training/sampler.py` and `src/evaluation/metrics.py` for the fixes.

---

## 1. Main model (K562, current architecture, no foundation-model embedding)

Script: `scripts/03_train.py --config configs/k562.yaml`
Protocol: sequential curriculum (localization → perturbation), memory
replay on, dual-evidence held out as zero-shot. `curriculum.protocol:
sequential` is the default — see Section 3 for why `all_at_once` is NOT the
default despite structurally avoiding forgetting.
Params: 315,015 (RNA/ATAC encoders + 2×2-layer relational GraphSAGE +
role-aware decoder; see `docs/figures/architecture_diagram_v2.png`)

| Tier | Role | AUPR | AUROC | EP | EPR |
|---|---|---:|---:|---:|---:|
| Localization | trained, then mildly forgotten | 0.573 | 0.515 | 0.547 | 1.02 |
| Perturbation | trained | 0.641 | 0.899 | 0.588 | 3.15 |
| Dual-evidence | **zero-shot** (never trained on) | **0.881** | **0.971** | 0.792 | 4.75 |

Source: `results/K562_results.json`.

---

## 2. What actually helped: architecture ablations (K562, current pipeline)

All full training budget (Stage1 30ep + Stage2 15ep) unless noted. Reported
metric is dual-evidence AUPR/AUROC (the zero-shot generalization test).

| Ablation | What it removes/changes | Dual AUPR | Dual AUROC |
|---|---|---:|---:|
| **Full model** | (nothing — reference) | **0.881** | **0.971** |
| `no_gnn` | no GNN message passing at all | 0.295 | 0.720 |
| `coexpr_only` | GNN sees only the co-expression graph | 0.814 | 0.953 |
| `tf_cand_only` | GNN sees only the TF-candidate graph | 0.785 | 0.933 |
| `rna_only` | ATAC features zeroed | 0.856 | 0.962 |

**Takeaways:**
- **GNN message passing is by far the most important component** — removing
  it entirely drops dual AUPR from 0.881 to 0.295. This has been true since
  before today's fixes and remains true after.
- **The ATAC fix (Section 4) made the TF-candidate graph genuinely useful
  for the first time.** Before the fix, `tf_cand_only` scored close to
  `no_gnn` (~0.26, since the graph's accessibility gate passed ~90% of
  genes and barely filtered anything). After the fix, `tf_cand_only` (0.785)
  is close to `coexpr_only` (0.814) — both individual graphs now carry
  substantial, comparable signal, and combining them (0.881) still helps but
  by a smaller margin than before, because neither one is a near-empty prior
  anymore.
- **ATAC now measurably helps** (Section 4): `rna_only` (0.856) is clearly
  below the full model (0.881) on every single tier, where before the fix
  there was almost no gap at all.

Ablations not yet rerun on the current (post-fix) pipeline —
`gated_fusion`, `concat_fusion`, `edge_mlp`, `with_replay`, `dual_only`,
`pert_only` — are intentionally omitted here rather than presented as
current. Their pre-fix numbers are in
`docs/archive/results_2026-09-07_pre_bugfix.md` for reference only.

---

## 3. Training protocol: sequential vs. all-at-once (K562, current pipeline)

An earlier (pre-bugfix) comparison found `all_at_once` (joint training on
both tiers at once, no stage boundary) beat the sequential curriculum on
every metric. That comparison used the buggy code, which specifically
crippled sequential training. Re-measured on the fixed pipeline, full
budget:

| Protocol | Loc AUPR | Loc AUROC | Pert AUPR | Pert AUROC | Dual AUPR | Dual AUROC |
|---|---:|---:|---:|---:|---:|---:|
| **Sequential (default)** | 0.573 | 0.515 | **0.641** | **0.899** | **0.881** | **0.971** |
| all_at_once | **0.937** | **0.929** | 0.351 | 0.741 | 0.791 | 0.944 |

`all_at_once` wins localization by a large margin and structurally avoids
forgetting, but loses on perturbation and on the headline zero-shot metric.
Working theory: joint training merges localization's ~818k positives with
perturbation's ~177k, diluting perturbation's signal roughly 4.6x; giving
perturbation a dedicated fine-tuning stage (as sequential does) transfers
better to the nested dual-evidence tier. **Sequential remains the default**;
`all_at_once` is kept available (`scripts/06_ablation.py --ablation
all_at_once`) as a legitimate, forgetting-free comparison point.

---

## 4. ATAC pipeline fix (root-cause diagnosis + before/after)

**Diagnosis.** The old pipeline linked peaks to genes via a binary ±100kb
window (no distance weighting), collapsed to a 2-dim `[mean, var]` that
turned out to be 96%+ correlated (effectively one number, not two), plus a
separate `openness` scalar that was a coarse, heavily-quantized
re-derivation of the same quantity (77 unique values across 22,943 genes,
64% correlated with the mean) — thresholded at `>0` for the TF-candidate
graph's accessibility gate, which passed ~90% of genes and barely filtered
anything.

**Fix.** Replaced with an exponential TSS-distance decay (regulatory
potential, MAESTRO/BETA-style), a proper 3-dim `[mean, var, detection]` of
RP-weighted gene activity, and a 4-dim locus-shape descriptor
(`mean_topK_RP`, `max_RP`, `mean_topK_distance`, `peak_count`) replacing the
collapsed scalar. The TF-candidate graph's accessibility gate now uses
`max_RP > 0.1` (~23kb proximal) instead of "any peak within 100kb."
`ATACEncoder` also got a real hidden layer (was a single Linear).

**Measured improvement in information content** (K562, real data):

| | Before | After |
|---|---|---|
| `openness` unique values / 22,943 genes | 77 | 20,694 / 6,875 / 20,205 (3 richest columns) |
| corr(openness, gene-activity mean) | 0.64 | 0.53 |
| Genes passing the TF-candidate graph's accessibility gate | ~90% (barely filters) | 84.4% at the new, more meaningful `max_RP>0.1` threshold |

**Measured improvement in what the model actually does with it**: see
Section 2 — `rna_only` now clearly underperforms the full model, and
`tf_cand_only` went from near-`no_gnn` to close-to-`coexpr_only`.

---

## 5. Foundation-model gene embeddings (Geneformer) — the single biggest lever found

**What we did.** Extracted Geneformer's (`ctheodoris/Geneformer`,
`Geneformer-V2-104M`, 768-dim) pretrained per-gene input-embedding table —
no fine-tuning, pure lookup — and added it as an optional additional input
channel (`use_fm: true`), fused additively into the RNA branch before the
GNN. Coverage: 77.0% of K562's 22,943 genes, 89.9% of Macrophage's, 89.6% of
MCF7's.

**Within-K562 result** (full budget):

| | Dual AUPR | Dual AUROC |
|---|---:|---:|
| Hand-crafted features only (Section 1 baseline) | 0.881 | 0.971 |
| `fm_only` (Geneformer alone, hand-crafted RNA zeroed) | **0.951** | **0.989** |
| `with_fm` (Geneformer + hand-crafted combined) | 0.949 | 0.988 |

The pretrained embedding *alone* slightly exceeds combining it with our
hand-crafted features — once it's present, our hand-built RNA/co-expression
features add essentially nothing. This is the single largest improvement
found in the whole project.

**Cross-cell-type transfer result — an important, literature-contradicting
finding.** Zero-shot transfer (train on K562 only, test with no retraining
at all) to two independent target cell types:

| Source → target | AUPR (no FM) | AUPR (with FM) |
|---|---:|---:|
| K562 → Macrophage | 0.607 | **0.447** |
| K562 → MCF7 | 0.758 | **0.597** |

The FM embedding — the single biggest in-domain win — **makes single-source
cross-cell-type transfer worse**, consistently, on two unrelated target cell
types. Working explanation: the small network reading the embedding
(`FMEncoder`) and everything downstream of it calibrates to whichever one
cell type it was trained on. See Section 6 for the fix.

---

## 6. Multi-cell-type generalization

**Reference ceiling** (train-from-scratch directly on the target cell
type's own data — "how good can a model get on this cell type at all"):

| Cell type | genes | TFs | localization edges | From-scratch AUPR | AUROC |
|---|---:|---:|---:|---:|---:|
| Macrophage | 16,202 | 22 | 112,383 | 0.841 | 0.915 |
| MCF7 | 16,730 | 250 | 1,209,558 | 0.933 | 0.921 |

**Joint multi-cell-type training** (new `scripts/12_joint_train.py`: one
shared model, one gradient step per cell type per epoch, K562 + MCF7
trained together, then evaluated zero-shot on Macrophage — a THIRD cell
type never seen in that run at all):

| Setup | Macrophage AUPR (zero-shot) | AUROC |
|---|---:|---:|
| K562 alone, no FM → Macrophage | 0.607 | 0.692 |
| K562 alone + FM → Macrophage | 0.447 | 0.648 |
| **K562 + MCF7 joint, no FM → Macrophage** | **0.743** | **0.835** |
| K562 + MCF7 joint + FM → Macrophage | 0.703 | 0.782 |

**Two honest findings, not one:**
1. **Joint training on more than one cell type is a big, unambiguous win for
   transfer, with or without the FM embedding** (no-FM: 0.607→0.743;
   with-FM: 0.447→0.703) — the strongest lever found for the transfer
   problem specifically, and reaches within 0.10 AUPR of Macrophage's own
   from-scratch ceiling (0.841) despite never training on Macrophage at all.
2. Joint training substantially shrinks the FM transfer penalty (gap
   narrows from -0.160 AUPR at single-cell-type to -0.040 AUPR jointly
   trained) but does **not fully reverse it** in this run — no-FM still
   edges out FM on this specific transfer number. In-domain performance on
   both training cell types also improved under joint training regardless
   of FM (e.g. no-FM joint K562 test AUPR 0.940 vs. 0.881 single-cell-type).

**This cross-cell-type, zero-shot-to-a-third-cell-type result is the
strongest, most complete finding in the project** and — as far as we can
tell — nobody else has run this kind of multi-cell-type experiment on
SC-MO-GRN-DB.

Sources: `results/joint_training/K562_MCF7_nofm_holdout_Macrophage.json`,
`results/transfer/*.json`.

---

## 7. Baselines (K562, same splits)

Script: `scripts/05_run_baselines.py`. AUPR/AUROC unaffected by the EPR fix
(only EPR itself changed); these are current.

| Method | Loc AUPR | Loc AUROC | Pert AUPR | Pert AUROC | Dual AUPR | Dual AUROC |
|---|---:|---:|---:|---:|---:|---:|
| **MEvD-GRN (Section 1)** | 0.573 | 0.515 | 0.641 | 0.899 | **0.881** | **0.971** |
| **MEvD-GRN + Geneformer (Section 5)** | — | — | — | — | **0.951** | **0.989** |
| GRNBoost2 | 0.535 | 0.499 | 0.206 | 0.550 | 0.173 | 0.525 |
| RegDiffusion | 0.546 | 0.515 | 0.185 | 0.493 | 0.165 | 0.490 |
| GMF-GAE | 0.526 | 0.514 | 0.218 | 0.565 | 0.160 | 0.516 |
| scMultiomeGRN (adapted) | *training in progress* | | | | | |

**scMultiomeGRN status**: this baseline is a reimplementation of Xu et al.
(NAR 2025), adapted to MEvD-GRN's TF→all-gene benchmark (the published
method is TF-TF only) — see `src/baselines/scmultiomegrn_wrapper.py`'s
docstring for the full list of deliberate adaptations (own features, not
MAESTRO RP/GRNBoost2/FIMO; own train/val/test split, not the paper's
tenfold+consensus protocol). Treat any comparison against it as
"architecture-inspired baseline on our benchmark," not "reproducing the
paper's reported numbers." Currently training on Ada (checkpointed every
improvement as of today — see `plan.md`); last known checkpoint: epoch 160,
localization val AUPR 0.853, still climbing. Will be updated here once it
finishes.

The RNA-only baselines (GRNBoost2, RegDiffusion, GMF-GAE) all cluster near
random on perturbation/dual (AUPR 0.16-0.22, AUROC 0.49-0.57) — MEvD-GRN's
multi-omic graph structure clearly matters relative to expression-only
methods.

---

## 8. Where to look for more detail

- `plan.md` — the living roadmap: why each change was made, what's still
  open, effort/risk estimates.
- `docs/figures/architecture_diagram_v2.png` — current architecture diagram.
- `docs/citations.md` — design rationale + literature citations per
  component.
- `docs/archive/` — pre-2026-09-11 results, presentations, and the original
  full theory spec, kept for historical reference. Numbers there predate
  the bugfixes in Section 0 and should not be treated as current.
- `paper/main.tex` — draft manuscript (needs updating with the numbers in
  this file — it still has the pre-bugfix headline numbers as of writing).
