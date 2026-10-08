# MEvD-GRN Results

Last updated: 2026-10-08 (dated corrections only: the leak banner below and
correction notes in Sections 2, 3 and 7; no number was changed in place).
Previous update 2026-09-18; first written 2026-09-11. This file replaces the
earlier results catalog,
which is preserved for the record at
`docs/archive/2026-09-07_results_pre_bugfix.md`. **Do not compare numbers
across the two files directly** — two real bugs (an EPR metric error and a
hard-negative train/test leakage bug, see Section 1) were fixed in between,
and several numbers moved substantially as a result. See `plan.md` for the
full narrative and rationale behind every change; this file is the numbers
reference.

> **Correction (2026-10-08): every MEvD-GRN number in this file was trained
> with a negative-sampling leak.** Until the 2026-09-30 fix (commit c26c03d,
> branch `leak-fix-and-benchmarks`), training drew random negatives from the
> whole 1M-pair pool, and that pool held 100% of every tier's val/test
> negatives. So each model was trained, as label-0 examples, on the very
> negatives it was later scored on. The baselines (scMultiomeGRN, GRNBoost2,
> RegDiffusion, GMF-GAE) never had this leak. Full write-up:
> `docs/experiments/leakfix_rerun.md`. Its 5-seed K562 rerun (Geneformer on,
> h384/l2), copied verbatim from that doc:
>
> | Run | n | loc AUPR | loc AUROC | pert AUPR | pert AUROC | dual AUPR | dual AUROC |
> |---|---|---|---|---|---|---|---|
> | all_at_once, **pre-fix** (paper headline) | 5 | 0.9676±0.0002 | 0.9609±0.0002 | 0.6955±0.0017 | 0.9096±0.0004 | 0.9517±0.0007 | 0.9880±0.0002 |
> | all_at_once, **fixed** | 5 | 0.9544±0.0007 | 0.9485±0.0007 | 0.5722±0.0114 | 0.8803±0.0039 | 0.9017±0.0040 | 0.9781±0.0009 |
> | full_curriculum, **fixed** | 5 | 0.7450±0.0040 | 0.7602±0.0032 | **0.8165±0.0006** | **0.9602±0.0001** | **0.9562±0.0008** | **0.9902±0.0002** |
> | scMultiomeGRN baseline (`results/baselines/scmultiomegrn_K562.json`; trains on train-split negatives only, so it was never leaky) | 1 | 0.9059 | 0.9141 | 0.6550 | 0.9185 | 0.8466 | 0.9705 |
> | GRNBoost2 baseline | 1 | 0.5347 | 0.4987 | 0.2060 | 0.5504 | 0.1728 | 0.5245 |
>
> - With the fix switched off (`exclude_eval_negatives: false`), the rerun
>   reproduces the old seed-42 numbers exactly, to 4 decimals, so every
>   difference comes from the fix alone.
> - After the fix, `all_at_once` (the recommendation in Sections 3 and 7)
>   loses to scMultiomeGRN on perturbation. The curriculum (`full_curriculum`,
>   i.e. the sequential protocol) wins perturbation and dual evidence, and
>   loses localization (0.745).
> - **A second leak is not yet fixed in these numbers either.** The
>   TF-candidate graph was built with every known positive excluded,
>   including val/test ones, so 0 val/test positives were graph edges against
>   about 2.8-3.0% of val/test negatives
>   (`docs/experiments/labelfree_graph_rerun.md` §1). The K562 rerun with a
>   label-free graph is running locally as of 2026-10-08; its numbers are
>   pending. The fixed rows above are therefore not final either.
> - Every number below is kept exactly as it was written, for the record.
>   Treat each MEvD-GRN number in this file as pre-fix. The external
>   benchmarks run since then use the negative-pool fix and label-free graphs
>   from the start: PBMC10k vs LINGER
>   (`docs/experiments/pbmc10k_linger_benchmark.md` §10-11) and BEAR-GRN
>   (`docs/experiments/bear_grn_benchmark.md` §10).

Cell type: K562 (main), plus Macrophage and MCF7 (multi-cell-type / transfer,
Section 5). Primary metric: AUPR. Also reported: AUROC, early precision
(EP), EPR (early precision ratio vs. random).

**Important caveat that applies to every table in this file**: every number
below is a **single run at a fixed seed (`seed: 42`, hard-coded project-wide)**.
There is no repeated-seed variance estimate, confidence interval, or
significance test anywhere in this codebase. Comparisons phrased as "beats"
or "within noise" (e.g. Section 3's model-size sweep, Section 7's baseline
comparison) are differences between n=1 runs, not statistically validated
claims — a difference described as "essentially noise-level" has not
actually been measured against a real noise floor. Treat every close call
in this file (anything within a few points of AUPR) as a plausible
direction, not a settled result, until repeated-seed runs exist.

> **Update (2026-09-18):** two qualifications to the caveat above.
> (1) There is now one exception. The recommended headline configuration
> (`all_at_once` + FM + h384/l2) has a real 5-seed mean±std (seeds 42-46,
> std ≤ 0.0017), reported in Section 7. The scMultiomeGRN baseline and every
> other table here are still n=1. (2) Before the 2026-09-18 fix, MEvD-GRN's
> training scripts never seeded model weight initialization, so runs made
> before then were "single run, unseeded" rather than truly "seed 42"
> (`plan.md` Section 9).

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
role-aware decoder; see `docs/figures/architecture_diagram_v3_1.png`)

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
| `gated_fusion` | legacy convex gated fusion + bilinear decoder | 0.883 | 0.970 |
| `concat_fusion` | legacy concat fusion + bilinear decoder | 0.890 | 0.971 |
| `edge_mlp` | nonlinear edge embedding instead of linear terms | 0.882 | 0.971 |
| `gated_relations` | learned per-layer relation weight (Section "gated combiner") | 0.876 | 0.969 |
| `pert_only` | train ONLY on perturbation (skip localization pretraining) | 0.441 | 0.757 |

> **Note (2026-10-08):** the `rna_only` row does not match its result file.
> `results/ablations/rna_only_K562.json` gives dual AUPR 0.8585 / AUROC
> 0.9633 (0.859 / 0.963), which is what `paper/main.tex` uses; the 0.856 /
> 0.962 above is left as written (probably an earlier run;
> `docs/mevd_vs_scmultiomegrn.md` §8, D18). Every row in this table is
> pre-fix (see the banner at the top), single-run, at the 315K-param size.

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
- **Integration style (role-aware vs. gated-fusion vs. concat-fusion vs.
  edge-MLP) is a genuine second-order effect**: all four sit within 0.01 AUPR
  of each other (0.882-0.890). The role-aware design is motivated by biology
  and interpretability, not by a raw accuracy gap over simpler fusion.
- **The gated relation combiner (Section "learned gated relation combiner")
  was neutral with only 2 relations** (0.876 vs. 0.881 sum) — as expected.
  **Update, 2026-09-18**: the real motif graph (28,204 PWM-backed edges) was
  finally used in a training run for the first time, via `motif_graph_only`
  (GNN restricted to *only* the motif relation, FM + h384/l2 size): AUPR
  0.724 / 0.818 / 0.955 (loc/pert/dual), respectable and reasonably close to
  the full sequential+FM+h384/l2 model's own 0.751/0.834/0.963 on the same
  size — the motif graph carries real, substantial signal on its own, not a
  near-empty relation. This makes it a legitimate 3rd relation for the gated
  combiner to weigh, unlike before. However, the `gated_relations` rerun at
  this size was mistakenly launched **without** `use_motif: true` in its
  config (an ablation-name override forces this for `motif_graph_only` but
  not for `gated_relations`), so it still only saw 2 relations — the actual
  "does the combiner do something interesting with 3 real relations"
  question remains open, rerun queued (job 2701095) and pending GPU
  availability (currently occupied by a scMultiomeGRN multi-seed run).
- **`pert_only` (0.441) is far below the full sequential curriculum (0.881)**:
  localization pretraining is doing real, load-bearing work, not just
  something later fine-tuning overwrites and could skip.

**Trained-on-dual, NOT zero-shot** (kept separate — comparing these to the
zero-shot rows above would be misleading):

| Ablation | What it does | Dual AUPR | Dual AUROC |
|---|---|---:|---:|
| `with_replay` | 3-stage curriculum, trains ON dual_evidence (frozen-encoder Stage 3) | 0.957 | 0.992 |
| `dual_only` | trains ONLY on dual_evidence's own (tiny) train split | 0.457 | 0.774 |

`with_replay`'s high score simply reflects that it was fine-tuned directly on
dual-evidence data before being scored on dual-evidence's test split — a much
easier task than the zero-shot claim the rest of this table makes, and not
directly comparable to it.

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
better to the nested dual-evidence tier. **Sequential remains the default
at this (307K-param, no-FM) size.**

**Update, 2026-09-18 — this call flips at the FM + bigger-model size.**
The perturbation-dilution tradeoff above assumed the tiny base model; once
FM embeddings and 10x the capacity are added (h384/l2), `all_at_once`
stops trading perturbation away and instead wins all three tiers outright,
including perturbation — see Section 7 for the full numbers and comparison
against scMultiomeGRN. **`all_at_once` + FM + h384/l2 is now the
recommended configuration**, not sequential. The mechanism above (positive
dilution) was real at 307K params; it evidently stops being the binding
constraint once the model has enough capacity to fit both tiers' signal
without one crowding out the other.

> **Correction (2026-10-08): the 2026-09-18 update above is wrong.** This
> file's own Section 7 table contradicts it. At FM + h384/l2, sequential
> gives perturbation AUPR 0.834 / dual 0.963, and `all_at_once` gives 0.696 /
> 0.952 (5-seed means 0.6955 / 0.9517). `all_at_once` wins only localization
> (0.968 vs 0.751); it does not win perturbation or dual evidence, so it does
> not "win all three tiers outright". The perturbation-dilution tradeoff did
> not go away at the bigger size. The leak fix (banner at the top) widened the
> gap: fixed 5-seed perturbation AUPR is 0.8165 for the curriculum and 0.5722
> for `all_at_once`, and dual AUPR 0.9562 vs 0.9017
> (`docs/experiments/leakfix_rerun.md`). The recommendation of `all_at_once`
> therefore rested on a misread comparison. That doc recommends making the
> curriculum the headline model again; the decision is pending (`plan.md`
> Section 11). The text and numbers above are unchanged.

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

## 5b. Model-size sweep — picking the final hidden_dim / depth

Once the input got richer (RP-weighted ATAC + FM embeddings, this section
and Section 4/5), `plan.md` Section 3 called for measuring capacity instead
of guessing it. Five `hidden_dim x n_gnn_layers` configs, K562 + FM, one per
GPU on Ada:

| config  | params    | hidden | layers | dual AUPR | dual AUROC | dual EP |
|---------|-----------|--------|--------|-----------|------------|---------|
| h128l3  | 562,567   | 128    | 3      | 0.9515    | 0.9889     | 0.8725  |
| h256l2  | 1,516,295 | 256    | 2      | 0.9572    | 0.9903     | 0.8809  |
| h256l3  | 2,042,631 | 256    | 3      | 0.9595    | 0.9908     | 0.8830  |
| h384l2  | 3,257,479 | 384    | 2      | 0.9628    | 0.9916     | 0.8888  |
| h384l3  | 4,440,199 | 384    | 3      | 0.9638    | 0.9917     | 0.8911  |

Monotonically better with size (confirming a richer input can now use more
capacity, unlike the original 306K-parameter model), but with sharply
diminishing returns — h384l2 to h384l3 is +36% params for +0.001 dual AUPR.
**Chosen size: hidden_dim=384, n_gnn_layers=2 (3.26M params)** — ties the
largest config within noise at well under its parameter count.

---

## 6. Multi-cell-type generalization

**Reference ceiling** (train-from-scratch directly on the target cell
type's own data — "how good can a model get on this cell type at all"):

| Cell type | genes | TFs | localization edges | From-scratch AUPR | AUROC |
|---|---:|---:|---:|---:|---:|
| Macrophage | 16,202 | 22 | 112,383 | 0.841 | 0.915 |
| MCF7 | 16,730 | 250 | 1,209,558 | 0.933 | 0.921 |

**Joint multi-cell-type training** (new `scripts/12_joint_train.py`: one
shared model, one gradient step per cell type per epoch, two cell types
trained together, evaluated zero-shot on the THIRD, never-seen one). We
tested all three choices of which cell type gets held out, each with and
without the FM embedding, for a full robustness check:

| Held out | Trained on | AUPR (zero-shot) | AUROC | vs. single-cell-type transfer |
|---|---|---:|---:|---|
| Macrophage | K562+MCF7, no FM | **0.743** | **0.835** | 0.607 (K562 alone) |
| Macrophage | K562+MCF7 + FM | 0.703 | 0.782 | 0.447 (K562 alone + FM) |
| MCF7 | K562+Macrophage, no FM | **0.837** | **0.802** | 0.758 (K562 alone) |
| MCF7 | K562+Macrophage + FM | 0.822 | 0.783 | 0.597 (K562 alone + FM) |
| K562 | Macrophage+MCF7, no FM | 0.704 | 0.691 | *(no single-source baseline)* |
| K562 | Macrophage+MCF7 + FM | **0.744** | **0.722** | *(no single-source baseline)* |

**Three honest findings:**
1. **Joint training on more than one cell type is a big, unambiguous win for
   transfer, with or without the FM embedding, for every choice of held-out
   cell type tested.** Where we have a single-cell-type-transfer baseline to
   compare against (Macrophage, MCF7 held out), joint training wins by
   0.05-0.26 AUPR every time. This is the strongest, most reliable lever
   found for the transfer problem specifically.
2. **Whether the FM embedding helps or hurts under joint training is
   holdout-dependent, not a fixed effect.** Holding out Macrophage or MCF7,
   FM is slightly worse than no-FM (by 0.02-0.04 AUPR) — consistent with the
   single-cell-type pattern in Section 5. Holding out K562, FM is slightly
   *better* (0.744 vs. 0.704) — the pattern reverses. All six numbers sit in
   a tight 0.70-0.84 band, so this is a small, second-order effect riding on
   top of the much larger, consistent joint-training win, not a contradiction
   of the headline result.
3. **The joint (no-FM) Macrophage-held-out and MCF7-held-out results both
   land within 0.10 AUPR of that cell type's own from-scratch ceiling**
   (0.743 vs. 0.841; 0.837 vs. 0.933) despite the model never training on
   that cell type at all — a genuinely strong generalization result by any
   standard.

**This cross-cell-type, zero-shot-to-a-third-cell-type result is the
strongest, most complete finding in the project** and — as far as we can
tell — nobody else has run this kind of multi-cell-type experiment on
SC-MO-GRN-DB.

Sources: `results/joint_training/*.json`, `results/transfer/*.json`.

---

## 7. Baselines (K562, same splits)

Script: `scripts/05_run_baselines.py`. AUPR/AUROC unaffected by the EPR fix
(only EPR itself changed); these are current.

| Method | Loc AUPR | Loc AUROC | Pert AUPR | Pert AUROC | Dual AUPR | Dual AUROC |
|---|---:|---:|---:|---:|---:|---:|
| MEvD-GRN, base (Section 1, 307K params, no FM) | 0.573 | 0.515 | 0.641 | 0.899 | 0.881 | 0.971 |
| MEvD-GRN, sequential + FM + h384/l2 (3.26M params) | 0.751 | 0.771 | 0.834 | **0.964** | 0.963 | **0.992** |
| **MEvD-GRN, `all_at_once` + FM + h384/l2, 5-seed mean±std** | **0.9676**±.0002 | **0.9609**±.0002 | **0.6955**±.0017 | 0.9096±.0004 | **0.9517**±.0007 | **0.9880**±.0002 |
| GRNBoost2 | 0.535 | 0.499 | 0.206 | 0.550 | 0.173 | 0.525 |
| RegDiffusion | 0.546 | 0.515 | 0.185 | 0.493 | 0.165 | 0.490 |
| GMF-GAE | 0.526 | 0.514 | 0.218 | 0.565 | 0.160 | 0.516 |
| scMultiomeGRN (adapted, single seed so far) | 0.906 | 0.914 | 0.655 | 0.919 | 0.847 | 0.971 |

> **Correction (2026-10-08):** the three MEvD-GRN rows above are pre-fix
> (banner at the top); the baseline rows are not. The "5 of 6 metrics" win in
> the next paragraph compared a leaky MEvD-GRN with a non-leaky
> scMultiomeGRN. With the fix (5 seeds, `docs/experiments/leakfix_rerun.md`),
> `all_at_once` + FM + h384/l2 scores loc 0.9544 / 0.9485, pert 0.5722 /
> 0.8803, dual 0.9017 / 0.9781 (AUPR / AUROC). Against scMultiomeGRN it now
> wins localization and dual evidence (4 of 6 metrics) and loses
> perturbation on both AUPR (0.572 vs 0.655) and AUROC (0.880 vs 0.919). The
> fixed curriculum wins perturbation and dual evidence by wide margins and
> loses localization. Those fixed numbers still carry the graph leak
> (`docs/experiments/labelfree_graph_rerun.md`; rerun pending).

**Head-to-head vs. the strongest baseline, resolved 2026-09-18, now with real
variance data**: a 5-seed rerun (seeds 42-46) of `all_at_once` + FM + h384/l2
landed with **remarkably low run-to-run variance** (std ≤ 0.0017 on every
metric) — this genuinely answers Limitation (1)/(3) below for this specific
configuration, not just asserts it. Comparing means against scMultiomeGRN's
current (still single-seed) point estimate: MEvD-GRN wins **5 of 6**
metrics by a margin that is tens to over a hundred times the measured
seed-to-seed std (localization AUPR +0.062 vs std 0.0002; dual-evidence AUPR
+0.105 vs std 0.0007; dual-evidence AUROC +0.017 vs std 0.0002; perturbation
AUPR +0.041 vs std 0.0017; localization AUROC +0.047 vs std 0.0002) — these
are not close calls. The one exception is **perturbation AUROC, which
MEvD-GRN now loses cleanly**: 0.9096±0.0004 vs scMultiomeGRN's 0.919, a
0.0094 gap that is also many times the measured std in the other direction
— a real, small, consistent loss, not the "near-tie" we described before
this rerun landed (that earlier read was an honest guess in the absence of
variance data, and turned out to be wrong in the specific direction of
underselling how real the gap is on both sides). scMultiomeGRN itself does
not yet have multi-seed data (2 additional seeds are running now, ETA
several more hours); until it does, comparisons against it specifically
should still be read with a grain of salt, even though MEvD-GRN's own side
of the comparison is now solid.

This is the first configuration found that doesn't trade tiers against each
other: earlier, sequential+FM+h384/l2 already won perturbation/dual but
lost localization (0.751 vs 0.906); switching only the curriculum protocol
(same model, same FM embedding) closes that gap **and improves the other
two tiers as well**, confirming Section 3's diagnosis that sequential
fine-tuning was actively damaging localization, not that the architecture
itself was weak there. `all_at_once` + FM + h384/l2 is now the
recommended configuration for the paper's headline numbers, superseding
the earlier "sequential is the default" call in Section 3 — that call was
correct for the small, no-FM model it was measured on, but doesn't hold
once FM + more capacity are in play.

> **Correction (2026-10-08): the paragraph above misreads the table.**
> Switching sequential + FM + h384/l2 to `all_at_once` (same model, same FM
> embedding) closes the localization gap (0.751 → 0.968), but it does **not**
> improve the other two tiers: perturbation AUPR falls from 0.834 to 0.696
> and dual evidence from 0.963 to 0.952 (rows 2 and 3 of the table above).
> `all_at_once` trades tiers too; it only wins localization. "Sequential is
> the default" therefore still held at this size, and the leak fix made the
> curriculum's lead larger (Section 3 correction;
> `docs/experiments/leakfix_rerun.md` item 4). The table does not support
> the `all_at_once` recommendation; which model becomes the paper's headline
> is pending (`plan.md` Section 11).

A secondary run, `with_replay` (3-stage curriculum, FM + h384/l2), was
also tested: localization 0.890 (better than sequential's 0.751, still
below `all_at_once`'s 0.967 and below scMultiomeGRN's 0.906), perturbation
0.731 (beats scMultiomeGRN). Its dual-evidence number (0.990 AUPR) is
**not zero-shot** — this protocol trains directly on dual_evidence in a
frozen-encoder Stage 3, so it isn't comparable to the other rows in this
table; it's excluded here for that reason (see Section 2's "trained-on-dual"
subtable for the established convention).

A third run, `fm_only` (Geneformer alone, hand-built RNA zeroed) at
h384/l2, essentially **ties** `with_fm` at the same size (loc 0.752 vs
0.751, pert 0.833 vs 0.834, dual 0.963 vs 0.963) — at this model size, the
hand-crafted RNA features add no measurable value once Geneformer and
enough capacity are both present.

**scMultiomeGRN status**: this baseline is a reimplementation of Xu et al.
(NAR 2025), adapted to MEvD-GRN's TF→all-gene benchmark (the published
method is TF-TF only) — see `src/baselines/scmultiomegrn_wrapper.py`'s
docstring for the full list of deliberate adaptations (own features, not
MAESTRO RP/GRNBoost2/FIMO; own train/val/test split, not the paper's
tenfold+consensus protocol). Treat any comparison against it as
"architecture-inspired baseline on our benchmark," not "reproducing the
paper's reported numbers."

**Final status (2026-09-18, completed 2026-09-12)**: an earlier attempt was
manually stopped after epoch 1145 (localization only, AUPR 0.901) because
the naive single-GPU run was too slow to reach perturbation within a
reasonable time — see below. Re-run as a full localization→perturbation
curriculum via `scripts/13_scmultiomegrn_full_ddp.py` (hand-rolled 4-GPU
data-parallel training, world_size=4, all 4 of Ada's 2080Tis), matching
MEvD-GRN's own `curriculum.main_curriculum_tiers` protocol: localization
and perturbation each trained for the full 2000-epoch budget (neither
early-stopped — both were still slowly improving at epoch 2000, e.g.
perturbation val AUPR kept climbing to 0.654), dual-evidence scored
**zero-shot** from the perturbation-tier model (dual-evidence is 100%
nested inside perturbation and 77.5% inside localization — training on it
directly would not be a comparable number to MEvD-GRN's own zero-shot
claim, so it wasn't). Total wall time: **5h 58m** end-to-end (data sync +
smoke test + both 2000-epoch tiers + final eval), source:
`results/logs/smg_ddp_full.log`, `results/baselines/scmultiomegrn_K562.json`.
On localization and perturbation specifically, this now beats MEvD-GRN's
own numbers (0.906 vs 0.573 localization AUPR; 0.655 vs 0.641 perturbation
AUPR) — plausible since scMultiomeGRN trains a dedicated model per tier
with no forgetting between stages, while MEvD-GRN's sequential curriculum
deliberately keeps a single model across both (Section 3). MEvD-GRN still
wins the tier that matters most for the paper's claim, zero-shot
dual-evidence (0.881/0.971 vs 0.847/0.971) — AUROC ties, AUPR favors
MEvD-GRN.

The RNA-only baselines (GRNBoost2, RegDiffusion, GMF-GAE) all cluster near
random on perturbation/dual (AUPR 0.16-0.22, AUROC 0.49-0.57) — MEvD-GRN's
multi-omic graph structure clearly matters relative to expression-only
methods.

---

## 8. Where to look for more detail

- `plan.md` — the living roadmap: why each change was made, what's still
  open (Section 10), effort/risk estimates.
- `docs/README.md` — index of all supporting documentation.
- `docs/figures/architecture_diagram_v3_1.png` — current architecture diagram.
- `docs/citations.md` — design rationale + literature citations per
  component.
- `docs/critical_review_independent.md` — the 2026-09-18 independent
  critical review and its final status update.
- `docs/archive/` — the pre-bugfix results catalog
  (`2026-09-07_results_pre_bugfix.md`), the original full theory spec, an
  early progress log and a pipeline walkthrough, kept for historical
  reference. Numbers there predate the bugfixes in Section 0 and should not
  be treated as current. Older weekly-update slide text is in
  `docs/weekly_updates/`, with the same caveat.
- `paper/main.tex` — manuscript. Its body was resynced to this file on
  2026-09-18, including the 5-seed headline numbers (`plan.md` Section 9).
  (2026-10-08: so every number in it is pre-fix too. The changes it needs
  are listed in `docs/paper_revision_plan.md`.)
- Added 2026-10-08:
  - `docs/experiments/leakfix_rerun.md`: the negative-sampling leak and its
    5-seed K562 rerun (banner at the top).
  - `docs/experiments/labelfree_graph_rerun.md`: the second (graph) leak and
    its K562 rerun (in progress).
  - `docs/experiments/pbmc10k_linger_benchmark.md`: PBMC10k Multiome vs
    LINGER, the main external benchmark (interim results in §10-11).
  - `docs/experiments/bear_grn_benchmark.md`: BEAR-GRN on SC-MO-GRN-DB
    datasets, with its hub and coverage artifacts (§10).
  - `docs/mevd_vs_scmultiomegrn.md`: comparison with scMultiomeGRN, novelty
    audit and the 18 code-vs-doc discrepancies.
