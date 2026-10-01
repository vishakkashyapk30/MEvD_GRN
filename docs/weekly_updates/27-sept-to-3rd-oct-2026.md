# Weekly update: 27 Sept – 3 Oct 2026

**Vishak Kashyap K · MeVD-GRN**
Status as of 2 Oct 2026, 04:00 IST. All work is on branch
`leak-fix-and-benchmarks`; `main` is unchanged since 29 Sept.

---

## 1. Summary

1. **We found and fixed a train/test leak, and it inflated the paper's
   headline numbers.** With the fix, the recommended model (`all_at_once`)
   loses to scMultiomeGRN on perturbation.
2. **The curriculum model is the stronger one, before and after the fix.**
   The paper's choice of `all_at_once` rested on a misread comparison. The
   curriculum should become the headline model again.
3. **We picked the main external benchmark after a survey of about 130
   papers:** 10x PBMC10k Multiome, scored against LINGER's Cistrome ChIP-seq
   evaluation.
   - On it, MeVD-GRN scores **AUROC 0.734 ± 0.005** against LINGER's
     published 0.714.
   - That holds with whole TFs held out of training and without external
     pretraining.
   - It is provisional until both methods are scored over identical
     candidate genes.
4. **BEAR-GRN (Nat Commun, Sept 2026) is the new baseline for our
   SC-MO-GRN-DB results,** and its pipeline is built. Its metrics are
   partly driven by hub genes: a trivial "number of targets" ranker beats
   every published method on K562.
5. **Ada was upgraded,** and the whole pipeline now runs on the new layout.

---

## 2. Negative-sampling leak (found and fixed 30 Sept)

**The bug:** the trainer drew random negatives from the *whole* negative
pool every epoch. That pool also supplies the val/test negatives.
- On K562 the localization stage asks for about 4M negatives per epoch from
  a 1M pool. So it trained on every val/test negative, as a label-0
  example, every epoch.
- The baselines never had this leak, so the comparison was unfair in our
  favour.

**The fix (commit c26c03d):**
- Every tier's val/test negatives are removed from the pool before
  training.
- `train_stage` refuses to run without this step.
- A config flag (`exclude_eval_negatives: false`) reproduces the old
  behaviour, for comparison.

**Rerun on Ada** (K562, Geneformer + h384/l2, 5 seeds; full detail in
`docs/experiments/leakfix_rerun.md`):

| Model | loc AUPR | pert AUPR | dual AUPR (zero-shot) | dual AUROC |
|---|---|---|---|---|
| `all_at_once`, before fix (paper headline) | 0.968 | 0.696 | 0.952 | 0.988 |
| `all_at_once`, fixed | 0.954 | 0.572 | 0.902 | 0.978 |
| **Curriculum, fixed** | 0.745 | **0.817** | **0.956** | **0.990** |
| scMultiomeGRN (never leaky) | 0.906 | 0.655 | 0.847 | 0.971 |

- **Sanity check:** with the fix switched off, the rerun reproduces the old
  numbers exactly, to 4 decimals. The drop comes from the fix alone.
- **The leak hit `all_at_once` much harder than the curriculum.** For
  `all_at_once`, perturbation fell by 0.12; for the curriculum, by 0.02.
  `all_at_once` turns hard negatives off, so all its negatives come from the
  leaky random pool.
- **The curriculum beat `all_at_once` on perturbation and dual evidence
  even before the fix.** results.md's claim that `all_at_once` "wins all
  three tiers" contradicts its own table.
- **The curriculum's remaining weakness is forgetting the localization
  tier** (0.745).

**A second, milder leak is still open:** the TF-candidate graph is built
with held-out positives excluded, so the graph depends on test labels.

---

## 3. Choosing the benchmark: literature survey

- **Knowledge base:** the 17 local papers and theses
  (`docs/reference/literature_knowledge_base.md`). It covers datasets,
  ground truths, protocols and the numbers to beat.
- **Web survey:** about 130 papers in three independent searches: high-impact
  journals, ML venues and preprints, and benchmarks and reviews.
  `docs/reference/multiome_grn_benchmark_consensus.md` and `lit_survey/`.
- **Consensus:** the **10x `pbmc_granulocyte_sorted_10k` Multiome**. It is
  truly paired RNA+ATAC from the same nuclei, and 17–22 papers use it,
  including SCENIC+, LINGER, scGLUE, KEGNI and scTFBridge.
- **Ground truth:** Cistrome blood ChIP-seq scored as LINGER does (Nat
  Biotechnol 2024), reused by scTFBridge, KEGNI and regX. The bar is
  **LINGER, AUROC 0.714**.
- **Is LINGER state of the art?** It is the accepted number on this
  benchmark and ranks first in BEAR-GRN. But other benchmarks rank
  differently: GRETA favours Dictys and SCENIC+, and Open Problems' geneRNIB
  favours GRNBoost2. No method does well against perturbation-based ground
  truth.
- **What reviewers will expect from supervised models** (InfoSEM, ICML
  2025; Stock et al. 2025):
  - TF-held-out splits
  - degree-matched negatives
  - trivial "hub" baselines
  Random edge splits let simple tricks match GNNs.
- **scMultiomeGRN's own benchmark was set aside.** Its labels are motif
  scans of the same ATAC data used as input, and its fetal-lung data is
  probably not truly paired. The pipeline is built but not run.

---

## 4. PBMC10k Multiome results (main benchmark)

**Protocol:**
- Evaluation is LINGER's own code (top-1000 Cistrome targets, AUC, AUPR
  ratio).
- 20 ChIP-seq datasets: 10 TFs × 4 cell types.
- LINGER's 9,543 cells; the cell-type counts match KEGNI's exactly.

**Our fair design, fixed before any results:**
- Train on CollecTRI with **all 10 evaluation TFs removed**.
- Degree-matched negatives; graphs built without labels.
- Model selection on validation only; 5 seeds.

All 220 runs are done. Detail: `docs/experiments/pbmc10k_linger_benchmark.md`
§10–11.

| Method | AUROC (19 datasets) | AUPR ratio |
|---|---|---|
| **MeVD-GRN** (Geneformer + h384, held-out TFs) | **0.734 ± 0.005** | 2.135 |
| LINGER (published) | 0.714 | **2.253** |
| MeVD-GRN, RNA only | 0.720 ± 0.030 | 2.178 |
| MeVD-GRN base (no Geneformer) | 0.670 ± 0.008 | 1.722 |
| MeVD-GRN base, RNA only | 0.544 ± 0.011 | 1.337 |
| MeVD-GRN, uniform negatives | 0.597 ± 0.014 | 1.698 |
| Trivial baselines (target count, gene ID, Pearson) | 0.50–0.59 | 1.2–1.6 |
| SCENIC+ / GENIE3 (published) | 0.548 / 0.539 | 1.29 / 1.17 |

**Takeaways:**
- **We beat LINGER on AUROC but not on AUPR ratio.** The win also has to
  hold once both methods are scored over the same candidate genes. That
  needs the LINGER re-run, which is being retried with more memory.
- **The signal is real, not just hub genes:** the trivial baselines reach
  only 0.50–0.59.
- **ATAC matters, but Geneformer largely substitutes for it.**
  - Without Geneformer, adding ATAC gains **+0.126 AUROC**.
  - With Geneformer, the gain is only +0.015, within noise.
- **Degree-matched negatives are essential:** 0.734 vs 0.597 with uniform
  negatives.
- **Open problem:** the model fails on target genes it never saw in
  training. With target genes held out, AUROC is about 0.45.

---

## 5. BEAR-GRN benchmark (new baseline for the SC-MO-GRN-DB results)

- **What BEAR-GRN is:** Karamveer…Uzun, *Nat Commun* (18 Sept 2026), from
  the same lab as SC-MO-GRN-DB. It benchmarks 9 multiome GRN methods on
  SC-MO-GRN-DB data against ChIP, knockout, their union, and their
  intersection. LINGER ranks first overall.
- **What we built:** BEAR-GRN's exact inputs and ground truths, plus a
  Python port of its R scoring.
  - The port reproduces their published LINGER and CellOracle scores
    (AUPRC exactly).
  - We train with TFs held out in 5 folds, the design the BEAR authors name
    in their peer-review file.
  - Detail: `docs/experiments/bear_grn_benchmark.md`.

**First K562 result (1 seed, a direction only):**

| K562 | MeVD-GRN | Ranking by target count | Best BEAR method |
|---|---|---|---|
| ChIP AUROC / AUPRC | 0.602 / 0.468 | **0.652 / 0.484** | 0.575 / 0.430 |
| Union AUROC / AUPRC | **0.634 / 0.347** | 0.626 / 0.344 | 0.570 / 0.344 |

**What the benchmark is measuring:**
- **A hub ranker beats every published method.** Ranking target genes by
  how many TFs regulate them beats all 9 published methods on K562 ChIP.
- **"Predict every measured pair" ties LINGER.** That baseline gets ChIP
  AUPRC 0.431, LINGER's score.
- **Knockout ground truth is unbeatable.** Every method is at or below
  random against the K562 knockouts.

So on K562, our margin over BEAR's methods is mostly the hub effect. A
hub-controlled check comes next. These artifacts are worth reporting in the
paper.

---

## 6. MeVD-GRN vs scMultiomeGRN: novelty and provenance study

`docs/mevd_vs_scmultiomegrn.md` covers a side-by-side comparison, where
each component came from, ranked novelty claims, risks, and 18 places where
the code and the docs/paper disagree.

**Claimable novelty** (none of it is a new architecture block):
1. **The evidence-tier protocol:** train on localization + perturbation,
   then test zero-shot on dual evidence.
2. **The Geneformer finding:** it boosts in-domain accuracy but hurts
   transfer to new cell types (K562→Macrophage 0.607→0.447).
3. **Joint training on two cell types with zero-shot transfer** to a third.
4. **The accessibility-gated decoder,** presented for interpretability only.

**Not novel:** Geneformer embeddings in a GRN GNN (scRegNet),
MAESTRO-style ATAC features, separate towers per modality, replay for GRNs
(LINGER), and hard-negative mining.

**Problems found in the paper:**
- The abstract centres on the curriculum, but the recommended model doesn't
  use one.
- The paper's descriptions of the RNA input, the ATAC encoder and the
  openness term don't match the code.
- The bilinear score is credited to DistMult; it is RESCAL.
- One citation (Duren et al.) has a DOI that doesn't resolve; it mixes two
  papers.
- No shipped config reproduces the headline.

---

## 7. Infrastructure and housekeeping

**Repo organisation (29 Sept):**
- Docs merged and renamed consistently.
- Weekly updates renamed by date.
- Superseded figures removed, archive dated.

**Ada after the upgrade** (details at the top of `docs/reference/ada.md`):
- The login gateway is now `ada-gw1`.
- Partition `u22`, with 2080 Ti nodes requested via `--constraint=2080ti`.
- `/share1` is visible only from the gateway, so jobs copy data to node
  scratch.
- Medium QoS caps each user at 4 GPUs, 40 CPUs, 125 GB of memory, and 20
  jobs.
- A new conda environment (torch 2.9.1, CUDA 12.6) and a clone of the repo
  are set up on Ada.

**Ada jobs this week:**

| Job | What | Status |
|---|---|---|
| 467 | K562 leak-fix rerun (12 runs) | done |
| 473 / 474 | PBMC10k build + 220 training runs | done |
| 475 | PBMC GRNBoost2 baseline | running |
| 819 → 2399 | PBMC LINGER re-run on our cells | 819 ran out of memory at 30 GB; resubmitted with 120 GB, queued |
| 587 → 2400 → 2401 | BEAR-GRN K562 + macrophage (5 seeds, 5 variants) | prep done; work and scoring queued behind LINGER |

---

## 8. Decisions needed

1. **Headline model:** switch the paper back to the curriculum. The rerun
   supports this.
2. **Story:** lead with PBMC10k against LINGER, and present SC-MO-GRN-DB and
   BEAR-GRN as the supporting results. Report BEAR's hub and coverage
   artifacts openly.
3. **Paper rewrite:** do it after the LINGER re-score and the remaining open
   checks. That covers the 18 code/paper discrepancies, the corrected
   numbers, and the reframed novelty.

## 9. Plan for next week

- **Finish the comparisons.**
  - Re-score PBMC on LINGER's own candidate genes once job 2399 finishes.
  - Add GRNBoost2.
  - Re-run SCENIC+ and scTFBridge if feasible.
- **Investigate model failures.**
  - The target-gene-disjoint failure on PBMC.
  - Localization forgetting in the curriculum.
- **Fix the remaining leak and test the margins.**
  - Fix the second leak (label-dependent TF-candidate graph), then rerun
    K562.
  - Run the hub-controlled check on BEAR-GRN.
  - Finish the BEAR-GRN 5-seed runs.
  - Run phase 2: iPSC and mouse.
- **Fix the paper text** (discrepancies and citations) and update
  `results.md`.
