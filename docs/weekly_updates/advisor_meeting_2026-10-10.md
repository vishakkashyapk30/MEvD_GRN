# Advisor meeting, 10 Oct 2026, 11:00: progress since 27 Sept

Status as of **2026-10-09 20:35 IST**. Prepared as raw material for the slides (due 01:00).
Everything below is a **zero-leakage** result: label-free input graphs, no val/test
negatives in the training pool, and (on BEAR-GRN and PBMC10k) the evaluated TFs never
seen in training. Old numbers from before 30 Sept are not used.

Architecture figure: `docs/figures/architecture_simple_2026-10-09.png` (SVG next to it).

---

## 1. Timeline: what finishes when

**What is running** (laptop RTX 4060, shared by two queues; Ada is reachable only from the
campus network or the IIIT VPN, and the laptop is currently on a phone hotspot):

| Job | What | ETA (laptop only) | Done before 01:00? |
|---|---|---|---|
| BEAR-GRN **M0** (K562 dev, seed 42; 5 folds) | baseline MeVD-GRN under BEAR's protocol | **~21:15** (3 of 5 folds done at 20:33, ~20 min/fold) | **yes** |
| BEAR-GRN **M1** (uniform negatives) | candidate 2 | ~22:55 | **probably** |
| BEAR-GRN **M2** (leak-free hub term) | candidate 3 | ~00:35 | borderline |
| BEAR-GRN **M3** (hub + motif features) | candidate 4 | ~02:15 | **no** |
| Model selection (M0-M3, K562 validation only) | picks the headline model | ~02:30, right after M3 | **no** |
| K562 zero-leak seed 46 (2 runs, ~1 h each) | 5th seed | ~02:30-04:30, after BEAR | no |
| PBMC `target_all` seeds 43-44 (supplementary) | confirms the held-out-target fix | ~05:00+ | no |

- After each BEAR-GRN training run there is a separate scoring step (short, not yet timed).
- **If you connect the laptop to campus Wi-Fi or the VPN before ~21:30,** M1-M3 can run
  on Ada's four GPUs in parallel instead (about 1-1.5 h each). That would complete the
  selection by about **23:30-00:00**, in time for the slides. Tell me and I'll start it.

**What will NOT exist by tomorrow**
- The BEAR-GRN headline grid (5 seeds x 9 datasets, then the compendium-only regime and
  the three ablations, about 150 runs of ~1-2 h each). **~4-6 days on the laptop alone;
  ~1.5-2 days if Ada's 4 GPUs are used.** No MeVD-GRN number across all 9 BEAR datasets
  can be shown tomorrow.
- LINGER re-scored on LINGER's own candidate genes (the PBMC comparison stays provisional;
  its environment needs one package fix and ~120 GB of memory on Ada).
- Jaccard stability on BEAR-GRN (needs a package that isn't installed).

---

## 2. Results (zero-leak only)

### 2a. K562, SC-MO-GRN-DB (three evidence tiers); the head-to-head with scMultiomeGRN

Same data, same splits, same protocol for all rows (train on localization + perturbation;
dual evidence is scored zero-shot). AUPR / AUROC. MeVD-GRN rows are mean over **4 seeds**
(42-45; seed 46 is queued). Baselines are single-seed runs of ours.

| Model | Localization | Perturbation | Dual evidence (zero-shot) |
|---|---|---|---|
| **MeVD-GRN, evidence curriculum** (Geneformer, 3.26M params) | 0.730 / 0.751 | **0.814 / 0.960** | **0.954 / 0.990** |
| MeVD-GRN, single-stage (`all_at_once`) | **0.953 / 0.947** | 0.582 / 0.883 | 0.902 / 0.979 |
| scMultiomeGRN (our re-run) | 0.906 / 0.914 | 0.655 / 0.919 | 0.847 / 0.971 |
| GRNBoost2 | 0.535 / 0.499 | 0.206 / 0.550 | 0.173 / 0.525 |
| RegDiffusion | 0.546 / 0.515 | 0.185 / 0.493 | 0.165 / 0.490 |
| GMF-GAE | 0.526 / 0.514 | 0.218 / 0.564 | 0.160 / 0.516 |

- Std over seeds is at most 0.001 for the curriculum and 0.012 for single-stage (perturbation).
- **Takeaway:** the curriculum beats scMultiomeGRN on perturbation (+0.16 AUPR) and on the
  zero-shot dual-evidence tier (+0.11), and loses on localization (-0.18; it forgets the
  first tier). Single-stage wins localization but loses perturbation.
- The second leak (label-dependent input graph) moved these numbers by at most ~0.01.
  The first leak (negatives) was the big one: it had inflated the old headline model by
  0.12 AUPR on perturbation.
- Caveats: the scMultiomeGRN row is our re-implementation (single seed), not their
  released code. This is not their own benchmark (fetal lung), which we did not run.

### 2b. PBMC10k Multiome vs LINGER (LINGER's Cistrome evaluation; 19 ChIP datasets)

Protocol: train on CollecTRI with all 10 evaluation TFs removed (TF-held-out), degree-matched
negatives, 5 seeds. Scored on all expressed genes (LINGER's published numbers use its own
candidate genes, so the comparison is provisional).

| Method | AUROC | AUPR ratio |
|---|---|---|
| **MeVD-GRN** (Geneformer + h384) | **0.734 ± 0.005** | 2.135 |
| MeVD-GRN, RNA only (no ATAC) | 0.720 ± 0.030 | 2.178 |
| MeVD-GRN, no Geneformer | 0.670 ± 0.008 | 1.722 |
| MeVD-GRN, no Geneformer and no ATAC | 0.544 ± 0.011 | 1.337 |
| MeVD-GRN, uniform instead of degree-matched negatives | 0.597 ± 0.014 | 1.698 |
| LINGER (published) | 0.714 | **2.253** |
| SCENIC+ / GENIE3 (published) | 0.548 / 0.539 | 1.29 / 1.17 |
| Simple baselines (target count, gene ID, Pearson) | 0.58 / 0.58 / 0.58 | 1.4-1.5 |
| **One TF-agnostic ranker: ATAC gene activity alone** | **0.814** | **3.680** |

- **ATAC matters when Geneformer is absent** (+0.126 AUROC) but adds only +0.015 with it (within noise).
- **Key caveat for the advisor:** a ranker that gives *every TF the same gene ranking*, using
  only ATAC accessibility, scores 0.814, above both MeVD-GRN and LINGER. The Cistrome
  metric rewards a per-gene prior (open chromatin), so the 0.734 vs 0.714 comparison
  says little about TF-specific regulation. We have a per-TF analysis in the runbook (§12.4).
- Reproduction caveat (open): a local rerun of the same configuration gives 0.690, not
  0.734, for reasons not yet found.
- Held-out target genes: the original split scored 0.409 (below random) because of how
  its negatives were built; the fixed split (`target_all`) scores 0.745 / 2.278 (1 seed).

### 2c. BEAR-GRN (Nat Commun 2026): baselines on all 9 datasets; MeVD-GRN pending

BEAR-GRN scores each method on its own ChIP ground truth. All numbers below are scored
with our Python port of BEAR's scoring (matches the paper's AUPRC exactly on the methods
checked). ChIP **AUPRC**; "trivial" baselines are ours.

| Dataset | LINGER | Best other published | Rank by target count (in-degree) | Predict every measured pair | Random |
|---|---|---|---|---|---|
| K562 | 0.430 | LINGER 0.430 | **0.484** | 0.431 | 0.331 |
| Macrophage S1 | 0.321 | TRIPOD 0.332 | **0.365** | 0.331 | 0.272 |
| Macrophage S2 | 0.324 | TRIPOD 0.336 | **0.374** | 0.342 | 0.271 |
| iPSC | 0.108 | TRIPOD 0.119 | **0.150** | 0.137 | 0.091 |
| Naive mESC | 0.160 | TRIPOD 0.192 | **0.308** | 0.269 | 0.153 |
| Mouse E7.5 rep 1 | 0.241 | LINGER 0.241 | **0.323** | 0.274 | 0.155 |
| Mouse E7.5 rep 2 | 0.225 | LINGER 0.225 | **0.306** | 0.268 | 0.153 |
| Mouse E8.5 rep 1 | 0.227 | LINGER 0.227 | **0.312** | 0.265 | 0.154 |
| Mouse E8.5 rep 2 | 0.224 | LINGER 0.224 | **0.313** | 0.263 | 0.155 |

- **On every dataset, ranking target genes by how many TFs regulate them (using no
  expression or ATAC data) beats every published method, LINGER included.** On AUROC,
  LINGER beats it only on the four mouse embryo sets (0.64-0.66 vs 0.60).
- "Predict every measured pair" alone gets within 0.001 of LINGER on K562.
- On the knockout and intersection ground truths every method (and every baseline) is
  at or below random on K562.
- **MeVD-GRN under BEAR's protocol:** the pre-registered pipeline (TF-held-out 5-fold, dense
  scoring, four pre-declared candidates M0-M3, model chosen on K562 validation only) is
  built and running. **No result yet.** The first one (M0, K562) lands at ~21:15.
  An earlier single run (older protocol, 1 seed) gave K562 ChIP 0.602 AUROC / 0.468 AUPRC:
  above LINGER (0.536 / 0.430), below the target-count baseline (0.652 / 0.484).

---

## 3. What we can and cannot claim tomorrow

**Supported**
1. **Found and fixed two train/test leaks.** The first inflated the old headline model by
   0.12 AUPR on perturbation. Zero-leak is now enforced in code (training refuses to
   run otherwise) and recorded in every result file.
2. **With zero leakage, the evidence-tier curriculum beats scMultiomeGRN (our re-run) on
   the perturbation and zero-shot dual-evidence tiers on K562**, with very small seed variance.
   It forgets localization.
3. **Benchmark audit:** on BEAR-GRN and on PBMC10k/LINGER, trivial per-gene rankers
   (target count; ATAC accessibility alone) match or beat the published methods. This
   matters for how any "we beat LINGER" claim must be made (against these baselines too).

**Not yet supported**
- That MeVD-GRN beats LINGER on BEAR-GRN. No result under the pre-registered protocol exists.
- A win on PBMC10k that reflects TF-specific signal (the AUPR ratio is below LINGER's, and
  the TF-agnostic ranker is higher than both).
- Anything on scMultiomeGRN's own benchmark.

---

## 4. Suggested slide outline

1. Goal and the question (diagram: `architecture_simple_2026-10-09.png`)
2. What changed since 27 Sept: benchmark survey (~130 papers), new standard datasets
3. The leak: how found, what it inflated, how it's now blocked
4. K562 zero-leak results vs scMultiomeGRN (table 2a)
5. PBMC10k vs LINGER, with the trivial-baseline caveat (table 2b)
6. BEAR-GRN: new baseline, what the benchmark rewards (table 2c)
7. Status of the BEAR-GRN runs and timeline (section 1)
8. Next steps: model selection, 5-seed grid on Ada, hub-controlled analysis, paper rewrite
