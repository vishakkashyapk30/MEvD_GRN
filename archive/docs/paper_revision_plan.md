# paper/main.tex revision plan

Written 2026-10-08. **This is a plan only: `paper/main.tex` has not been
edited.** The user has not yet decided the headline model (gate G1 below).

- **Line numbers** refer to `paper/main.tex` as of 2026-10-08 (1,309 lines).
  Each item also quotes a short anchor phrase, so it can be found again after
  earlier edits shift the lines. Within a pass, edit from the bottom up.
- **Sources of truth for numbers.** Copy numbers from these docs; never
  retype them from memory.

| Topic | Source |
|---|---|
| K562, negative-sampling fix (5 seeds) | `docs/experiments/leakfix_rerun.md`, Results table |
| K562, label-free graph (second fix) | `docs/experiments/labelfree_graph_rerun.md` §1 (the leak), §5 (results, **pending**) |
| PBMC10k vs LINGER | `docs/experiments/pbmc10k_linger_benchmark.md` §2 (numbers to beat), §3-4a (protocol), §10-11 (results) |
| BEAR-GRN | `docs/experiments/bear_grn_benchmark.md` §5 (numbers to beat), §8 (port validation), §10 (results) |
| Novelty, risks, discrepancies | `docs/mevd_vs_scmultiomegrn.md` §4, §6, §7, §8 |
| Older K562 numbers (all pre-fix) | `results.md` (dated corrections of 2026-10-08 at the top and in Sections 2, 3, 7) |

**Pending-number tags** used below. A number with a tag must not go into the
paper until the tagged work lands.

| Tag | What is pending | Where it will appear |
|---|---|---|
| **[P-LF]** | K562 rerun with the label-free TF-candidate graph and the negative-pool fix (both leaks fixed). Running locally since 2026-10-08. Until it lands, the `leakfix_rerun.md` numbers are interim: they still carry the graph leak. | `labelfree_graph_rerun.md` §5 |
| **[P-LINGER]** | LINGER re-run on our cells (job 2399) and the re-score of every PBMC method on LINGER's own target set (`linger_tg`). Every PBMC number so far is on the `expressed` gene universe. | `pbmc10k_linger_benchmark.md` (new section) |
| **[P-GB2]** | PBMC GRNBoost2 (job 475; 1 of 19 datasets scored at the last compile). | same |
| **[P-BEAR]** | BEAR-GRN 5-seed runs and scoring on Ada (K562, Macrophage S1/S2), the sparse-matched N over all 8 methods, the indegree-residual control (planned, not built), and phase 2 (iPS, mouse). | `bear_grn_benchmark.md` §10 |
| **[P-SMG]** | scMultiomeGRN with more than one seed (still n = 1). | `results/baselines/` |
| **[P-EP]** | EP/EPR of the fixed runs. Neither `leakfix_rerun.md` nor `labelfree_graph_rerun.md` reports them. | to aggregate from the result JSONs |
| **[NOT RUN]** | No rerun is scheduled: the 315K-param ablation grid, the 3-stage `with_replay` run, Macrophage/MCF7 transfer and joint training, `rna_only` at h384 + FM, and K562 hub controls. Each is pre-fix and n = 1. | `plan.md` §10 items 1, 4, 11 |

Ada has been inaccessible since 2026-10-08 (`plan.md` §11.5), so [P-LINGER],
[P-GB2] and [P-BEAR] cannot progress until access returns.

---

## Phase 0. Gates (decide before editing)

- [ ] **G1. Headline model** (the user's decision). Option (a), recommended
  by `leakfix_rerun.md` item 3: the curriculum, i.e. `full_curriculum` =
  sequential localization (30 epochs, lr 1e-3) → perturbation (15 epochs, lr
  3e-4), replay of 10% of earlier-tier positives at weight 0.1, half of the
  negatives hard, Geneformer + h384/l2, 3,257,479 params. (Settings as read
  from `configs/default.yaml` → `configs/k562.yaml` → `configs/k562_fm.yaml`
  → `configs/sweep/k562_fm_h384_l2.yaml` on 2026-10-08; re-check before
  writing them into the paper, because another session is editing configs.)
  Option (b): keep `all_at_once`. Phases 1 and 2 give both variants where
  they differ.
- [ ] **G2. Wait for [P-LF]** before writing any K562 number, or write the
  `leakfix_rerun.md` numbers explicitly labelled as interim.
- [ ] **G3. Wait for [P-LINGER]** before claiming a win over LINGER.
  `docs/weekly_updates/27-sept-to-3rd-oct-2026.md` §8 suggests rewriting
  after both G2 and G3.

---

## Phase 1. Headline and framing

Reframe the novelty as in `docs/mevd_vs_scmultiomegrn.md` §4 and §7.
Claimable: (1) the evidence-tier supervision and evaluation protocol; (2) the
Geneformer finding (an in-domain gain that hurts single-source transfer); (3)
multi-cell-type joint training with transfer to a third cell type; (4) the
accessibility-gated decoder, as an interpretable prior with no accuracy gain.
**Not claimable as novel:** FM embeddings in a GRN GNN (scRegNet), MAESTRO-style
RP features, modality-specific towers (GraFRank / scMultiomeGRN),
GraphSAGE / R-GCN message passing, hard-negative mining in general, replay or
continual learning for GRNs (LINGER), and evaluation against an intersection
gold standard in itself (scMTNI, McCalla *et al.*).

- [ ] **1.1 Title**, L30 (`\icmltitlerunning{...Curriculum Learning over
  Confidence-Tiered...}`) and L35-36 (`\icmltitle{...via Curriculum Learning
  over Confidence-Tiered Evidence}`).
  - (a): can stay. The curriculum is then part of the headline model.
  - (b): drop "via Curriculum Learning" (D5). Lead with the evidence-tier
    protocol instead, e.g. "...Supervised by Confidence-Tiered Experimental
    Evidence".
- [ ] **1.2 Abstract**, L53-85. Rewrite it whole:
  - L61-67 ("treats this confidence structure as a training curriculum"): keep
    under (a); under (b) describe the tiers as the split and held-out
    evaluation, not as training stages (D5).
  - L67-69 ("315,015-parameter model reaches AUPR 0.881 / AUROC 0.971"):
    pre-fix. Drop it, or replace it with the headline's fixed dual-evidence
    number [P-LF].
  - L69-75 ("outperforms ... scMultiomeGRN on AUPR across all three evidence
    tiers ... (0.950 vs 0.847)"): **false after the fix.** Fixed
    `all_at_once` loses perturbation to scMultiomeGRN (0.572 vs 0.655); the
    fixed curriculum wins perturbation and dual evidence and loses
    localization (0.745 vs 0.906). Replace it with whichever holds for the G1
    choice [P-LF]. (Also E2: 0.950 is not the 5-seed mean, 0.9517.)
  - Add the main external result: PBMC10k Multiome, LINGER's Cistrome
    evaluation, evaluation TFs held out of training, 5 seeds: AUROC 0.7343 ±
    0.0046 vs LINGER's published 0.7143, AUPR ratio 2.135 vs 2.2526
    [P-LINGER]. State both metrics, since MeVD-GRN is lower on AUPR ratio.
  - Add one sentence on BEAR-GRN: on SC-MO-GRN-DB data its metrics are
    largely matched by a coverage baseline and a target-hub ranker, and
    MeVD-GRN is compared against those as well as the published methods
    [P-BEAR].
  - L75-79 (FM transfer 0.607 → 0.447; joint training): keep only with
    "single run, localization labels, before the leak fixes" [NOT RUN], or
    move it to the body.
  - L80-84: remove "no result in this paper yet has a repeated-seed variance
    estimate" (D6). Mention the two leaks found and fixed.
- [ ] **1.3 Editorial status note**, L87-107 (LaTeX comment). Replace it with
  a short note pointing to this file, or delete it before any submission.
- [ ] **1.4 Introduction, central contribution**, L164-175 ("Our central
  methodological contribution is to instead treat the tier hierarchy
  explicitly, as an evidence-tier curriculum ... depends on model capacity").
  - Delete "the better choice depends on model capacity": both before and
    after the fix, the curriculum beats `all_at_once` on perturbation and dual
    evidence at h384 + FM (`leakfix_rerun.md` item 4).
  - Claim novelty only for "the specific use of the experimental evidence
    types of one resource as training stages and held-out targets", and only
    "to our knowledge" (risk R8).
  - Call the dual-evidence test "held-out intersection-tier edges", not
    generalisation to a new label type (R4).
- [ ] **1.5 Contributions**, L177-195. Replace (1)-(4) with:
  1. the evidence-tier protocol: hierarchy-consistent global split,
     training on the constituent tiers, scoring held-out intersection-tier
     edges. Under (a), add the evidence-type curriculum, with its
     localization forgetting as a stated cost;
  2. the external benchmarks: PBMC10k vs LINGER with whole TFs held out,
     degree-matched negatives and trivial baselines; BEAR-GRN with its hub
     and coverage artifacts reported;
  3. the Geneformer finding, framed as a failure mode with its caveats;
  4. joint multi-cell-type training with transfer to a third cell type, with
     its caveats.
  The current (1) role-aware decoder becomes a design choice "motivated by
  biology and interpretability, with no measured accuracy gain" (R10). The
  current (2) "provable guarantee against a model being evaluated on an edge
  it was trained on" covers *positives* only; say so (see 2.2). The current
  (4) "reaching a single configuration that outperforms the strongest
  available baseline across all three evidence tiers" is false; delete it.
- [ ] **1.6 Related Work, scMultiomeGRN**, L218-225: apply D16 and use the
  wording of `mevd_vs_scmultiomegrn.md` §7.
- [ ] **1.7 Related Work, add the missing precedents** (new paragraphs after
  L225; references in Phase 6):
  - scRegNet, the frozen FM embeddings + GNN precedent;
  - GENELink and GNNLink, supervised TF→target link prediction;
  - LINGER (lifelong learning from external bulk data; also the PBMC bar);
  - KEGNI, scTFBridge and regX, other users of the PBMC10k benchmark;
  - BEAR-GRN;
  - scMTNI and McCalla *et al.*, intersection gold standards used for
    evaluation;
  - PECA and CellOracle, accessibility × expression;
  - InfoSEM, Stock *et al.* and Yılmaz *et al.*, hub/degree bias in
    supervised link prediction.
- [ ] **1.8 Related Work, replay claim**, L236-240 ("we ... find that in our
  setting replay alone is insufficient without an additional consolidation
  stage"): this rests on pre-fix 315K-param numbers (Table `tab:forgetting`).
  Rewrite it from the fixed curriculum numbers [P-LF], or drop it.
- [ ] **1.9 Positioning**, L242-246: rewrite along `mevd_vs_scmultiomegrn.md`
  §1 (directed TF→all-genes vs undirected TF-TF; experimental vs motif-derived
  labels; label-free structural graphs vs the training adjacency), and drop
  the "differs from all of the above" framing.
- [ ] **1.10 Conclusion**, L1186-1212. Rewrite it to the new contributions.
  - L1190-1194 ("outperforms the strongest available supervised, multi-omic
    baseline on AUPR across all three evidence tiers") is false after the
    fix.
  - L1202-1205 ("no variance estimates") contradicts the 5-seed results.
  - L1206-1212 ("no prior work has used this way"): soften to "as far as we
    can tell". BEAR-GRN benchmarks unsupervised methods on the same resource
    but does not train jointly.

---

## Phase 2. Replace every pre-fix number (D9)

Every MeVD-GRN number in `main.tex` was trained with the negative-sampling
leak (100% of val/test negatives in the training pool) and with the
label-dependent TF-candidate graph. In that graph, 0 val/test positives were
edges, against about 2.8-3.0% of val/test negatives
(`labelfree_graph_rerun.md` §1). The baselines had neither leak.

- [ ] **2.1 Data table**, L323 (`Random negative pool & 1{,}000{,}000`):
  split it into "negative pool 1,000,000; training pool after removing the
  300,000 val/test negatives: 700,000" (`leakfix_rerun.md`;
  `labelfree_graph_rerun.md` §3).
- [ ] **2.2 Split section**, L406-430.
  - Proposition L410-418 is correct for positives only. Add a remark that it
    says nothing about negatives or input graphs.
  - L424-430 ("We use 5 sampled negatives per positive for training and 5
    for evaluation, drawn from a pre-sampled pool of 1,000,000 ...")
    describes the leak. Replace it with the fixed procedure (val/test
    negatives are drawn from the pool and then removed from the training
    pool; `train_stage` refuses to run otherwise).
  - Add one sentence on the graph fix (see D4).
- [ ] **2.3 Forgetting table and text**, L611-653 (`tab:forgetting`: 0.509
  and 0.719). Both are pre-fix, 315K params.
  - (a): replace with the fixed curriculum's localization AUROC after Stage 2
    (interim 0.7602 ± 0.0032 from `leakfix_rerun.md`; final [P-LF]).
  - The 3-stage row has no fixed rerun [NOT RUN]: drop it or label it
    pre-fix.
- [ ] **2.4 Training-curve figure**, L664-673 (`fig7_training_curves.png`,
  shipped 2-stage model): regenerate it from a fixed run [P-LF], or label it
  pre-fix.
- [ ] **2.5 Results caveat**, L759-767: rewrite it. K562 headline: 5 seeds,
  both fixes [P-LF]. PBMC: 5 seeds. Baselines and ablations: n = 1.
- [ ] **2.6 Main results text**, L769-785 ("Switching *only* the training
  protocol ... removes that damage *and* improves the other two tiers").
  False even before the fix (sequential 0.834 / 0.963 vs `all_at_once` 0.696
  / 0.952). Rewrite it for the G1 choice.
- [ ] **2.7 Table `tab:main`**, L787-820.
  - Base-model rows L802-805 (0.573 / 0.641 / 0.881 ...): pre-fix, no rerun
    [NOT RUN]. Drop them or mark them pre-fix.
  - Recommended rows L807-810: replace with the G1 headline, both fixes
    [P-LF]. Interim values: `leakfix_rerun.md` (curriculum 0.7450 / 0.7602,
    0.8165 / 0.9602, 0.9562 / 0.9902; `all_at_once` 0.9544 / 0.9485, 0.5722 /
    0.8803, 0.9017 / 0.9781). EP/EPR [P-EP].
  - Caption L790-793 and footnote L816-819 ("resolved under
    `all_at_once`"): rewrite.
- [ ] **2.8 Figures** `fig1_mevdgrn_per_tier` (L822-827) and
  `fig5_catastrophic_forgetting_auroc` (L829-839): regenerate from fixed
  runs [P-LF].
- [ ] **2.9 Head-to-head text**, L841-864 ("MEvD-GRN wins on AUPR across all
  three tiers by margins of tens to over a hundred times ... including the
  harder ... perturbation tier"): false after the fix. Rewrite it from the
  [P-LF] numbers.
  - (a): the interim curriculum wins perturbation and dual evidence against
    scMultiomeGRN and loses localization.
  - (b): the interim `all_at_once` loses perturbation on both AUPR and AUROC.
  - Either way, scMultiomeGRN is n = 1 [P-SMG], and it is an adapted
    reimplementation (D17, R9): use the wording "outperforms an adapted
    scMultiomeGRN on our benchmark".
- [ ] **2.10 Table `tab:baselines`**, L866-887. Replace the MeVD-GRN row
  L878 [P-LF]. Under (a), bold localization moves to scMultiomeGRN (0.906).
  The baseline rows are unchanged; they were never leaky.
- [ ] **2.11 Alternative-schedule section**, L896-913 ("only switching to
  joint training at the larger size wins all three tiers outright"), and
  table `tab:head2head-alt`, L915-939.
  - Rewrite the text: at both sizes `all_at_once` wins only localization.
  - Table: either replace the 3.26M rows (L930-931) with the fixed rows
    [P-LF], or keep the whole table as an explicitly pre-fix history table
    in the appendix.
  - Regenerate `fig9_mevdgrn_vs_scmultiomegrn` (L941-946).
- [ ] **2.12 Ablation section**, L948-1003 (`tab:ablation`, heatmap), and
  appendix `tab:full-ablation`, L1231-1299. Every row is pre-fix, 315K
  params, n = 1 [NOT RUN]. Either rerun under both fixes or label the tables
  "development-time ablations, trained before the two leak fixes; directional
  only". Also fix E1 (the `with_fm` row).
- [ ] **2.13 Discussion**, L1005-1082. Keep only claims that survive 2.12.
  - L1014-1027 (ATAC "measurably below on every tier", 0.859 vs 0.872):
    n = 1, pre-fix, small model, and never measured at the headline size
    (R7). Add the PBMC result: ATAC adds +0.126 AUROC without Geneformer and
    +0.015 with it, within noise (`pbmc10k_linger_benchmark.md` §11).
  - L1045-1061: E3, and rewrite it per 2.6.
  - L1063-1082: E1.
- [ ] **2.14 Limitation (1)**, L1091-1111: rewrite it with the fixed
  numbers. Fixed `all_at_once` loses perturbation on both metrics, so the
  "five of the six method-tier cells" text is wrong. Keep the "scMultiomeGRN
  still n = 1" caveat [P-SMG].
- [ ] **2.15 Limitation (6)** L1163-1171, abstract L75-79 and Data L288-293
  (transfer): all transfer and joint-training numbers are pre-fix, n = 1,
  h128, localization only [NOT RUN]. Also E4.
- [ ] **2.16 Add a limitation on both leaks.** Describe what each did, how it
  was found and fixed, and the sanity check: with the fix off, the rerun
  reproduces the old numbers exactly (`leakfix_rerun.md`). This replaces the
  current Limitation (3), L1131-1141, which covers only the `all_at_once`
  split fix.

---

## Phase 3. New external benchmarks

- [ ] **3.1 New section: PBMC10k Multiome vs LINGER (main external
  benchmark).** Insert it after Baselines (after L754, before
  `\section{Results}` at L756), with `\label{sec:pbmc}`, and reference it in
  the abstract and introduction. Content (all from
  `pbmc10k_linger_benchmark.md`):
  - Data: 10x `pbmc_granulocyte_sorted_10k`, truly paired RNA+ATAC from the
    same nuclei; LINGER's 9,543 labelled cells in 4 cell types; no
    mitochondrial filter (amendment A1b); gene universe of 25,477 genes
    (§7).
  - Ground truth and metric: LINGER's `bm_trans` code (top-1000 Cistrome
    targets per dataset; AUC and AUPR ratio); 20 ChIP datasets, 10 TFs, 4
    cell types; the headline is the mean over the 19 datasets of LINGER
    Table S7.
  - Supervision, pre-registered before any result (§4): CollecTRI labels
    with all 10 evaluation TFs removed as regulators; TF-disjoint splits;
    degree-matched negatives; label-free graphs; negative-pool fix on;
    selection on validation only; 5 seeds. Training is single-tier, so the G1
    choice does not apply here; say so.
  - Results table from §11: `fm_h384` 0.7343 ± 0.0046 / 2.135;
    `fm_h384_rnaonly` 0.7197 ± 0.0301 / 2.178; `base` 0.6701 ± 0.0083 /
    1.722; `base_rnaonly` 0.5438 ± 0.0111 / 1.337; `fm_h384_uniformneg`
    0.5973 ± 0.0136 / 1.698; DoRothEA A-B 0.6578 ± 0.0175 / 1.933; random
    split 0.6686 ± 0.0147 / 1.987; trivial baselines (degree, gene-ID,
    Pearson) 0.50-0.59.
  - Published rows (Table S7): LINGER 0.7143 / 2.2526; SCENIC+ 0.5481 /
    1.2905; PCC 0.5408 / 1.2287; GENIE3 0.5387 / 1.1686; PIDC 0.5292 /
    1.1717. KEGNI 0.699 (its paper).
  - Mark every MeVD-GRN and baseline number as `expressed`-space until
    [P-LINGER], then report both spaces (amendment A3). GRNBoost2
    [P-GB2]. SCENIC+ and scTFBridge are "published, not re-run" (§9).
  - Required honesty points:
    - MeVD-GRN is lower than LINGER on AUPR ratio.
    - The target-disjoint regime fails (0.455 / 0.458).
    - With Geneformer, ATAC's contribution is within noise.
    - Do not write "without external pretraining": MeVD-GRN uses
      CollecTRI supervision and Geneformer, which was pretrained on external
      single-cell data. LINGER uses external bulk data. Name what each
      method uses instead.
  - Also add the PBMC data to Section 4 (Data), the metrics to Section 6
    (Evaluation, L675-696), and LINGER, the published Table S7 methods and
    the trivial baselines to Section 7 (Baselines, L698-754).
- [ ] **3.2 New section: BEAR-GRN (secondary benchmark on SC-MO-GRN-DB
  data).** Insert it after the K562 results (after L946, before
  `\section{Ablation Study}` at L948), with `\label{sec:bear}`. Content
  (from `bear_grn_benchmark.md`):
  - What BEAR-GRN is: 9 methods × 9 datasets; ChIP, KO, union and
    intersection ground truths; same lab as SC-MO-GRN-DB.
  - Our setup: BEAR's own input matrices and ground truths; a Python port of
    BEAR's R scoring, verified on the released GRNs (AUPRC exact in all 10
    checked cases, §8); TF-disjoint 5-fold cross-fitting, the design the BEAR
    authors name in their peer-review file (§6.1); degree-matched negatives;
    label-free graphs.
  - The **hub and coverage artifacts** (§10.1), as a finding in their own
    right. A constant score for every measured pair gets K562 ChIP AUPRC
    0.431, equal to the best published method (LINGER 0.430). A TF-disjoint
    target in-degree ranker beats every BEAR method on K562 ChIP (AUROC
    0.652, AUPRC 0.484). Nothing is above random AUPRC on the K562 KO and
    intersection ground truths.
  - The rule for MeVD-GRN: credit only AUPRC above the coverage baseline,
    require beating the in-degree ranker, and report the sparse-matched
    variant (§6.6, §10.1).
  - MeVD-GRN results [P-BEAR]. The seed-42 numbers in §10.2 are "a
    direction, not a result" and must not be reported as results. Note what
    they suggest: on ChIP below the in-degree ranker, on Union marginally
    above it.
  - Protocol note: BEAR's L1 regime was pre-registered as `all_at_once`
    (§6.2). Keep it whatever G1 decides, and say so.
- [ ] **3.3 Introduction**: add one paragraph on the evaluation strategy
  (PBMC10k as the community benchmark, SC-MO-GRN-DB K562 for the tier
  protocol, BEAR-GRN as the SC-MO-GRN-DB community baseline) and on the
  controls a supervised model needs (TF-held-out splits, degree-matched
  negatives, trivial hub baselines; InfoSEM, Stock *et al.*).

---

## Phase 4. The 18 code-vs-documentation discrepancies

One item per discrepancy in `docs/mevd_vs_scmultiomegrn.md` §8. Where a
discrepancy lives outside `main.tex`, the item says what `main.tex` needs
and where the outside fix stands.

- [ ] **D1. RNA input.** §Modality encoders, L450-453 ("concatenated with a
  50-dimensional co-expression signature $s_g$ ..."). The code uses
  `rna_in_dim: 3`, so $x^{\text{rna}}_g \in \mathbb{R}^3$ only (mean,
  variance, detection). Move the signature sentence to where it is used:
  graph construction (GNN towers, L470-475) and the decoder's
  $\mathrm{coexpr}(i,j)$ term (L500-501). $W_1 \in \mathbb{R}^{64\times3}$ at
  L463 is already correct.
- [ ] **D2. ATAC input and encoder.**
  - L454-455 ($x^{\text{atac}}_g \in \mathbb{R}^2$, "mean and variance"):
    change to $\mathbb{R}^3$, [mean, variance, detection] of RP-weighted gene
    activity.
  - Equation L460 (one Linear + LN + GELU): use the same 2-layer MLP form as
    L459.
  - L463-465 ($W_a \in \mathbb{R}^{128\times2}$): change to 3 → 64 → $d$,
    which agrees with Table `tab:params` L530-531.
  - Write $d$ = 128 (base) or 384 (headline) throughout L463-465 and L497
    ($a_j \in [0,1]^{128}$).
  - Related stale text in §K562 data footprint, L298-301 ("multiplying the
    cell×peak matrix by a peak→gene TSS-proximity incidence matrix (window
    ±100kb)"): change to the RP weighting, i.e. an exponential TSS-distance
    decay (10 kb) within ±100 kb (`mevd_vs_scmultiomegrn.md` §2,
    "Gene/TF features").
- [ ] **D3. Openness term.** L455-456 ("a scalar locus openness feature"),
  L495 ($w_o \cdot \mathrm{open}_j$) and L501-502 ("the locus-openness
  scalar"). Change to the 4-d locus descriptor (mean and max RP of the
  top-10 peaks, mean signed distance, peak count) through `Linear(4,1)`:
  $w_o^\top \mathrm{open}_j$, with $\mathrm{open}_j \in \mathbb{R}^4$.
- [ ] **D4. "Labels cannot leak" through the graphs.** §Dual relational GNN
  towers, L467-470 ("built once from structure only (never from labeled
  positive edges, to avoid leaking label information ...)").
  - The paper's graph excluded every known positive, so graph membership
    implied label 0: 0 val/test positives were edges, against 2.8-3.0% of
    val/test negatives (`labelfree_graph_rerun.md` §1).
  - After [P-LF], describe the label-free construction
    (`prior_exclude_positives: false`; no positive is added or removed), and
    put the leak in the limitations (2.16).
  - The Discussion's graph ablations (L1014-1027, `tf_cand_only`) used the
    leaky graph; label them (2.12).
  - Outside the paper: `docs/citations.md` §8 still says "so labels cannot
    leak" (open; `plan.md` §10 item 12).
- [ ] **D5. Curriculum, hard negatives and replay described as the headline's
  training.** Abstract L61-66; Introduction L164-175, L177-195; Training
  L569-609.
  - (a): these become accurate. State the headline's exact settings (G1) and
    that it is the headline.
  - (b): state that `all_at_once` trains one 45-epoch stage at lr 1e-3 on
    localization ∪ perturbation, with hard negatives disabled
    (`build_all_at_once_stage`) and no replay. Stop presenting the
    curriculum, hard negatives and replay as part of the headline, and apply
    1.1.
- [ ] **D6. Variance claim in the abstract.** L83-84 ("no result in this paper
  yet has a repeated-seed variance estimate"). Delete it, and state the real
  status: K562 headline 5 seeds [P-LF], PBMC 5 seeds, BEAR 5 seeds
  [P-BEAR], scMultiomeGRN and the ablations n = 1.
- [ ] **D7. "Three flags".** §Architectural ablation hooks, L558-564. List
  every switch: `integration` (role-aware / gated / concat), `use_atac`,
  `use_gnn`, `graph_mode` (which prior graphs the GNN sees), `combine_mode`
  (sum / gated), `use_fm`, `use_motif`, `use_edge_mlp`, plus the training
  switch `curriculum.protocol`.
- [ ] **D8. No shipped config reproduces the headline.** Appendix
  "Environment and Reproducibility", L1301-1307. Add the exact command:
  - (a): `scripts/06_ablation.py --ablation full_curriculum --seed N` with
    `configs/sweep/k562_fm_h384_l2_labelfree.yaml`
    (`labelfree_graph_rerun.md` §4), or `configs/sweep/k562_fm_h384_l2.yaml`
    for the negative-fix-only numbers (`leakfix_rerun.md`).
  - (b): `--ablation all_at_once` with the same configs.
  - Outside the paper: making a shipped default reproduce it is `plan.md`
    §10 item 1.
- [ ] **D9. Every reported number predates the fix.** This is all of
  Phase 2.
- [ ] **D10. Bilinear score attribution.** Decoder, L494 and L499-500
  ("\text{raw} is an asymmetric bilinear TF→target compatibility score").
  - Cite RESCAL (`nickel2011rescal`, new; Phase 6) here and at L561-562
    ("a classical asymmetric bilinear decoder").
  - Do not cite DistMult for it: DistMult's diagonal $W$ makes the score
    symmetric.
  - Outside the paper: fixed in `docs/citations.md` §4 on 2026-10-08.
- [ ] **D11. Duren *et al.* DOI.** Fixed in `paper/references.bib` on
  2026-10-08: the `duren2018integrative` entry now holds the 2017 PECA paper
  its title names (PNAS 114(25):E4914-E4923, doi 10.1073/pnas.1704553114),
  and the bib carries an `@comment` explaining the change.
  - In `main.tex` (L141): optionally rename the key to `duren2017peca`,
    changing the bib key at the same time.
  - Also cite PECA in the decoder paragraph (L488-490), as the biological
    basis for gating (`mevd_vs_scmultiomegrn.md` C12).
  - Outside the paper: also fixed in `docs/citations.md` §2 and §6.
- [ ] **D12. Gated relation combiner, motivation and result.** Table row
  L983, Discussion L1034-1043, Limitation (5) L1154-1161.
  - Wherever the combiner is discussed, cite HAN (`wang2019han`, new) as the
    closest precedent; `main.tex` does not currently cite the Relational
    Graph Transformer.
  - Add the h384 + FM result: near-uniform weights (0.43-0.58) and the same
    accuracy (dual AUPR 0.962 vs 0.963), pre-fix and n = 1
    (`mevd_vs_scmultiomegrn.md` C14).
  - Outside the paper: `docs/citations.md` §17 is still open.
- [ ] **D13. Motif relation.**
  - The figure caption at L547-550 ("up to three shared structural graphs
    ... and a real TF-motif graph") also contradicts the figure itself:
    `architecture_diagram_v3_1` does not draw the motif relation
    (`docs/README.md`, figures).
  - At L1039-1041 and L1154-1161, say that the motif scan is restricted to
    TF-candidate pairs, so the relation is a *subset* of the TF-candidate
    relation, and that it is off in every reported run.
  - If `motif_graph_only` is mentioned (`results.md` §2: 0.724 / 0.818 /
    0.955), label it sequential, pre-fix, n = 1.
  - Outside the paper: `docs/citations.md` §19 and `plan.md` §2/§4 are
    still open.
- [ ] **D14. Replay "on by default", half the negatives hard.** Training
  L602-609 (hard negatives) and L611-616 (replay).
  - (a): true for the headline; say "in the headline configuration".
  - (b): mark both as unused by the headline.
  - Outside the paper: `docs/citations.md` §10 and §12 are still open.
- [ ] **D15. K562 RNA and ATAC are not cell-paired.** §K562 data footprint,
  L296-304.
  - Add that the RNA (DS019, 953 cells) and ATAC (DS025 multiome, 1,126
    cells) come from different datasets and are not cell-paired, and that
    the model uses only gene-level summaries (`configs/k562.yaml`).
  - Never call the K562 input "paired". (Abstract L56-57 is about the field
    and is fine.)
  - Contrast with PBMC10k (3.1), which is truly paired.
  - Outside the paper: the `docs/reference/literature_knowledge_base.md`
    header is still open.
- [ ] **D16. scMultiomeGRN called "a supervised graph neural network".**
  Related Work L218-222. Rewrite with `mevd_vs_scmultiomegrn.md` §7: it
  adopts GraFRank's modality-specific neighbour aggregation and cross-modal
  attention, predicts undirected TF-TF links, and is trained
  semi-supervised on motif-derived labels from its own scATAC.
- [ ] **D17. scMultiomeGRN's paper vs its code.** Baselines L729-744.
  - State which version our adapter implements where the paper and the
    official code differ: decoder σ(hhᵀ) in Eq. 11 vs the code's
    symmetrised concatenation-MLP; FIMO thresholds. Check
    `src/baselines/scmultiomegrn_wrapper.py` before writing.
  - List the adaptations (R9): our pseudobulk features instead of MAESTRO
    per-cell RP and GRNBoost2 vectors; neighbours capped at 500 per TF; a
    200k training-edge cap per epoch; a symmetric decoder on a directed
    task; single seed; not converged at 2,000 epochs.
- [ ] **D18. results.md `rna_only` 0.856 / 0.962 vs the JSON's 0.859 /
  0.963.** No `main.tex` change: L977, L1022 and L1265 already use the
  JSON's 0.859 / 0.963. `results.md` §2 was annotated on 2026-10-08.

---

## Phase 5. Further errors found while writing this plan

These are not in the list of 18.

- [ ] **E1. The `with_fm` row mixes two runs.** Table `tab:ablation` L986,
  text L1067-1068 and appendix L1290-1292 label it "3.26M params".
  - Its localization and perturbation (0.751 / 0.771, 0.834 / 0.964) come
    from `results/ablations/with_fm_h384l2_K562.json`.
  - Its dual-evidence AUPR / AUROC / EPR (0.949 / 0.988 / 5.2) match the
    earlier `results/ablations/with_fm_K562.json` (0.9488 / 0.9884 /
    5.2152).
  - The h384/l2 file has dual 0.9628 / 0.9916 / 5.333, i.e. the 0.963 /
    0.992 already in `tab:head2head-alt` L930.
  - Use one run consistently (checked against both JSONs on 2026-10-08).
    All of these are pre-fix.
- [ ] **E2. Abstract L75 "0.950 vs 0.847".** The 5-seed mean is 0.9517 (L810);
  0.950 is an older single-run value (`plan.md` §3b). It is superseded by
  Phase 2 anyway.
- [ ] **E3. Self-contradictory sentence**, L1048-1051: "`all_at_once`
  dominates on localization (0.938 vs 0.568) and dual-evidence zero-shot
  generalization (0.798 vs 0.872 -- actually *worse* here ...)". Rewrite it:
  at 315K params `all_at_once` wins only localization.
- [ ] **E4. The transfer study is called future work**, L288-293 ("reserved
  for the zero-shot cross-cell-type transfer study outlined as future
  work"). The Macrophage/MCF7 transfer and joint-training results exist
  (`results.md` §6; Limitation (6) L1163-1171). Make the two passages agree,
  with the 2.15 caveats.
- [ ] **E5. Stale environment text**, L1305-1306: "All experiments were run on
  a single consumer GPU (RTX 4060); the full model has under 400K
  parameters". The headline has 3,257,479 params. The size sweep, the
  leak-fix 5-seed rerun, the scMultiomeGRN baseline and the PBMC grid ran on
  Ada (RTX 2080 Ti). The BEAR seed-42 run and the label-free rerun ran on the
  RTX 4060.

---

## Phase 6. Citations

**Fixes:**

- [ ] **C1. Duren *et al.***: done in the bib (D11). Optional key rename.
- [ ] **C2. RESCAL for the bilinear score**: add the entry and cite it (D10).
- [ ] **C3. HAN for the relation combiner**, if the combiner stays in the
  paper (D12).

**New references.** Every DOI below was checked against the Crossref API
on 2026-10-08, and the arXiv IDs against the arXiv API. Volume, issue and
pages are as Crossref returns them. Bib keys are suggestions.

| Key (suggested) | Reference | DOI / ID | Use in the paper |
|---|---|---|---|
| `nickel2011rescal` | Nickel M, Tresp V, Kriegel H-P. A three-way model for collective learning on multi-relational data. ICML 2011 | no DOI (none in Crossref); https://icml.cc/2011/papers/438_icmlpaper.pdf (title checked from the PDF) | decoder (D10) |
| `yuan2025linger` | Yuan Q, Duren Z. Inferring gene regulatory networks from single-cell multiome data using atlas-scale external data. *Nat Biotechnol* 43(2):247-257 (2025; online 12 Apr 2024) | 10.1038/s41587-024-02182-7 | PBMC bar, related work (replay), BEAR |
| `li2025kegni` | Li P, Li L, Nan J, Chen J, Sun J, Cao Y. KEGNI: knowledge graph enhanced framework for gene regulatory network inference. *Genome Biol* 26:294 (2025) | 10.1186/s13059-025-03780-7 | PBMC (0.699) |
| `karamveer2026beargrn` | Karamveer K, Moeller E, Valensi H, Manful E-E, Uzun Y. BEAR-GRN: Systematic assessment of single-cell multi-omics-based gene regulatory network inference methods. *Nat Commun* (2026, published 18 Sep 2026) | 10.1038/s41467-026-77838-w (Crossref has no volume or article number yet) | BEAR section |
| `cui2025infosem` | Cui T, Xu S-J, Moskalev A, Li S, Mansi T, Prakash M, Liao R. InfoSEM: a deep generative model with informative priors for gene regulatory network inference. ICML 2025 | arXiv:2503.04483 (no Crossref DOI) | hub-bias controls |
| `stock2025hidden` | Stock M, Ratajczak F, Bertin P, Hoermanseder E, Bengio Y, Hartford J, *et al.* Hidden sampling biases inflate performance in gene regulatory network inference. bioRxiv (posted 23 Dec 2025) | 10.64898/2025.12.19.695616 | degree-matched negatives |
| `yilmaz2025bias` | Yılmaz S, Yorgancioglu K, Koyutürk M. Bias-aware training and evaluation of link prediction algorithms in network biology. *PNAS* 122(24):e2416646122 (2025) | 10.1073/pnas.2416646122 | degree baselines |
| `mullerdott2023collectri` | Müller-Dott S, Tsirvouli E, Vazquez M, Ramirez Flores RO, Badia-i-Mompel P, Fallegger R, *et al.* Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities (CollecTRI). *Nucleic Acids Res* 51(20):10934-10949 (2023) | 10.1093/nar/gkad841 | PBMC training labels |
| `garciaalonso2019dorothea` | Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. Benchmark and integration of resources for the estimation of human transcription factor activities (DoRothEA). *Genome Res* 29(8):1363-1375 (2019) | 10.1101/gr.240663.118 | PBMC secondary labels |
| `zheng2019cistrome` or `taing2024cistrome` | Zheng R, Wan C, Mei S, *et al.* Cistrome Data Browser: expanded datasets and new tools for gene regulatory analysis. *Nucleic Acids Res* 47(D1):D729-D735 (2019). / Taing L, Dandawate A, L'Yi S, Gehlenborg N, Brown M, Meyer CA. Cistrome Data Browser: integrated search, analysis and visualization of chromatin data. *Nucleic Acids Res* 52(D1):D61-D66 (2024) | 10.1093/nar/gky1094 / 10.1093/nar/gkad1069 | PBMC ground truth. Cite the version LINGER cites (not checked yet) |
| `kommu2025scregnet` | Kommu S, Wang Y, Wang Y, Wang X. Prediction of gene regulatory connections with joint single-cell foundation models and graph-based learning (scRegNet). *Bioinformatics* 41(Supplement_1):i619-i627 (2025) | 10.1093/bioinformatics/btaf217 | FM + GNN precedent |
| `zhang2023scmtni` | Zhang S, Pyne S, Pietrzak S, *et al.* Inference of cell type-specific gene regulatory networks on cell lineages from single cell omic datasets (scMTNI). *Nat Commun* 14:3064 (2023) | 10.1038/s41467-023-38637-9 | intersection gold standards |
| `mccalla2023` | McCalla SG, Fotuhi Siahpirani A, Li J, *et al.* Identifying strengths and weaknesses of methods for computational network inference from single-cell RNA-seq data. *G3* 13(3):jkad004 (2023) | 10.1093/g3journal/jkad004 | intersection gold standards |
| `wang2025sctfbridge` | Wang F-a, Yi C, Chen J, He R, Liu J, Li Y, *et al.* scTFBridge: a disentangled deep generative model informed by TF-motif binding for gene regulation inference in single-cell multi-omics. *Nat Commun* 16:9166 (2025) | 10.1038/s41467-025-64227-y | PBMC (published STAT1 number) |
| `xi2025regx` | Xi X, Li J, Jia J, Meng Q, Li C, Wang X, *et al.* A mechanism-informed deep neural network enables prioritization of regulators that drive cell state transitions (regX). *Nat Commun* 16:1284 (2025) | 10.1038/s41467-025-56475-9 | PBMC benchmark users |
| `sankar2021grafrank` | Sankar A, Liu Y, Yu J, Shah N. Graph neural networks for friend ranking in large-scale social platforms (GraFRank). WWW 2021, pp. 2535-2546 | 10.1145/3442381.3450120 | related work (D16), encoders |
| `li2022deeptfni` | Li H, Sun Y, Hong H, *et al.* Inferring transcription factor regulatory networks from single-cell ATAC-seq data based on graph neural networks (DeepTFni). *Nat Mach Intell* 4(4):389-400 (2022) | 10.1038/s42256-022-00469-5 | scMultiomeGRN's label construction |
| `kamimoto2023celloracle` | Kamimoto K, Stringa B, Hoffmann CM, Jindal K, Solnica-Krezel L, Morris SA. Dissecting cell identity via network inference and in silico gene perturbation (CellOracle). *Nature* 614:742-751 (2023) | 10.1038/s41586-022-05688-9 | accessibility masks; BEAR method |
| `bravogonzalezblas2023scenicplus` | Bravo González-Blas C, De Winter S, Hulselmans G, *et al.* SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks. *Nat Methods* 20(9):1355-1367 (2023) | 10.1038/s41592-023-01938-4 | PBMC and BEAR baseline |
| `chan2017pidc` | Chan TE, Stumpf MPH, Babtie AC. Gene regulatory network inference from single-cell data using multivariate information measures (PIDC). *Cell Syst* 5(3):251-267.e3 (2017) | 10.1016/j.cels.2017.08.014 | PBMC published baseline (Table S7) |
| `wang2020maestro` | Wang C, Sun D, Huang X, *et al.* Integrative analyses of single-cell transcriptome and regulome using MAESTRO. *Genome Biol* 21:198 (2020) | 10.1186/s13059-020-02116-x | RP features (shared with scMultiomeGRN) |
| `chen2022genelink` | Chen G, Liu Z-P. Graph attention network for link prediction of gene regulations from single-cell RNA-sequencing data (GENELink). *Bioinformatics* 38(19):4522-4529 (2022) | 10.1093/bioinformatics/btac559 | supervised TF→target precedent |
| `pratapa2020beeline` | Pratapa A, Jalihal AP, Law JN, Bharadwaj A, Murali TM. Benchmarking algorithms for gene regulatory network inference from single-cell transcriptomic data (BEELINE). *Nat Methods* 17(2):147-154 (2020) | 10.1038/s41592-019-0690-6 | EPR convention (L691-692) |
| `wang2019han` | Wang X, Ji H, Shi C, *et al.* Heterogeneous graph attention network (HAN). WWW 2019, pp. 2022-2032 | 10.1145/3308558.3313562 | relation combiner (D12) |
| `dauphin2017glu` | Dauphin YN, Fan A, Auli M, Grangier D. Language modeling with gated convolutional networks (GLU). ICML 2017 | arXiv:1612.08083 | sigmoid gate in the decoder |
| `saito2015prc` | Saito T, Rehmsmeier M. The precision-recall plot is more informative than the ROC plot when evaluating binary classifiers on imbalanced datasets. *PLoS ONE* 10(3):e0118432 (2015) | 10.1371/journal.pone.0118432 | AUPR as primary metric (L678-681) |
| `kipf2016vgae` | Kipf TN, Welling M. Variational graph auto-encoders. arXiv (2016) | arXiv:1611.07308 | link-prediction framing (L125-131) |
| `qi2026universal` | Qi J, Li H, Cui Y, Zheng Y, Huang J. Towards universal gene regulatory network inference: unlocking generalizable regulatory knowledge in single-cell foundation models. arXiv (2026) | arXiv:2605.08128 | FM transfer finding (counterpoint) |
| `hegde2026grnformer` | Hegde A, Cheng J. GRNFormer: accurate gene regulatory network inference using graph transformer. *Bioinformatics* 42(4):btag144 (2026) | 10.1093/bioinformatics/btag144 | cross-cell-type test precedent |
| `fleck2023pando` | Fleck JS, Jansen SMJ, Wollny D, *et al.* Inferring and perturbing cell fate regulomes in human brain organoids (Pando). *Nature* 621:365-372 (2023; online 2022) | 10.1038/s41586-022-05279-8 | BEAR method |
| `kartha2022figr` | Kartha VK, Duarte FM, Hu Y, *et al.* Functional inference of gene regulation using single-cell multi-omics (FigR). *Cell Genomics* 2(9):100166 (2022) | 10.1016/j.xgen.2022.100166 | BEAR method |
| `zhang2022directnet` | Zhang L, Zhang J, Nie Q. DIRECT-NET: an efficient method to discover cis-regulatory elements and construct regulatory networks from single-cell multiomics data. *Sci Adv* 8(22):eabl7393 (2022) | 10.1126/sciadv.abl7393 | BEAR method |
| `jiang2022tripod` | Jiang Y, Harigaya Y, Zhang Z, Zhang H, Zang C, Zhang NR. Nonparametric single-cell multiomic characterization of trio relationships between transcription factors, target genes, and cis-regulatory regions (TRIPOD). *Cell Syst* 13(9):737-751.e4 (2022) | 10.1016/j.cels.2022.08.004 | BEAR method |
| `kamal2023granie` | Kamal A, Arnold C, Claringbould A, *et al.* GRaNIE and GRaNPA: inference and evaluation of enhancer-mediated gene regulatory networks. *Mol Syst Biol* 19(6):e11627 (2023) | 10.15252/msb.202311627 | BEAR method |
| `badiaimompel2024greta` (optional) | Badia-i-Mompel P, Casals-Franch R, Wessels L, Müller-Dott S, Trimbour R, Yang Y, *et al.* Comparison and evaluation of methods to infer gene regulatory networks from multimodal single-cell data (GRETA). bioRxiv (2024) | 10.1101/2024.12.20.629764 | "benchmarks rank methods differently" |
| (no key; footnote or data statement) | 10x Genomics, "PBMC from a healthy donor, granulocytes removed through cell sorting (10k)", Cell Ranger ARC 2.0.0 | no DOI; URL in `pbmc10k_linger_benchmark.md` §1 | PBMC data |

The existing entries `valensi2026scmogrndb`, `xu2025scmultiomegrn`,
`theodoris2023geneformer`, `moerman2019grnboost2` and `huynhthu2010genie3`
stay. Also cite the scMultiomeGRN code (Zenodo, doi 10.5281/zenodo.14848389;
a DataCite DOI, not in Crossref) where the adapter is described.

---

## Phase 7. Finishing

- [ ] **F1. Regenerate the figures** flagged in the old editorial note
  (fig1, fig4, fig5, fig6, fig9) from fixed runs, and add one PBMC figure
  (per-dataset AUROC vs LINGER) and one BEAR figure (methods vs the coverage
  and in-degree baselines). Figure scripts live in `scripts/`, the other
  session's lane, so coordinate first.
- [ ] **F2. Check every number against its source doc** in the table at the
  top, and remove any number still carrying a pending tag.
- [ ] **F3. Build the paper** (`paper/README.md`) and check that every new
  `\cite` resolves.
- [ ] **F4. Update `results.md` and `plan.md` §10** to say the paper was
  resynced, with the date.
