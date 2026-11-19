# The Marginal Dominance Hypothesis for GRN inference

Draft framing for the paper, 2026-11-14. It is built from every MeVD-GRN experiment since July 2026: the
SC-MO-GRN-DB tiers, PBMC10k/Cistrome, BEAR-GRN and the scMultiomeGRN benchmark. It borrows its structure from
Jiang et al., *What Makes a Representation Good for Single-Cell Perturbation Prediction?* (ICML 2026,
arXiv 2605.19343).

Only numbers that come from result files are quoted. The archived zero-shot transfer numbers (0.607 → 0.447) are
**not** used: their no-FM and FM result files are byte-identical, and the runs predate the leak fixes.

---

## 1. The analogy

Jiang et al. argue that in single-cell perturbation data, perturbation-invariant background programs dominate
expression, while perturbation-specific signal is sparse. Methods that do not separate the two fail in one of two
ways:
- causal representation learning lets the invariant part leak into the perturbation latent;
- foundation models suppress the perturbation signal.

Their fix (PerturbedVAE) splits the latent into an invariant block and a responsive block. It aligns the invariant
block with matched **control cells**, so that the responsive block only has to explain what differs.

GRN inference has the same structure. The roles map as follows:

| Perturbation Suppression Hypothesis | GRN analogue (this document) |
|---|---|
| Invariant background programs `z_ι` dominate expression | **Node marginals** dominate edge labels: how "big" a regulator is, and how accessible, expressed, measured or hub-like a target is |
| Sparse perturbation-specific signal `z_ν` | Sparse **pair-specific** regulation: does *this* TF regulate *this* gene, beyond both genes' marginals? |
| CRL entangles the invariant part into `z_ν` | Supervised GNNs and link predictors write hubness into edge scores (entanglement) |
| Foundation models suppress perturbation signal | Pretrained gene embeddings (Geneformer) encode gene identity, family and vocabulary, which is context-invariant and suppresses context-specific edges |
| Label-free baselines are confounded by background | Co-expression and co-accessibility across cells are dominated by cell-state programs, so label-free statistics are at chance |
| Control cells as contrastive anchors | **Marginal-matched non-edges** (same TF, same target degree; or 2x2 tetrads) as anchors |
| PerturbedVAE: invariant / responsive split + alignment | A model with an explicit marginal head plus an interaction head trained by **conditional likelihood**, which cancels the marginals exactly (section 5) |

## 2. Statement

Let `Y_ijc ∈ {0,1}` be the ground-truth label for regulator `i`, target `j`, context `c`. Write its log-odds, and
any method's score, as a two-way decomposition:

```
logit P(Y_ij = 1) = α_i + β_j + γ_ij          (ground truth)
s_ij              = a_i + b_j + g_ij          (method; additive on a logit or rank scale, g double-centred)
```

- `α_i` is the regulator marginal: ChIP depth, motif promiscuity, TF out-degree.
- `β_j` is the target marginal: accessibility, expression, being measured, being in the Geneformer vocabulary,
  in-degree.
- `γ_ij` is the pair-specific part. Its rows and columns are centred, so it carries no node-level information.

**Marginal Dominance Hypothesis (MDH).**
1. In the ground truths used by current GRN benchmarks, the node marginals `α + β` carry most of the predictable
   label signal. The label matrix is close to rank one in log-odds.
2. Supervised training on random non-edges, and generic pretrained gene representations, both encode the marginals
   first. Marginals are the cheapest and most transferable signal, and they are shared across TFs and contexts.
   The pair-specific term `g` is therefore either entangled with the marginals or suppressed.
3. Pooled edge-level metrics (AUROC, AUPR, EPR over all pairs) rank methods by how well `a + b` matches `α + β`,
   not by how well `g` matches `γ`.

**Consequence.** Methods that do not explicitly separate marginals from interactions fail systematically on any
evaluation that controls the marginals: stratified or matched tests, held-out nodes, context-differential edges.
Their pooled-metric wins are matched by marginal-only baselines.

## 3. The evidence we already have, read through MDH

| Observation (source) | What MDH says it is |
|---|---|
| Logistic regression on the two TFs' training link counts: 0.931 AUPR, vs scMultiomeGRN 0.919 and MeVD-GRN 0.935 (deck §6) | Marginal-only model ≈ the best method on a pooled metric (claim 1, 3) |
| Floor model + Geneformer 0.931; ATAC tower, GNN and graphs each change AUPR by ≤ 0.001 (deck §9) | The representation learns marginals; the pair-specific machinery adds nothing (claim 2) |
| Caliper-matched (marginals fixed): MeVD 0.560, SVD 0.567, scMultiomeGRN 0.524, all heuristics 0.48–0.57 (deck §10) | Once `α + β` is controlled, little `γ` is recovered by anyone |
| Held-out TFs: one unseen 0.753 vs the link-count rule 0.779; both unseen ≤ 0.575 for every model; caliper 0.487 (deck §11) | Performance tracks whether a node's marginal is knowable (prediction P4) |
| BEAR-GRN: TF-disjoint target in-degree ranker beats all 9 published methods on ChIP AUPRC 9/9, and MeVD on 7/9 (deck §12) | The target marginal `β` dominates ChIP |
| Ranking by a TF's GT out-degree gives AUROC 0.87–0.93 on BEAR ChIP (selfsup analysis) | The regulator marginal `α` dominates ChIP too |
| "All measured pairs" equals the best published ChIP AUPRC (K562 0.431 vs LINGER 0.430) (bear_grn §10.1) | Being measured is a target marginal |
| PBMC Cistrome: TF-agnostic ATAC-activity ranker 0.814 > MeVD 0.734 > LINGER 0.714 (pbmc §12.1) | Cistrome RP labels ≈ target accessibility marginal (structurally the same RP transform) |
| PBMC own-TF row 0.690 vs another TF's row 0.685; TF-swap gaps ≤ 0.033 for every method audited (pbmc §12.4) | `g ≈ 0`: per-TF scores are the shared marginal ranking |
| Geneformer token presence (0/1) alone: AUROC 0.6955; no-token genes are 0.1% positive (pbmc §12.1–12.2) | A vocabulary marginal inside the foundation model |
| ATAC adds +0.126 AUROC without Geneformer, +0.015 with it (pbmc §11) | ATAC and Geneformer encode the same target marginal, so they are redundant |
| Target-disjoint split 0.455 (below chance) → 0.745 after showing never-labelled genes as negatives; a TF-blind gene prior reproduces both numbers (0.409 / 0.770) (pbmc §12.2–12.3) | The model's output is a per-target marginal whose sign is set by the negative design |
| Degree-matched training negatives: +0.137 on PBMC; uniform negatives better on BEAR inner-val (+0.019) (pbmc §11, bear §14) | The negative design selects which marginal is learned, and the "best" design flips with whichever marginal the benchmark rewards (P5) |
| Every label-free statistic is at 0.50–0.54 on scMultiomeGRN (deck §4) | Covariation across cells is driven by cell-state background, the literal PSH mechanism |
| Within in-degree strata, MeVD keeps ChIP 0.530–0.573 (human), KO 0.516–0.657; Intersection 0.28 (bear §15, §10.3) | A small `γ` exists, larger for functional (KO) ground truth; Intersection is anti-marginal (P6, P7) |
| ChIP and KO targets overlap 1.3%; KO targets are not enriched among measured genes (selfsup §1.3, bear §10.1) | Functional ground truths are less marginal-dominated than binding ones |

## 4. Implicit assumptions MDH says the field gets wrong

| # | Assumption usually treated as obvious | Why it fails | Our evidence |
|---|---|---|---|
| A1 | Pooled edge-level AUROC/AUPR measures edge-level regulation | It mostly measures node marginals | link-count LR ≈ best; in-degree > all BEAR methods |
| A2 | Random non-edges are valid negatives | They differ from positives mostly in marginals, so the task becomes degree prediction (Stock et al. 2025 make this point for evaluation) | caliper tests |
| A3 | Ground truths are pair-specific | ChIP ≈ accessibility × TF ChIP depth; motif-in-open-promoter ≈ motif promiscuity × openness; Cistrome RP top-1000 ≈ peak density near the gene | TF out-degree 0.87–0.93; ATAC ranker 0.814 |
| A4 | A TF-disjoint split is sufficient for honest evaluation | Target marginals transfer across TFs | in-degree ranker on BEAR |
| A5 | Foundation-model gene embeddings add regulatory knowledge | They add identity, family and vocabulary marginals, the same in every context | floor+GF ≈ link-count LR; token-presence 0.6955 |
| A6 | More modalities give more regulatory signal | ATAC enters as a target marginal, redundant with the FM | +0.126 → +0.015 |
| A7 | Co-variation across cells contains regulation | Cell-state programs dominate it | label-free stats at chance |
| A8 | Dense and sparse outputs are comparable under AUPR | Coverage is itself a marginal | "all measured pairs" = LINGER |
| A9 | The best negative-sampling scheme is a tuning choice | It decides which marginal the model learns | PBMC vs BEAR sign flip |

## 5. Method that targets the hypothesis directly

Working name **MID-GRN (Marginal–Interaction Decomposition)**. It can be built as a MeVD-GRN revision; most of the
input pipeline is reused.

1. **Marginal head `M_ij = α̂_i + β̂_j`.**
   - Empirical-Bayes shrunk training degrees of `i` and `j`.
   - For unseen nodes, the degree is predicted from node features: ATAC gene activity, expression, detection, the
     Geneformer embedding.
   - By construction this is at least as good as the in-degree, link-count and ATAC rankers. They become special
     cases instead of rivals.
2. **Interaction head `g_ij`, from pair-specific evidence only.**
   - Cis evidence: TF `i`'s motif in peaks linked to gene `j` (peak–gene co-accessibility).
   - Background-removed covariation (the GRN counterpart of PSH's invariant block): TF activity (chromVAR deviation
     or expression) against target expression after regressing both on the top cell-state programs or metacell
     factors.
   - Geneformer only through `e_i ⊙ e_j`, after projecting out the directions that predict degree.
   - `g` is double-centred inside each mini-batch, so it cannot carry node-level information.
3. **Loss that cancels marginals exactly.**
   - Use a **tetrad conditional logit** (Graham 2017, *Econometrica* 85(4); the Rasch / β-model conditional
     likelihood, Andersen 1970).
   - For a 2x2 tetrad `{i,k} × {j,l}` with one alternating pattern observed (`Y_ij = Y_kl = 1`,
     `Y_il = Y_kj = 0`, or the reverse), the conditional probability depends only on `g_ij + g_kl − g_il − g_kj`.
     Every `α` and `β` cancels.
   - A simpler variant: softmax over matched sets (one TF, targets in the same in-degree stratum). This cancels
     `α_i` exactly and `β_j` approximately.
   - This plays the role of PerturbedVAE's control-cell alignment, with marginal-matched non-edges as the controls.
4. **Score `s = M + λ g`**, with `λ` chosen on validation **within-stratum** AUROC, not pooled AUPR.
5. **Theory hook**, the counterpart of Jiang's Theorem 1:
   - Under the two-way logistic model, the conditional likelihood is free of the node effects, so `γ` is
     identified up to double-centring.
   - Pooled BCE with random negatives has gradients dominated by the marginal terms; `g` is not identified and
     absorbs marginal error.
   - This is a short, standard proposition; check the exact conditions in Graham (2017) before stating it.

**Why this can beat the baselines.**
- The marginal head equals or beats every marginal baseline on pooled metrics.
- Any `g` with within-stratum AUROC > 0.5 then adds on top of it.
- MeVD-GRN already shows within-stratum ChIP 0.53–0.57 and KO up to 0.66. It loses to in-degree only because its
  marginal is a noisy copy of the in-degree ranker (bear §10.3).

**Honesty rule.** The headline is the within-stratum (double-centred) metric. Pooled wins are reported as a
consequence, never as the claim.

## 6. Field-testable predictions, with what would falsify them

| ID | Prediction | Test | Falsified if | Cost |
|---|---|---|---|---|
| P1 | Ground truths are near-additive in log-odds | Fit the two-way logistic model on training folds (scMultiomeGRN, BEAR 9 × 4 GTs, PBMC 19) and compare it with the best method | On ≥ half the datasets, a published method beats it by > 0.05 AUROC | CPU, hours |
| P2 | Published methods' scores are marginal-dominated | Double-centred / within-stratum AUROC and TF-swap gap for all 9 BEAR released networks (on disk in `~/.cache/bear/data/INFERRED.GRNS`) and LINGER on PBMC | Any method ≥ 0.60 within-stratum on most datasets | CPU, a day |
| P3 | Pretrained gene embeddings encode marginals, not interactions | Linear probes from Geneformer (and scGPT/UCE if cheap) to log degree / accessibility (expect high R²) vs to the pair residual `γ` on held-out TFs (expect ≈ 0) | Probe of `γ` comparable to the probe of degree | CPU/GPU, a day |
| P4 | Generalisation follows marginal knowability | Across splits and datasets, model AUROC vs the two-way model's AUROC on the same split | Correlation < 0.6 | Existing scores |
| P5 | The negative design trades pooled against within-stratum performance | Uniform vs degree-matched vs tetrad training on BEAR and scMultiomeGRN; report both metrics | The pooled-best design is also within-stratum-best | GPU, days |
| P6 | Pair-specific signal is larger in functional and context-differential GTs | `γ`-share (P1 residual) for KO vs ChIP; on BEAR CellTypeExclusive and across the 9 lung cell types, does a model's Δscore across contexts track ΔGT? FM-based models should not | ChIP `γ`-share ≥ KO `γ`-share | CPU |
| P7 | Intersection (binding ∧ function) is anti-marginal | Across BEAR methods, correlation of ChIP AUROC with Intersection AUROC should be negative | Positive correlation | CPU, minutes |
| P8 | Removing marginals in training recovers more `γ` | MID-GRN within-stratum AUROC vs MeVD-GRN and published methods, same splits | No gain over MeVD's 0.53–0.57 | GPU, 2–3 weeks |

## 7. Order of work (cheapest first; P1, P2, P7 need no training)

1. **Week 1 (laptop CPU, BEAR files already local).**
   - Run P1, P2 and P7 on BEAR.
   - Composition check on the existing BEAR outputs: in-degree + `λ` × (MeVD within-stratum residual), scored with
     the port. If this beats in-degree on pooled ChIP, the method recipe already has a working prototype.
2. **Week 1–2.**
   - P1 and P4 on scMultiomeGRN (saved score matrices on Ada).
   - P3 probes (Geneformer cache).
   - P6 on the 9 lung cell types and BEAR CellTypeExclusive.
3. **Weeks 2–5.** Implement the tetrad / matched-set interaction head; run P5 and P8 on BEAR and scMultiomeGRN.
4. **Decision gate (mid December).**
   - If P2 holds across published methods and P8 shows a within-stratum gain, write a hypothesis + method paper.
   - If not, write a hypothesis paper: the predictions tested across all published methods are the contribution,
     and the conditional-likelihood estimator is a diagnostic.

## 8. Novelty boundary and risks

- **Stock et al. (bioRxiv, Dec 2025)** show degree bias in GRN evaluation negatives. MDH goes further in four ways:
  1. It places the bias in the ground truths themselves, in pretrained representations and in label-free
     statistics, not only in negative sampling.
  2. It predicts generalisation behaviour (P4) and the trade-off from negative design (P5).
  3. It explains why ChIP and KO differ (P6, P7).
  4. It removes the marginals in **training** with an estimator that has an identification argument.

  Read Stock et al.'s full text to confirm they do not use conditional likelihood or double-centring.
- **Risk: `γ` may be tiny in every current ground truth.** Then no method can win within strata. The paper would
  then say that current GRN benchmarks cannot measure pair-specific regulation, which P1 and P6 support directly.
- **Risk: tetrad sampling may be sparse.** Tetrads are rare in very sparse KO ground truths. The matched-set
  variant is the fallback.
- **Rule:** before quoting any number here, check it against the result JSON. The archive showed write-ups
  drifting from the artifacts several times.

## References

- Jiang W. et al. What Makes a Representation Good for Single-Cell Perturbation Prediction? ICML 2026.
  arXiv:2605.19343.
- Stock M. et al. Hidden sampling biases inflate performance in gene regulatory network inference. bioRxiv (2025).
  doi:10.64898/2025.12.19.695616.
- Graham B. S. An econometric model of network formation with degree heterogeneity. *Econometrica* 85(4) (2017).
- Andersen E. B. Asymptotic properties of conditional maximum-likelihood estimators. *JRSS-B* 32 (1970).
- Chatterjee S., Diaconis P., Sly A. Random graphs with a given degree sequence. *Ann. Appl. Probab.* (2011)
  (β-model).
- Park Y., Marcotte E. M. Flaws in evaluation schemes for pair-input computational predictions. *Nat. Methods* 9
  (2012).
- BEAR-GRN, Nat. Commun. (2026); scMultiomeGRN, NAR (2025); LINGER (Cistrome protocol).
