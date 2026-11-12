# MEvD-GRN — Critical Pre-Submission Review (Reviewer 2)

Scope: independent read of the codebase, docs, and paper draft, with claims
checked against actual code and data where feasible (not just against the
prose in `plan.md`/`results.md`). File paths and numbers below are cited
directly from the repository as of this review; several were verified by
re-running the relevant code against `data/processed/K562` and
`data/splits/` rather than trusting the docstrings.

---

## 1. Summary / Verdict

MEvD-GRN is a competent, honestly-documented applied ML project with real
engineering care (the team found and fixed two genuine bugs, built a
hierarchy-aware split protocol, and reports negative/mixed results rather
than hiding them). That honesty is the project's biggest asset and should be
preserved. But the work is **not currently submittable** to either target
venue in its present state, for reasons that go beyond polish:

- The paper draft (`paper/main.tex`) reports pre-bugfix numbers that the
  project's own docs say are wrong (Section 8 below) — this alone blocks
  submission regardless of anything else.
- The single most important scientific claim in the paper — a "formal
  non-leakage guarantee" (Proposition 1) licensing the zero-shot
  dual-evidence evaluation — **does not actually hold by construction** for
  the training protocol (`all_at_once`) that `plan.md`/`results.md` now
  recommend as the paper's headline configuration. It happens to hold for
  the current K562 data by a coincidence of implementation details, not by
  the mechanism the paper describes (Section 2 — this is the top finding of
  this review).
- Across the whole project timeline, the "never trained on" dual-evidence
  tier has been used repeatedly as the deciding metric for which
  architecture, curriculum protocol, and model size to keep. That is
  model/protocol selection on the test set in every meaningful sense except
  literal gradient exposure (Section 3).
- Every number in the project is a single run, single seed. Differences of
  0.01–0.05 AUPR are repeatedly treated as decisive (e.g., "beats
  scMultiomeGRN," "ties within noise") with no variance estimate anywhere
  (Section 4).
- The ML methodological novelty is thin for ICLR/ICML/NeurIPS: the
  architecture is a recombination of well-known components (GraphSAGE,
  R-GCN-style relation summing, a gating decoder, replay-based continual
  learning). The paper's own ablations show the "novel" role-aware decoder
  is a second-order effect (±0.02–0.03 AUPR vs. plain concat fusion). This
  reads as a solid empirical/applied bioinformatics contribution, not a
  methods paper (Section 6).

**Recommendation**: target a bioinformatics/computational-biology venue
(Bioinformatics, Genome Biology, PLOS Computational Biology, or a strong
journal track) rather than ICLR/ICML/NeurIPS, and do not submit anywhere
until Sections 2, 3, 4, and 8 are addressed.

---

## 2. TOP FINDING — The zero-shot guarantee is coincidental, not structural, for the currently-recommended headline configuration

### 2.1 What the paper claims

`paper/main.tex` (Algorithm 1, Proposition 1, `sec:split`) describes a
"global, hierarchy-consistent" edge-splitting protocol: every unique
(TF, gene) pair, across the union of all evidence tiers, is assigned to
exactly one split (train/val/test) once, and each tier's own split is
derived from that single assignment. This is implemented correctly as
`create_global_edge_splits` in `src/data/graph_builder.py:335-405`, is run
once in `scripts/02_preprocess.py:123-142`, saved to `data/splits/`, and
audited by `verify_no_cross_tier_leakage` (`graph_builder.py:408-424`),
which prints `LEAK`/`OK` per tier. This machinery is real and does what it
says — I re-ran it against `data/processed/K562` and confirmed
`create_global_edge_splits` reproduces the on-disk splits exactly.

### 2.2 What actually trains the model under `all_at_once`

`plan.md` Section 3b / `results.md` Section 3 and 7 (dated 2026-09-18) now
recommend `all_at_once + FM + h384/l2` as **"the recommended configuration
for the paper's headline numbers"** — the numbers that beat scMultiomeGRN on
all three tiers (localization 0.967, perturbation 0.690, dual-evidence
0.950 AUPR, `results.md` Section 7).

But `all_at_once` does **not** use the global, leak-checked splits at all
for its training positives. `src/training/curriculum.py:31-57`
(`build_all_at_once_stage`) does this instead:

```python
main_tiers = set(cfg["curriculum"].get("main_curriculum_tiers", splits_per_tier))
merge_tiers = [t for t in splits_per_tier if t in main_tiers]
pos = torch.cat([data.evidence[t] for t in merge_tiers], dim=1)   # FULL evidence, not train-split-only
pos = torch.unique(pos, dim=1)
...
merged = create_edge_splits(pos, data.negative_pool, ...)          # a FRESH, independent random split
```

`data.evidence[t]` is the tier's entire positive set (train+val+test
pooled), not `splits_per_tier[t]["train"]["pos"]`. `create_edge_splits`
(`graph_builder.py:293-332`) is a **different, older, per-tier-only**
splitting function whose own docstring says:

> "For anything that spans multiple evidence tiers, use
> `create_global_edge_splits` instead — splitting each tier independently
> lets a held-out edge in one tier leak into another tier's training set
> whenever tiers are nested... Measured on K562: 85.8% of dual_evidence's
> 'test' edges were already present in localization-train or
> perturbation-train [under this function]."

`build_all_at_once_stage` uses exactly this leak-prone function, on exactly
the nested multi-tier data it warns against, to produce the training
positives for the paper's proposed headline protocol. Since dual_evidence
is (per the project's own numbers) 100% nested in perturbation, and
perturbation is one of the two merged tiers, most of dual_evidence's
positive edges are re-shuffled into a training pool with no relationship to
the global split whatsoever. The same applies to `scripts/06_ablation.py`'s
`--ablation all_at_once` path (`scripts/06_ablation.py:172-176`), which
calls the identical function — so every `all_at_once` number reported
anywhere in `results.md` (Section 3's 0.937/0.351/0.791 and Section 7's
0.967/0.690/0.950) rests on this code path.

### 2.3 What I actually measured

I reproduced `build_all_at_once_stage`'s split logic directly against
`data/processed/K562/evidence_*.pt` and `data/splits/K562_*_splits.pt`:

```
merged union pos (localization ∪ perturbation)      : 1,377,658
all_at_once merged train_pos                          : 964,360
dual-evidence held-out (val+test) eval set            : 7,943
  -> overlap with all_at_once train_pos               : 0 / 7,943   (0.0%)
localization TEST split                               : 175,237
  -> overlap with all_at_once train_pos               : 0 / 175,237 (0.0%)
perturbation TEST split                                : 34,457
  -> overlap with all_at_once train_pos               : 0 / 34,457  (0.0%)
```

I then compared `merged["train"]["pos"]` directly against the union of the
official, on-disk, leak-checked train splits for localization+perturbation:
**the two sets are byte-for-byte identical** (0 symmetric difference).

**So: for K562 as currently configured, there is no actual leakage in the
`all_at_once` numbers.** But this is not because the code guarantees it —
it is because three independent conditions happen to coincide:

1. Both `create_edge_splits` (called from `build_all_at_once_stage`) and
   `create_global_edge_splits` (called from `scripts/02_preprocess.py`) are
   seeded with the exact same value (`cfg["data"]["seed"] = 42` in both
   places), using `torch.Generator().manual_seed(seed)` — same PRNG, same
   seed, same `n` → same permutation.
2. Because dual_evidence is (currently, for K562) **100%** nested inside
   perturbation, `localization ∪ perturbation ∪ dual_evidence` (the
   population `create_global_edge_splits` splits over) is set-equal to
   `localization ∪ perturbation` (the population `build_all_at_once_stage`
   splits over) — so the two functions are, by luck, splitting the *same*
   underlying set.
3. `torch.unique(pos, dim=1)`'s column sort order happens to match Python's
   `sorted()` over the same tuples, so the two functions enumerate that
   identical set in the identical order before permuting it.

None of this is asserted, tested, or even mentioned anywhere in the code or
docs. It is an accidental invariant, not a designed one, and it is fragile
in ways that matter directly to this project's own stated future work:

- **It requires the exact same seed at both call sites.** Nothing in the
  code enforces this; `cfg["data"]["seed"]` is read independently by
  `scripts/02_preprocess.py` (at data-prep time) and by
  `build_all_at_once_stage` (at train time, from whatever config is passed
  to `scripts/03_train.py`/`06_ablation.py`). A config refactor, a
  `--seed` CLI override, or simply preprocessing once and training many
  times with different configs would silently break the coincidence.
- **It requires ~100% nesting.** The project's own numbers show this does
  *not* hold for the loc↔pert pair itself (only 9.0% overlap) — nesting
  fractions are data-dependent and cell-type-dependent. `plan.md`/`results.md`
  Section 6 report from-scratch and joint-training runs on Macrophage and
  MCF7; if either cell type's dual/localization/perturbation nesting is not
  ~100%, or if `dual_evidence_mode: reconstruct` is used instead of `as_is`
  for a future cell type, this invariant breaks silently and
  `all_at_once` training would directly leak held-out test/val positives
  as label=1 training examples — with no error, warning, or leak-check to
  catch it, because `verify_no_cross_tier_leakage` is **never called** on
  the output of `build_all_at_once_stage`.
- **It requires two independently-implemented sort routines
  (`torch.unique(dim=1)` and Python `sorted()`) to keep agreeing** across
  PyTorch versions — an implementation detail neither routine documents as
  a stable contract.

This is precisely the failure mode the project spent real effort finding
and fixing once already (`results.md` Section 0's hard-negative leakage
bug — also a cross-tier train/eval contamination issue, also caught only by
manual investigation after suspicious numbers, not by an automated check).
The same class of bug is sitting, currently dormant, in the exact protocol
now slated to produce the paper's headline table.

### 2.4 Fix (concrete, small)

`build_all_at_once_stage` should not re-split anything. It should build its
training pool directly from the **already-correct, already-leak-checked**
per-tier train splits:

```python
train_pos = torch.cat([splits_per_tier[t]["train"]["pos"] for t in merge_tiers], dim=1)
train_pos = torch.unique(train_pos, dim=1)
```

This makes the "no edge trained under one tier can be held out under
another" guarantee hold **by construction**, matching what Proposition 1 in
the paper actually claims, instead of holding by a seed/nesting/sort-order
coincidence that nobody checked. It should be followed by actually running
`verify_no_cross_tier_leakage`-style assertions against the merged training
pool before trusting any `all_at_once` number again, and that check should
be re-run for Macrophage/MCF7 before their `all_at_once` (or joint-training)
numbers are trusted for a paper table. This is the single highest-priority
fix in this review — it directly undermines the paper's stated main
methodological contribution (Section 5.2's "we formalize an evidence-tier
curriculum with a provable guarantee against a model being evaluated on an
edge it was trained on under a different tier").

*(Separately, `scripts/12_joint_train.py`, the multi-cell-type joint
training script, is fine — it trains each cell type only on its own
on-disk `localization` split loaded via `load_splits`, and evaluates the
held-out cell type on that cell type's full evidence, which is a distinct
gene/edge universe from the training cell types. That script does not
share `build_all_at_once_stage`'s problem.)*

---

## 3. The "zero-shot" tier has been used as a de facto validation/selection metric all project long

Separately from Section 2's literal leakage risk, there is a methodological
problem the project should confront honestly: **dual-evidence AUPR has been
the deciding number for nearly every architecture and protocol choice made
in this project**, even though the entire point of holding it out is that
it is supposed to be a clean, unbiased generalization estimate for the
*final* model.

Evidence from `plan.md`/`results.md`, in order:

- Section 1.2: the sequential-vs-all_at_once training-schedule decision was
  made "because [all_at_once] does *worse* on... the score that matters most
  for the paper" (dual-evidence AUPR).
- Section 2: the ATAC pipeline rewrite is validated by "the actual 'does
  ATAC help' test," measured on dual-evidence AUPR (0.856 vs 0.881).
- Section 3 (model-size sweep): five configs are ranked and the final size
  chosen by "dual AUPR" (0.9515 → 0.9638).
- Section 3b/7 (2026-09-18): the curriculum protocol recommendation is
  reversed a *second* time — "the tradeoff does not hold at this size" —
  again decided by comparing dual-evidence AUPR across configurations.
- Section 5: Geneformer is adopted specifically because it improves
  dual-evidence AUPR (0.881 → 0.951).
- Section 6: the joint-training recipe, and even *which* direction the FM
  embedding helps or hurts, is adjudicated by dual/held-out-cell-type AUPR.

Every one of these is a legitimate thing to *measure*, but collectively
they mean dual-evidence AUPR has functioned as the project's validation
metric for at least six major decisions (protocol, ATAC featurization,
model size, use of FM, curriculum choice a second time, joint-training
recipe). A benchmark that is used to pick between six-plus candidate
configurations and is *also* reported as "the zero-shot generalization
number" for the winning configuration is no longer an unbiased estimate of
generalization to that tier — it is closer to a validation set that the
final reported number is implicitly optimistic with respect to (classic
"researcher degrees of freedom" / multiple-comparisons inflation, distinct
from literal gradient leakage).

**Recommended fix**: designate a true, untouched final-check split — e.g.,
carve out a small fixed slice of dual-evidence (or a 4th tier, or a
different cell type entirely) that is *never* looked at until the paper's
very last table, and use localization/perturbation's own validation splits
(which are legitimately trained-tier validation data) for all of the
architecture/protocol decisions instead. Failing that, at minimum the paper
should state plainly that dual-evidence AUPR was used iteratively for model
selection across the project and that the reported number should be read
with that caveat — reviewers at a top ML venue will ask this question
directly.

---

## 4. Statistical Rigor

- **Every number in `results.md` and the paper is a single run, single
  seed** (`seed: 42` in `configs/default.yaml`, never varied). There is no
  code path anywhere in `scripts/` that repeats a run across seeds, and no
  variance/CI is reported for any table.
- Decisions are routinely made on differences that are plausibly within
  seed noise: the model-size sweep concludes "diminishing returns" from
  h384l2→h384l3 being "+36% params for +0.001 dual AUPR" (`plan.md` Section
  3) — a difference indistinguishable from run-to-run noise with n=1 per
  config. The paper's own head-to-head table calls a 0.011 AUROC gap on
  perturbation "very close" without ever establishing what the noise floor
  actually is.
- `results.md` Section 7 reports MEvD-GRN beating scMultiomeGRN by margins
  as small as 0.011 AUROC and calls it a win outright. Without repeated
  seeds (even 3–5 would help enormously here, since training is fast —
  16.6–30 minutes per the paper/`results.md`), none of these claims can be
  distinguished from noise, and a competent reviewer will say so.
- `compute_early_precision`/`compute_epr` (`src/evaluation/metrics.py:22-41`)
  are deterministic given scores/labels, but the negative *sampling* that
  feeds into the eval sets (`create_negative_pool`, `create_edge_splits`)
  is itself a random draw — the EP/EPR numbers therefore also carry
  unquantified sampling variance from the negative draw that is never
  reported.

**Fix**: rerun the headline configuration (and the two or three baselines
it's compared against) across 3–5 seeds and report mean ± std (or at least
range) for AUPR/AUROC on all three tiers before any submission. This is
cheap given the reported per-run training times.

---

## 5. Baseline Fairness

- **scMultiomeGRN's training budget is inconsistent between documents.**
  The paper (`paper/main.tex`, Section "Baselines" and Limitation (4)) says
  it is trained for "100 epochs per tier (an exploratory budget; its
  published default is 2,000 epochs)" and explicitly lists this as an open
  limitation not yet resolved. `results.md` Section 7, dated later
  (2026-09-18), says a full 2000-epoch × 2-tier run was since completed via
  4-GPU DDP, and that comparison is what backs the "MEvD-GRN beats
  scMultiomeGRN on all three tiers" headline claim. So the paper draft and
  the actual current comparison are not just stale in wording — they
  describe two *different experiments* against two *different trained
  versions* of the same baseline, and the paper's caveat ("not yet at the
  baseline's full budget... we intend to re-run before treating [it] as
  final") has been silently resolved in `results.md` but not reflected
  anywhere in the manuscript.
- **The scMultiomeGRN baseline itself flags, in its own results.md
  writeup, that it "hadn't converged"**: "both were still slowly improving
  at epoch 2000, e.g. perturbation val AUPR kept climbing to 0.654." A
  baseline that the authors' own log says is still improving at the budget
  cutoff is not a fully-tuned baseline; declaring victory over it on
  perturbation/localization AUPR margins of a few points is not solid
  ground. `learning_rate=1e-5` (`src/baselines/scmultiomegrn_wrapper.py:379`)
  is also unusually low and undocumented as tuned — there is no evidence of
  a learning-rate sweep for this baseline anywhere, in contrast to the
  multi-config sweep MEvD-GRN itself received.
- **GRNBoost2, RegDiffusion, and GMF-GAE are unsupervised methods being
  compared against supervised link-prediction AUPR.** None of the three
  ever sees the localization/perturbation/dual-evidence positive/negative
  labels during fitting — GRNBoost2 regresses expression on TF expression,
  RegDiffusion learns a diffusion prior over expression alone, and GMF-GAE
  is an unsupervised graph autoencoder trained with random negative
  sampling unrelated to the evidence tiers (`src/baselines/gmfgrn_wrapper.py:94-103`
  samples negatives via `torch.randint`, not from `data.negative_pool` or
  any tier-aware split). Their near-random AUPR (0.16–0.22) on
  perturbation/dual is close to a foregone conclusion given they never see
  supervision at all — comparing them is standard practice as a *sanity
  floor*, but the paper's framing ("The RNA-only classical and diffusion
  baselines... trail on every tier, consistent with none of them exploiting
  chromatin accessibility or the evidence-tier structure") undersells that
  this is an apples-to-oranges comparison (supervised vs. unsupervised),
  not purely a multi-omic-vs-RNA-only comparison. scMultiomeGRN is the only
  fair, like-for-like (supervised, multi-omic) comparison point, and it is
  the one where the margins are smallest and the baseline is
  under-converged.
- **The "official" GMFGRN baseline was never actually run.**
  `run_gmfgrn()` (`src/baselines/gmfgrn_wrapper.py:31-46`) raises
  `NotImplementedError` unless the official repo is manually cloned and
  wired up — this was evidently never done, and `run_gmf_gae()` (a
  self-authored graph-autoencoder stand-in) is used instead. This is
  disclosed in the code's own docstring, which is good practice, but the
  paper (`paper/main.tex`, Section "Baselines") describes "GMF-GAE... as a
  graph-autoencoder baseline... standing in for the GMFGRN family," which
  is honest, but a reviewer at a bioinformatics venue may reasonably ask
  why the actual published method (with a public repo) was not used, given
  it is cited as a direct comparison point.
- **No hyperparameter search effort is reported for any baseline**,
  whereas MEvD-GRN itself received a 5-point model-size sweep and multiple
  curriculum-protocol comparisons. This asymmetry (extensive tuning effort
  on the proposed method, essentially none on the baselines beyond what
  ships in each wrapper's defaults) is a classic source of favorable-by-
  construction comparisons and should be flagged explicitly in the paper
  if not corrected.

---

## 6. Novelty & Venue Fit

Read charitably, MEvD-GRN's contributions are:
1. An evidence-tier curriculum over ChIP-seq/knockout/dual-evidence tiers
   (a domain-specific training-order heuristic, not a new algorithm — it's
   an application of curriculum learning (Bengio et al. 2009) and replay
   (Rolnick et al. 2019), both cited correctly as prior work, not claimed
   as new).
2. A "role-aware" decoder where ATAC gates RNA multiplicatively rather than
   being concatenated. This is a reasonable, biologically-motivated design
   choice, but the project's *own* ablation table
   (`results.md` Section 2, and `paper/main.tex` Table `tab:ablation`)
   shows it is a second-order effect: role-aware / gated-fusion /
   concat-fusion sit within ~0.01–0.03 AUPR of each other. The paper says
   this outright: "the role-aware decoder's advantage is currently more
   interpretive/biologically-motivated... than a large raw-metric win." A
   design whose authors themselves report it doesn't move the needle much
   is not going to read as a methodological advance to an ICLR/ICML
   audience.
3. A global, hierarchy-consistent edge-splitting protocol
   (Proposition 1). This is good engineering hygiene and a genuinely useful
   idea for anyone else working with nested-evidence GRN benchmarks, but it
   is a dataset-handling protocol, not a new learning algorithm or
   theoretical result — Proposition 1 itself is a one-line consequence of
   assigning each item to one partition (the "proof," such as it is, is
   definitional). This belongs in a Methods section of a bioinformatics
   paper, not as a proposition in an ICML submission. It is also, per
   Section 2 above, currently *not actually enforced* for the recommended
   headline training protocol — undercutting exactly the claim it's meant
   to support.
4. A learned gated combiner over 2 relations summed via a per-layer GNN
   combiner, an optional gated softmax over relations
   (`combine_mode: gated`), and a motif-based third relation that exists in
   code (`src/data/motif_scan.py`) but has never actually been run/ablated
   (`use_motif: false` by default everywhere; no
   `results/ablations/*motif*` file exists) — so the mechanism specifically
   designed to make the "novel" gated-relation idea interesting
   ("[it] doesn't move the numbers yet... will matter once the motif graph
   exists," `results.md` Section 2) is still unrealized at review time.

None of this is a knock against the work's value — it is a well-executed,
honestly-reported empirical study on a real, underused resource
(SC-MO-GRN-DB's tier structure, and the multi-cell-type joint-training
result in particular is a genuinely interesting empirical finding). But
**the ML novelty bar for ICLR/ICML/NeurIPS is not met**: there is no new
architecture class, no theoretical contribution beyond a bookkeeping
proposition, and the ablations the authors ran themselves show the
architectural choices are not where the performance comes from (the GNN
backbone and the FM embedding are — both are off-the-shelf components,
GraphSAGE and Geneformer respectively, used as intended). This paper's
comparative advantage is the domain framing, the dataset use, and the
honest empirical rigor about a real biological problem — exactly the
profile of a strong bioinformatics/computational-biology paper (Genome
Biology, Bioinformatics, PLOS Computational Biology), not a top-tier ML
venue.

---

## 7. Ablation and Sweep Rigor

- The model-size sweep (`plan.md` Section 3, 5 configs) is a genuinely
  good practice — "picked the smallest model within noise of the best"
  is exactly the right instinct — but with n=1 run per config (Section 4
  above), "within noise" is asserted, not measured; there is no repeated-
  seed data to know what the noise band actually is.
- The curriculum-protocol conclusion **flips twice** over the project's
  timeline (sequential wins → all_at_once wins, `plan.md` Section 3b/7),
  each time attributed to a plausible mechanistic story (perturbation
  dilution, then "enough capacity to avoid the dilution"). Both stories are
  reasonable post-hoc explanations, but neither is tested directly (e.g. no
  controlled experiment varying capacity alone at fixed FM/protocol to
  isolate whether capacity, FM, or both explain the flip) — this reads as
  narrative-fitting after the fact rather than a pre-registered ablation
  design. Combined with Section 3's point about dual-evidence AUPR being
  the deciding metric for this exact flip, this is one of the weaker parts
  of the empirical story as currently told.
- The FM-transfer "sign flip by holdout choice" (`results.md` Section 6,
  "holding out K562 instead... FM actually helps") is reported honestly as
  a "small, second-order, holdout-dependent wobble," which is good
  practice, but with only 3 held-out-cell-type conditions × 1 seed each,
  there is no way to know if this is a real interaction effect or noise.
- The motif-graph relation (Section 6 above) and the gated relation
  combiner it's meant to make useful are both implemented but essentially
  untested against real data — the "gated_relations" ablation
  (`results/ablations/gated_relations_K562.json`, dual AUPR 0.876 vs 0.881)
  only has two, already-similar-strength relations to weigh, so it cannot
  currently demonstrate anything the mechanism is actually for. Presenting
  this as a finished, "done" line item (`plan.md` Section 4, "Done") is a
  bit generous — it's an infrastructure piece that has not yet been
  exercised under conditions where it could show a real effect.

---

## 8. Reproducibility

- `requirements.txt` pins loose version ranges (`torch>=2.1`,
  `torch-geometric>=2.4`, etc.) with no lockfile, no `environment.yml`,
  and no Docker/conda-lock artifact — reasonable for a research repo but
  insufficient for exact reproducibility, especially since Section 2's
  entire "coincidental non-leak" finding depends on the *exact* sort
  behavior of `torch.unique(dim=1)`, which is not a documented, versioned
  API contract.
- Data dependencies (SC-MO-GRN-DB raw networks, single-cell matrices,
  hg38 genome FASTA, JASPAR PFM file, Geneformer's pretrained embedding
  table) are fetched via cell-type-specific globs and external downloads
  (`scripts/01_download_data.sh`, `scripts/09_download_transfer_data.sh`,
  `scripts/11_extract_fm_embeddings.py`) rather than bundled or hashed —
  an outside researcher would need to independently source and verify each
  of these matches what was used here; there's no checksum/version pinning
  visible for the JASPAR release or the exact Geneformer checkpoint beyond
  a docstring mention (`ctheodoris/Geneformer`, `Geneformer-V2-104M`).
- `plan.md` explicitly documents cluster paths specific to one author's
  account (`/share1/vishakkashyap.k/mevd_grn`) and a specific compute
  environment (Ada, IIIT-H, 4×2080Ti). This is fine as project notes but
  should not leak into anything claimed as a reproducibility artifact.
- Seed is fixed (`seed: 42`) and used consistently for the parts of the
  pipeline that matter for the split-safety argument (Section 2), which is
  good — but as noted, nothing enforces that preprocessing and training
  configs actually share this value if someone runs them independently in
  the future.
- The official GMFGRN and full RegDiffusion/GRNBoost2 dependencies are
  commented out as optional in `requirements.txt` ("install on demand"),
  and GMFGRN's official path is an unfinished stub
  (Section 5) — a researcher trying to reproduce the baseline table exactly
  as described in `docs/citations.md` Section 14 cannot, without first
  writing the missing integration code themselves.

**Net assessment**: a motivated outside researcher could reproduce the
MEvD-GRN main-model numbers with moderate effort (the training scripts
are clear and config-driven), but could not reproduce the full baseline
comparison table without additional implementation work, and could not
verify the split-safety guarantee's robustness without independently
discovering the coincidence documented in Section 2.

---

## 9. Writing / Paper-Draft Issues

This is the most immediately actionable problem, and should block
submission on its own until fixed:

- **`paper/main.tex` reports pre-bugfix numbers throughout its Results,
  Ablation, and Limitations sections**, and says so explicitly in its own
  "EDITORIAL STATUS NOTE" (lines 87-100): "The body sections below this
  note... still contain PRE-BUGFIX numbers and have not yet had a full
  rewrite pass — treat results.md as the source of truth until that pass
  is done." The Abstract (lines 53-84) already reflects post-bugfix
  numbers (0.881/0.971, 0.951/0.989), but the body directly contradicts it
  a few pages later: Table `tab:main` reports dual-evidence AUPR
  **0.558**/AUROC **0.838** (the pre-bugfix, since-retracted numbers from
  `results.md`'s own Section 0 "two bugs fixed" writeup), and Table
  `tab:forgetting` reports localization AUROC **0.326** post-Stage-2 — the
  exact below-random number `results.md` Section 0 identifies as the
  *symptom of the leakage bug*, not the corrected result. A reviewer
  reading this paper today would be reading numbers the authors' own repo
  says are wrong.
- **The EPR definition described in the paper is the bugged one.**
  `paper/main.tex` Section "Evaluation Protocol" defines EPR as "EP
  normalized by the random-guessing base rate over the full bipartite
  candidate space ($225\times22{,}942=5{,}161{,}950$ pairs)" and reports
  EPR values of 12.5, 80.2, 272.5 in Table `tab:main` — this is exactly the
  metric bug `results.md` Section 0 says was fixed ("early precision was...
  normalized against a genome-wide TF×gene candidate count, inflating EPR
  into the hundreds... Fixed to normalize against the actual eval set").
  The actual code (`src/evaluation/metrics.py:36-41`,
  `compute_epr(early_precision, n_positives, n_total)` called with
  `n_total_candidates=None` i.e. the eval split's own size,
  `src/training/trainer.py:205-209`) does **not** implement what the paper
  describes. This is a direct, checkable contradiction between the paper's
  stated methodology and the shipped code.
- **The headline number has changed at least four times** across the
  project's lifetime with no single passage anywhere reconciling them:
  paper abstract says 0.881 (then "further additions push this... to
  0.951"); `results.md` Section 7's current recommended config reports
  0.950 (all_at_once+FM+h384/l2) alongside 0.963 (sequential+FM+h384/l2)
  as apparently-competing "headline" rows in the same table. Before
  submission, the paper needs exactly one unambiguous headline
  configuration and number, clearly reconciled with `results.md`.
- The paper's own Limitations Section (4) says scMultiomeGRN was trained at
  a reduced 100-epoch budget "and we intend to re-run it at full budget
  before treating Table `tab:baselines` as final" — this re-run has since
  happened per `results.md`, but the paper was never updated to reflect it
  (Section 5 above).
- Author list is a placeholder ("MEvD-GRN Project," `mevd-grn@example.org`)
  — fine for a working draft, not for a submission checklist item to
  forget.

---

## 10. Architecture / Domain Critique (secondary points)

- **RNA/ATAC node features are extremely low-dimensional** (3-dim RNA
  mean/var/detection + a 50-dim SVD co-expression signature; ATAC now 3+4
  dims post-fix). `plan.md` Section 3 itself argues this is *why* the
  model was kept small (306K params) — a reasonable diagnosis, but it also
  means the pre-FM "hand-built feature" story that much of Sections 1-4 of
  `results.md` is built on was, by the team's own later admission (Section
  5), the weak link the whole time; a large fraction of the ablation table
  in the current paper draft is characterizing a featurization the team
  has since mostly abandoned in favor of Geneformer. The paper needs to be
  honest that the "RNA-only ablation adds nothing" (`rna_only` vs. full
  model) finding is now superseded by the FM finding that *hand-built*
  features add nothing once Geneformer is present (`results.md` Section
  5) — these are two different, easily-conflated claims and the current
  draft only reflects the first, now-stale one.
- **K562 is a single, immortalized cancer cell line with only 953 RNA /
  1,126 ATAC cells.** The truncated-SVD co-expression signature and the
  gene-activity RP scores are both computed from these small, noisy
  single-cell matrices; no cross-validation or bootstrap of the
  feature-computation step itself is reported, so it's unclear how much of
  the reported signal is stable vs. dataset-specific noise from ~1,000
  cells. The multi-cell-type results (Macrophage, MCF7) partially address
  the single-cell-type generalizability concern, which is good, but all
  three cell types are still evaluated with n=1 run each (Section 4).
- **Localization ChIP-seq edges (1.4M raw, 1.17M aligned) are being
  treated as ground truth "positives"** at face value; ChIP-seq binding
  calls are well known in the field to have substantial false-positive
  rates and cell-line/antibody/batch dependence. The paper's own framing
  ("high recall... but low precision") acknowledges this qualitatively but
  doesn't quantify or discuss how much label noise in the localization
  tier could be inflating or deflating any of the reported metrics —
  reasonable for a first paper on this dataset, but worth a sentence of
  acknowledgment for a bioinformatics-venue reviewer.
- **`negative_pool.pt` and `evidence_*.pt` are cell-type/gene-universe
  specific**, so cross-cell-type "zero-shot" comparisons (Section 6 of
  `results.md`) are not literally comparing the same candidate space across
  cell types (different gene universes, different TF counts: 22 TFs for
  Macrophage vs 250 for MCF7 vs 225 for K562) — worth a caveat that
  "zero-shot to Macrophage" and "zero-shot to MCF7" are not measuring
  exactly the same kind of task difficulty.

---

## 11. Recommended Fixes, Ranked by Severity

1. **(Blocking, correctness)** Fix `build_all_at_once_stage`
   (`src/training/curriculum.py:31-57`) to build its training pool from
   `splits_per_tier[t]["train"]["pos"]` instead of re-splitting
   `data.evidence[t]` from scratch, and add an explicit
   `verify_no_cross_tier_leakage`-style assertion on the result. Re-run
   this check for Macrophage and MCF7 before trusting any `all_at_once` or
   joint-training number for those cell types. (Section 2)
2. **(Blocking, integrity)** Rewrite `paper/main.tex` body sections
   (Results, Ablation Study, Limitations, Evaluation Protocol's EPR
   definition) to match `results.md`'s current numbers and the actual EPR
   implementation, and settle on exactly one headline configuration/number
   before any submission. Remove the "remove before real submission"
   editorial note only once this is actually done. (Sections 8, 9)
3. **(High, rigor)** Re-run the headline configuration and its strongest
   baseline (scMultiomeGRN) across ≥3-5 seeds; report mean±std/range for
   AUPR/AUROC on all three tiers. Do not report point-estimate margins
   below ~2× the measured seed-to-seed spread as a "win." (Section 4)
4. **(High, methodology)** Explicitly address the "dual-evidence AUPR used
   as de facto validation metric across ~6 major decisions" concern
   (Section 3) in the paper — either restructure model/protocol selection
   to use a genuinely untouched final check, or add an explicit limitations
   paragraph acknowledging the iterative-selection process and its effect
   on how the final zero-shot number should be interpreted.
5. **(Medium, fairness)** Give scMultiomeGRN (and ideally GRNBoost2/
   RegDiffusion/GMF-GAE) a comparable tuning effort to what MEvD-GRN
   received (a small LR/epoch-budget sweep), and report scMultiomeGRN as
   converged (or explicitly flag it as not) before using it for a
   head-to-head "who wins" table. Clarify in the paper that
   GRNBoost2/RegDiffusion/GMF-GAE are unsupervised sanity floors, not
   like-for-like competitors. (Section 5)
6. **(Medium, venue fit)** Reframe the contribution honestly for the
   target venue: if pursuing ICLR/ICML/NeurIPS, either find and foreground
   a genuine methodological/theoretical contribution beyond the split
   protocol and role-aware gate (both already shown to be second-order
   empirically), or retarget to a bioinformatics venue where the dataset
   framing, curriculum idea, and multi-cell-type transfer result are a
   strong, novel-enough empirical contribution on their own. (Section 6)
7. **(Low-medium, completeness)** Either finish and report the motif-graph
   ablation (the mechanism the gated-relation combiner needs to be
   interesting) or remove/soften claims that this workstream is "done."
   (Section 7)
8. **(Low, reproducibility)** Pin exact dependency versions (or provide a
   lockfile/Docker image), and document/checksum the exact external
   resources (JASPAR release, Geneformer checkpoint, genome FASTA build)
   used, separate from the personal-cluster-path notes in `plan.md`.
   (Section 8)
