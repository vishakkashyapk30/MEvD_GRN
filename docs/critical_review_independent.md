# Independent Critical Review — MEvD-GRN

*Written 2026-09-18 and moved from the repo root to `docs/` on 2026-09-29.
File and line references below are as of 2026-09-18. The items still open
are tracked in `plan.md` Section 10.*

Reviewer stance: skeptical, code-first. Every claim below was checked against
the actual source in `src/`, the committed configs/results, and — where
possible — reproduced or falsified with a throwaway script against the real
data on disk. Line/file references are exact as of the current working tree
(`main`, with `plan.md`/`results.md` locally modified but uncommitted).

## Bottom line

This is a genuine, hard-working empirical project with two real, well-diagnosed
bugs already fixed and an honest habit of reporting negative results. But the
newest and most important claims in the project — the ones that would anchor a
paper's abstract — are **not backed by any artifact in the repository**, and
the manuscript draft still reports numbers the authors themselves know are
wrong. As it stands, this is not close to ICLR/ICML-main-track, and it is not
yet submission-ready anywhere. It could become a solid applied
bioinformatics-venue paper (PLOS One / Bioinformatics-tier) once the items in
Sections 1–3 are resolved and the manuscript is actually re-synced to
`results.md`.

---

## 1. [CRITICAL] The paper's newest headline result has no supporting artifact anywhere in the repo

`plan.md` §3b ("Resolved 2026-09-18") and `results.md` §7 both state, as the
project's central new finding, that **`all_at_once` + FM + h384/l2 (3.26M
params) beats scMultiomeGRN on AUPR across all three tiers** — localization
0.967 vs 0.906, perturbation 0.690 vs 0.655, dual-evidence 0.950 vs 0.847 — and
this configuration is explicitly recommended as *"the recommended
configuration for the paper's headline numbers."*

I could not find a single piece of evidence for this run anywhere in the repo:

- **No config**: no committed YAML combines `protocol: all_at_once` with
  `hidden_dim: 384` / `n_gnn_layers: 2` / `use_fm: true`. The five
  `configs/sweep/*.yaml` files (which do set hidden_dim/n_gnn_layers/use_fm)
  all inherit `protocol: sequential` from `configs/default.yaml:102` and never
  override it. Producing the claimed run would require
  `scripts/06_ablation.py --ablation all_at_once --config configs/sweep/k562_fm_h384_l2.yaml`,
  which is nowhere documented as having been run.
- **No results JSON**: `results/ablations/all_at_once_K562.json` (the only
  `all_at_once` artifact that exists) reports localization AUPR **0.9367**,
  perturbation **0.3477**, dual **0.7876** — these are the *old*, small
  (307K-param, no-FM) numbers from `results.md` §3's original table, not the
  0.967/0.690/0.950 numbers claimed for the new configuration. I grepped every
  JSON under `results/` for `0.967`/`0.966`/`0.950`; none contains them.
- **No checkpoint**: every file under `results/checkpoints/K562/*.pt` is
  1.2–1.3 MB, consistent with the ~307K/315K-param model (fp32, ~4
  bytes/param), and all are dated Sept 11. Nothing dated Sept 18, nothing
  sized for a 3.26M-param model, exists anywhere on disk (checked with `find
  ... -newer`, `ls -la`, and `git status --porcelain -uall`).
- `git status` shows only `plan.md` and `results.md` as modified — no new
  config, log, checkpoint, or results file was added alongside these claims.

The same is true, even more starkly, of the "full 3-tier scMultiomeGRN" run
that `results.md` §7 also treats as settled ("Total wall time: 5h 58m...
source: `results/logs/smg_ddp_full.log`, `results/baselines/scmultiomegrn_K562.json`"):

- `results/logs/smg_ddp_full.log` **does not exist**.
- `results/baselines/scmultiomegrn_K562.json` contains only the *old*,
  explicitly-labeled-incomplete run: `"note": "Manually stopped after epoch
  1145 ... Perturbation tier was never reached"`, with only a `localization`
  entry (AUPR 0.9009). There is no `perturbation` or `dual_evidence` entry in
  this file, despite `results.md` quoting perturbation AUPR 0.655 and dual
  AUPR 0.847 from it.
- `results/checkpoints/scmultiomegrn_K562/` contains only
  `localization_best.pt`; there is no `perturbation_best.pt`, even though a
  "full curriculum" run would need to save one (the wrapper's own
  `_fit_tier` checkpoints every tier it fits, per
  `src/baselines/scmultiomegrn_wrapper.py:413-421`).

**Why this matters**: these two results are not incidental — they are the
justification for switching the paper's recommended default configuration
(§3b explicitly supersedes the earlier, well-supported "sequential is the
default" conclusion) and for the claim that MEvD-GRN now beats the strongest
baseline outright. A reviewer who tried to reproduce either number from what's
checked into the repo today would find nothing to run. Until a config, log,
checkpoint, and results JSON exist for both runs, these numbers should be
treated as unverified, not as "resolved."

> **UPDATE (2026-09-18, after this review was written): half-fixed, and the
> other half is worse than this section originally described.** The runs
> themselves were real — done on a remote cluster, fetched and reported
> live via SSH — and the missing-artifact gap this section describes is
> real too, just not for the reason it assumed. `results/ablations/all_at_once_fm_h384l2_K562.json`
> (0.967/0.690/0.950), `with_replay_fm_h384l2_K562.json`,
> `fm_only_h384l2_K562.json`, and the full 3-tier `results/baselines/scmultiomegrn_K562.json`
> now exist **on this machine's disk**, synced from the remote cluster. But
> `.gitignore` has a blanket `results/` entry — `git ls-files results/`
> returns zero files, and always has. None of this was ever going to be
> visible to `git clone`, with or without syncing; syncing only fixed the
> "does this file exist anywhere" half of the problem, not the "is it part
> of the repository" half. The right fix is a `.gitignore` carve-out for
> the small results/ablations/*.json and results/baselines/*.json text
> files specifically (checkpoints and logs can stay ignored — large and
> regenerable), followed by actually committing them. Not done as of this
> update; flagging rather than doing it unilaterally since it's a
> repo-convention change. Missing config note above no longer applies —
> `configs/sweep/k562_fm_h384_l2.yaml` was run with `--ablation all_at_once`
> directly (the ablation script overrides the protocol via its `--ablation`
> flag rather than needing a config-level override). This section's
> diagnosis was correct at the time it
> was written; leaving it unedited above as the record of what was found.

## 2. [CRITICAL] The manuscript itself still reports numbers the project has retracted

`paper/main.tex`'s abstract (lines 68-75) correctly reports the post-bugfix
headline (AUPR 0.881/AUROC 0.971 dual-evidence), but the body tables
disagree with it, in the same document:

- Table at line 722-726: dual-evidence **AUPR 0.558 / AUROC 0.838**, EPR
  272.5, localization AUROC 0.326.
- These are exactly the pre-bugfix, explicitly-flagged-as-wrong numbers.
  `results.md` §0 states in so many words that the EPR metric was inflated
  "up to ~450" by a normalization bug, and that pre-fix dual AUPR was 0.558
  (post-fix: 0.870/0.881). Localization AUROC 0.326 is the *bug's* number
  (post-fix: ~0.51-0.52).
- The paper has a self-aware "EDITORIAL STATUS NOTE" (lines 85-99) admitting
  the body is stale and instructing the reader to "treat results.md as the
  source of truth until that pass is done" — so this is known, not hidden.
  But given how much has been layered on top since that note was written
  (ATAC fix, Geneformer, multi-cell-type transfer, model-size sweep, motif
  scan, and now the Section 1 claims above), the gap between the draft and
  the current state of the project has only grown, not shrunk. A document
  that has an abstract and a body arguing different headline numbers is not
  in a state that could be submitted, workshop or not.

## 3. [HIGH] Every number in the project is n=1 — no seeds, no variance, no CI anywhere

I grepped the entire codebase (`src/`, `scripts/`) for any variance-related
machinery: `bootstrap`, `confidence interval`, `std_dev`, repeated-seed
looping. The only `.std(...)` calls found are feature standardization inside
model wrappers (`scmultiomegrn_wrapper.py:36`, `regdiffusion_wrapper.py:70`),
not result variance. `configs/default.yaml:23` hard-codes `seed: 42`, and it
is the *only* seed value used anywhere for MEvD-GRN training, the
scMultiomeGRN baseline (`scripts/05_run_baselines.py:68`, `scripts/13_...py:247`),
and every ablation. `MEvDTrainer.__init__` (`src/training/trainer.py:32`)
seeds a single `torch.Generator` from this same config value; nothing in
`scripts/06_ablation.py` or `scripts/03_train.py` loops over multiple seeds.

Concrete consequences:
- Every AUPR/AUROC gap discussed in `results.md` — including the headline
  "MEvD-GRN beats scMultiomeGRN," every ablation gap, and every sweep-point
  ranking — is a difference between two single training runs. There is no way
  to know, from anything in this repo, whether e.g. dual-evidence AUPR
  0.881 vs. scMultiomeGRN's 0.847 (a 0.034 gap) is a real, stable effect or
  within the run-to-run noise of AdamW + dropout + random negative sampling
  on this size of network.
- The model-size sweep's own conclusion ("h384l2 to h384l3 is +36% params for
  +0.001 dual AUPR — essentially noise-level," `plan.md` line 243) *assumes*
  a noise floor without ever measuring one. Nothing in the repo runs the same
  config twice to establish what "noise-level" actually is for this pipeline;
  the claim is asserted, not verified.
- This is fixable relatively cheaply (the model trains in minutes on this
  data — `scripts/03_train.py` timing in `results/K562_results.json` shows
  full curricula finishing quickly) — a 3-5 seed repeat of the headline
  comparison would cost little and would substantially change how much weight
  the paper's central claim can bear. As it stands, no comparison in this
  project rises above "one run beat another run."

## 4. [MEDIUM-HIGH] Baseline comparisons are not apples-to-apples in several specific, checkable ways

`src/baselines/scmultiomegrn_wrapper.py`'s own docstring (lines 1-16) lists
real, deliberate deviations from Xu et al. (NAR 2025): own RNA/ATAC features
instead of the paper's MAESTRO RP score/GRNBoost2/FIMO pipeline, own
train/val/test split instead of the paper's tenfold + consensus-network
protocol, TF→all-gene instead of the published TF→TF-only task, and edges
decoded on demand instead of a materialized N×N matrix. The docstring is
honest about these, and re-scoping the published method to a different task
(TF→all-gene) is a legitimate adaptation given MEvD-GRN's own benchmark — but
it does mean every "vs. scMultiomeGRN" comparison in this project is, as
`results.md` §7 itself says, "architecture-inspired baseline on our
benchmark," not a reproduction of the paper's reported numbers. That caveat
is stated once in `results.md` but is easy to lose sight of given how many
tables present the comparison as a head-to-head "ours vs. theirs" result.

Two further asymmetries that are *not* flagged anywhere in `results.md` or
`plan.md`:

- **Training budget is wildly unequal and not compute-matched.** MEvD-GRN's
  full curriculum is 30 + 15 = 45 epochs total
  (`configs/default.yaml:99-101`), with early stopping (patience-based). The
  scMultiomeGRN baseline is explicitly run for the "full 2000-epoch budget"
  per tier (`results.md` §7: "neither early-stopped — both were still slowly
  improving at epoch 2000"). That's ~45x more epochs per tier for the
  baseline. This cuts both ways rhetorically (MEvD-GRN wins the headline
  metric despite far less training; but scMultiomeGRN's localization/
  perturbation wins could just as easily be a training-budget artifact rather
  than an architectural one) — but the point is that no version of this
  comparison controls for compute, and the paper doesn't acknowledge the gap
  either way.
- **Three of the five baselines are unsupervised, RNA-only methods being
  scored against a supervised, multi-omic link-prediction benchmark.**
  GRNBoost2 (`src/baselines/grnboost2_wrapper.py`) is a per-target
  gradient-boosted regressor on TF expression with **no access to ATAC and no
  access to any TF→gene ground-truth label at any point** — it never sees the
  evidence tiers used to construct the benchmark. The same is true of
  RegDiffusion and GMF-GAE (dense-matrix / embedding scorers, per
  `src/baselines/common.py`). Comparing them against MEvD-GRN and
  scMultiomeGRN — both of which are trained directly on the tier being
  evaluated — is a legitimate "does supervision + multi-omics help at all"
  sanity check, and it's standard practice in the GRN literature to include
  classic unsupervised methods as a floor. But `results.md` §7's summary line
  ("MEvD-GRN's multi-omic graph structure clearly matters relative to
  expression-only methods") reads as an architecture claim when the more
  accurate framing is "supervised methods with labels beat unsupervised
  methods without labels," which is a much weaker and less interesting
  statement. The paper should be explicit that these three are a different
  category of comparison, not peers of scMultiomeGRN.

## 5. [MEDIUM] The `all_at_once` re-split path is a fragile coincidence, not a checked invariant

`src/training/curriculum.py:31-57` (`build_all_at_once_stage`) merges the
**full, unsplit** evidence tensors (`data.evidence[t]`, i.e. train+val+test
together) for the main curriculum tiers and re-splits them with the old,
non-hierarchy-aware `create_edge_splits` (the function whose own docstring at
`src/data/graph_builder.py:301-305` warns it should not be used "for anything
that spans multiple evidence tiers"). This looked, on first read, like a
serious train/test leak for `all_at_once`: it re-derives train/val/test from
scratch, independently of the disk-persisted, hierarchy-safe splits in
`data/splits/*_splits.pt` that `scripts/03_train.py`'s final evaluation block
(lines 97-105) uses to score the model afterward.

I verified this empirically rather than assuming it (script below, run
against the real K562 processed data):

```
union pos (loc+pert, full evidence): 1,377,658
merged train pos (via create_edge_splits, seed=42): 964,360
localization held-out(val+test)=350,107  overlap_with_all_at_once_train=0  (0.0%)
perturbation held-out(val+test)=69,277   overlap_with_all_at_once_train=0  (0.0%)
dual_evidence held-out(val+test)=7,943   overlap_with_all_at_once_train=0  (0.0%)
```

The overlap is genuinely zero — **but only because** `create_edge_splits` and
`create_global_edge_splits` both sort edges into the same canonical order and
both seed `torch.Generator().manual_seed(42)` with matching train/val ratios,
so the two independently-written split functions happen to reproduce
*bit-identical* partitions when applied to the same edge universe. That's an
emergent property of two unrelated functions sharing a seed and an ordering
convention, not a designed guarantee — nothing in the code asserts it, tests
it, or would catch it silently breaking (e.g. if `main_curriculum_tiers` were
changed to a subset that changes the edge universe's size or ordering
relative to what `create_global_edge_splits` saw, or if this pattern were
reused for Macrophage/MCF7 where the edge-count arithmetic could differ).
`verify_no_cross_tier_leakage` (`src/data/graph_builder.py:408`) — the
function that *would* catch this — is only ever called from
`scripts/02_preprocess.py:142`, never from the `all_at_once` code path. I'd
flag this as a latent risk that happens not to have fired yet, and recommend
either (a) making `build_all_at_once_stage` derive its train pool by
filtering the *existing* per-tier train splits rather than re-splitting the
full evidence union, or (b) adding an explicit overlap assertion at
training time.

## 6. [MEDIUM] Several claimed architectural contributions show no measurable ablation gap

`results.md` §2's own ablation table is commendably honest about this, but the
paper (per `paper/main.tex`'s introduction, which foregrounds the role-aware
decoder as a headline design contribution) doesn't carry that honesty
through:

- **Role-aware decoder vs. simpler fusion**: `gated_fusion` (0.883),
  `concat_fusion` (0.890), `edge_mlp` (0.882) vs. the role-aware default
  (0.881) — all four sit within **0.01 dual AUPR** of each other, i.e.
  statistically indistinguishable even before accounting for the n=1 problem
  in Section 3. `results.md` calls this "a genuine second-order effect...
  motivated by biology and interpretability, not by a raw accuracy gap" —
  which is an honest characterization, but it means the paper's most
  biologically-motivated design choice (ATAC-as-accessibility-gate) is not
  actually supported by the evidence as an accuracy-improving component. It
  may still be worth keeping for interpretability, but it should not be
  framed as a performance contribution.
- **Learned relation gating (`combine_mode: "gated"`)** is explicitly neutral
  (0.876 vs 0.881 sum, `results.md` §2) and is defended as "the mechanism the
  real motif graph will actually make useful" — but see Section 7 below: the
  motif graph exists in the data pipeline and is never turned on in any
  shipped config, so this mechanism currently has no component to be useful
  for. It is being carried in the architecture (and the architecture
  diagrams) as a contribution whose payoff is perpetually "next."

## 7. [MEDIUM] A committed, coded component is unused by every reported configuration

`data/processed/K562/motif_edges.pt` exists, `src/data/motif_scan.py` (272
lines) implements real PWM/FIMO-style motif scanning, `src/models/mevd_grn.py`
has a working `use_motif` flag and a third relation slot, and
`scripts/06_ablation.py`'s `ABLATIONS` list includes `motif_graph_only`. But:

- `use_motif: false` in every single committed config (`configs/default.yaml:76`
  and I grepped every other config — none override it).
- There is no `results/ablations/motif_graph_only_*.json` — the ablation is
  registered in code but was never run (or was run and the output was never
  saved).
- `docs/citations.md` and `plan.md` present the motif-scanning workstream as
  substantially "closed" (a real, non-accessibility-derived third relation),
  and the gated-relation combiner (Section 6) is explicitly justified by
  reference to this graph, but no reported number anywhere in `results.md`
  reflects a model that actually used it.

This is exactly the "vestigial component presented as load-bearing" pattern
the task asked me to check for: the architecture diagrams
(`docs/figures/architecture_diagram_v3_1/v3_2`) and the design-rationale doc
both describe a 3-relation GNN, but every number in `results.md` — including
every headline and sweep number — comes from the 2-relation configuration.

## 8. [MEDIUM] Overfitting risk at 3.26M params / ~8K positives is asserted, not checked

`plan.md` §3 explicitly worries about this before scaling up ("some of our
tiers have very few positive examples ... a much bigger model risks
overfitting those few examples") and then resolves the worry purely by
observing that validation AUPR increases monotonically with model size across
the sweep. That is not a check for overfitting — it's consistent with
overfitting the *validation set* too, especially with a fixed seed, a fixed
split, and no independent re-validation. Nothing in the repo:
- plots train vs. val loss/AUPR curves over epochs for the sweep points (no
  learning-curve artifacts exist for the sweep, only final-epoch metrics in
  the table),
- re-validates the winning size on a fresh split or fresh seed,
- reports train-set AUPR alongside val/test AUPR to check for a widening gap
  at the larger sizes.
`dual_evidence` has only 7,943 positive edges total (confirmed directly from
`evidence_dual_evidence.pt`/the JSON artifacts above) being evaluated zero-shot
against a 3.26M-parameter model trained on tiers that heavily nest it. Given
Section 3's finding that n=1 noise isn't even characterized, "the gains
flatten out, so 384/2 is safely within noise of 384/3" is a plausible-sounding
but unverified claim.

## 9. Novelty assessment

The "evidence-confidence curriculum" — training on cheaper, larger-but-noisier
labels first and evaluating zero-shot on a smaller, higher-confidence held-out
tier — is a specific, reasonable application of ideas that are well
established under other names: curriculum learning, multi-fidelity
supervision, and weak-to-strong generalization. None of the individual
pieces (relational GraphSAGE, gated modality fusion, distance-weighted
regulatory potential, foundation-model embedding lookup as an auxiliary
feature) are new in themselves; MAESTRO/BETA-style regulatory potential and
Geneformer lookup embeddings are both used essentially as published. The
project's actual novel contribution is empirical and dataset-level: applying
this specific combination to SC-MO-GRN-DB's tiered evidence structure, and —
per `plan.md` §6 — apparently being the first to run a multi-cell-type
joint-training + zero-shot-third-cell-type-transfer experiment on this
dataset. That multi-cell-type result (Section 6 of `results.md`) is, on its
own evidence (setting aside the n=1 caveat), the most genuinely interesting
and least "textbook-recombination" finding in the project — it's the part I'd
lead with if reframing this for a venue.

A knowledgeable ICLR/ICML reviewer would very likely read the method
description and conclude "solid, careful engineering and honest bug-hunting,
incremental combination of known techniques, not a new idea" — which by
itself would likely be a borderline-to-negative signal for a methods-focused
main-track venue, independent of the issues in Sections 1-3.

## 10. Smaller documentation/consistency issues

- **Parameter-count mismatch**: `results.md` and `results.md`'s own §1 report
  315,015 params for the base model; `paper/main.tex` line 69 agrees
  (315,015); but `paper/main.tex`'s ablation table caption (line 824) says
  "306,691 params" for the same "shipped default" model. Three different
  numbers (306,691 / 306K-ish prose in `plan.md` / 315,015) circulate for
  what should be one fixed quantity.
- **Stale, contradicted docstring**: `src/training/curriculum.py:39-41`
  still says all_at_once "already beats the sequential 2-stage curriculum on
  every metric in prior ablations (see results.md)" — this is the
  *pre-bugfix* conclusion that `results.md` §3 explicitly documents as having
  been reversed after the leakage fix. Anyone reading the code comment in
  isolation (rather than the current `results.md`) would be actively misled
  about the project's own current conclusion.
- `configs/default.yaml`'s inline comment block (lines 79-96) is, by
  contrast, correctly updated with the reversed conclusion — so the
  contradiction is specifically in `curriculum.py`'s docstring, not
  everywhere.

## 11. What's actually solid

To be fair, and briefly: the two original bugs (EPR normalization,
hard-negative val/test leakage in `src/training/sampler.py`) were real,
non-obvious, correctly diagnosed, and the fix is legible in the code
(`get_hard_negatives`'s docstring at `src/training/sampler.py:20-45` explains
the mechanism precisely and the fix is minimal and correct). The
hierarchy-safe global split (`create_global_edge_splits`,
`src/data/graph_builder.py:335-405`) is a genuinely careful piece of
engineering — it correctly handles the fact that dual-evidence is deeply
nested inside the other tiers, includes a self-check
(`verify_no_cross_tier_leakage`) that is actually invoked in the main
preprocessing path, and the structural prior graph construction
(`build_prior_graph`) correctly excludes known positives from message-passing
edges to prevent that flavor of leakage. The project's narrative style —
writing up negative/mixed results (the FM cross-cell-type transfer penalty,
the all_at_once vs. sequential reversal) as first-class findings rather than
hiding them — is a genuine strength and above the median for a project at
this stage.

## Overall verdict

Not publication-ready at any venue right now, for a reason that is fixable
faster than the other issues: **the two most recent and most important
claims in the project cannot be reproduced from anything in the repository.**
Until (a) a config/checkpoint/log/results-JSON exists for the `all_at_once` +
FM + h384/l2 run and the full scMultiomeGRN 3-tier run, (b) the headline
comparisons are repeated across at least a handful of seeds so "beats
scMultiomeGRN" and "within noise" mean something quantitative, and (c)
`paper/main.tex`'s body is resynced to `results.md`, this is not something I
would send to any venue, workshop included. Once those are fixed, the
strongest, most defensible piece of the story is the multi-cell-type
joint-training / zero-shot-transfer result (Section 6 of `results.md`) — that
is a more novel and more interesting claim than the within-K562 architecture
ablations, and I'd suggest building the paper's contribution framing around
it rather than around the marginal (and currently unverifiable) baseline
comparisons.

---

> **FINAL STATUS UPDATE (2026-09-18, end of day)**, added by a second
> session working the same review in parallel -- all three blocking items
> in this verdict have since been addressed:
>
> (a) The `all_at_once`+FM+h384/l2 and full 3-tier scMultiomeGRN runs are
> now real, committed artifacts: `results/` was entirely gitignored (a
> deeper version of this section's finding -- `git ls-files results/`
> returned zero files, full stop, not just "these two runs are missing"),
> fixed with a `.gitignore` carve-out for the small JSON result files and
> a commit adding them, including several other results this project's
> docs cited that had never even been synced from the remote cluster.
>
> (b) The headline comparison now has real 5-seed variance data: std
> $\leq0.0017$ on every tier/metric for the recommended configuration.
> This resolved in both directions -- it confirmed 5 of 6 metrics beat
> scMultiomeGRN by tens to 100x+ the measured noise (not close calls), and
> it revealed that perturbation AUROC is a real, small, consistent *loss*
> for MEvD-GRN, which earlier single-run reporting had been too generous
> to (called a "near-tie"). scMultiomeGRN's own multi-seed rerun is still
> in progress as of this update.
>
> (c) `paper/main.tex`'s body has been fully resynced to `results.md`:
> current numbers throughout, the EPR description fixed to match the
> already-corrected code, explicit baseline-fairness caveats added, and
> Limitations expanded to cover the statistical-rigor and metric-reuse
> concerns raised elsewhere in this review. Verified it compiles cleanly.
>
> Separately, this session also found and fixed a bug neither review
> caught: MEvD-GRN's own training scripts never seeded model weight
> initialization, so "single run, seed 42" was optimistic -- it was really
> "single run, unseeded." Fixed and used for the 5-seed rerun above.
>
> The `all_at_once` split-safety concern (this document's own §5, and the
> other session's §2) is also fixed now -- built from the persisted
> per-tier splits with an explicit leak assertion, not a coincidental
> re-split.
>
> Not yet done: the motif-graph relation now has a real single-relation
> result (`motif_graph_only`, respectably close to the full model), but
> `gated_relations` was accidentally rerun without the real motif graph
> enabled (a config oversight), so whether the learned combiner does
> anything interesting with 3 real relations is still open. The
> dual-evidence-as-selection-metric concern (§ elsewhere in this review)
> remains a flagged limitation, not something resolved by carving out a
> new untouched split -- deliberately left as a call for the project's
> owner, not something to unilaterally implement.
