# Docs index

The primary, actively maintained documents are at the repo root:

| File | What it is |
|---|---|
| `README.md` | Setup, data download, how to run, repo layout, design notes |
| `plan.md` | Living roadmap: what changed and why, and what is still open (Section 10). Code comments cite it as "plan.md Section N", so its section numbers are stable |
| `results.md` | Current numbers (the source of truth). Each table names its `results/*.json` source |
| `paper/main.tex` | Manuscript, resynced to `results.md` on 2026-09-18. Build instructions are in `paper/README.md` |

Where documents disagree, `results.md` wins. Anything under `archive/`, and
the older weekly-update text under `weekly_updates/`, has pre-bugfix numbers
that should not be compared with current results.

## Layout of `docs/`

```
docs/
├── README.md                          this index
├── citations.md                       design rationale + literature citations, one section per model component
├── critical_review_independent.md     2026-09-18 independent critical review + its final status update
├── figures/                           architecture diagrams (PNG + HTML source)
├── reference/                         external reference material
├── experiments/                       runbooks + findings for competitor-dataset benchmarks
├── weekly_updates/                    weekly-update decks (PDF) + markdown text for some of them
└── archive/                           superseded documents, kept for the record
```

### `citations.md`
One section per design choice (dataset and tiers, ATAC featurization, encoders,
role-aware decoder, GraphSAGE, prior graphs, curriculum, replay, splits,
negatives, metrics, baselines, edge MLP, learned gated relation combiner,
Geneformer embeddings, motif scanning, multi-cell-type joint training). Each
says what we do, why, and which paper it comes from.

### `critical_review_independent.md`
A skeptical, code-first review done on 2026-09-18 to check publication
readiness. The original findings are kept unedited, with dated updates added
afterwards. Most blocking items are resolved (see its final status block and
`plan.md` Section 9). The items still open are tracked in `plan.md`
Section 10. (`scripts/06_ablation.py` refers to this file by name.)

### `figures/`
| File | Status |
|---|---|
| `architecture_diagram_v3_1.png` / `.html` | **Current** (2026-09-12). A top-to-bottom data-flow redesign with a clearer step-by-step role-aware decoder panel, in Poppins with a validated colorblind-safe palette. Its subtitle gives 315,015 base params / 430,471 with Geneformer (the base count was corrected from a stale 306,691 on 2026-09-18). The paper uses the copy at `paper/figures/architecture_diagram_v3_1.png`. It does **not** draw the PWM motif relation, which was left out because no training run had used it at the time. |
| `architecture_diagram_v3_2.png` / `.html` | Alternate, more didactic view of the **same** architecture. It is split into 5 labeled blocks (RNA processing / ATAC processing / prior graphs and ground-truth networks side by side / encoders + GraphSAGE / role-aware decoder), which helps when the audience needs the prior-graph vs. ground-truth-network distinction spelled out. |
| `architecture_diagram.png` / `.html` | The **original** (pre-2026-09-11) architecture, kept for comparison. Its name is unchanged because the update 06 and 08 texts in `weekly_updates/` point to it. |

The `.html` files are the diagram sources, and each PNG is a rendering of its
HTML. The interim `architecture_diagram_v2` and `architecture_diagram_v3` were
deleted on 2026-09-29 as superseded. v3 drew the motif relation as a real
component, while v2 showed it only as "planned". Both can be recovered from
git history, e.g.
`git show af8185c:docs/figures/architecture_diagram_v3.html`.

### `reference/`
External reference material: the SC-MO-GRN-DB dataset catalog
(`SC_MO_GRN_DB.md`, `dataset_info.md`), the Ada HPC cluster guide (`ada.md`),
the scMultiomeGRN baseline paper (`scMultiomeGRN.pdf`), and
`literature_knowledge_base.md`. That last file covers the 17 competitor papers
and theses in `~/research/`. It has a section per paper, a catalogue of which
datasets each paper trains and evaluates on (and whether paired RNA+ATAC
exists), the reported numbers MeVD-GRN must beat, and protocol caveats.

### `experiments/`
One runbook per competitor benchmark. It records the dataset, the competitor's
exact protocol, the target numbers, the files that implement it, and the Ada
launch commands.
- `scmultiomegrn_generalization.md`: MeVD-GRN on scMultiomeGRN's fetal-lung
  benchmark (GSM4508936). Built and smoke-tested on 2026-09-30; not yet run on
  Ada.

### `weekly_updates/`
Files are named `<start>_to_<end>_update_<NN>.<ext>` with ISO dates, so they
sort chronologically. The dates come from the original filenames. The title
slides of 04 and 07 show slightly different ranges (10-21 July and 23 Aug - 5 Sept).

| Update | Period | Files |
|---|---|---|
| 01 | 2026-05-14 to 2026-06-25 | `.pdf` |
| 02 | 2026-06-26 to 2026-07-02 | `.pdf` |
| 03 | 2026-07-03 to 2026-07-09 | `.pdf` |
| 04 | 2026-07-10 to 2026-07-23 | `.pdf` |
| 05 | 2026-07-24 to 2026-08-01 | `.pdf` |
| 06 | 2026-08-02 to 2026-08-22 | `.pdf` + `.md` (slide text, formerly `docs/archive/presentation.md`) |
| 07 | 2026-08-22 to 2026-09-05 | `.pdf` + `.md` (slide text, formerly `docs/archive/presentation_5th_sept_2026.md`) |
| 08 | 2026-09-06 to 2026-09-12 | `.pdf` + `.md` (long-form write-up, formerly `presentation_this_week.md`) |

The PDFs are the Canva decks as delivered. The 06 and 07 `.md` files carry
pre-bugfix numbers. The 08 `.md` has numbers as of 2026-09-12, except its
scMultiomeGRN rows, which were later updated to the full 3-tier run. It still
lists motif scanning and the model-size sweep as not done, and both have since
finished (`plan.md` Sections 3 and 10).

### `archive/`
Superseded documents, named `<date>_<description>.md`. Each one starts with a
banner saying what superseded it.

| File | What it is |
|---|---|
| `2026-07-02_progress_log_session_1.md` | Handoff log from the first implementation session: first prototype numbers (incl. ESC), environment notes, and the baseline-strategy notes (supervision and modality brackets) |
| `2026-07-20_mevd_grn_plan_v2_original_spec.md` | The original full theory and implementation spec (v2) the code was first built from |
| `2026-09-07_results_pre_bugfix.md` | Results catalog from before the 2026-09-11 EPR and hard-negative-leakage fixes. **Do not compare with `results.md`** |
| `2026-09-07_toy_pipeline_walkthrough.md` | Shape-by-shape walkthrough of the data flow, from before the 2026-09-11 ATAC rewrite |

The local-only `status.md` TODO tracker, which was gitignored, was folded into
`plan.md` Section 10 on 2026-09-29.
