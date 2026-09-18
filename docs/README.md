# Docs index

Primary, actively-maintained documents live at the repo root:
`README.md` (setup/usage), `plan.md` (living roadmap), `results.md`
(current numbers).

This directory holds supporting material:

- `citations.md` — design rationale + literature citations for every model
  component.
- `figures/` — architecture diagrams. `architecture_diagram_v3_1.png` is
  current (2026-09-12: v3's top-to-bottom data-flow redesign, with the PWM
  motif prior removed since no training run has actually used it yet, and
  a clearer step-by-step role-aware decoder panel). `architecture_diagram_v3_2.png`
  is an alternate, more didactic view of the SAME architecture, restructured
  into 5 explicit bounded blocks (RNA processing / ATAC processing / prior
  graphs & ground-truth networks side by side / encoders+GraphSAGE / role-aware
  decoder) — useful when the audience needs the prior-graphs-vs-ground-truth-networks
  distinction spelled out explicitly. `architecture_diagram_v3.png`,
  `architecture_diagram_v2.png` and `architecture_diagram.png` are earlier
  iterations, kept for comparison.
- `reference/` — external reference material: the SC-MO-GRN-DB dataset
  catalog (`SC_MO_GRN_DB.md`, `dataset_info.md`), the Ada HPC cluster guide
  (`ada.md`), and the scMultiomeGRN baseline paper PDF.
- `archive/` — superseded documents, kept for historical record. Notably
  `results_2026-09-07_pre_bugfix.md`, whose numbers predate two bugfixes
  (see `results.md` Section 0) and should not be compared directly against
  current results. Also the original full theory/spec
  (`MEvD_GRN_plan_v2_superseded.md`), old presentation decks, and an early
  progress log.
- `weekly_updates/` — weekly status-update PDFs.
