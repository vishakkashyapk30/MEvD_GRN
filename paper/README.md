# MEvD-GRN paper (ICML-style LaTeX)

`main.tex` is a rigorous, self-contained write-up of the project (Sections 1-12
+ appendices), typeset with the official ICML 2024 style package
(`icml2024.sty`/`.bst`, `algorithm(ic).sty`, `fancyhdr.sty` — vendored here
verbatim from `media.icml.cc`). It is a standalone scientific paper, not a
conference submission.

## Build

```bash
cd paper
pdflatex -interaction=nonstopmode main.tex
bibtex main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

Requires a standard TeX Live install (`pdflatex`, `bibtex`). `main.pdf`
(15 pages as of the 2026-09-18 build) is committed so the paper can be read
without a LaTeX toolchain.

## Layout
- `main.tex` — paper source
- `references.bib` — bibliography (SC-MO-GRN-DB, scMultiomeGRN, GRNBoost2,
  RegDiffusion, GMFGRN, GraphSAGE, curriculum learning / catastrophic
  forgetting literature)
- `figures/` — the 8 figures `main.tex` includes, copied from
  `results/figures/` (regenerate the sources with
  `scripts/07_compile_results.py` and `scripts/10_presentation_figures.py`,
  then re-copy). `architecture_diagram_v3_1.png` is a copy of
  `docs/figures/architecture_diagram_v3_1.png`, and its HTML source lives
  there.
- `icml2024.sty`, `icml2024.bst`, `algorithm.sty`, `algorithmic.sty`,
  `fancyhdr.sty` — vendored ICML 2024 style package

## Relationship to the weekly-update decks
The slide-deck walkthroughs of the project are in `docs/weekly_updates/`
(PDF decks, some with markdown slide text). The old root-level
`presentation.md` that this section used to describe is now
`docs/weekly_updates/2026-08-02_to_2026-08-22_update_06.md`, and its numbers
are pre-bugfix. `paper/main.tex` is the peer-review-style write-up of the
current results in the root `results.md`. It is aimed at a reader who wants
the full methodological rigor (formal split-consistency guarantee, related
work, ablation discussion) rather than a slide-by-slide narrative.
