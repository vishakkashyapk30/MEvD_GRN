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
(12 pages) is committed so the paper can be read without a LaTeX toolchain.

## Layout
- `main.tex` — paper source
- `references.bib` — bibliography (SC-MO-GRN-DB, scMultiomeGRN, GRNBoost2,
  RegDiffusion, GMFGRN, GraphSAGE, curriculum learning / catastrophic
  forgetting literature)
- `figures/` — 8 figures copied from `results/figures/` (regenerate the
  sources with `scripts/07_compile_results.py` and
  `scripts/10_presentation_figures.py`, then re-copy)
- `icml2024.sty`, `icml2024.bst`, `algorithm.sty`, `algorithmic.sty`,
  `fancyhdr.sty` — vendored ICML 2024 style package

## Relationship to `presentation.md`
`presentation.md` (repo root) is the slide-deck-oriented walkthrough of the
same project, written for making a PowerPoint. `paper/main.tex` is the
peer-review-style write-up of the same underlying results, aimed at a reader
who wants the full methodological rigor (formal split-consistency guarantee,
related work, ablation discussion) rather than a slide-by-slide narrative.
