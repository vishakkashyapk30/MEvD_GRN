#!/usr/bin/env python
"""Build docs/meetings/presentation_15th_nov_22nd_nov.md (advisor feedback of 19 Nov 2026) from compiled results:
  docs/results/bear/matched.csv (+ _per_dataset.csv)        scripts/bear/27_bear_matched_table.py
  docs/results/scmgrn/final/*.json                           scripts/scmgrn/25_scmgrn_baseline_table.py
Rows with missing results print "pending". Run: python docs/meetings/build_deck_11_22.py
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
REPO = HERE.parents[1]
OUT = HERE / "presentation_15th_nov_22nd_nov.md"
_s = importlib.util.spec_from_file_location("d14", HERE / "build_deck_11_14.py")
d14 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(d14)
_t = importlib.util.spec_from_file_location("t27", REPO / "scripts/bear/27_bear_matched_table.py")
t27 = importlib.util.module_from_spec(_t)
_t.loader.exec_module(t27)

BM = REPO / "docs/results/bear/matched.csv"
B = pd.read_csv(BM) if BM.exists() else pd.DataFrame(columns=["gt", "method", "kind"])
BKM = REPO / "docs/results/bear/matched_kocov.csv"
BK = pd.read_csv(BKM) if BKM.exists() else pd.DataFrame(columns=["gt", "method", "kind"])
BCOLS = [("AUROC_inspace", "AUROC"), ("AUROC_perTF", "AUROC per TF"), ("AUROC_degmatched", "degree-matched"),
         ("AUROC_caliper", "caliper-matched"), ("AUROC_adjusted", "covariate-adjusted"),
         ("AUROC_aiyappa", "degree-corrected")]


def f3(x):
    return "pending" if x is None or x != x else f"{x:.3f}"


def bget(gt, m, col, D=None):
    D = B if D is None else D
    r = D[(D["gt"] == gt) & (D["method"] == m)]
    return float(r[col].iloc[0]) if len(r) and col in r else None


def bn(gt, m, D=None):
    D = B if D is None else D
    r = D[(D["gt"] == gt) & (D["method"] == m)]
    return int(r["n_datasets"].iloc[0]) if len(r) else 0


def btable(gt, items, cols=BCOLS, D=None):
    """items: [(method, label)]; bold = best, italics = second best per column among the shown rows."""
    vals = {m: [bget(gt, m, c, D) for c, _ in cols] for m, _ in items}
    best = []
    for j in range(len(cols)):
        col = sorted({round(v[j], 3) for v in vals.values() if v[j] is not None and v[j] == v[j]}, reverse=True)
        best.append(col[:2])
    out = ["| Method | Datasets | " + " | ".join(h for _, h in cols) + " |", "|---|---|" + "---|" * len(cols)]
    for m, lab in items:
        cells = []
        for j, x in enumerate(vals[m]):
            t = f3(x)
            if x is not None and x == x and best[j]:
                if round(x, 3) == best[j][0]:
                    t = f"**{t}**"
                elif len(best[j]) > 1 and round(x, 3) == best[j][1]:
                    t = f"*{t}*"
            cells.append(t)
        out.append(f"| {lab} | {bn(gt, m, D) or '–'} | " + " | ".join(cells) + " |")
    return "\n".join(out)


RELEASED = [("LINGER", "LINGER"), ("TRIPOD", "TRIPOD"), ("CellOracle", "CellOracle"), ("SCENIC+", "SCENIC+"),
            ("FigR", "FigR"), ("GRaNIE", "GRaNIE"), ("Pando-GLM", "Pando-GLM"), ("Pando-XGB", "Pando-XGB")]
BASE = [("baseline_indegree", "In-degree (training TFs' links)"), ("baseline_pearson", "Pearson"),
        ("baseline_grnboost2", "GRNBoost2"), ("baseline_coverage", "All measured pairs")]
MEVD = [("mevd_fm_h384_hubmotif", "MeVD-GRN"), ("mevd_fm_h384_hubmotif_L2", "MeVD-GRN, L2 (no cell-type labels)")]
ATAC_BEAR = [("mevd_fm_h384_hubmotif", "MeVD-GRN, full data"),
             ("mevd_fm_h384_hubmotif_atacnone", "MeVD-GRN, ATAC removed (all four routes)"),
             ("mevd_fm_h384_hubmotif_atacshuffle", "MeVD-GRN, ATAC shuffled across genes"),
             ("LINGER", "LINGER, released network (BEAR-GRN)"),
             ("linger_rerun_full", "LINGER re-run, full data"),
             ("linger_rerun_peakshuffle", "LINGER re-run, ATAC peaks shuffled"),
             ("linger_rerun_random", "LINGER re-run, ATAC fully randomised")]


def nparams(model, nf):
    cov = nf * (nf + 1) // 2
    return {"lr": nf + 1, "lda": 2 * nf + cov, "qda": 2 * (nf + cov), "gnb": 4 * nf, "gmm": 2 * 4 * (nf + cov + 1),
            "rbf": 301, "mlp": nf * 32 + 32 + 32 + 1, "hgb": "200 trees"}[model]


SIMPLE = []
for mdl in ["lr", "lda", "qda", "gnb", "gmm", "rbf", "mlp", "hgb"]:
    for fs, nf in (("all", 20), ("data", 19)):
        SIMPLE.append((f"simple_{mdl}_{fs}", f"{t27.MODEL[mdl]}, {'all 20 features' if fs == 'all' else 'cell data only'}"
                       f" ({nparams(mdl, nf)}{'' if mdl == 'hgb' else ' parameters'})"))
SIMPLE_KO = [(f"simple_{m}_{fs}_perturbation", f"{t27.MODEL[m]}, {'all 20 features' if fs == 'all' else 'cell data only'}, "
              "trained on knockout labels only") for m in ["lr", "lda", "gnb", "hgb"] for fs in ("data", "all")]
LABELFREE = [("lf_ggm_ridge", "Gaussian graphical model (ridge neighbourhood regression)"),
             ("lf_regdiffusion", "RegDiffusion (diffusion model with a learned gene x gene adjacency)"),
             ("lf_deepsem", "DeepSEM (sparse structural-equation VAE)"),
             ("baseline_pearson", "Pearson"), ("baseline_grnboost2", "GRNBoost2"),
             ("SCENIC+", "SCENIC+ (released)"), ("LINGER", "LINGER (released)")]

SMG_ATAC = [(d14.MEVD, "MeVD-GRN, full data"), ("mevd_fm_h384_obs_atacnone", "MeVD-GRN, every ATAC input removed"),
            ("mevd_fm_h384_obs_atacshuffle", "MeVD-GRN, every ATAC input shuffled across genes"),
            (d14.SMG, "scMultiomeGRN, full data"), ("official_scmultiomegrn_scrna_only", "scMultiomeGRN, ATAC removed (RNA only)"),
            ("official_scmultiomegrn_atac_shuffle", "scMultiomeGRN, ATAC shuffled across TFs")]
LADDER = [("lr_degree_only", "Logistic regression on the two link counts (4 parameters)"),
          ("ladder_beta", "beta-model (one parameter per TF)"),
          ("ladder_lsm2", "Latent-space model, 2 dimensions"), ("ladder_lsm4", "Latent-space model, 4 dimensions"),
          ("ladder_lsm8", "Latent-space model, 8 dimensions"),
          ("ladder_dclsm2", "Degree-corrected latent-space model, 2 dimensions"),
          ("ladder_dclsm8", "Degree-corrected latent-space model, 8 dimensions"),
          ("ladder_lda", "LDA on cell-data pair features (Gaussian)"), ("ladder_qda", "QDA on cell-data pair features"),
          ("ladder_gnb", "Gaussian naive Bayes on cell-data pair features"),
          ("ladder_gmm", "Gaussian mixture classifier on cell-data pair features")]


def best_of(gt, items, col):
    xs = [(bget(gt, m, col), lab) for m, lab in items if bget(gt, m, col) is not None]
    xs = [x for x in xs if x[0] == x[0]]
    return max(xs) if xs else (None, "–")


def main():
    S = []
    all_methods = RELEASED + BASE + MEVD + SIMPLE + LABELFREE[:3]
    bc = best_of("ChIP", all_methods, "AUROC_caliper")
    ba = best_of("ChIP", all_methods, "AUROC_adjusted")
    S.append(f"""# Degree-controlled evaluation, ATAC value, simple models and non-GNN models

Advisor meeting, follow-up on the feedback of 19 November 2026.

## 1. Summary

- **Citations (item 1).** Degree-matched test negatives are published (Aiyappa et al., ICML 2025; Cappelletti et al.,
  Bioinformatics Advances 2024; BioGraphBench-PPI, Frontiers in Bioinformatics 2026). The caliper test is the
  covariate-adjusted ROC of biostatistics (Janes & Pepe, Biometrika 2009) with degree as the covariate. No
  link-prediction paper uses a caliper on degree; the statistical method is established.
- **BEAR-GRN under degree control (item 2).** On ChIP, the best caliper-matched AUROC of any method is
  {f3(bc[0])} ({bc[1]}); the best covariate-adjusted AUROC is {f3(ba[0])} ({ba[1]}). MeVD-GRN:
  {f3(bget('ChIP', 'mevd_fm_h384_hubmotif', 'AUROC_caliper'))} caliper-matched,
  {f3(bget('ChIP', 'mevd_fm_h384_hubmotif', 'AUROC_adjusted'))} covariate-adjusted.
- **ATAC value (item 3), MeVD-GRN on BEAR-GRN, ChIP AUROC:** full {f3(bget('ChIP', 'mevd_fm_h384_hubmotif', 'AUROC_inspace'))},
  ATAC removed {f3(bget('ChIP', 'mevd_fm_h384_hubmotif_atacnone', 'AUROC_inspace'))}, ATAC shuffled
  {f3(bget('ChIP', 'mevd_fm_h384_hubmotif_atacshuffle', 'AUROC_inspace'))}.
- **Simple models (item 4) and non-GNN models (item 5):** slides 7-10.""")

    S.append("""## 2. Sources for the degree-controlled tests (item 1)

| Test | Closest published method | Published, or own |
|---|---|---|
| Degree-matched test (each test link paired with a non-link of similar endpoint degrees) | Degree-corrected negatives, Aiyappa et al., ICML 2025; degree-aware evaluation negatives, Cappelletti et al., Bioinformatics Advances 2024; "degree-matched" negatives preserving the endpoint log-degree-bin distribution, BioGraphBench-PPI, Frontiers in Bioinformatics 2026 | Published; ours matches each pair instead of the distribution |
| Caliper-matched test (non-link within 0.1 log-degree on both endpoints, else the link is dropped) | Matched case-control accuracy studies and the covariate-adjusted ROC (Janes & Pepe, Biometrics 2008, Biometrika 2009; Pepe, Fan & Seymour, Academic Radiology 2013); caliper matching (Rosenbaum & Rubin 1985; Austin 2011); dropped units (Stuart 2010) | Published statistical method, applied with degree as the covariate |
| AUROC within in-degree groups (BEAR-GRN) | Covariate-adjusted AUC = AUROC within covariate strata weighted by positives (Janes & Pepe 2009; Pepe et al. 2013); popularity-percentile evaluation (Bellogin et al., Information Retrieval Journal 2017) | Published statistical method |
| Degree-corrected test on BEAR-GRN (new) | Aiyappa et al., ICML 2025, directly | Published |

- No paper was found that applies a caliper on node degree to link-prediction negatives; the honest claim is the
  transfer of standard matching methods to a covariate (degree) that recent work shows confounds the evaluation.
- Full survey (33 sources, verification log): `docs/reference/lit_survey/degree_controlled_evaluation.md`.""")

    S.append(f"""## 3. BEAR-GRN under degree control: ChIP-seq ground truth (item 2)

- Covariate = target in-degree among the training TFs of the fold that holds the scored TF out (never the scored TF's
  own labels). Every method is scored on identical pairs; non-output pairs score 0.
- Degree-matched and caliper-matched pair each link (t, g) with a non-link (t, g') of the same TF, so the TF's own
  number of targets cannot help either. Covariate-adjusted = AUROC within (TF, in-degree group) cells, weighted by
  links. Degree-corrected = non-links drawn with the links' endpoint-degree distributions (Aiyappa et al.).
- The in-degree baseline is the covariate, so it scores about 0.5 on the caliper-matched and covariate-adjusted tests.

{btable('ChIP', RELEASED + BASE + MEVD)}""")

    S.append(f"""## 4. BEAR-GRN under degree control: knockout ground truth (item 2)

{btable('KO', RELEASED + BASE + MEVD)}""")

    S.append(f"""## 5. How much ATAC adds on BEAR-GRN (item 3)

- MeVD-GRN uses ATAC in four places: the ATAC features, promoter openness, the peak-derived motif features and the
  accessibility filter of the TF-candidate graph. "Removed" cuts all four; "shuffled" gives every gene another gene's
  ATAC values (same distribution, broken gene mapping). Seeds 42-44.
- LINGER is re-run on Macrophage_S1 with LINGER's own bulk reference: full data, ATAC peaks shuffled (genomic location
  broken, cell structure kept) and ATAC fully randomised (no ATAC information left).

ChIP-seq ground truth:

{btable('ChIP', ATAC_BEAR)}""")

    S.append(f"""## 6. How much ATAC adds on the scMultiomeGRN benchmark (item 3)

- The ground-truth network of this benchmark is itself built from the ATAC peaks (motif in an open promoter), so ATAC
  could only look more useful here than it is; BEAR-GRN (previous slide) is the clean test.

{d14.table(d14.OFF, SMG_ATAC, d14.OFF_COLS, d14.OFF_HEAD)}""")

    S.append(f"""## 7. Simple and Gaussian models on BEAR-GRN, ChIP (item 4)

- Same TF-held-out folds as MeVD-GRN (seed 42), so no label of a scored TF is used. 20 interpretable features per
  (TF, gene) pair: target in-degree among training TFs (left out for the pair's own link), Pearson correlation,
  GRNBoost2 importance, RNA statistics of TF and gene, the gene's ATAC and promoter-openness summaries, 3 motif scores,
  Geneformer cosine. Training negatives are non-links of the same TF, so the models rank genes within a TF.
- Their scores are not calibrated across held-out TFs, so the per-TF and matched columns are the fair comparison.

{btable('ChIP', [('baseline_indegree', 'In-degree'), ('mevd_fm_h384_hubmotif', 'MeVD-GRN')] + SIMPLE)}""")

    S.append(f"""## 8. Simple and Gaussian models on BEAR-GRN, knockout ground truth (item 4)

- Training labels as for MeVD-GRN: binding (ChIP-derived) and knockout labels pooled. Binding labels dominate, and
  the models then rank knockout targets below chance within a TF. The second block trains the same models on the
  knockout labels only (6 datasets with knockout labels; same folds).

{btable('KO', [('baseline_indegree', 'In-degree'), ('baseline_pearson', 'Pearson'), ('mevd_fm_h384_hubmotif', 'MeVD-GRN')] + SIMPLE)}

Trained on knockout labels only:

{btable('KO', [('mevd_fm_h384_hubmotif', 'MeVD-GRN'), ('simple_deg_all_perturbation', 'In-degree rule, knockout labels only (no fit)')] + SIMPLE_KO)}

- The matched tests control the in-degree counted over all label files. Trained on knockout labels only, the models
  pick up "genes that many knockouts affect", a hub of the knockout network that this covariate does not control;
  the in-degree rule on knockout labels alone measures how much of the gain that is.

Same knockout ground truth, matched tests re-drawn with the covariate = target in-degree over the knockout labels
of the training TFs (5 datasets):

{btable('KO', [('simple_deg_all_perturbation', 'In-degree rule, knockout labels only (= covariate)'),
               ('baseline_indegree', 'In-degree, all labels'), ('baseline_pearson', 'Pearson'),
               ('mevd_fm_h384_hubmotif', 'MeVD-GRN'),
               ('simple_lr_data', 'Logistic regression, cell data only, all labels'),
               ('simple_lda_data', 'LDA, cell data only, all labels')] + SIMPLE_KO, D=BK)}""")

    S.append(f"""## 9. Parameter ladder on the scMultiomeGRN benchmark (item 4)

- Network models fitted on the training links only: the beta-model (one parameter per TF; Chatterjee, Diaconis & Sly
  2011), the latent-space model (Hoff, Raftery & Handcock 2002) and its degree-corrected form (Krivitsky et al. 2009).
  Gaussian classifiers use MAESTRO and GRNBoost2 pair features only (32 numbers, no training-graph input).

Official split:

{d14.table(d14.OFF, [(d14.MEVD, 'MeVD-GRN')] + LADDER, d14.OFF_COLS, d14.OFF_HEAD)}

Held-out-TF split:

{d14.table(d14.TF, [(d14.MEVD_TF, 'MeVD-GRN')] + LADDER, d14.TF_COLS, d14.TF_HEAD)}""")

    S.append(f"""## 10. Label-free, non-GNN models on BEAR-GRN (item 5)

- No labels at all, one model per dataset scores every TF, as the released BEAR-GRN methods do.

ChIP-seq ground truth:

{btable('ChIP', LABELFREE)}

Knockout ground truth:

{btable('KO', LABELFREE)}""")

    S.append("""## 11. Conclusions and next steps

- To be written from the final numbers.""")
    OUT.write_text("\n\n".join(S) + "\n")
    print(f"wrote {OUT} ({len(S)} sections)")


if __name__ == "__main__":
    main()
