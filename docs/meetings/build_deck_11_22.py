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
BPD = REPO / "docs/results/bear/matched_per_dataset.csv"
MS1 = (pd.read_csv(BPD).query("dataset == 'Macrophage_S1'").assign(n_datasets=1) if BPD.exists()
       else pd.DataFrame(columns=["gt", "method", "kind"]))
BLS = REPO / "docs/results/bear/matched_linger_per_seed.csv"
LPS = pd.read_csv(BLS) if BLS.exists() else pd.DataFrame(columns=["gt", "method", "dataset", "seed"])
LMETA = (pd.read_csv(REPO / "docs/results/bear/linger_rerun_meta.csv")
         if (REPO / "docs/results/bear/linger_rerun_meta.csv").exists() else pd.DataFrame())
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
ATAC_BEAR = [("mevd_fm_h384_hubmotif_seeds42to44", "MeVD-GRN, full data"),
             ("mevd_fm_h384_hubmotif_atacnone", "MeVD-GRN, ATAC removed (all four routes)"),
             ("mevd_fm_h384_hubmotif_atacshuffle", "MeVD-GRN, ATAC shuffled across genes")]
LINGER_COND = [("LINGER", "Released (BEAR-GRN)"), ("linger_rerun_full", "Re-run, full data"),
               ("linger_rerun_peakshuffle", "Re-run, ATAC peaks shuffled"),
               ("linger_rerun_random", "Re-run, ATAC fully randomised")]
LINGER_DS = ["Macrophage_S1", "Macrophage_S2", "K562", "iPS"]
NEG = [("baseline_indegree", "In-degree (training TFs' links)"),
       ("mevd_fm_h384_hubmotif_seeds42to44", "MeVD-GRN, degree-matched negatives (headline)"),
       ("mevd_fm_h384_hubmotif_sametf", "MeVD-GRN, same-TF negatives"),
       ("simple_lr_data", "Logistic regression, cell data only (20 parameters; same-TF negatives)"),
       ("simple_lda_data", "LDA, cell data only (same-TF negatives)")]


def linger_rows(gt="ChIP"):
    """(dataset, seed) pairs that have every re-run condition, in LINGER_DS order."""
    L = LPS[(LPS["gt"] == gt) & LPS["method"].str.startswith("linger_rerun_")]
    out = []
    for ds in LINGER_DS:
        x = L[L["dataset"] == ds]
        for sd in sorted(x["seed"].unique()):
            if set(x[x["seed"] == sd]["method"]) >= {m for m, _ in LINGER_COND[1:]}:
                out.append((ds, int(sd)))
    return out


def lget(ds, seed, m, col, gt="ChIP"):
    L = LPS[(LPS["gt"] == gt) & (LPS["method"] == m) & (LPS["dataset"] == ds)]
    if m != "LINGER":
        L = L[L["seed"] == seed]
    return float(L[col].mean()) if len(L) else None


def lmeta(ds, seed):
    if LMETA.empty:
        return "–"
    x = LMETA[(LMETA["dataset"] == ds) & (LMETA["seed"] == seed)]
    if not len(x):
        return "–"
    k, mc, c = int(x["k_pseudo_labels"].iloc[0]), int(x["metacells"].iloc[0]), int(x["cells"].iloc[0])
    return f"{c} cells, one population, {mc} metacells" if k == 1 else f"{c} cells, k-means k = {k}, {mc} metacells"


def linger_table(col, gt="ChIP"):
    rows = linger_rows(gt)
    head = "| Dataset | Seed | Pseudo-labels | " + " | ".join(l for _, l in LINGER_COND) + " |"
    out = [head, "|" + "---|" * (3 + len(LINGER_COND))]
    acc = {m: [] for m, _ in LINGER_COND}
    for ds, sd in rows:
        v = [lget(ds, sd, m, col, gt) for m, _ in LINGER_COND]
        for (m, _), x in zip(LINGER_COND, v):
            if x is not None:
                acc[m].append(x)
        out.append(f"| {ds} | {sd} | {lmeta(ds, sd)} | " + " | ".join(f3(x) for x in v) + " |")
    if rows:
        out.append(f"| **Mean over the {len(rows)} runs** | | | " +
                   " | ".join(f3(sum(a) / len(a)) if a else "–" for a in acc.values()) + " |")
    return "\n".join(out)


def lmean(m, col, gt="ChIP"):
    v = [lget(ds, sd, m, col, gt) for ds, sd in linger_rows(gt)]
    v = [x for x in v if x is not None]
    return sum(v) / len(v) if v else None


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
  {f3(bget('ChIP', 'mevd_fm_h384_hubmotif', 'AUROC_adjusted'))} covariate-adjusted, level with the in-degree rule. On
  KO, the in-degree counted over knockout labels alone scores
  {f3(bget('KO', 'simple_deg_all_perturbation', 'AUROC_caliper'))} caliper-matched; with it as the covariate MeVD-GRN
  falls to {f3(bget('KO', 'mevd_fm_h384_hubmotif', 'AUROC_caliper', BK))}.
- **ATAC value (item 3):** no model gains from ATAC. MeVD-GRN on BEAR-GRN ChIP, caliper-matched AUROC: full
  {f3(bget('ChIP', 'mevd_fm_h384_hubmotif', 'AUROC_caliper'))}, ATAC removed
  {f3(bget('ChIP', 'mevd_fm_h384_hubmotif_atacnone', 'AUROC_caliper'))}, ATAC shuffled
  {f3(bget('ChIP', 'mevd_fm_h384_hubmotif_atacshuffle', 'AUROC_caliper'))} (9 datasets, 3 seeds). LINGER re-run on
  Macrophage_S1: full {f3(bget('ChIP', 'linger_rerun_full', 'AUROC_caliper', MS1))}, peaks shuffled
  {f3(bget('ChIP', 'linger_rerun_peakshuffle', 'AUROC_caliper', MS1))}, fully randomised
  {f3(bget('ChIP', 'linger_rerun_random', 'AUROC_caliper', MS1))}. MeVD-GRN and scMultiomeGRN are also unchanged on
  the scMultiomeGRN benchmark.
- **Simple models (item 4):** on BEAR-GRN ChIP a 20-parameter logistic regression on cell data reaches
  {f3(bget('ChIP', 'simple_lr_data', 'AUROC_caliper'))} caliper-matched AUROC (MeVD-GRN
  {f3(bget('ChIP', 'mevd_fm_h384_hubmotif', 'AUROC_caliper'))}); on scMultiomeGRN the beta-model matches MeVD-GRN on the
  official test. Slides 8, 9 and 11.
- **Label-free, non-GNN models (item 5):** DeepSEM, RegDiffusion and a ridge Gaussian graphical model score at the
  level of the released methods. Slide 12.""")

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
  ATAC values (same distribution, broken gene mapping). All three rows: the same seeds 42-44 on all 9 datasets.

MeVD-GRN, ChIP-seq ground truth:

{btable('ChIP', ATAC_BEAR)}

MeVD-GRN, knockout ground truth:

{btable('KO', ATAC_BEAR)}""")

    lr = linger_rows()
    lds = sorted({d for d, _ in lr}, key=LINGER_DS.index)
    S.append(f"""## 6. LINGER re-run under the ATAC conditions (item 3)

- LINGER (LingerGRN 1.110) re-run on BEAR-GRN's own input matrices with LINGER's bulk reference: full data, ATAC
  peaks shuffled (each peak takes another peak's cell profile: genomic location broken, cell structure kept) and ATAC
  fully randomised (every value permuted across peaks and cells: no ATAC information left).
- Datasets: the human BEAR-GRN datasets ({', '.join(lds)}); LINGER's bulk reference is human only.
- Pseudo-labels: LINGER builds floor(sqrt(N)) + 1 metacells per cell label and its training needs >= 50. On iPS a
  fixed random 2,500 of its 5,711 cells give 51 metacells from one population, so LINGER runs there as published,
  with no split (all 5,711 cells need about twice the memory and time). The other datasets are one population of
  412-529 cells (21-24 metacells), so their cells are split by k-means on RNA PCA (no labels used).
- Seeds: LINGER fixes its training seed at 42. Seed s > 0 changes the k-means split, the metacell draw, the ATAC
  permutation and LINGER's training seed (42 + s).

Caliper-matched AUROC, ChIP-seq ground truth:

{linger_table('AUROC_caliper')}

AUROC (no hub control), ChIP-seq ground truth:

{linger_table('AUROC_inspace')}

Covariate-adjusted AUROC, ChIP-seq ground truth:

{linger_table('AUROC_adjusted')}""")

    S.append(f"""## 7. How much ATAC adds on the scMultiomeGRN benchmark (item 3)

- The ground-truth network of this benchmark is itself built from the ATAC peaks (motif in an open promoter), so ATAC
  could only look more useful here than it is; BEAR-GRN (previous slide) is the clean test.

{d14.table(d14.OFF, SMG_ATAC, d14.OFF_COLS, d14.OFF_HEAD)}""")

    S.append(f"""## 8. Simple and Gaussian models on BEAR-GRN, ChIP (item 4)

- Same TF-held-out folds as MeVD-GRN (seed 42), so no label of a scored TF is used. 20 interpretable features per
  (TF, gene) pair: target in-degree among training TFs (left out for the pair's own link), Pearson correlation,
  GRNBoost2 importance, RNA statistics of TF and gene, the gene's ATAC and promoter-openness summaries, 3 motif scores,
  Geneformer cosine. Training negatives are non-links of the same TF, so the models rank genes within a TF.
- Their scores are not calibrated across held-out TFs, so the per-TF and matched columns are the fair comparison.

{btable('ChIP', [('baseline_indegree', 'In-degree'), ('mevd_fm_h384_hubmotif', 'MeVD-GRN')] + SIMPLE)}""")

    S.append(f"""## 9. Simple and Gaussian models on BEAR-GRN, knockout ground truth (item 4)

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

    S.append(f"""## 10. MeVD-GRN retrained with the simple models' negatives

- Only the training negatives change. Headline MeVD-GRN: a negative's TF is drawn by out-degree and its target by
  in-degree, independently, so a TF's negatives match its links only on average and targets follow the in-degree.
  Same-TF: every TF gets non-links in proportion to its links, targets drawn uniformly among its non-links (the rule
  of the simple models on slides 8-9). On dense datasets a hub TF has few non-links and gets all of them.
- Same 9 datasets, TF-held-out folds, model and seeds 42-44 as the headline rows shown.
- Result: on ChIP, same-TF negatives raise MeVD-GRN's caliper-matched and covariate-adjusted AUROC by about 0.02, to
  just below the logistic regression; on KO nothing changes, and with the knockout-label covariate both versions stay
  at chance.
- Result: on ChIP, same-TF negatives raise MeVD-GRN's caliper-matched and covariate-adjusted AUROC by about 0.02, to
  just below the logistic regression; on KO nothing changes, and with the knockout-label covariate both versions stay
  at chance.

ChIP-seq ground truth:

{btable('ChIP', NEG)}

Knockout ground truth:

{btable('KO', NEG)}

Knockout ground truth, covariate = target in-degree over the knockout labels of the training TFs (5 datasets):

{btable('KO', [NEG[1], NEG[2], ('simple_lda_data', 'LDA, cell data only, all labels'), ('simple_lda_data_perturbation', 'LDA, cell data only, knockout labels only')], D=BK)}""")

    S.append(f"""## 11. Parameter ladder on the scMultiomeGRN benchmark (item 4)

- Network models fitted on the training links only: the beta-model (one parameter per TF; Chatterjee, Diaconis & Sly
  2011), the latent-space model (Hoff, Raftery & Handcock 2002) and its degree-corrected form (Krivitsky et al. 2009).
  Gaussian classifiers use MAESTRO and GRNBoost2 pair features only (32 numbers, no training-graph input).

Official split:

{d14.table(d14.OFF, [(d14.MEVD, 'MeVD-GRN')] + LADDER, d14.OFF_COLS, d14.OFF_HEAD)}

Held-out-TF split:

{d14.table(d14.TF, [(d14.MEVD_TF, 'MeVD-GRN')] + LADDER, d14.TF_COLS, d14.TF_HEAD)}""")

    S.append(f"""## 12. Label-free, non-GNN models on BEAR-GRN (item 5)

- No labels at all, one model per dataset scores every TF, as the released BEAR-GRN methods do.

ChIP-seq ground truth:

{btable('ChIP', LABELFREE)}

Knockout ground truth:

{btable('KO', LABELFREE)}""")

    M = "mevd_fm_h384_hubmotif"
    o = lambda m, k: f3(d14.v(d14.OFF, m, k))
    t = lambda m, k: f3(d14.v(d14.TF, m, k))
    rel_c = best_of("ChIP", RELEASED, "AUROC_caliper")
    S.append(f"""## 13. Conclusions and next steps

**Hub-controlled evaluation (items 1-2)**
- Every hub-controlled test has a published counterpart: covariate-adjusted AUROC (Janes & Pepe 2009), matched
  estimation of it (caliper; Austin 2011), degree-corrected negatives (Aiyappa et al., ICML 2025). Only the 0.1 caliper
  width on log-degree is a new choice.
- BEAR-GRN ChIP: no released method exceeds {f3(rel_c[0])} caliper-matched AUROC ({rel_c[1]}). MeVD-GRN reaches
  {f3(bget('ChIP', M, 'AUROC_caliper'))} caliper-matched and {f3(bget('ChIP', M, 'AUROC_adjusted'))} covariate-adjusted,
  level with the in-degree rule ({f3(bget('ChIP', 'baseline_indegree', 'AUROC_caliper'))},
  {f3(bget('ChIP', 'baseline_indegree', 'AUROC_adjusted'))}).
- BEAR-GRN KO: the in-degree rule counted over knockout labels alone scores
  {f3(bget('KO', 'simple_deg_all_perturbation', 'AUROC_caliper'))} caliper-matched, above every method. With that
  in-degree as the covariate, MeVD-GRN falls to {f3(bget('KO', M, 'AUROC_caliper', BK))}; the degree-corrected test
  gives the same picture ({f3(bget('KO', M, 'AUROC_aiyappa'))}).

**Simple models (item 4)**
- BEAR-GRN ChIP: logistic regression on 19 cell-data features (20 parameters) reaches
  {f3(bget('ChIP', 'simple_lr_data', 'AUROC_caliper'))} caliper-matched and
  {f3(bget('ChIP', 'simple_lr_data', 'AUROC_adjusted'))} covariate-adjusted, above MeVD-GRN. Tested: retraining MeVD-GRN with negatives drawn from the same TF as each link raises ChIP
  caliper-matched AUROC from {f3(bget('ChIP', 'mevd_fm_h384_hubmotif_seeds42to44', 'AUROC_caliper'))} to {f3(bget('ChIP', 'mevd_fm_h384_hubmotif_sametf', 'AUROC_caliper'))} and covariate-adjusted from {f3(bget('ChIP', 'mevd_fm_h384_hubmotif_seeds42to44', 'AUROC_adjusted'))} to {f3(bget('ChIP', 'mevd_fm_h384_hubmotif_sametf', 'AUROC_adjusted'))} (9 datasets, seeds
  42-44), about half of the gap to the logistic regression; the rest remains. On KO it changes nothing
  ({f3(bget('KO', 'mevd_fm_h384_hubmotif_sametf', 'AUROC_caliper'))} vs
  {f3(bget('KO', 'mevd_fm_h384_hubmotif_seeds42to44', 'AUROC_caliper'))}; {f3(bget('KO', 'mevd_fm_h384_hubmotif_sametf', 'AUROC_caliper', BK))}
  with the knockout-label covariate).
- BEAR-GRN KO, knockout-label covariate: knockout-trained cell-data LDA keeps
  {f3(bget('KO', 'simple_lda_data_perturbation', 'AUROC_caliper', BK))} caliper-matched, the largest signal left
  beyond hubness on KO (Pearson {f3(bget('KO', 'baseline_pearson', 'AUROC_caliper', BK))}).
- scMultiomeGRN: the beta-model (one parameter per TF) matches MeVD-GRN on the official test
  ({o('ladder_beta', 'test_aupr')} vs {o(d14.MEVD, 'test_aupr')} AUPR); a 4-dimensional latent-space model gives the best
  caliper-matched AUROC ({o('ladder_lsm4', 'test_degcaliper_auroc')} vs {o(d14.MEVD, 'test_degcaliper_auroc')}). On pairs
  with one unseen TF every network model scores {t('ladder_lsm2', 'test_c2_auroc')}-{t('ladder_beta', 'test_c2_auroc')}
  vs MeVD-GRN {t(d14.MEVD_TF, 'test_c2_auroc')}.

**ATAC (item 3)**
- MeVD-GRN gains nothing from ATAC on either benchmark: BEAR-GRN ChIP caliper-matched
  {f3(bget('ChIP', M, 'AUROC_caliper'))} full, {f3(bget('ChIP', M + '_atacnone', 'AUROC_caliper'))} removed,
  {f3(bget('ChIP', M + '_atacshuffle', 'AUROC_caliper'))} shuffled; scMultiomeGRN AUPR {o(d14.MEVD, 'test_aupr')} /
  {o('mevd_fm_h384_obs_atacnone', 'test_aupr')} / {o('mevd_fm_h384_obs_atacshuffle', 'test_aupr')}.
- scMultiomeGRN (official code) is unchanged with its ATAC features shuffled ({o(d14.SMG, 'test_aupr')} vs
  {o('official_scmultiomegrn_atac_shuffle', 'test_aupr')} AUPR), although that benchmark's ground truth is built from
  ATAC peaks.
- LINGER re-run (Macrophage_S1 and Macrophage_S2 in slide 6; iPS, K562 and a second Macrophage_S1 seed still
  running at the time of writing), Macrophage_S1 caliper-matched AUROC: full {f3(bget('ChIP', 'linger_rerun_full', 'AUROC_caliper', MS1))},
  peaks shuffled {f3(bget('ChIP', 'linger_rerun_peakshuffle', 'AUROC_caliper', MS1))}, fully randomised
  {f3(bget('ChIP', 'linger_rerun_random', 'AUROC_caliper', MS1))}: no loss without ATAC (one dataset, one seed).

**Label-free, non-GNN models (item 5)**
- DeepSEM, RegDiffusion and the ridge Gaussian graphical model score
  {f3(bget('ChIP', 'lf_ggm_ridge', 'AUROC_caliper'))}-{f3(bget('ChIP', 'lf_deepsem', 'AUROC_caliper'))}
  caliper-matched on ChIP, the level of the released BEAR-GRN methods.

**Next steps (for discussion)**
- Close the remaining ChIP gap between MeVD-GRN with same-TF negatives and the 20-parameter logistic regression.
- Report BEAR-GRN with the covariate-adjusted and degree-corrected AUROC next to the standard metrics, and KO with the
  knockout-label covariate.
- Decide the role of ATAC in MeVD-GRN: no test shows a contribution.""")
    OUT.write_text("\n\n".join(S) + "\n")
    print(f"wrote {OUT} ({len(S)} sections)")


if __name__ == "__main__":
    main()
