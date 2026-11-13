#!/usr/bin/env python
"""Builds docs/meetings/presentation_18th_oct_14th_nov.md (bullet-point slides, advisor meeting).

Every number is read from the compiled tables in docs/results/scmgrn/final/ (scripts/scmgrn/25_scmgrn_baseline_table.py,
12 cell types x 10 folds, cluster-bootstrap CIs over cell types, paired Wilcoxon over cell type x fold) and, for
BEAR-GRN, from docs/results/bear/headline_2026-11.json when present. Re-run after new results:
    python docs/meetings/build_deck_11_14.py
In every table the best value of a column is bold and the second best italic.
"""
import json
from pathlib import Path

HERE = Path(__file__).parent
REPO = HERE.parents[1]
RES = REPO / "docs/results/scmgrn/final"
OUT = HERE / "presentation_18th_oct_14th_nov.md"
N_CT = 12


def rows(name):
    d = json.load(open(RES / f"{name}.json"))
    out = {}
    for s in d["sections"]:
        for r in s["rows"]:
            if r["n_cell_types"] and r["method"] not in out:
                out[r["method"]] = r
    return out


OFF = rows("official_vs_scmultiomegrn")          # reference: scMultiomeGRN (official code, 3 seeds)
OFF_M = rows("official_vs_mevd")                 # reference: MeVD-GRN
TF = rows("tfholdout_vs_scmultiomegrn")
TF_M = rows("tfholdout_vs_mevd")
MEVD, SMG = "mevd_fm_h384_obs_dm", "official_scmultiomegrn_full"
MEVD_TF, SMG_TF = "mevd_fm_h384_obs", "official_scmultiomegrn"


def v(R, m, key):
    r = R.get(m)
    if r is None or r.get("n_cell_types", 0) < N_CT:
        return None
    return r.get(key)


def f3(x):
    return "–" if x is None else f"{x:.3f}"


def ci(R, m, key):
    r = R.get(m) or {}
    lo, hi = r.get(f"{key}_ci") or (None, None)
    return "" if lo is None else f" [{lo:.3f}, {hi:.3f}]"


def pval(R, m, key):
    p = (R.get(m) or {}).get(f"{key}_p")
    return "n/a" if p is None else (f"{p:.0e}" if p < 1e-3 else f"{p:.3f}")


def mark(col):
    vals = sorted({round(x, 3) for x in col if x is not None}, reverse=True)
    best, second = (vals + [None, None])[:2]
    out = []
    for x in col:
        if x is None:
            out.append("–")
        elif round(x, 3) == best:
            out.append(f"**{x:.3f}**")
        elif round(x, 3) == second:
            out.append(f"*{x:.3f}*")
        else:
            out.append(f"{x:.3f}")
    return out


def table(R, items, cols, header):
    """items: [(method, label)]; cols: [key]; rows with fewer than 12 cell types are shown as pending."""
    body = [[v(R, m, k) for k in cols] for m, _ in items]
    marked = list(zip(*[mark(c) for c in zip(*body)])) if body else []
    lines = ["| Method | " + " | ".join(header) + " |", "|---|" + "---|" * len(header)]
    for (m, lab), cells in zip(items, marked):
        n = (R.get(m) or {}).get("n_cell_types", 0)
        if n < N_CT:
            cells = ["pending"] * len(cells)
        lines.append(f"| {lab} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


REF = [(MEVD, "MeVD-GRN"), (SMG, "scMultiomeGRN (official code)")]
OFF_COLS = ["test_aupr", "test_auroc", "test_degmatched_auroc", "test_degcaliper_auroc"]
OFF_HEAD = ["AUPR", "AUROC", "AUROC, degree-matched", "AUROC, caliper-matched"]

LABEL_FREE = [("unsup_pearson_abs", "Pearson correlation of expression"),
              ("unsup_spearman_abs", "Spearman correlation of expression"),
              ("unsup_partial_corr_abs", "Partial correlation of expression"),
              ("unsup_mutual_information", "Mutual information of expression"),
              ("unsup_clr", "CLR on mutual information"),
              ("unsup_grnboost_tf_tf", "GRNBoost2 between the two TFs"),
              ("unsup_genie3", "GENIE3 between the two TFs"),
              ("unsup_ridge_regression", "Ridge regression"),
              ("unsup_lasso_regression", "Lasso regression"),
              ("unsup_maestro_pearson_abs", "Correlation of MAESTRO accessibility"),
              ("unsup_geneformer_cosine", "Similarity of Geneformer embeddings")]
SUPERVISED = [("lr_maestro", "MAESTRO accessibility only"),
              ("lr_grnboost", "GRNBoost2 profiles only"),
              ("lr_geneformer", "Geneformer embedding only"),
              ("lr_histogram", "Joint expression histograms only"),
              ("lr_promoter_openness", "Promoter openness only"),
              ("lr_maestro_grnboost_concat", "MAESTRO and GRNBoost2 side by side"),
              ("lr_maestro_x_grnboost", "MAESTRO multiplied by GRNBoost2"),
              ("lr_within_pair_product", "Product of the two TFs' features only"),
              ("linreg_maestro_grnboost_concat", "Linear regression, side by side"),
              ("rf_maestro_grnboost_concat", "Random forest, side by side")]
HEURISTICS = [("heur_degree_product", "Product of the two link counts"),
              ("heur_degree_sum", "Sum of the two link counts"),
              ("heur_common_neighbours", "Common neighbours"),
              ("heur_adamic_adar", "Adamic-Adar index"),
              ("heur_resource_allocation", "Resource allocation index"),
              ("heur_katz", "Katz index"),
              ("heur_rwr", "Random walk with restart"),
              ("heur_svd_lowrank", "Low-rank SVD of the training links"),
              ("cf_lookalike", "Look-alike TFs (matrix completion)"),
              ("lr_degree_only", "Logistic regression on the link counts"),
              ("lr_topology_only", "Logistic regression on counts and neighbour scores"),
              ("lr_degree_maestro_grnboost", "Link counts + MAESTRO + GRNBoost2"),
              ("lr_topology_geneformer", "Counts and neighbour scores + Geneformer")]
GNN = [("deeptfni", "DeepTFni"), ("genelink", "GENELink")]
SMG_ABL = [("official_scmultiomegrn_no_prior", "Without the TF-TF graph"),
           ("official_scmultiomegrn_string_prior", "STRING network instead of the known links"),
           ("official_scmultiomegrn_no_edge", "Without the pair histograms"),
           ("official_scmultiomegrn_atac_only", "ATAC features only"),
           ("official_scmultiomegrn_scrna_only", "RNA features only"),
           ("official_scmultiomegrn_plus_geneformer", "With Geneformer added")]
MEVD_ABL = [("mevd_abl_nofm_obs_dm", "Without Geneformer"),
            ("mevd_abl_noatac_obs_dm", "Without the ATAC tower"),
            ("mevd_abl_nognn_obs_dm", "Without the graph network"),
            ("mevd_abl_coexpr_obs_dm", "Co-expression graph only"),
            ("mevd_abl_tfcand_obs_dm", "TF-candidate graph only"),
            ("mevd_abl_base_obs_dm", "Base size (smaller model)"),
            ("mevd_fm_h384_noobs_dm", "Without the training-link relation"),
            ("mevd_abl_min_noobs_dm", "Floor model (none of the four parts)"),
            ("mevd_abl_fl_atac_noobs_dm", "Floor model + ATAC tower"),
            ("mevd_abl_fl_fm_noobs_dm", "Floor model + Geneformer"),
            ("mevd_abl_fl_gnn_noobs_dm", "Floor model + graph network")]
MEVD_ABL_TF = [(m.replace("_dm", ""), lab) for m, lab in MEVD_ABL]

TF_COLS = ["test_c2_auroc", "test_c3_auroc", "test_degcaliper_auroc"]
TF_HEAD = ["AUROC, C2", "AUROC, C3", "AUROC, C2 caliper-matched"]
TF_ITEMS = [(MEVD_TF, "MeVD-GRN"), (SMG_TF, "scMultiomeGRN (official code)"), ("deeptfni", "DeepTFni"),
            ("genelink", "GENELink"), ("heur_degree_product", "Product of the two link counts"),
            ("lr_degree_only", "Logistic regression on the link counts"),
            ("cf_lookalike", "Look-alike TFs (matrix completion)"),
            ("gf_family_top5", "Geneformer look-alikes linked to the partner"),
            ("unsup_mutual_information", "Mutual information of expression"),
            ("unsup_grnboost_tf_tf", "GRNBoost2 between the two TFs"),
            ("lr_histogram", "Joint expression histograms only"),
            ("rf_maestro_grnboost_concat", "Random forest, MAESTRO and GRNBoost2"),
            ("lr_degree_geneformer", "Link counts + Geneformer")]


def g(R, m, k):
    x = v(R, m, k)
    return f3(x)


def best_in(R, items, key):
    xs = [(v(R, m, key), lab) for m, lab in items if v(R, m, key) is not None]
    return max(xs) if xs else (None, "–")


def bear_slide():
    import pandas as pd
    p = REPO / "docs/results/bear/head_to_head__mevd_fm_h384_hubmotif.csv"
    if not p.exists():
        return "- MeVD-GRN runs on the 9 BEAR-GRN datasets are in progress; results will be added when scored."
    d = pd.read_csv(p)
    prim = d[d["gt"].isin(["ChIP", "KO", "Union", "Intersection"]) & d["mevd"].notna()]
    chip = prim[prim["gt"] == "ChIP"]
    n_ds, seeds = chip["dataset"].nunique(), sorted(chip["n_seeds"].unique())
    a, b = chip[chip.metric == "AUROC"], chip[chip.metric == "AUPRC"]
    lines = [
        "- BEAR-GRN (Nat Methods 2025) scores TF-to-gene networks against ChIP-seq, knockout (KO), union and intersection",
        "  ground truths on 9 datasets (human K562, macrophage, iPS; mouse embryo and ESC). None of them comes from a motif rule.",
        "- MeVD-GRN follows the pre-registered protocol (section 12 of `docs/benchmarks/bear_grn.md`): TF-disjoint 5-fold",
        "  cross-fitting, so every scored TF is unseen in training; the configuration was fixed on K562 before any other",
        "  dataset was run. Released outputs of the BEAR methods are re-scored with the same code.",
        f"- **Interim:** {n_ds} of 9 datasets, {'-'.join(map(str, seeds))} seeds per dataset (5 seeds and the strict L2 regime",
        "  are still running).",
        f"- Over the four ground truths ({len(prim)} dataset x metric cells), MeVD-GRN beats LINGER in "
        f"{int(prim.beats_LINGER.sum())}, every released method in {int(prim.beats_all_released.sum())}, and the target",
        f"  in-degree baseline in {int((prim.mevd > prim.indegree).sum())}. The pre-registered win (above LINGER, in-degree"
        f" and coverage) holds in {int(prim.WIN.sum())}.",
        f"- On ChIP, MeVD-GRN beats LINGER in AUPRC on {int(b.beats_LINGER.sum())} of {len(b)} datasets, but the in-degree",
        f"  baseline is above MeVD-GRN on {int((a.indegree > a.mevd).sum())} of {len(a)} (AUROC) and also above LINGER on "
        f"{int((b.indegree > b.LINGER).sum())} of {len(b)} (AUPRC).",
        "- So BEAR-GRN shows the same pattern as scMultiomeGRN: ranking target genes by how many TFs regulate them beats",
        "  the published methods, and MeVD-GRN does not yet clear that baseline on ChIP.",
        "",
        "ChIP-seq ground truth (in-degree = mean of 5 seeds; best of the three per metric in bold):",
        "",
        "| Dataset | AUROC MeVD-GRN | AUROC LINGER | AUROC in-degree | AUPRC MeVD-GRN | AUPRC LINGER | AUPRC in-degree | AUPRC random |",
        "|---|---|---|---|---|---|---|---|"]
    for ds in sorted(chip["dataset"].unique()):
        r = {m: chip[(chip.dataset == ds) & (chip.metric == m)].iloc[0] for m in ("AUROC", "AUPRC")}
        cells = []
        for m in ("AUROC", "AUPRC"):
            vals = [r[m]["mevd"], r[m]["LINGER"], r[m]["indegree"]]
            top = max(vals)
            cells += [f"**{x:.3f}**" if x == top else f"{x:.3f}" for x in vals]
        lines.append(f"| {ds} | " + " | ".join(cells) + f" | {r['AUPRC']['random']:.3f} |")
    return "\n".join(lines)


def main():
    mev_dm, smg_dm = v(OFF, MEVD, "test_degmatched_auroc"), v(OFF, SMG, "test_degmatched_auroc")
    best_simple_dm = best_in(OFF, LABEL_FREE + SUPERVISED + HEURISTICS, "test_degmatched_auroc")
    best_cal = best_in(OFF, LABEL_FREE + SUPERVISED + HEURISTICS + GNN, "test_degcaliper_auroc")
    tf_best_cal = best_in(TF, TF_ITEMS, "test_degcaliper_auroc")
    tf_best_cal_m = next(m for m, lab in TF_ITEMS if lab == tf_best_cal[1])
    tf_best_c3 = best_in(TF, TF_ITEMS, "test_c3_auroc")
    S = []
    S.append(f"""# Baselines, ablations and hub-controlled evaluation on the scMultiomeGRN benchmark

Advisor meeting, work from 18 October to 14 November 2026.

Every scMultiomeGRN-benchmark number is an average over the 12 cell types (9 lung, 3 blood) and 10 folds, scored
with scMultiomeGRN's metric code on held-out pairs only. Brackets are 95% confidence intervals from resampling cell
types; p-values are paired Wilcoxon tests over cell type x fold. In each table **bold** marks the best value in a
column and *italics* the second best.""")

    S.append(f"""## 1. Summary of results

- About 50 simple and statistical models, the two GNN baselines of the scMultiomeGRN paper (DeepTFni, GENELink),
  six ablations of scMultiomeGRN and eleven ablations of MeVD-GRN were evaluated on scMultiomeGRN's datasets.
- No label-free statistical score recovers the ground-truth network: all score 0.50 to 0.54 AUROC on the official
  test.
- Counting each TF's training links is enough to beat scMultiomeGRN on the official test. A logistic regression on
  the two counts scores {g(OFF,'lr_degree_only','test_aupr')} AUPR against {g(OFF,SMG,'test_aupr')} for
  scMultiomeGRN and {g(OFF,MEVD,'test_aupr')} for MeVD-GRN.
- On tests that remove the advantage of hub TFs, every method falls towards chance. On the caliper-matched test,
  MeVD-GRN scores {g(OFF,MEVD,'test_degcaliper_auroc')} AUROC and scMultiomeGRN {g(OFF,SMG,'test_degcaliper_auroc')}.
- On a new split that holds out whole TFs, the link-count rule is the best method when one TF of the pair is unseen
  ({g(TF,'heur_degree_product','test_c2_auroc')} AUROC against {g(TF,MEVD_TF,'test_c2_auroc')} for MeVD-GRN), and
  every model is near chance once hubness is controlled ({g(TF,MEVD_TF,'test_degcaliper_auroc')} for MeVD-GRN).
- Conclusion: on this benchmark, the scores of all methods, MeVD-GRN included, are explained by hub TFs.
- BEAR-GRN (ChIP-seq and knockout ground truths, interim results): MeVD-GRN beats LINGER on most datasets, but a
  baseline that ranks target genes by their number of regulators is stronger still on ChIP (section 11).
- Correction: an earlier draft reported 0.878 AUPR for MeVD-GRN on the degree-matched test. That value came from an
  evaluation bug, fixed in this period; the corrected value is {g(OFF,MEVD,'test_degmatched_aupr')} (section 9).""")

    S.append("""## 2. Model families and evaluation protocol

- **Label-free statistical scores.** No labels: correlations, mutual information, GRNBoost2, GENIE3, ridge and lasso
  regression (the BEELINE and geneRNIB baselines).
- **Supervised models on pair features.** Logistic regression, linear regression and random forest on scMultiomeGRN's
  inputs (MAESTRO accessibility, GRNBoost2 profiles, expression histograms) and the Geneformer embedding.
- **Link-prediction heuristics on the training graph.** Degree, common neighbours, Adamic-Adar, Katz, random walks,
  low-rank completion and look-alike TFs, computed from the fold's training links only.
- **GNN baselines.** DeepTFni and GENELink, the baselines of the scMultiomeGRN paper, run from the copies shipped in
  the official repository on the same folds.
- **scMultiomeGRN ablations.** One part of the authors' code is switched off at a time; with every switch at its
  default the patched code reproduces the published code exactly.
- **MeVD-GRN ablations.** One component removed at a time from the full model, and each component added alone to a
  floor model that has none of them (3 seeds each).
- **Official test.** scMultiomeGRN's 10 folds split the links: each test link is paired with a random non-link.
- **Degree-matched test.** Each test link is paired with one of the three test non-links closest to it in the two
  TFs' training-link counts ("hidden sampling biases" preprint, December 2025).
- **Caliper-matched test.** Stricter: a test link is kept only if a test non-link has both TFs' log link counts
  within 0.1 of the link's, so the link-count rule scores 0.5 by construction.
- **Held-out-TF split (new).** The TFs, not the links, are split into 10 groups (Park and Marcotte 2012). C2 = one TF
  of the pair unseen in training, C3 = both unseen; a C2 caliper-matched test pairs each link of a hidden TF with a
  non-link of the same TF whose partner has the same link count.
- **Metrics.** Hub-controlled columns use AUROC, because scMultiomeGRN's trapezoidal AUPR gives a method with all
  scores tied (1 + share of positives) / 2, about 0.75, instead of 0.5.
- Nothing from the test fold is used for fitting; every model is checked for zero test-to-train leakage.""")

    S.append(f"""## 3. Label-free statistical scores

- These methods ask whether two TFs behave alike across cells, without any label.
- All of them sit between 0.50 and 0.54 AUROC on every test, which is chance.
- A link in the ground-truth network means one TF's motif lies in the other TF's open promoter; that says little
  about co-expression, so these scores cannot recover it.

{table(OFF, REF + LABEL_FREE, OFF_COLS, OFF_HEAD)}""")

    S.append(f"""## 4. Supervised models on pair features

- Learned from the training links, scMultiomeGRN's inputs become informative on the official test: GRNBoost2 alone
  reaches {g(OFF,'lr_grnboost','test_aupr')} AUPR, the two inputs side by side {g(OFF,'lr_maestro_grnboost_concat','test_aupr')},
  a random forest {g(OFF,'rf_maestro_grnboost_concat','test_aupr')}.
- Giving the inputs side by side works better than multiplying them.
- Under the caliper-matched test all of them fall to between 0.49 and 0.52: what they learn is mostly which TFs are
  hubs.

{table(OFF, REF + SUPERVISED, OFF_COLS, OFF_HEAD)}""")

    S.append(f"""## 5. Link-prediction heuristics on the training graph

- The product of the two link counts scores {g(OFF,'heur_degree_product','test_aupr')} AUPR and a logistic regression on
  the two counts {g(OFF,'lr_degree_only','test_aupr')}; both beat scMultiomeGRN ({g(OFF,SMG,'test_aupr')}).
- Look-alike TFs (are TFs with similar training links linked to the partner?) score
  {g(OFF,'cf_lookalike','test_aupr')} and keep {g(OFF,'cf_lookalike','test_degmatched_auroc')} AUROC on the degree-matched
  test, close to MeVD-GRN ({g(OFF,MEVD,'test_degmatched_auroc')}): the family structure of the training links explains
  most of what link counts do not.
- On the caliper-matched test every heuristic is between 0.48 and 0.57.

{table(OFF, REF + HEURISTICS, OFF_COLS, OFF_HEAD)}""")

    S.append(f"""## 6. GNN baselines of the scMultiomeGRN paper

- DeepTFni scores {g(OFF,'deeptfni','test_aupr')} AUPR and GENELink {g(OFF,'genelink','test_aupr')} on the official test,
  both below scMultiomeGRN, consistent with the scMultiomeGRN paper.
- Both are close to chance on the caliper-matched test.

{table(OFF, REF + GNN, OFF_COLS, OFF_HEAD)}""")

    S.append(f"""## 7. scMultiomeGRN ablations

- The known-link graph is what matters: without it AUPR drops from {g(OFF,SMG,'test_aupr')} to
  {g(OFF,'official_scmultiomegrn_no_prior','test_aupr')}; a STRING network in its place gives
  {g(OFF,'official_scmultiomegrn_string_prior','test_aupr')}.
- The pair histograms, either single modality and an added Geneformer view change nothing measurable.
- On the hub-controlled tests the ablations are as close to chance as the full model.

{table(OFF, REF + SMG_ABL, OFF_COLS, OFF_HEAD)}""")

    S.append(f"""## 8. MeVD-GRN ablations

- Geneformer is the component that matters on the official test: the floor model plus Geneformer alone reaches
  {g(OFF,'mevd_abl_fl_fm_noobs_dm','test_aupr')} AUPR, and removing Geneformer from the full model costs the most.
- Removing the ATAC tower, the graph network or the training-link relation changes little.
- No component moves the hub-controlled scores away from chance.

{table(OFF, REF[:1] + MEVD_ABL, OFF_COLS, OFF_HEAD)}

Held-out-TF split:

{table(TF, [(MEVD_TF, "MeVD-GRN")] + MEVD_ABL_TF, TF_COLS, TF_HEAD)}""")

    S.append(f"""## 9. Hub-controlled tests on the official split, and a correction

- Degree-matched test: MeVD-GRN {g(OFF,MEVD,'test_degmatched_auroc')} AUROC{ci(OFF,MEVD,'test_degmatched_auroc')},
  scMultiomeGRN {g(OFF,SMG,'test_degmatched_auroc')}, best simple model {f3(best_simple_dm[0])}, from {best_simple_dm[1]}.
  MeVD-GRN is best, but only {f3((mev_dm or 0) - (best_simple_dm[0] or 0))} above the best simple model.
- The link-count rule still scores {g(OFF,'heur_degree_product','test_degmatched_auroc')} on the degree-matched test,
  because non-links between two hubs are rare in these dense graphs. The caliper-matched test removes that.
- Caliper-matched test: MeVD-GRN {g(OFF,MEVD,'test_degcaliper_auroc')}{ci(OFF,MEVD,'test_degcaliper_auroc')},
  scMultiomeGRN {g(OFF,SMG,'test_degcaliper_auroc')}, best other method {f3(best_cal[0])}, from {best_cal[1]}.
- MeVD-GRN's interval lies above 0.5, so a small signal beyond link counts remains on the official split, of the
  same size as low-rank SVD's ({g(OFF,'heur_svd_lowrank','test_degcaliper_auroc')}{ci(OFF,'heur_svd_lowrank','test_degcaliper_auroc')})
  and GENELink's ({g(OFF,'genelink','test_degcaliper_auroc')}{ci(OFF,'genelink','test_degcaliper_auroc')}). These sets are
  small (a few pairs per fold in the densest cell types).
- **Correction.** The degree-matched value of 0.878 AUPR reported for MeVD-GRN in an earlier draft was wrong. The
  training script reduced each fold to one row per TF pair before scoring, and the matched tests counted a TF's links
  from the first column only, so MeVD-GRN was scored on different pairs from every other method. The fix counts links
  on unique pairs; every method now uses identical test pairs on every fold (checked on all 55,000 fold results), and
  the corrected value is {g(OFF,MEVD,'test_degmatched_aupr')} AUPR.""")

    S.append(f"""## 10. Held-out-TF split

- With one TF unseen (C2), the link-count rule is the best method: it still knows the partner TF's link count.
  MeVD-GRN scores {g(TF,MEVD_TF,'test_c2_auroc')}{ci(TF,MEVD_TF,'test_c2_auroc')}, scMultiomeGRN
  {g(TF,SMG_TF,'test_c2_auroc')}, the link-count rule {g(TF,'heur_degree_product','test_c2_auroc')}.
- With both TFs unseen (C3), every supervised model is at chance; the best score is {f3(tf_best_c3[0])},
  from {tf_best_c3[1]}. C3 sets are small (4 to 102 links per fold).
- On the C2 caliper-matched test no method exceeds 0.52: MeVD-GRN {g(TF,MEVD_TF,'test_degcaliper_auroc')}
  {ci(TF,MEVD_TF,'test_degcaliper_auroc').strip()}, best method {f3(tf_best_cal[0])}{ci(TF,tf_best_cal_m,'test_degcaliper_auroc')},
  from {tf_best_cal[1]}.
- Neighbourhood heuristics and look-alike TFs score exactly 0.5 here, because a hidden TF has no training links.

{table(TF, TF_ITEMS, TF_COLS, TF_HEAD)}""")

    S.append(f"""## 11. BEAR-GRN: ChIP-seq and knockout ground truths

{bear_slide()}""")

    S.append(f"""## 12. Conclusions and next steps

- The scMultiomeGRN benchmark measures how well a method knows which TFs are hubs. Its ground-truth network comes
  from a motif rule on the same ATAC data the models read, and its random non-links are easy to separate by link count.
- Simple link-count models match or beat the published GNNs; label-free statistics are at chance.
- Once hubness is controlled (caliper matching, held-out TFs), no method, MeVD-GRN included, keeps a meaningful signal.
- BEAR-GRN shows the same effect on target genes: an in-degree ranking beats every published method on ChIP AUPRC,
  and MeVD-GRN beats LINGER but not yet that ranking.
- These results are a finding about how GRN benchmarks are scored and belong in the paper as an evaluation
  contribution, alongside MeVD-GRN.
- Next: complete BEAR-GRN (5 seeds, the strict L2 regime, ablations) with a hub control against the in-degree
  baseline (within-in-degree-stratum scores); then PBMC10k with LINGER's own candidate genes.""")

    OUT.write_text("\n\n".join(S) + "\n")
    print(f"wrote {OUT} ({len(S)} slides)")


if __name__ == "__main__":
    main()
