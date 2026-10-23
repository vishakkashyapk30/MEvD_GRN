#!/usr/bin/env python
"""Builds docs/presentation/presentation_18th_oct_24th_oct.md (bullet-point slides).

Numbers come from docs/experiments/scmgrn_results/baselines/vs_smg_lung_hema.json (scripts/25_scmgrn_baseline_table.py,
12 cell types, held-out test fold). In every table the best value of a column is bold and the second best italic.
"""
import json
from pathlib import Path

HERE = Path(__file__).parent
OUT = HERE / "presentation_18th_oct_24th_oct.md"
RES = json.load(open(HERE.parent / "experiments/scmgrn_results/baselines/vs_smg_lung_hema.json"))
ROW = {r["method"]: r for s in RES["sections"] for r in s["rows"]}


def mark(vals, nd=3):
    shown = [None if v is None else f"{v:.{nd}f}" for v in vals]
    uniq = sorted({float(s) for s in shown if s is not None}, reverse=True)
    best = uniq[0] if uniq else None
    second = uniq[1] if len(uniq) > 1 else None
    out = []
    for s in shown:
        if s is None:
            out.append("–")
        elif float(s) == best:
            out.append(f"**{s}**")
        elif float(s) == second:
            out.append(f"*{s}*")
        else:
            out.append(s)
    return out


COLS = [("test_aupr", "AUPR"), ("test_auroc", "AUROC"),
        ("test_degmatched_aupr", "AUPR, degree-matched"), ("test_degmatched_auroc", "AUROC, degree-matched")]


def table(items, cols=COLS):
    """items: list of (method key, label). Columns are marked down the column."""
    vals = {k: mark([ROW[m].get(k) for m, _ in items]) for k, _ in cols}
    lines = ["| Method | " + " | ".join(h for _, h in cols) + " |", "|---|" + "---|" * len(cols)]
    for i, (m, label) in enumerate(items):
        lines.append(f"| {label} | " + " | ".join(vals[k][i] for k, _ in cols) + " |")
    return "\n".join(lines)


def v(m, k="test_aupr"):
    return f"{ROW[m][k]:.3f}"


REF = [("mevd_fm_h384_obs_dm", "MeVD-GRN"), ("official_scmultiomegrn_full", "scMultiomeGRN (unchanged)")]

S = []
add = S.append
add("""# Simple baselines and scMultiomeGRN ablations on the scMultiomeGRN benchmark

Advisor meeting, week of 18 to 24 October 2026.

Every number is an average over the 12 cell types (9 lung, 3 blood) and the 10 official folds, scored on the held-out test links with scMultiomeGRN's own metric code. In each table, **bold** marks the best value in a column and *italics* the second best.""")

add(f"""## 1. This week in one slide

- We tried 38 simple and statistical models on scMultiomeGRN's datasets and splits, and ran six ablations of scMultiomeGRN itself with the authors' code.
- None of the classic network-inference scores from the cell data finds the answer key. Correlation, mutual information, GRNBoost2, GENIE3 and regression all score between 0.50 and 0.53 AUPR, which is chance on this test.
- Counting each TF's known training links is enough to beat scMultiomeGRN. A logistic regression on the two counts scores {v('lr_degree_only')} against {v('official_scmultiomegrn_full')}, and MeVD-GRN scores {v('mevd_fm_h384_obs_dm')}.
- scMultiomeGRN depends on the known links it passes messages over. Without them it drops to {v('official_scmultiomegrn_no_prior')}, and with a STRING network in their place to {v('official_scmultiomegrn_string_prior')}. Adding Geneformer to it changes nothing.
- We also scored every model on a harder test where each true link is paired with a non-link whose TFs have similar numbers of links.
- On that test, scMultiomeGRN and every simple model that counts links fall to between 0.70 and 0.76 AUPR, close to the counting rule at {v('heur_degree_product', 'test_degmatched_aupr')}. MeVD-GRN keeps {v('mevd_fm_h384_obs_dm', 'test_degmatched_aupr')}.""")

add("""## 2. What we tried and how we scored it

- **Statistical scores from the cell data.** These use no labels at all. They are the methods from the BEELINE and geneRNIB benchmarks: correlations, mutual information, GRNBoost2, GENIE3, and ridge and lasso regression.
- **Simple supervised models on features.** These are logistic regressions on scMultiomeGRN's own inputs, MAESTRO accessibility and GRNBoost2 profiles, plus the Geneformer embedding. We give the inputs one at a time, side by side, and multiplied together.
- **Counting rules on the known links.** These are the classic link-prediction heuristics: degree, common neighbours, Adamic-Adar, Katz and random walks. All of them use only the fold's training links.
- **scMultiomeGRN ablations.** We changed one part of the authors' code at a time and checked that with no change it gives exactly their numbers.
- **A degree-matched test.** This follows a December 2025 preprint, which showed that hub TFs inflate scores in this kind of benchmark. Each test link is paired with the test non-link closest to it in how many training links its two TFs have, so knowing the hubs no longer helps much.
- Nothing from the test fold is used for fitting, and the regularisation of every logistic regression is chosen on the training pairs only.""")

add("""## 3. Statistical scores from the cell data

- None of these methods sees a label. They ask whether two TFs behave alike across cells.
- All of them land near 0.5, which is what a random guess gets on a test with as many links as non-links.
- The reason is the answer key itself. A link means that one TF's motif sits in the other TF's open promoter, and that says little about whether the two TFs are co-expressed.
- GRNBoost2 and GENIE3, the standard tools for gene networks, do no better than plain correlation here.
- This matches the scMultiomeGRN paper, where GENIE3 also came last among its baselines.

""" + table(REF + [
    ("unsup_pearson_abs", "Pearson correlation of expression"),
    ("unsup_spearman_abs", "Spearman correlation of expression"),
    ("unsup_partial_corr_abs", "Partial correlation of expression"),
    ("unsup_mutual_information", "Mutual information of expression"),
    ("unsup_clr", "CLR on mutual information"),
    ("unsup_grnboost_tf_tf", "GRNBoost2 between the two TFs"),
    ("unsup_genie3", "GENIE3 between the two TFs"),
    ("unsup_ridge_regression", "Ridge regression"),
    ("unsup_lasso_regression", "Lasso regression"),
    ("unsup_maestro_pearson_abs", "Correlation of MAESTRO accessibility"),
    ("unsup_geneformer_cosine", "Similarity of Geneformer embeddings"),
]))

add(f"""## 4. Simple models on scMultiomeGRN's inputs

- Once a logistic regression can learn from the training links, the same inputs become useful. MAESTRO alone reaches {v('lr_maestro')} AUPR and GRNBoost2 alone {v('lr_grnboost')}.
- Giving the two inputs side by side works better than multiplying them. Side by side reaches {v('lr_maestro_grnboost_concat')}, the per-TF product {v('lr_maestro_x_grnboost')}, and the product of the two TFs' features only {v('lr_within_pair_product')}.
- Geneformer on its own reaches {v('lr_geneformer')}, below the cell data, even though it knows nothing about these cells.
- A plain linear regression on the 0/1 label gives the same as the logistic regression. A random forest is the best of this group at {v('rf_maestro_grnboost_concat')}.
- All of these models stay well below scMultiomeGRN, so the cell data alone does not explain its score.

""" + table(REF + [
    ("lr_maestro", "MAESTRO accessibility only"),
    ("lr_grnboost", "GRNBoost2 profiles only"),
    ("lr_geneformer", "Geneformer embedding only"),
    ("lr_histogram", "Joint expression histograms only"),
    ("lr_promoter_openness", "Promoter openness only"),
    ("lr_maestro_grnboost_concat", "MAESTRO and GRNBoost2 side by side"),
    ("lr_maestro_x_grnboost", "MAESTRO multiplied by GRNBoost2"),
    ("lr_within_pair_product", "Product of the two TFs' features only"),
    ("linreg_maestro_grnboost_concat", "Linear regression, side by side"),
    ("rf_maestro_grnboost_concat", "Random forest, side by side"),
]))

add(f"""## 5. Counting the known links

- These rules use only the links in the fold's training part, never the cell data.
- The simplest one wins. The product of the two TFs' link counts scores {v('heur_degree_product')}, and a logistic regression on the two counts {v('lr_degree_only')}. Both beat scMultiomeGRN.
- Rules based on shared neighbours do worse. Common neighbours scores {v('heur_common_neighbours')}, and the Jaccard index ranks pairs almost backwards.
- Adding the cell data or Geneformer to the counts does not help. Every combination stays between {v('lr_degree_maestro_grnboost')} and {v('lr_topology_geneformer')}, no better than the counts alone.
- This is the pattern the recent literature warns about: on this benchmark, knowing which TFs are hubs is worth more than any feature.

""" + table(REF + [
    ("heur_degree_product", "Product of the two link counts"),
    ("heur_degree_sum", "Sum of the two link counts"),
    ("heur_common_neighbours", "Common neighbours"),
    ("heur_adamic_adar", "Adamic-Adar index"),
    ("heur_resource_allocation", "Resource allocation index"),
    ("heur_katz", "Katz index"),
    ("heur_rwr", "Random walk with restart"),
    ("heur_jaccard", "Jaccard index"),
    ("lr_degree_only", "Logistic regression on the link counts"),
    ("lr_topology_only", "Logistic regression on counts and neighbour scores"),
    ("lr_degree_maestro_grnboost", "Link counts + MAESTRO + GRNBoost2"),
    ("lr_degree_geneformer", "Link counts + Geneformer"),
    ("lr_topology_geneformer", "Counts and neighbour scores + Geneformer"),
]))

add(f"""## 6. Taking scMultiomeGRN apart

- We patched the authors' code so that each part can be switched off. With every switch at its default, the patched code gives exactly the published code's numbers.
- The known-link graph is what matters. Passing messages over no links at all drops AUPR from {v('official_scmultiomegrn_full')} to {v('official_scmultiomegrn_no_prior')}. The drop is small on lung but large on blood, where it falls to 0.585.
- A STRING protein network in place of the known links is worse still, at {v('official_scmultiomegrn_string_prior')}. STRING contains only 2 to 6% of the answer-key links.
- The pair histograms, the ATAC view and the RNA view can each be removed without a measurable change.
- Adding the Geneformer embedding as a third view also changes nothing, at {v('official_scmultiomegrn_plus_geneformer')}. On blood almost all of scMultiomeGRN's score comes from the known links, and none of its extra inputs adds anything measurable.

""" + table(REF + [
    ("official_scmultiomegrn_no_prior", "Without the TF-TF graph"),
    ("official_scmultiomegrn_string_prior", "STRING network instead of the known links"),
    ("official_scmultiomegrn_no_edge", "Without the pair histograms"),
    ("official_scmultiomegrn_atac_only", "ATAC features only"),
    ("official_scmultiomegrn_scrna_only", "RNA features only"),
    ("official_scmultiomegrn_plus_geneformer", "With Geneformer added"),
]))

add(f"""## 7. The degree-matched test

- On the official test, a non-link is any random pair, so it usually involves TFs with few links. A model can do well just by knowing the hubs.
- On the degree-matched test, each test link is paired with one of the three test non-links closest to it in the two TFs' training-link counts. The pairs come from the official test set, so nothing new is revealed to any model.
- The counting rule falls from {v('heur_degree_product')} to {v('heur_degree_product', 'test_degmatched_aupr')}. In these dense graphs non-links between two hubs are rare, so it cannot reach exactly 0.5.
- scMultiomeGRN falls to {v('official_scmultiomegrn_full', 'test_degmatched_aupr')} and the best simple model to {v('lr_topology_only', 'test_degmatched_aupr')}. Both end up close to the counting rule.
- MeVD-GRN keeps {v('mevd_fm_h384_obs_dm', 'test_degmatched_aupr')}. This is the clearest sign so far that it uses something beyond hubness.
- The test can only reuse the official test non-links, so hard non-links are scarce, and the next step is a split with more hard non-links, such as holding out whole TFs.

""" + table(REF + [
    ("heur_degree_product", "Product of the two link counts"),
    ("lr_degree_only", "Logistic regression on the link counts"),
    ("lr_topology_only", "Logistic regression on counts and neighbour scores"),
    ("heur_rwr", "Random walk with restart"),
    ("rf_maestro_grnboost_concat", "Random forest on MAESTRO and GRNBoost2"),
    ("official_scmultiomegrn_no_prior", "scMultiomeGRN without the TF-TF graph"),
]))

add(f"""## 8. What we take from this and what comes next

- On the official test, this benchmark mostly rewards knowing the hub TFs. A two-number logistic regression beats the published model, as recent papers found in other GRN benchmarks.
- Label-free network-inference methods cannot recover a motif-based answer key, so they are not a useful baseline for this task.
- scMultiomeGRN's features add almost nothing on top of its known-link graph, and Geneformer does not help it.
- MeVD-GRN is the only model that keeps most of its score once hubs are controlled for, at {v('mevd_fm_h384_obs_dm', 'test_degmatched_aupr')} against {v('official_scmultiomegrn_full', 'test_degmatched_aupr')} for scMultiomeGRN.
- Next, we want to confirm that gap on a split that holds out whole TFs, where no model can count a TF's links.
- After that, we want to check which part of MeVD-GRN carries the gap, using the ablations we already have under the degree-matched test.""")

OUT.write_text("\n\n".join(S) + "\n")
print("written", OUT)
