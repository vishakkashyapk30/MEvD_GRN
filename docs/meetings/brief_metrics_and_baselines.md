# Brief: the math behind the metrics, tests and baselines

Companion to `presentation_18th_oct_14th_nov.md`. Every definition below is what the code actually computes
(file names in brackets). Toy numbers were computed with the same code paths.

Notation: a **link** is a pair in the ground-truth network (label 1), a **non-link** a pair that is not (label 0).
A model gives every pair a **score**; higher = "more likely a link". On scMultiomeGRN the pairs are TF-TF and
undirected; on BEAR-GRN they are TF -> target gene.

---

## 1. Metrics

### 1.1 AUROC (area under the ROC curve)
**Meaning:** pick one link and one non-link at random; AUROC = probability the link gets the higher score
(a tie counts 1/2). 0.5 = random, 1.0 = perfect.

**Toy:** 4 links scored 0.9, 0.8, 0.4, 0.3; 4 non-links scored 0.7, 0.35, 0.2, 0.1.
There are 4 x 4 = 16 (link, non-link) pairs; the link wins 13 of them -> AUROC = 13/16 = **0.8125**.

AUROC does not depend on the share of positives, so it is comparable across test sets.

### 1.2 AUPR (area under the precision-recall curve), as scMultiomeGRN computes it
Walk down the ranking from the top score. At each cut-off:
precision = links above the cut / pairs above the cut; recall = links above the cut / all links.
AUPR = area under precision (y) vs recall (x).

scMultiomeGRN's `metric_fn` computes the area with the **trapezoidal rule** (`sklearn.metrics.auc(recall,
precision)`), i.e. it joins the points with straight lines (`src/benchmarks/scmgrn_protocol.py: metric_fn`).

**The tie problem (why hub-controlled columns use AUROC):** if a method gives every pair the same score, the curve
has only two points: (recall 0, precision 1) and (recall 1, precision p), p = share of links. The trapezoid between
them has area **(1 + p) / 2**. With p = 0.5 (the official test is balanced) that is **0.75 for a method that knows
nothing**, while its AUROC is 0.5. Linear interpolation in precision-recall space overestimates the area; this is
a known result (Davis and Goadrich, ICML 2006).

**Toy:** the 8 scores of 1.1 give trapezoidal AUPR 0.835. Setting all 8 scores to 0.5 gives AUROC 0.500 but
trapezoidal AUPR **0.750**.

### 1.3 AP (average precision)
AP = sum over each link of (precision at that link's rank) / number of links, i.e. step-wise area without
straight-line interpolation (`sklearn.metrics.average_precision_score`). All-tied scores give AP = p = 0.5.
Toy of 1.1: AP = 0.854. The tables report it next to AUPR.

### 1.4 Random level
AUROC: 0.5 always. AUPR / AP: the share of positives in the test set (0.5 on balanced scMultiomeGRN tests, but
note 1.2: trapezoidal AUPR of a constant score is 0.75). Tables print it as "Random AUPR".

---

## 2. Test sets on the scMultiomeGRN benchmark

### 2.1 Official test (scMultiomeGRN paper; re-implemented byte-for-byte)
The links of a cell type are split into 10 folds (KFold, seed 666). Fold k: test = fold k's links **plus the same
number of random non-links**; the other links are training (80%) and validation (20%).
Problem: a random non-link usually involves low-degree TFs ("degree" = number of training links of a TF), while a
link usually involves a hub. So "is either TF a hub?" separates them without any biology.

### 2.2 Degree-matched test (own variant; idea from Yu et al. 2010 and the Stock et al. 2025 preprint)
Same test links. For each test link (u, v), compute the two TFs' log training degrees, sorted:
lo = log(1 + smaller degree), hi = log(1 + larger degree). Among non-links not used in training or validation,
take the **3 nearest** by distance |lo' - lo| + |hi' - hi| and pick one at random as its partner.
(`degree_matched_test`, deterministic per fold so every method sees the same pairs.)

**Toy:** link between a TF with degree 10 and one with degree 2: lo = log(3) = 1.099, hi = log(11) = 2.398.
Candidate non-links: degrees (2, 11): distance 0 + 0.087 = 0.087; (1, 10): 0.406 + 0 = 0.406;
(3, 12): 0.287 + 0.167 = 0.454; (5, 6): far. The 3 nearest are the first three; one is drawn at random.
Limitation: in dense graphs there may be no close non-link, so the 3 "nearest" can still differ in degree. That is
why the degree product still scores 0.653 AUROC here instead of 0.5.

### 2.3 Caliper-matched test (own construction; caliper idea from propensity-score matching, Austin 2011)
Stricter version of 2.2: a test link is kept **only if** some non-link has **both** |lo' - lo| <= 0.1 and
|hi' - hi| <= 0.1; one such non-link is drawn at random. Links without a match are **dropped**.
Every kept link has a non-link partner with (almost) the same degrees, so any rule that only looks at degrees scores
~0.5 AUROC by construction. (`degree_caliper_test`.)

**Toy:** in 2.2, only the (2, 11) candidate is inside the caliper (differences 0 and 0.087), so it is the partner.
If it did not exist, the link would be dropped.

### 2.4 Held-out-TF split, C2 and C3 (Park and Marcotte, Nature Methods 2012)
TFs, not links, are split into 10 groups (seed 666). Fold k: group k = test TFs, group k+1 = validation TFs, the rest
training TFs. Training uses only pairs of two training TFs.
- **C2** = pairs with one test TF and one training TF (the test TF has **no** training links).
- **C3** = pairs of two test TFs.
Each class gets as many random non-links as links. Pairs of a test TF with a validation TF are used nowhere.
(`make_tf_holdout_splits`.) A test TF's degree is 0 in training, so a model cannot count its links.

### 2.5 C2 caliper-matched test (own construction)
For a C2 link (hidden TF t, training TF r), the partner non-link is (t, w) with the **same hidden TF t** and a
training TF w whose log degree is within 0.1 of r's. Links without such a w are dropped.
(`tfholdout_degree_matched_test`.) Both ends are matched: t is the same TF, r and w have the same degree.

**Toy:** link t-r with deg(r) = 20 (log 21 = 3.045). Non-links of t: w1 with degree 21 (3.091, inside 0.1) and w2
with degree 5 (1.792, outside). Partner = w1.

---

## 3. Baselines on the scMultiomeGRN benchmark

All use only the fold's **training** links (or no links at all); every score matrix is made symmetric.
(`scripts/scmgrn/22_scmgrn_simple_baselines.py` and `26_scmgrn_cf_probe.py`, branch `baselines`.)

**Toy training graph** for 3.1: TFs A, B, C, D, E with links A-B, A-C, A-D, B-C, D-E.
Degrees: A 3, B 2, C 2, D 2, E 1. Score the non-linked pair **B-D**.

### 3.1 Link-prediction heuristics (no cell data)
| Name | Formula for pair (u, v) | Toy B-D |
|---|---|---|
| Degree product (preferential attachment) | (d_u + 1)(d_v + 1) | 3 x 3 = 9 |
| Degree sum | d_u + d_v | 4 |
| Common neighbours | number of TFs linked to both | 1 (A) |
| Jaccard | common / (number linked to either) | 1 / 3 (A out of {A, C, E}) |
| Adamic-Adar | sum over common neighbours w of 1 / log(d_w) | 1 / log 3 = 0.910 |
| Resource allocation | sum over common neighbours w of 1 / d_w | 1 / 3 |
| Katz | sum over L = 1..4 of 0.05^L x (walks of length L from u to v) | walks: L2 = 1, L3 = 1, L4 = 5 -> 0.0027 |
| Random walk with restart | R = 0.15 (I - 0.85 P)^-1, P = row-normalised adjacency; score R_uv + R_vu | (a random walker from u that restarts at u w.p. 0.15: how often it is at v) |
| Low-rank SVD | rank-16 reconstruction of the training adjacency (U_k S_k V_k^T)_uv | "matrix completion" |
| Low-rank spectral | rank-16 reconstruction of D^-1/2 A D^-1/2 | degree-normalised version |

Adamic-Adar and resource allocation down-weight common neighbours that are hubs.

### 3.2 Look-alike TFs (collaborative filtering)
C_uw = cosine similarity of TF u's and TF w's **training-link rows** (diagonal 0).
score(u, v) = sum_w C_uw x A_wv + the same with u and v swapped: "are TFs that look like u linked to v?"
**Toy:** B's row is (A, C) and C's row is (A, B): similarity 0.5 with C and D, 0.408 with A.
For B-D: TFs linked to D are A and E, so score(B -> D) = C_BA = 0.408; from D's side 0.5; total 0.908.
Why it matters: the ground truth is a motif scan, so TFs of one family (similar motifs) have near-identical rows.
"Normalised" divides by sum_w C_uw. "Geneformer look-alikes": the same with the 5 most similar TFs by Geneformer
embedding instead of link rows.

### 3.3 Label-free statistical scores (no labels; BEELINE / geneRNIB family)
From the cells' expression (or MAESTRO accessibility) of the two TFs:
- |Pearson| and |Spearman| correlation across cells.
- |Partial correlation|: P = inverse covariance (Ledoit-Wolf shrinkage) of all TFs; rho_uv = -P_uv / sqrt(P_uu P_vv).
  It is the correlation left after removing all other TFs.
- Mutual information: each TF's expression cut into 8 equal-width bins; MI = sum p(a, b) log[p(a, b) / (p(a) p(b))].
- CLR: z-score each MI value against row u and row v (negative z set to 0), score = sqrt(z_u^2 + z_v^2).
- GRNBoost2 between the two TFs: max(importance u -> v, v -> u) from GRNBoost2's gene-boosting regression.
- GENIE3: random forest predicting each TF from all other TFs; importance, symmetrised.
- Ridge (alpha 1) and lasso (alpha 0.01) regression of each TF on all others: |coefficient|, symmetrised.
- Cosine similarity of GRNBoost2 profiles / Geneformer embeddings.

### 3.4 Supervised models on pair features
Per-TF vectors (MAESTRO, GRNBoost2 profile, Geneformer) are first reduced to 32 numbers by SVD (no labels).
**Pair feature** = [x_u * x_v (element-wise), |x_u - x_v|]: symmetric, so (u, v) and (v, u) get the same score.
- "X only": logistic regression on the pair features of one input.
- "Side by side" (given independently): concatenate the MAESTRO and GRNBoost2 pair features.
- "MAESTRO multiplied by GRNBoost2": per TF, m_u * g_u first, then pair features.
- "Product of the two TFs' features only": only the x_u * x_v parts.
- Logistic regression: L2 strength chosen from {0.001, 0.01, 0.1, 1, 10} by 3-fold CV on training pairs only.
  Linear regression = ridge on the 0/1 label. Random forest = 300 trees.
**Toy:** x_A = (1, 0), x_B = (1, 1): features = [(1, 0), (0, 1)] -> logistic regression gives a weight per feature.

### 3.5 Link counts + features
Logistic regression on [log(1 + d_u) + log(1 + d_v), |difference|, product] (plus the 3.1 scores for "neighbour
scores"), optionally with the 3.4 features. Leak guard: a training link's own features are computed with that link
removed from the graph (5 chunks), because a test link is never in the training graph.

---

## 4. BEAR-GRN (TF -> target gene)

Scored with a line-by-line Python port of BEAR-GRN's R code, checked against R to 7 digits
(`src/benchmarks/bear_metrics.py`).

### 4.1 Universe
Universe = (TFs in the ground truth) x (target genes in the ground truth). Pairs a method does not output get
score 0.

### 4.2 AUROC (BEAR-GRN definition)
Computed **only on the edges the method outputs** (inside the universe), with the larger class down-sampled so links
and non-links are 1:1. A method that outputs few edges is judged only on those.

### 4.3 AUPRC (BEAR-GRN definition)
Computed on the **whole universe**: all links + a random sample of 10 x as many non-links (all non-links if fewer),
unoutput pairs scored 0, area by PRROC's continuous interpolation. Repeated over 20 random draws; mean reported.
**Random AUPRC** = the same with random scores = about the share of links (K562 ChIP 0.331).

**Toy:** 3 TFs x 4 genes = 12 pairs, 4 links (3 to gene g1, 1 to g2); g4 is not measured in the data.
- Random: 4 / 12 = 0.333.
- "All measured pairs" (score 1 for the 9 pairs with g1-g3, 0 for g4): **0.444**, above random with no ranking at
  all, because unmeasured genes hold no links.
- A sparse method that outputs only T1-g1 (0.9) and T1-g3 (0.8): 0.529.
- In-degree (g1 = 3, g2 = 1, g3 = g4 = 0): 0.949.

### 4.4 Baselines run on BEAR-GRN (13 in total, all 9 datasets)
- **9 released methods** (CellOracle, DIRECT-NET, FigR, GRaNIE, LINGER, Pando-GLM, Pando-XGB, SCENIC+, TRIPOD):
  their published networks, re-scored. DIRECT-NET has no scores, so no AUROC / AUPRC.
- **In-degree** (`scripts/bear/23_bear_baselines.py`): same TF-disjoint 5 folds as MeVD-GRN. For a held-out TF t
  and gene g: score = number of **training** TFs that have g as a target. Every held-out TF gets the same gene
  ranking. Uses labels (of other TFs), no cell data. Averaged over 5 fold seeds.
  **Toy:** training TFs T2, T3, T4 target {g1, g2}, {g1}, {g1, g3}: in-degree g1 = 3, g2 = 1, g3 = 1, g4 = 0, so
  the held-out T1 is predicted g1 first.
- **Pearson**: |correlation| of TF and gene expression across cells (CP10k + log1p). No labels.
- **GRNBoost2**: importance of the TF for the gene. No labels.
- **All measured pairs** ("coverage"): every TF x measured-gene pair scores 1. No ranking; AUROC = 0.5 by ties.

### 4.5 Pre-registered "WIN" (per dataset x ground truth x metric)
MeVD-GRN (mean over seeds) > LINGER **and** > in-degree **and** > all measured pairs **and** > random.

### 4.6 Hub control: AUROC within in-degree groups (own construction, pre-registered)
Genes are grouped by in-degree: group 0 = in-degree 0, then 20 quantile bins of the non-zero in-degrees. Inside each
group, each score is replaced by its percentile within the group; the pooled percentiles are scored with BEAR's
AUROC. In-degree itself is constant inside a group, so it scores 0.5; whatever is left above 0.5 is ranking
information **beyond** "this gene is a hub". (`src/benchmarks/bear_hub.py`.)

**Toy:** 8 edges. Hub genes (in-degree 8): links 0.9, 0.7, 0.6, non-link 0.8. Non-hub genes (in-degree 1): link 0.5,
non-links 0.4, 0.3, 0.2. Raw AUROC = 0.81, in-degree alone = 0.75. Within groups: among hubs the method ranks the
links poorly (1 of 3 above the non-link), among non-hubs perfectly; within-group AUROC = **0.625**.

---

## 5. Statistics in the tables

- **Per cell type:** mean over the 10 folds (and over seeds).
- **Headline number:** mean over the 12 cell types.
- **95% CI: cluster bootstrap over cell types.** Draw 12 cell types with replacement, average, repeat 2000 times;
  CI = 2.5th and 97.5th percentiles. Cell types, not folds, are resampled because folds of one cell type share data.
- **diff (p):** mean of (method - reference) over the cell type x fold pairs (120 pairs); p from a paired
  **Wilcoxon signed-rank test** on those differences (rank |differences|, compare the sums of ranks of the positive and
  negative ones). Only fold-paired, so both methods are compared on the same test pairs.

---

## 6. Likely questions, short answers

- **"Is the degree-matched / caliper test from a paper?"** The idea (negatives matched in degree) is published for
  protein interactions (Yu et al., Bioinformatics 2010) and as a GRN preprint (Stock et al., bioRxiv 2025). The
  nearest-3 and the 0.1 caliper are own choices; caliper matching itself is standard in causal inference.
- **"Why AUROC in the hub-controlled columns?"** scMultiomeGRN's trapezoidal AUPR gives 0.75 to a constant score
  (section 1.2); AUROC gives 0.5.
- **"Why does the degree product still score 0.653 on the degree-matched test?"** Non-links between two hubs are rare
  in these dense graphs, so the 3 nearest non-links are often still less hub-like. The caliper test fixes this by
  dropping unmatched links: there the degree product scores 0.506.
- **"Why does 'all measured pairs' beat LINGER on BEAR-GRN AUPRC?"** BEAR-GRN scores every unoutput pair as 0 over
  the whole universe (4.3), so a dense output earns AUPRC and a sparse one (LINGER) is penalised; AUROC (4.2) does
  not have this issue.
- **"Is the in-degree baseline fair to compare with?"** It uses exactly the labels MeVD-GRN trains on (training TFs
  of the same fold), and no cell data, so beating it is the minimum for a supervised model to show it uses the data.
- **"Any test-to-train leakage?"** Every split is built so that test pairs (and, in the held-out-TF split, test TFs)
  never appear in training; scripts assert this per fold.
