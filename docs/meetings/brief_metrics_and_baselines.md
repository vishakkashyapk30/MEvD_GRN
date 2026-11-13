# Brief: the math behind the metrics, tests and baselines

Companion to `presentation_18th_oct_14th_nov.md`. Every definition below is what the code actually computes
(file names in brackets). Toy numbers were computed with the same code paths.

## Contents

1. [Notation](#notation)
2. [Metrics](#s1)
   - [1.1 AUROC](#s1-1)
   - [1.2 AUPR, as scMultiomeGRN computes it](#s1-2)
   - [1.3 AP (average precision)](#s1-3)
   - [1.4 Random level](#s1-4)
3. [Test sets on the scMultiomeGRN benchmark](#s2)
   - [2.1 Official test](#s2-1)
   - [2.2 Degree-matched test](#s2-2)
   - [2.3 Caliper-matched test](#s2-3)
   - [2.4 Held-out-TF split](#s2-4)
   - [2.5 Caliper-matched test on pairs with one unseen TF](#s2-5)
4. [Baselines on the scMultiomeGRN benchmark](#s3)
   - [3.1 Link-prediction heuristics](#s3-1)
   - [3.2 Look-alike TFs](#s3-2)
   - [3.3 Label-free statistical scores](#s3-3)
   - [3.4 Supervised models on pair features](#s3-4)
   - [3.5 Link counts + features](#s3-5)
5. [BEAR-GRN](#s4)
   - [4.1 Universe](#s4-1)
   - [4.2 AUROC (BEAR-GRN definition)](#s4-2)
   - [4.3 AUPRC (BEAR-GRN definition)](#s4-3)
   - [4.4 Baselines run on BEAR-GRN](#s4-4)
   - [4.5 Pre-registered "WIN"](#s4-5)
   - [4.6 Hub control: AUROC within in-degree groups](#s4-6)
6. [Statistics in the tables](#s5)
7. [Likely questions, short answers](#s6)

---

<a id="notation"></a>
## Notation

- A **link** is a pair in the ground-truth network (label $y = 1$); a **non-link** is a pair that is not ($y = 0$).
- A model gives every pair a **score** $s$; higher means "more likely a link".
- On scMultiomeGRN the pairs are TF–TF and undirected; on BEAR-GRN they are TF $\to$ target gene.
- $A$ = training adjacency matrix ($A_{uv} = 1$ if $u$–$v$ is a training link), $d_u = \sum_w A_{uw}$ = degree of
  TF $u$ (its number of training links).
- $\ell_u = \log(1 + d_u)$ = log degree, used by the matched tests.

---

<a id="s1"></a>
## 1. Metrics

<a id="s1-1"></a>
### 1.1 AUROC (area under the ROC curve)

Pick one link and one non-link at random. AUROC is the probability that the link gets the higher score:

$$
\text{AUROC} = \frac{1}{n_+ \, n_-} \sum_{i \in \text{links}} \sum_{j \in \text{non-links}}
\Big( \mathbb{1}[s_i > s_j] + \tfrac{1}{2}\,\mathbb{1}[s_i = s_j] \Big)
$$

0.5 = random, 1.0 = perfect. It does not depend on the share of positives, so it is comparable across test sets.

**Toy:** links scored $0.9, 0.8, 0.4, 0.3$; non-links scored $0.7, 0.35, 0.2, 0.1$.
There are $4 \times 4 = 16$ (link, non-link) pairs; the link wins 13 of them, so AUROC $= 13/16 = \mathbf{0.8125}$.

<a id="s1-2"></a>
### 1.2 AUPR (area under the precision-recall curve), as scMultiomeGRN computes it

Walk down the ranking from the top score. At each cut-off $t$:

$$
\text{precision}(t) = \frac{\#\{\text{links with } s \ge t\}}{\#\{\text{pairs with } s \ge t\}}, \qquad
\text{recall}(t) = \frac{\#\{\text{links with } s \ge t\}}{\#\{\text{links}\}}
$$

AUPR is the area under precision ($y$-axis) against recall ($x$-axis). scMultiomeGRN's `metric_fn` uses the
**trapezoidal rule** (`sklearn.metrics.auc(recall, precision)`), i.e. it joins consecutive points
$(r_k, p_k)$ with straight lines:

$$
\text{AUPR}_{\text{trap}} = \sum_k (r_{k+1} - r_k)\,\frac{p_k + p_{k+1}}{2}
$$

(`src/benchmarks/scmgrn_protocol.py: metric_fn`)

**The tie problem (why hub-controlled columns use AUROC).** If a method gives every pair the same score, the curve
has only two points: $(r, p) = (0, 1)$ and $(1, \pi)$, where $\pi$ = share of links. The trapezoid between them is

$$
\text{AUPR}_{\text{trap}} = \frac{1 + \pi}{2}
$$

With $\pi = 0.5$ (the official test is balanced) this is $\mathbf{0.75}$ **for a method that knows nothing**, while
its AUROC is 0.5. Linear interpolation in precision-recall space overestimates the area; this is a known result
(Davis and Goadrich, ICML 2006).

**Toy:** the 8 scores of 1.1 give $\text{AUPR}_{\text{trap}} = 0.835$. Setting all 8 scores to 0.5 gives
AUROC $= 0.500$ but $\text{AUPR}_{\text{trap}} = \mathbf{0.750}$.

<a id="s1-3"></a>
### 1.3 AP (average precision)

Step-wise area, no straight-line interpolation (`sklearn.metrics.average_precision_score`):

$$
\text{AP} = \sum_k (r_k - r_{k-1})\, p_k
\;=\; \frac{1}{n_+} \sum_{i \in \text{links}} \text{precision at the rank of } i \quad (\text{no ties})
$$

All-tied scores give $\text{AP} = \pi = 0.5$. Toy of 1.1: AP $= 0.854$. The tables report it next to AUPR.

<a id="s1-4"></a>
### 1.4 Random level

- AUROC: 0.5 always.
- AUPR / AP: the share of positives $\pi$ in the test set (0.5 on the balanced scMultiomeGRN tests), but see 1.2:
  trapezoidal AUPR of a constant score is $(1 + \pi)/2 = 0.75$. Tables print $\pi$ as "Random AUPR".

---

<a id="s2"></a>
## 2. Test sets on the scMultiomeGRN benchmark

<a id="s2-1"></a>
### 2.1 Official test (scMultiomeGRN paper; re-implemented byte-for-byte)

The links of a cell type are split into 10 folds (KFold, seed 666). Fold $k$: test = fold $k$'s links **plus the
same number of random non-links**; the other links are training (80%) and validation (20%).

Problem: a random non-link usually involves low-degree TFs, while a link usually involves a hub. So "is either TF
a hub?" separates them without any biology.

<a id="s2-2"></a>
### 2.2 Degree-matched test (own variant; idea from Yu et al. 2010 and the Stock et al. 2025 preprint)

Same test links. For a test link $(u, v)$ take the sorted log degrees

$$
\text{lo} = \min(\ell_u, \ell_v), \qquad \text{hi} = \max(\ell_u, \ell_v)
$$

and for every candidate non-link $(u', v')$ not used in training or validation

$$
\text{dist} = \lvert \text{lo}' - \text{lo} \rvert + \lvert \text{hi}' - \text{hi} \rvert .
$$

Take the **3 nearest** candidates and draw one at random as the link's partner (`degree_matched_test`;
deterministic per fold, so every method sees the same pairs).

**Toy:** link between TFs of degree 10 and 2: $\text{lo} = \log 3 = 1.099$, $\text{hi} = \log 11 = 2.398$.

| Candidate non-link (degrees) | lo' | hi' | dist |
|---|---|---|---|
| (2, 11) | 1.099 | 2.485 | 0 + 0.087 = **0.087** |
| (1, 10) | 0.693 | 2.398 | 0.406 + 0 = **0.406** |
| (3, 12) | 1.386 | 2.565 | 0.287 + 0.167 = **0.454** |
| (5, 6) | 1.792 | 1.946 | 0.693 + 0.452 = 1.145 |

The first three are the nearest; one is drawn at random.

Limitation: in dense graphs there may be no close non-link, so the 3 "nearest" can still differ in degree. That is
why the degree product still scores 0.653 AUROC here instead of 0.5.

<a id="s2-3"></a>
### 2.3 Caliper-matched test (own construction; caliper idea from propensity-score matching, Austin 2011)

Stricter version of 2.2: a test link is kept **only if** some candidate non-link satisfies

$$
\lvert \text{lo}' - \text{lo} \rvert \le 0.1 \quad \text{and} \quad \lvert \text{hi}' - \text{hi} \rvert \le 0.1 ;
$$

one such non-link is drawn at random. Links without a match are **dropped**. Every kept link has a partner with
(almost) the same degrees, so any rule that only looks at degrees scores about 0.5 AUROC by construction
(`degree_caliper_test`).

**Toy:** in 2.2 only the (2, 11) candidate is inside the caliper (differences 0 and 0.087), so it is the partner.
If it did not exist, the link would be dropped.

<a id="s2-4"></a>
### 2.4 Held-out-TF split: pairs with one unseen TF, pairs of two unseen TFs (Park and Marcotte, Nature Methods 2012)

TFs, not links, are split into 10 groups (seed 666). Fold $k$: group $k$ = test TFs, group $k+1$ = validation TFs,
the rest = training TFs. Training uses only pairs of two training TFs.

- **Pairs with one unseen TF** (the paper calls them C2): one test TF and one training TF. The test TF has **no**
  training links ($d = 0$), so a model cannot count its links.
- **Pairs of two unseen TFs** (C3): both TFs of the pair are test TFs.

Each class gets as many random non-links as links. Pairs of a test TF with a validation TF are used nowhere
(`make_tf_holdout_splits`).

<a id="s2-5"></a>
### 2.5 Caliper-matched test on pairs with one unseen TF (own construction)

For a link $(t, r)$ with hidden TF $t$ and training TF $r$, the partner non-link is $(t, w)$: the **same hidden TF**
$t$ and a training TF $w$ with

$$
\lvert \ell_w - \ell_r \rvert \le 0.1 .
$$

Links without such a $w$ are dropped (`tfholdout_degree_matched_test`). Both ends are matched: $t$ is the same TF,
$r$ and $w$ have the same degree.

**Toy:** link $t$–$r$ with $d_r = 20$ ($\ell_r = \log 21 = 3.045$). Non-links of $t$: $w_1$ with $d = 21$
($\ell = 3.091$, difference 0.046, inside) and $w_2$ with $d = 5$ ($\ell = 1.792$, outside). Partner $= w_1$.

---

<a id="s3"></a>
## 3. Baselines on the scMultiomeGRN benchmark

All use only the fold's **training** links (or no links at all); every score matrix is made symmetric
(`scripts/scmgrn/22_scmgrn_simple_baselines.py` and `26_scmgrn_cf_probe.py`, branch `baselines`).

**Toy training graph** used in 3.1 and 3.2: TFs A, B, C, D, E with links A–B, A–C, A–D, B–C, D–E.

```
    B ─── C
     \   /
       A ─── D ─── E
```

Degrees: $d_A = 3$, $d_B = d_C = d_D = 2$, $d_E = 1$. Score the non-linked pair **B–D**.
$N(u)$ = set of TFs linked to $u$: $N(B) = \{A, C\}$, $N(D) = \{A, E\}$, common neighbours $= \{A\}$.

<a id="s3-1"></a>
### 3.1 Link-prediction heuristics (no cell data)

| Name | Formula for pair $(u, v)$ | Toy B–D |
|---|---|---|
| Degree product (preferential attachment) | $(d_u + 1)(d_v + 1)$ | $3 \times 3 = 9$ |
| Degree sum | $d_u + d_v$ | $4$ |
| Common neighbours | $\lvert N(u) \cap N(v) \rvert$ | $1$ |
| Jaccard | $\dfrac{\lvert N(u) \cap N(v) \rvert}{\lvert N(u) \cup N(v) \rvert}$ | $1/3$ |
| Adamic-Adar | $\sum_{w \in N(u) \cap N(v)} \dfrac{1}{\log d_w}$ | $1/\log 3 = 0.910$ |
| Resource allocation | $\sum_{w \in N(u) \cap N(v)} \dfrac{1}{d_w}$ | $1/3$ |
| Katz | $\sum_{L=1}^{4} 0.05^L \, (A^L)_{uv}$ | walks: $L{=}2$: 1, $L{=}3$: 1, $L{=}4$: 5 $\Rightarrow 0.0027$ |
| Random walk with restart | $R = 0.15\,(I - 0.85\,P)^{-1}$, $P = D^{-1} A$; score $R_{uv} + R_{vu}$ | how often a walker from $u$ (restarting at $u$ w.p. 0.15) is at $v$ |
| Low-rank SVD | $(U_{16} \Sigma_{16} V_{16}^\top)_{uv}$ of $A$ | matrix completion |
| Low-rank spectral | rank-16 reconstruction of $D^{-1/2} A D^{-1/2}$ | degree-normalised version |

$(A^L)_{uv}$ counts walks of length $L$ from $u$ to $v$. Adamic-Adar and resource allocation down-weight common
neighbours that are hubs.

<a id="s3-2"></a>
### 3.2 Look-alike TFs (collaborative filtering)

Similarity of two TFs = cosine of their training-link rows (diagonal set to 0):

$$
C_{uw} = \frac{A_u \cdot A_w}{\lVert A_u \rVert \, \lVert A_w \rVert}, \qquad
\text{score}(u, v) = \sum_w C_{uw} A_{wv} \;+\; \sum_w C_{vw} A_{wu}
$$

"Are TFs that look like $u$ linked to $v$?" (and the same from $v$'s side).

**Toy (B–D):** B's row $\{A, C\}$; its similarities are $C_{BC} = C_{BD} = 0.5$, $C_{BA} = 0.408$, $C_{BE} = 0$.
TFs linked to D are A and E, so $\sum_w C_{Bw} A_{wD} = 0.408$. From D's side, TFs linked to B are A and C:
$C_{DA} = 0$, $C_{DC} = 0.5$, giving 0.5. Score $= 0.908$.

Why it matters: the ground truth is a motif scan, so TFs of one family (similar motifs) have near-identical rows.
"Normalised" divides by $\sum_w C_{uw}$. "Geneformer look-alikes": the same idea with the 5 most similar TFs by
Geneformer embedding instead of link rows.

<a id="s3-3"></a>
### 3.3 Label-free statistical scores (no labels; BEELINE / geneRNIB family)

$x_u$ = TF $u$'s expression across cells (or MAESTRO accessibility for that row).

| Score | Definition |
|---|---|
| Pearson / Spearman | $\lvert \text{corr}(x_u, x_v) \rvert$ (Spearman on ranks) |
| Partial correlation | $P = \Sigma^{-1}$ (Ledoit-Wolf shrinkage), $\rho_{uv} = \dfrac{-P_{uv}}{\sqrt{P_{uu} P_{vv}}}$; correlation left after removing all other TFs |
| Mutual information | 8 equal-width bins per TF; $\text{MI} = \sum_{a, b} p(a, b) \log \dfrac{p(a, b)}{p(a)\,p(b)}$ |
| CLR | $z_u = \max\!\big(0, \frac{\text{MI}_{uv} - \mu_u}{\sigma_u}\big)$ over row $u$ (same for $v$); score $= \sqrt{z_u^2 + z_v^2}$ |
| GRNBoost2 between the two TFs | $\max(\text{imp}_{u \to v}, \text{imp}_{v \to u})$ from gradient-boosting regression |
| GENIE3 | random forest predicting each TF from all others; importance, symmetrised |
| Ridge / lasso | regress each TF on all others ($\alpha = 1$ / $0.01$); $\lvert \beta \rvert$, symmetrised |
| Cosine of GRNBoost2 profiles / Geneformer | $\cos(g_u, g_v)$, $\cos(f_u, f_v)$ |

<a id="s3-4"></a>
### 3.4 Supervised models on pair features

Per-TF vectors (MAESTRO, GRNBoost2 profile, Geneformer) are first reduced to 32 numbers by SVD (no labels). The
**pair feature** is symmetric, so $(u, v)$ and $(v, u)$ get the same score:

$$
\phi(u, v) = \big[\, x_u \odot x_v ,\; \lvert x_u - x_v \rvert \,\big] \qquad (\odot = \text{element-wise product})
$$

- "X only": logistic regression on $\phi$ of one input.
- "Side by side" (given independently): $[\phi_{\text{MAESTRO}}, \phi_{\text{GRNBoost2}}]$.
- "MAESTRO multiplied by GRNBoost2": per TF first, $z_u = m_u \odot g_u$, then $\phi$ on $z$.
- "Product of the two TFs' features only": only the $x_u \odot x_v$ parts.
- Logistic regression: L2 strength chosen from $\{0.001, 0.01, 0.1, 1, 10\}$ by 3-fold CV on training pairs only.
  Linear regression = ridge on the 0/1 label. Random forest = 300 trees.

**Toy:** $x_A = (1, 0)$, $x_B = (1, 1)$: $\phi(A, B) = [(1, 0), (0, 1)]$; logistic regression learns one weight per
entry.

<a id="s3-5"></a>
### 3.5 Link counts + features

Logistic regression on

$$
\big[\, \ell_u + \ell_v ,\; \lvert \ell_u - \ell_v \rvert ,\; \ell_u \ell_v \,\big]
$$

(plus the 3.1 scores for "neighbour scores"), optionally with the 3.4 features. Leak guard: a training link's own
features are computed with that link removed from the graph (5 chunks), because a test link is never in the training
graph.

---

<a id="s4"></a>
## 4. BEAR-GRN (TF $\to$ target gene)

Scored with a line-by-line Python port of BEAR-GRN's R code, checked against R to 7 digits
(`src/benchmarks/bear_metrics.py`).

<a id="s4-1"></a>
### 4.1 Universe

$$
\mathcal{U} = \{\text{TFs in the ground truth}\} \times \{\text{target genes in the ground truth}\}
$$

Pairs in $\mathcal{U}$ that a method does not output get score 0.

<a id="s4-2"></a>
### 4.2 AUROC (BEAR-GRN definition)

Computed **only on the edges the method outputs** (inside $\mathcal{U}$), with the larger class down-sampled so links
and non-links are 1:1. A method that outputs few edges is judged only on those.

<a id="s4-3"></a>
### 4.3 AUPRC (BEAR-GRN definition)

Computed on the **whole universe**: all links plus a random sample of $\min(10\,n_+,\; n_-)$ non-links, unoutput pairs
scored 0, area by PRROC's continuous interpolation. Repeated over 20 random draws; the mean is reported.
**Random AUPRC** = the same with random scores $\approx$ the share of links (K562 ChIP: 0.331).

**Toy:** 3 TFs $\times$ 4 genes $= 12$ pairs, 4 links (3 to gene g1, 1 to g2); g4 is not measured in the data.

| Method | AUPRC |
|---|---|
| Random | $4/12 = 0.333$ |
| All measured pairs (score 1 for the 9 pairs with g1–g3, 0 for g4) | **0.444**: above random with no ranking, because unmeasured genes hold no links |
| Sparse method: only T1–g1 (0.9) and T1–g3 (0.8) | 0.529 |
| In-degree (g1 = 3, g2 = 1, g3 = g4 = 0) | 0.949 |

<a id="s4-4"></a>
### 4.4 Baselines run on BEAR-GRN (13 in total, all 9 datasets)

- **9 released methods** (CellOracle, DIRECT-NET, FigR, GRaNIE, LINGER, Pando-GLM, Pando-XGB, SCENIC+, TRIPOD):
  their published networks, re-scored. DIRECT-NET has no scores, so no AUROC / AUPRC.
- **In-degree** (`scripts/bear/23_bear_baselines.py`): same TF-disjoint 5 folds as MeVD-GRN. For a held-out TF $t$
  and gene $g$:

  $$
  \text{score}(t, g) = \#\{\, t' \in \text{training TFs} : t' \to g \text{ in the ground truth} \,\}
  $$

  Every held-out TF gets the same gene ranking. Uses labels (of other TFs), no cell data. Averaged over 5 fold
  seeds.
  **Toy:** training TFs T2, T3, T4 target $\{g1, g2\}$, $\{g1\}$, $\{g1, g3\}$: in-degree g1 = 3, g2 = 1, g3 = 1,
  g4 = 0, so the held-out T1 is predicted g1 first.
- **Pearson**: $\lvert \text{corr}(x_t, x_g) \rvert$ across cells (CP10k + log1p). No labels.
- **GRNBoost2**: importance of the TF for the gene. No labels.
- **All measured pairs** ("coverage"): every TF $\times$ measured-gene pair scores 1. No ranking; AUROC $= 0.5$ by
  ties.

<a id="s4-5"></a>
### 4.5 Pre-registered "WIN" (per dataset $\times$ ground truth $\times$ metric)

$$
\text{WIN} \iff m_{\text{MeVD}} > m_{\text{LINGER}} \;\wedge\; m_{\text{MeVD}} > m_{\text{in-degree}}
\;\wedge\; m_{\text{MeVD}} > m_{\text{all measured pairs}} \;\wedge\; m_{\text{MeVD}} > m_{\text{random}}
$$

with $m$ = AUROC or AUPRC and $m_{\text{MeVD}}$ = mean over seeds.

<a id="s4-6"></a>
### 4.6 Hub control: AUROC within in-degree groups (own construction, pre-registered)

1. Group genes by in-degree: group 0 = in-degree 0, then 20 quantile bins of the non-zero in-degrees.
2. Inside each group $G$, replace each score by its percentile within the group:
   $\tilde{s}_i = \text{rank}_G(s_i) / \lvert G \rvert$.
3. Score the pooled $\tilde{s}$ with BEAR's AUROC.

In-degree is constant inside a group, so it scores 0.5; whatever is left above 0.5 is ranking information **beyond**
"this gene is a hub" (`src/benchmarks/bear_hub.py`).

**Toy:** 8 edges.

| Group | Links (scores) | Non-links (scores) | Within-group order |
|---|---|---|---|
| Hub genes (in-degree 8) | 0.9, 0.7, 0.6 | 0.8 | 1 of 3 links above the non-link |
| Non-hub genes (in-degree 1) | 0.5 | 0.4, 0.3, 0.2 | the link above all non-links |

Raw AUROC $= 0.81$ and in-degree alone $= 0.75$, but within groups AUROC $= \mathbf{0.625}$.

---

<a id="s5"></a>
## 5. Statistics in the tables

- **Per cell type** $c$: $\bar m_c$ = mean over the 10 folds (and over seeds).
- **Headline number:** $\bar m = \frac{1}{12} \sum_c \bar m_c$.
- **95% CI, cluster bootstrap over cell types:** for $b = 1, \dots, 2000$ draw 12 cell types with replacement and
  take $\bar m^{(b)}$ = their mean; CI = 2.5th and 97.5th percentiles of $\bar m^{(b)}$. Cell types, not folds, are
  resampled because folds of one cell type share data.
- **diff (p):** $\Delta_{c,k} = m_{c,k}^{\text{method}} - m_{c,k}^{\text{ref}}$ for every cell type $c$ and fold $k$
  (120 pairs); the table shows the mean $\Delta$ and the p-value of a **Wilcoxon signed-rank test** on the
  $\Delta_{c,k}$ (rank the $\lvert \Delta \rvert$, compare the rank sums of the positive and negative ones). Both
  methods are compared on the same test pairs.

---

<a id="s6"></a>
## 6. Likely questions, short answers

- **"Is the degree-matched / caliper test from a paper?"** The idea (negatives matched in degree) is published for
  protein interactions (Yu et al., Bioinformatics 2010) and as a GRN preprint (Stock et al., bioRxiv 2025). The
  nearest-3 rule and the 0.1 caliper are own choices; caliper matching itself is standard in causal inference.
- **"Why AUROC in the hub-controlled columns?"** scMultiomeGRN's trapezoidal AUPR gives $(1 + \pi)/2 = 0.75$ to a
  constant score ([1.2](#s1-2)); AUROC gives 0.5.
- **"Why does the degree product still score 0.653 on the degree-matched test?"** Non-links between two hubs are rare
  in these dense graphs, so the 3 nearest non-links are often still less hub-like. The caliper test fixes this by
  dropping unmatched links: there the degree product scores 0.506.
- **"Why does 'all measured pairs' beat LINGER on BEAR-GRN AUPRC?"** BEAR-GRN scores every unoutput pair as 0 over
  the whole universe ([4.3](#s4-3)), so a dense output earns AUPRC and a sparse one (LINGER) is penalised; AUROC
  ([4.2](#s4-2)) does not have this issue.
- **"Is the in-degree baseline fair to compare with?"** It uses exactly the labels MeVD-GRN trains on (training TFs
  of the same fold) and no cell data, so beating it is the minimum for a supervised model to show it uses the data.
- **"Any test-to-train leakage?"** Every split is built so that test pairs (and, in the held-out-TF split, test TFs)
  never appear in training; the scripts assert this per fold.
