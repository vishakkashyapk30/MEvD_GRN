# Brief: the math behind the metrics, tests and baselines

Companion to `presentation_18th_oct_14th_nov.md`. Every definition below is what the code actually computes
(file names in brackets). Toy numbers were computed with the same code paths.

## Contents

1. [Notation](#notation)
2. [Glossary of terms used in the deck](#glossary)
   - [Training graph](#g-training-graph) · [Degree, hub](#g-degree) · [In-degree](#g-in-degree)
   - [Link-prediction heuristic](#g-heuristic) · [Label-free vs supervised](#g-label-free)
   - [MeVD-GRN parts: Geneformer embedding](#g-geneformer) · [ATAC tower](#g-atac-tower) ·
     [Graph network](#g-graph-network) · [Training-link relation](#g-training-link-relation) ·
     [Floor model](#g-floor-model)
   - [Degree-matched and caliper-matched, in one paragraph](#g-matched)
   - [Held-out-TF split](#g-held-out-tf) · [TF-disjoint cross-fitting (BEAR-GRN)](#g-cross-fitting)
   - [L1 and L2 regimes (BEAR-GRN)](#g-l1-l2) · [Universe, all measured pairs](#g-universe) ·
     [Hub control](#g-hub-control)
3. [Metrics](#s1)
   - [1.1 AUROC](#s1-1) · [1.2 AUPR, as scMultiomeGRN computes it](#s1-2) · [1.3 AP](#s1-3) ·
     [1.4 Random level](#s1-4)
4. [Test sets on the scMultiomeGRN benchmark](#s2)
   - [2.1 Official test](#s2-1) · [2.2 Degree-matched test](#s2-2) · [2.3 Caliper-matched test](#s2-3) ·
     [2.4 Held-out-TF split](#s2-4) · [2.5 Caliper-matched test on pairs with one unseen TF](#s2-5)
5. [Baselines on the scMultiomeGRN benchmark](#s3)
   - [3.1 Link-prediction heuristics](#s3-1):
     [degree product](#m-degree-product) · [degree sum](#m-degree-sum) ·
     [common neighbours](#m-common-neighbours) · [Jaccard](#m-jaccard) · [Adamic-Adar](#m-adamic-adar) ·
     [resource allocation](#m-resource-allocation) · [Katz](#m-katz) · [random walk with restart](#m-rwr) ·
     [low-rank SVD](#m-svd) · [low-rank spectral](#m-spectral)
   - [3.2 Look-alike TFs](#s3-2): [look-alike TFs](#m-lookalike) · [normalised](#m-lookalike-norm) ·
     [Geneformer look-alikes](#m-gf-lookalike)
   - [3.3 Label-free statistical scores](#s3-3):
     [Pearson / Spearman](#m-pearson) · [partial correlation](#m-partial) · [mutual information](#m-mi) ·
     [CLR](#m-clr) · [GRNBoost2](#m-grnboost2) · [GENIE3](#m-genie3) · [ridge / lasso](#m-ridge-lasso) ·
     [cosine of profiles / embeddings](#m-cosine)
   - [3.4 Supervised models on pair features](#s3-4):
     inputs: [MAESTRO accessibility](#m-in-maestro) ([RP formula](#m-maestro-rp)) · [GRNBoost2 profiles](#m-in-grnboost) ·
     [Geneformer](#m-in-geneformer) · [joint expression histograms](#m-in-histogram) · [promoter openness](#m-in-openness);
     [pair feature](#m-pair-feature); variants: [single input](#m-single) · [side by side](#m-side-by-side) ·
     [multiplied](#m-multiplied) · [product only](#m-product-only) · [linear regression](#m-linreg) ·
     [random forest](#m-rf); [classifiers](#m-models); [reading the results](#m-reading)
   - [3.5 Link counts + features](#s3-5)
6. [BEAR-GRN](#s4)
   - [4.1 Universe](#s4-1) · [4.2 AUROC](#s4-2) · [4.3 AUPRC](#s4-3)
   - [4.4 Baselines](#s4-4): [released methods](#m-bear-released) · [in-degree](#m-bear-indegree) ·
     [Pearson](#m-bear-pearson) · [GRNBoost2](#m-bear-grnboost2) · [all measured pairs](#m-bear-coverage)
   - [4.5 Pre-registered "WIN"](#s4-5) · [4.6 Hub control](#s4-6)
7. [Statistics in the tables](#s5)
8. [Likely questions, short answers](#s6)

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

<a id="glossary"></a>
## Glossary of terms used in the deck

<a id="g-training-graph"></a>
**Training graph.** For one fold, the TF-TF network made only of that fold's **training links**: nodes = the cell
type's TFs, edges = known links in the training part. Test links (and validation links) are removed, so a model
that reads the training graph never sees the answer for a test pair. Example: if the ground truth has 1,000 links,
fold 1 tests 100 of them; the training graph holds about 720 (the rest minus 20% for validation).

<a id="g-degree"></a>
**Degree, hub.** A TF's degree = its number of links in the training graph. A **hub** is a TF with many links.
On this benchmark hubs dominate: most true links involve a hub, most random non-links do not.

<a id="g-in-degree"></a>
**In-degree (BEAR-GRN).** For a target gene: how many TFs regulate it in the (training) ground truth. Directed
version of degree, counted on the target side.

<a id="g-heuristic"></a>
**Link-prediction heuristic.** A fixed formula on the training graph that scores an unlinked pair, with **no
learning and no cell data** (e.g. "how many partners do the two TFs share?"). Why it is run: it measures how much
of the benchmark can be solved from the **shape of the known network alone**. If such a formula matches a deep
model, the benchmark rewards network structure (hubs, shared partners), not biology from the cells.

<a id="g-label-free"></a>
**Label-free vs supervised.** Label-free = no ground-truth links used at all (e.g. correlation of expression).
Supervised = fitted on the fold's training links (e.g. logistic regression, MeVD-GRN, scMultiomeGRN).

<a id="g-geneformer"></a>
**Geneformer embedding.** A 768-number vector per gene from Geneformer, a model pre-trained on tens of millions of single-cell transcriptomes
elsewhere; it is not trained on these cells or labels.

<a id="g-atac-tower"></a>
**ATAC tower.** The branch of MeVD-GRN that reads each gene's chromatin accessibility (MAESTRO activity summary:
mean, variance, detection rate; plus promoter-openness descriptors) and turns it into a vector.

<a id="g-graph-network"></a>
**Graph network.** MeVD-GRN's graph neural network: each gene's vector is updated with messages from its
neighbours in label-free graphs built from the cells, a **co-expression graph** (genes with similar expression
across cells) and a **TF-candidate graph** (TF to genes that are accessible / co-expressed with it).

<a id="g-training-link-relation"></a>
**Training-link relation.** An optional extra edge type for the graph network: the **training graph** itself.
With it, MeVD-GRN also passes messages along the fold's known training links, which is exactly what
scMultiomeGRN's encoder does (same information for both). "Without the training-link relation" = the training
links are used only as labels to learn from, not as an input graph. Test links are never in it.

<a id="g-floor-model"></a>
**Floor model.** MeVD-GRN with all four parts above removed (no Geneformer, no ATAC tower, no graph network, no
training-link relation), smaller size. Each TF = 3 RNA numbers (mean, variance, detection rate) + a co-expression
term. It shows what is left with almost nothing.

<a id="g-matched"></a>
**Degree-matched and caliper-matched, in one paragraph.** On the official test, a true link usually involves a hub
and a random non-link usually does not, so "is it a hub?" already wins. Both tests replace the random non-links by
non-links whose two TFs have the **same numbers of training links** as the true link's TFs, so that hub-counting no
longer helps. **Degree-matched**: the partner is one of the 3 closest non-links (may still differ somewhat).
**Caliper-matched**: the partner must be within 0.1 on the log scale for both TFs, otherwise the link is dropped;
a pure counting rule then scores 0.5. Details: [2.2](#s2-2), [2.3](#s2-3).

<a id="g-held-out-tf"></a>
**Held-out-TF split.** Instead of hiding some links, hide whole TFs: about one tenth of the TFs have **no** links in
training. Test pairs have one hidden TF or two. Details: [2.4](#s2-4).

<a id="g-cross-fitting"></a>
**TF-disjoint cross-fitting (BEAR-GRN).** The ground-truth TFs are split into 5 groups; a model is trained on 4
groups' labels and scores the 5th group's TFs, rotating 5 times. Every scored TF is unseen in training.

<a id="g-l1-l2"></a>
**L1 and L2 regimes (BEAR-GRN).** Which labels MeVD-GRN may train on.
- **L1 (headline, all results in the deck):** labels of **other TFs in the same cell type** (the BEAR-GRN ground
  truth with the scored TFs removed, by cross-fitting).
- **L2 (strict, still running):** **no labels from this cell type at all**: training labels come only from a generic
  ChIP-seq compendium (SC-MO-GRN-DB's non-cell-type-specific ChIP network, RN005 human / RN006 mouse), with every TF
  that appears in the dataset's ground truth removed; one model scores all ground-truth TFs. This is the closest
  analogue to LINGER, which also learns only from external data.

<a id="g-universe"></a>
**Universe, all measured pairs (BEAR-GRN).** The universe = every (ground-truth TF, ground-truth target) pair; pairs
a method does not output count as score 0. "All measured pairs" outputs every pair whose gene is measured, all with
the same score. Details: [4.1](#s4-1), [4.3](#s4-3).

<a id="g-hub-control"></a>
**Hub control (BEAR-GRN).** Compare scores only among target genes with similar in-degree, so "this gene is a hub"
cannot help. Details: [4.6](#s4-6).

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
| <a id="m-degree-product"></a>Degree product (preferential attachment) | $(d_u + 1)(d_v + 1)$ | $3 \times 3 = 9$ |
| <a id="m-degree-sum"></a>Degree sum | $d_u + d_v$ | $4$ |
| <a id="m-common-neighbours"></a>Common neighbours | $\lvert N(u) \cap N(v) \rvert$ | $1$ |
| <a id="m-jaccard"></a>Jaccard | $\dfrac{\lvert N(u) \cap N(v) \rvert}{\lvert N(u) \cup N(v) \rvert}$ | $1/3$ |
| <a id="m-adamic-adar"></a>Adamic-Adar | $\sum_{w \in N(u) \cap N(v)} \dfrac{1}{\log d_w}$ | $1/\log 3 = 0.910$ |
| <a id="m-resource-allocation"></a>Resource allocation | $\sum_{w \in N(u) \cap N(v)} \dfrac{1}{d_w}$ | $1/3$ |
| <a id="m-katz"></a>Katz | $\sum_{L=1}^{4} 0.05^L \, (A^L)_{uv}$ | walks: $L{=}2$: 1, $L{=}3$: 1, $L{=}4$: 5 $\Rightarrow 0.0027$ |
| <a id="m-rwr"></a>Random walk with restart | $R = 0.15\,(I - 0.85\,P)^{-1}$, $P = D^{-1} A$; score $R_{uv} + R_{vu}$ | how often a walker from $u$ (restarting at $u$ w.p. 0.15) is at $v$ |
| <a id="m-svd"></a>Low-rank SVD | $(U_{16} \Sigma_{16} V_{16}^\top)_{uv}$ of $A$ | matrix completion |
| <a id="m-spectral"></a>Low-rank spectral | rank-16 reconstruction of $D^{-1/2} A D^{-1/2}$ | degree-normalised version |

$(A^L)_{uv}$ counts walks of length $L$ from $u$ to $v$. Adamic-Adar and resource allocation down-weight common
neighbours that are hubs.

<a id="s3-2"></a>
### 3.2 Look-alike TFs (collaborative filtering)

<a id="m-lookalike"></a>Similarity of two TFs = cosine of their training-link rows (diagonal set to 0):

$$
C_{uw} = \frac{A_u \cdot A_w}{\lVert A_u \rVert \, \lVert A_w \rVert}, \qquad
\text{score}(u, v) = \sum_w C_{uw} A_{wv} \;+\; \sum_w C_{vw} A_{wu}
$$

"Are TFs that look like $u$ linked to $v$?" (and the same from $v$'s side).

**Toy (B–D):** B's row $\{A, C\}$; its similarities are $C_{BC} = C_{BD} = 0.5$, $C_{BA} = 0.408$, $C_{BE} = 0$.
TFs linked to D are A and E, so $\sum_w C_{Bw} A_{wD} = 0.408$. From D's side, TFs linked to B are A and C:
$C_{DA} = 0$, $C_{DC} = 0.5$, giving 0.5. Score $= 0.908$.

Why it matters: the ground truth is a motif scan, so TFs of one family (similar motifs) have near-identical rows.
<a id="m-lookalike-norm"></a>"Normalised" divides by $\sum_w C_{uw}$. <a id="m-gf-lookalike"></a>"Geneformer look-alikes": the same idea with the 5 most similar TFs by
Geneformer embedding instead of link rows.

<a id="s3-3"></a>
### 3.3 Label-free statistical scores (no labels; BEELINE / geneRNIB family)

$x_u$ = TF $u$'s expression across cells (or MAESTRO accessibility for that row).

| Score | Definition |
|---|---|
| <a id="m-pearson"></a>Pearson / Spearman | $\lvert \text{corr}(x_u, x_v) \rvert$ (Spearman on ranks) |
| <a id="m-partial"></a>Partial correlation | $P = \Sigma^{-1}$ (Ledoit-Wolf shrinkage), $\rho_{uv} = \dfrac{-P_{uv}}{\sqrt{P_{uu} P_{vv}}}$; correlation left after removing all other TFs |
| <a id="m-mi"></a>Mutual information | 8 equal-width bins per TF; $\text{MI} = \sum_{a, b} p(a, b) \log \dfrac{p(a, b)}{p(a)\,p(b)}$ |
| <a id="m-clr"></a>CLR | $z_u = \max\!\big(0, \frac{\text{MI}_{uv} - \mu_u}{\sigma_u}\big)$ over row $u$ (same for $v$); score $= \sqrt{z_u^2 + z_v^2}$ |
| <a id="m-grnboost2"></a>GRNBoost2 between the two TFs | $\max(\text{imp}_{u \to v}, \text{imp}_{v \to u})$ from gradient-boosting regression |
| <a id="m-genie3"></a>GENIE3 | random forest predicting each TF from all others; importance, symmetrised |
| <a id="m-ridge-lasso"></a>Ridge / lasso | regress each TF on all others ($\alpha = 1$ / $0.01$); $\lvert \beta \rvert$, symmetrised |
| <a id="m-cosine"></a>Cosine of GRNBoost2 profiles / Geneformer | $\cos(g_u, g_v)$, $\cos(f_u, f_v)$ |

<a id="s3-4"></a>
### 3.4 Supervised models on pair features

These are the rows of the deck's "Supervised models on pair features" slide. Each one: (1) describe every TF by a
vector built from the cells, (2) turn two TFs' vectors into one **pair feature**, (3) fit a simple classifier on the
fold's training pairs (training links = 1, the remaining training non-links = 0), (4) score every test pair.
No neural network, no training graph as input.

<a id="m-inputs"></a>
#### 3.4.1 The inputs (what describes a TF)

| Input (deck name) | Per TF | Size before reduction | What it measures |
|---|---|---|---|
| <a id="m-in-maestro"></a>MAESTRO accessibility | RP score in every cell (formula below) | one number per cell | how open the chromatin around the TF's gene is, cell by cell |
| <a id="m-in-grnboost"></a>GRNBoost2 profiles | GRNBoost2 importance of this TF for every gene | one number per gene | which genes this TF's expression helps predict (its expression "targets") |
| <a id="m-in-geneformer"></a>Geneformer embedding | pre-trained Geneformer vector | 768 | general knowledge about the gene from other datasets |
| <a id="m-in-histogram"></a>Joint expression histograms | (a pair feature directly, see below) | 136 per pair | how the two TFs' expression values co-occur across cells |
| <a id="m-in-openness"></a>Promoter openness | 4 summary numbers | 4 | mean and max RP weight of the nearest peaks, their mean distance, peak count |

MAESTRO, GRNBoost2 and Geneformer vectors are compressed to **32 numbers per TF** by SVD (centre the TF x feature
matrix, keep the top 32 components, standardise each). No labels are used for this.

<a id="m-maestro-rp"></a>**MAESTRO regulatory potential** (scMultiomeGRN's ATAC input; MAESTRO "Enhanced" model,
decay $d_0 = 1$ kb, re-implemented in `src/benchmarks/scmgrn_features.py`). For gene $g$ in cell $c$:

$$
\text{RP}(g, c) = \sum_{\text{peaks } p} w_{g,p}\,\text{reads}(p, c), \qquad
w_{g,p} = \begin{cases}
1 / (\text{exon length of } g \text{ in kb}) & \text{peak in an exon of } g \\
2^{-\lvert d \rvert / d_0} & \text{peak within } \pm 2 \text{ kb of the TSS, or outside every gene body within 15 kb} \\
0 & \text{intron, another gene's body, or farther than 15 kb}
\end{cases}
$$

with $d$ = distance from the peak centre to the TSS. Toy: a promoter peak 0.5 kb away with 2 reads gives
$2^{-0.5} \times 2 = 1.41$; a peak 3 kb away with 1 read gives $2^{-3} = 0.125$; an exon peak (exons 2 kb) with 1
read gives $0.5$; a peak 20 kb away gives 0. RP $= 2.04$.

**Joint expression histogram** (scMultiomeGRN's edge feature). For TFs $u, v$: put every cell into a $16 \times 16$
grid by $(\log_{10}(x_u + 0.01), \log_{10}(x_v + 0.01))$, count cells per square $H$, convert to
$e = (\log_{10}(H / n_\text{cells} + 10^{-4}) + 4)/4$, average with its transpose so the order of $u, v$ does not
matter, and keep the 136 unique values. Already a pair feature, so no SVD and no $\phi$.

<a id="m-pair-feature"></a>
#### 3.4.2 From two TF vectors to one pair feature

$$
\phi(u, v) = \big[\, x_u \odot x_v ,\; \lvert x_u - x_v \rvert \,\big] \qquad (\odot = \text{element-wise product})
$$

Both halves are unchanged when $u$ and $v$ are swapped, so the score is symmetric. The product is large when both TFs
are high on the same component; the difference is small when the two TFs look alike.

<a id="m-variants"></a>
#### 3.4.3 The variants (one toy, all rows of the table)

Toy, with 2 numbers per TF instead of 32: MAESTRO $m_u = (1, 2)$, $m_v = (3, 1)$; GRNBoost2 $g_u = (0, 1)$,
$g_v = (2, 2)$.

| Deck row | Feature given to the classifier | Toy feature | Length |
|---|---|---|---|
| <a id="m-single"></a>MAESTRO accessibility only | $\phi(m_u, m_v)$ | $[(3, 2), (2, 1)]$ | 4 |
| GRNBoost2 profiles only | $\phi(g_u, g_v)$ | $[(0, 2), (2, 1)]$ | 4 |
| Geneformer embedding only | $\phi(f_u, f_v)$ | (same recipe, 32 + 32 numbers) | 64 |
| Joint expression histograms only | the 136 histogram values | – | 136 |
| Promoter openness only | $[\phi(o_u, o_v),\; o_u + o_v]$ | – | 12 |
| <a id="m-side-by-side"></a>MAESTRO and GRNBoost2 side by side ("given independently") | $[\phi(m_u, m_v),\; \phi(g_u, g_v)]$ | $[(3, 2), (2, 1), (0, 2), (2, 1)]$ | 8 |
| <a id="m-multiplied"></a>MAESTRO multiplied by GRNBoost2 | first $z_u = m_u \odot g_u$, then $\phi(z_u, z_v)$ | $z_u = (0, 2)$, $z_v = (6, 2)$ $\Rightarrow [(0, 4), (6, 0)]$ | 4 |
| <a id="m-product-only"></a>Product of the two TFs' features only | $[m_u \odot m_v,\; g_u \odot g_v]$ (no difference terms) | $[(3, 2), (0, 2)]$ | 4 |
| <a id="m-linreg"></a>Linear regression, side by side | side-by-side feature, ridge regression on the 0/1 label | as side by side | 8 |
| <a id="m-rf"></a>Random forest, side by side | side-by-side feature, random forest | as side by side | 8 |

<a id="m-models"></a>
#### 3.4.4 The classifiers

- **Logistic regression** (every row unless stated otherwise):
  $P(\text{link} \mid u, v) = \sigma(w^\top \phi + b)$, $\sigma(z) = 1 / (1 + e^{-z})$, fitted with an L2 penalty
  $\lambda \lVert w \rVert^2$. $\lambda$ is chosen from 5 values by 3-fold cross-validation on the training pairs only
  (scored by average precision). Features are standardised with training-pair statistics.
  Toy: side-by-side feature $(3, 2, 2, 1, 0, 2, 2, 1)$, all weights 0.1, $b = -1$: $\sigma(1.3 - 1) = 0.57$.
- **Linear regression**: ridge regression ($\alpha = 1$) predicting the 0/1 label directly; the prediction is the score.
- **Random forest**: 300 trees, at least 2 training pairs per leaf; score = the trees' average link probability. The one
  non-linear model of the slide.

<a id="m-reading"></a>
#### 3.4.5 Reading the results

- **One input at a time:** every input carries some signal on the official test (AUPR 0.59 to 0.73), GRNBoost2
  profiles most.
- **Side by side beats multiplied (0.788 vs 0.721):** with side by side the classifier can weigh MAESTRO and
  GRNBoost2 separately; multiplying first mixes them, and a near-zero value in one input erases the other
  (toy: $g_u$'s first entry 0 wipes out $m_u$'s first entry).
- **Product only is near chance (0.523):** dropping the $\lvert x_u - x_v \rvert$ terms removes most of the signal.
- **Linear regression = logistic regression (0.788 both); random forest is best (0.836):** some non-linear structure
  helps, still below scMultiomeGRN (0.919) and the link-count rules.
- **Caliper-matched test: all between 0.49 and 0.52.** The likely reading: what the classifiers learn from these profiles is mostly which TFs are hubs, so once both TFs of every test pair have matched link counts, nothing is
  left.

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

- <a id="m-bear-released"></a>**9 released methods** (CellOracle, DIRECT-NET, FigR, GRaNIE, LINGER, Pando-GLM, Pando-XGB, SCENIC+, TRIPOD):
  their published networks, re-scored. DIRECT-NET has no scores, so no AUROC / AUPRC.
- <a id="m-bear-indegree"></a>**In-degree** (`scripts/bear/23_bear_baselines.py`): same TF-disjoint 5 folds as MeVD-GRN. For a held-out TF $t$
  and gene $g$:

  $$
  \text{score}(t, g) = \#\{\, t' \in \text{training TFs} : t' \to g \text{ in the ground truth} \,\}
  $$

  Every held-out TF gets the same gene ranking. Uses labels (of other TFs), no cell data. Averaged over 5 fold
  seeds.
  **Toy:** training TFs T2, T3, T4 target $\{g1, g2\}$, $\{g1\}$, $\{g1, g3\}$: in-degree g1 = 3, g2 = 1, g3 = 1,
  g4 = 0, so the held-out T1 is predicted g1 first.
- <a id="m-bear-pearson"></a>**Pearson**: $\lvert \text{corr}(x_t, x_g) \rvert$ across cells (CP10k + log1p). No labels.
- <a id="m-bear-grnboost2"></a>**GRNBoost2**: importance of the TF for the gene. No labels.
- <a id="m-bear-coverage"></a>**All measured pairs** ("coverage"): every TF $\times$ measured-gene pair scores 1. No ranking; AUROC $= 0.5$ by
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
