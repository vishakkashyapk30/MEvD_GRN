# Preparation guide for the 18–24 October deck

This guide goes through `presentation_18th_oct_24th_oct.md` slide by slide. For each slide it gives the maths, explains the idea in plain words, and lists questions your advisor may ask, with answers. Open it in the VS Code markdown preview (Ctrl+Shift+V) so the equations render.

---

## Part 0. Background you need for every slide

### 0.1 The task

Each cell type gives one graph. The nodes are transcription factors (TFs), and an undirected edge $u - v$ means "u and v are linked" in the answer key. We write the answer key as a symmetric 0/1 matrix $Y$, with $Y_{uv} = 1$ for a link.

A model gives every pair a score $s(u,v)$, with $s(u,v) = s(v,u)$. A higher score means "more likely linked".

### 0.2 Where the answer key comes from

This is the most important fact in the deck. Almost every result follows from it.

$$Y_{uv} = 1 \iff \big(\text{motif of } u \text{ in an open promoter of } v\big) \ \text{or}\ \big(\text{motif of } v \text{ in an open promoter of } u\big)$$

- "Motif" means a FIMO hit for u's HOCOMOCO v11 motif with $p \le 10^{-6}$.
- "Open promoter" means a peak within ±2 kb of v's start site that is open in more than 10% of the cells.

So a link is a fact about **DNA sequence and open chromatin**. It is **not** a measurement of co-expression or of physical binding.

**Why this creates hubs.** Whether $u - v$ is a link depends mostly on two per-TF properties:

- does u have a common motif (many hits genome-wide)?
- is v's promoter open?

A TF with an open promoter, or a promiscuous motif, collects many links. To a first approximation,

$$P(Y_{uv}=1) \approx g(a_u, a_v)$$

where $a_u$ is a per-TF "propensity". This is the structure of a **configuration model** in network science. In such a graph, the degree of a node almost determines its link probability. That is why counting links (slide 5) works so well.

### 0.3 The splits

The official scMultiomeGRN splits use seed 666 and 10 folds. In each fold:

- **Test:** 1/10 of the links, plus the **same number** of random non-links. The test set is therefore balanced, 50% positives.
- **Validation:** 20% of the remaining links, which is 18% of all links, plus non-links.
- **Training:** everything else, 72% of the links. The training file also lists every remaining non-link (`rest_all_train`).

The split holds out **links, not TFs**. Every TF appears in training with most of its links. This is called the **transductive** setting.

### 0.4 The metrics

**AUROC.** Pick one random test link and one random test non-link. AUROC is the probability that the link gets the higher score:

$$\text{AUROC} = P\big(s(\text{link}) > s(\text{non-link})\big)$$

- Random guessing gives 0.5.
- A perfect model gives 1.
- A model that ranks pairs **backwards** gives a value below 0.5. This happens with the Jaccard index on slide 5.

**AUPR.** Sort the pairs by score and move a threshold down the list. At each threshold:

$$\text{precision} = \frac{TP}{TP+FP}, \qquad \text{recall} = \frac{TP}{TP+FN}$$

AUPR is the area under the precision–recall curve. The official code computes it with the trapezoid rule, `auc(recall, precision)`, which differs slightly from "average precision".

**Key fact:** for a random score, precision equals the fraction of positives at every threshold. Our test set is 50% positives, so **random AUPR = 0.5**. On a balanced test, AUPR and AUROC behave alike. That is why the two columns move together in every table.

**How the numbers are averaged.** Each value is averaged over the 10 folds of a cell type, then over the 12 cell types. MeVD-GRN results are also averaged over its 5 seeds.

### 0.5 Significance test (the p-values in the full tables)

To compare method X with a reference R, take the paired differences over every (cell type, fold):

$$d_i = \text{AUPR}_X(i) - \text{AUPR}_R(i), \qquad i = 1,\dots,120$$

The **Wilcoxon signed-rank test** asks whether these differences are centred on zero:

1. Rank the $|d_i|$ from smallest to largest.
2. Add up the ranks of the positive $d_i$.
3. Compare that sum with what you would expect if positive and negative differences were equally likely.

It does not assume the differences are normally distributed. "Paired" means each fold of X is compared with the same fold of R, which removes the variation between folds and cell types.

| Comparison (test AUPR) | Difference | p |
|---|---|---|
| MeVD-GRN vs scMultiomeGRN | +0.017 | 4e-15 |
| Logistic regression on link counts vs scMultiomeGRN | +0.014 | 2e-13 |
| MeVD-GRN vs logistic regression on link counts | +0.003 | 9e-5 |

### 0.6 Notation for slides 5 and 7

$A$ is the adjacency matrix of the **training links only**:

- $A_{uv} = 1$ if $u - v$ is a training link;
- $A$ is symmetric with a zero diagonal.

From $A$:

- The **degree** of u is $d_u = \sum_v A_{uv}$, the number of training links of u.
- $N(u) = \{w : A_{uw}=1\}$ is the set of neighbours of u.
- $D = \text{diag}(d_1, \dots, d_n)$ is the diagonal matrix of degrees.

---

## Slide 1. "This week in one slide"

This slide is the summary. You need the one-line reason behind each bullet.

| Bullet | The reason in one line |
|---|---|
| Statistics from the cell data are at chance | The answer key is a motif rule, not co-expression (0.2) |
| Counting links beats scMultiomeGRN | The answer key is hub-structured (0.2), and the split holds out links, not TFs (0.3) |
| scMultiomeGRN depends on its known-link graph | Its message passing runs over training links, which carry the degree information |
| Degree-matched test | It removes the advantage of knowing the hubs |
| MeVD-GRN keeps 0.878 | It ranks hard non-links (pairs of hubs) lower than the other models do |

**Likely questions**

- *"38 models: list them."*
  - 10 counting rules: degree product, degree sum, common neighbours, Jaccard, Adamic–Adar, resource allocation, Katz, random walk with restart, low-rank SVD, spectral reconstruction.
  - 12 label-free scores: Pearson, Spearman, partial correlation, mutual information, CLR, GRNBoost2, GENIE3, ridge, lasso, MAESTRO correlation, GRNBoost2-profile cosine, Geneformer cosine.
  - 10 supervised feature models: logistic regression on MAESTRO, GRNBoost2, Geneformer, histograms and promoter openness; side by side; multiplied per TF; within-pair product; linear regression; random forest.
  - 6 models on link counts plus features.
- *"Is using the training links cheating?"* No. Test links never enter the training graph, and scMultiomeGRN itself passes messages over exactly these training links. It is the normal transductive link-prediction setup. What it exposes is a **weakness of the benchmark**, not leakage.

---

## Slide 2. "What we tried and how we scored it"

This slide is about method families. The maths is on the following slides, and the degree-matched test is on slide 7.

**Likely questions**

- *"How did you choose the regularisation without the test set?"*
  - Each logistic regression tries $C \in \{10^{-3}, 10^{-2}, 10^{-1}, 1, 10\}$, where $C = 1/\lambda$.
  - It picks the value by 3-fold cross-validation **inside the training pairs**, scored by average precision.
  - The test fold is never touched.
- *"How do you know none of these leaked?"*
  - Each script asserts that no test or validation pair appears among the training pairs.
  - Counting rules are computed on $A$, which contains training links only.
  - Label-free scores never see any label at all.

---

## Slide 3. Statistical scores from the cell data

All of these turn the cell data into a TF × TF score matrix **without any labels**. The input is TF expression across cells, $x_u \in \mathbb{R}^{C}$ (C cells), after $\log(1+\text{count})$.

### Pearson correlation

$$r_{uv} = \frac{\sum_c (x_{uc}-\bar x_u)(x_{vc}-\bar x_v)}{\sqrt{\sum_c (x_{uc}-\bar x_u)^2}\sqrt{\sum_c (x_{vc}-\bar x_v)^2}}, \qquad s(u,v) = |r_{uv}|$$

This measures a straight-line relationship, and $r \in [-1, 1]$. We use $|r|$ because a link could go with either a positive or a negative relationship.

### Spearman correlation

Spearman is Pearson applied to the **ranks** of the values instead of the values themselves. It captures any relationship that always goes in the same direction, even if not straight, and is robust to outliers.

### Partial correlation

Two TFs can correlate only because a third TF drives both. Partial correlation removes the effect of all other TFs.

1. Estimate the covariance matrix $\Sigma$ of the TFs.
2. Invert it to get the **precision matrix** $\Theta = \Sigma^{-1}$.
3. Compute

$$\rho_{uv\cdot\text{rest}} = -\frac{\Theta_{uv}}{\sqrt{\Theta_{uu}\Theta_{vv}}}$$

With more TFs than reliable samples, $\Sigma$ is badly estimated. So we use **Ledoit–Wolf shrinkage**, $\hat\Sigma = (1-\delta)S + \delta\,\mu I$. It pulls the sample covariance S toward a scaled identity, choosing $\delta$ from the data. This is the idea behind PPCOR in BEELINE.

### Mutual information (MI)

MI measures **any** dependence, not only linear:

$$I(u;v) = \sum_{a,b} p(a,b)\,\log\frac{p(a,b)}{p(a)\,p(b)}$$

- We cut each TF's values into 8 equal-width bins.
- $p(a,b)$ is the fraction of cells with u in bin a and v in bin b.
- $I = 0$ exactly when u and v are independent.

### CLR (context likelihood of relatedness)

Some TFs have high MI with everything, for example highly variable genes. CLR asks whether the u–v value is unusual **for u** and **for v**:

$$z_u(v) = \max\Big(0, \frac{I_{uv} - \mu_u}{\sigma_u}\Big), \qquad \text{CLR}_{uv} = \sqrt{z_u(v)^2 + z_v(u)^2}$$

Here $\mu_u$ and $\sigma_u$ are the mean and standard deviation of u's MI with all TFs.

### GENIE3 (random forests)

For each TF v, train a random forest to predict v's expression from all the other TFs:

$$x_v \approx f_v(x_{-v})$$

The score of u for v is u's **feature importance** in $f_v$, which is how much splitting on u reduces the prediction error.

- We symmetrise with $W + W^\top$.
- We use 300 trees, on 3,000 randomly sampled cells.

### GRNBoost2

GRNBoost2 is the same idea, but with gradient boosting (trees added one at a time, each fixing the previous errors) instead of a random forest. We took scMultiomeGRN's own GRNBoost2 matrix (TF × all genes), read off the TF → TF entries, and used $\max(u\to v,\ v\to u)$.

### Ridge and lasso regression

For each TF v, fit a linear model on standardised data and use $|w_u|$ as the score.

**Ridge** (L2 penalty, $\lambda = 1$):
$$\min_w \ \lVert x_v - X_{-v} w\rVert_2^2 + \lambda\lVert w\rVert_2^2$$

**Lasso** (L1 penalty, $\alpha = 0.01$):
$$\min_w \ \tfrac{1}{2n}\lVert x_v - X_{-v} w\rVert_2^2 + \alpha\lVert w\rVert_1$$

The L1 penalty sets most coefficients to exactly zero, which gives a sparse network. This is the idea behind TIGRESS and the Inferelator.

### MAESTRO correlation and Geneformer cosine

**MAESTRO correlation** is Pearson correlation applied to the TFs' ATAC "regulatory potential" across cells.

**Geneformer cosine** compares the two TFs' pretrained embeddings $e_u$ and $e_v$:

$$\cos(e_u, e_v) = \frac{e_u \cdot e_v}{\lVert e_u\rVert\,\lVert e_v\rVert}$$

### Why they all land at 0.50–0.53

- They all measure **whether two TFs vary together across cells**.
- The answer key (0.2) asks **whether one TF's motif sits in the other's open promoter**.
- These are different biological questions, so the scores carry almost no information about the key.
- This is not a bug in the methods. They answer a different question.

**Likely questions**

- *"The literature says simple statistics beat deep learning. Why not here?"*
  - Those papers (geneRNIB; Ahlmann-Eltze et al. 2025) evaluate against answer keys from perturbations or ChIP-seq, where co-expression is part of the signal.
  - Here the answer key is a motif scan, so statistics on the cell data cannot recover it.
  - The literature's point still holds in a different form: the **simplest supervised** model, counting links, beats the deep model (slide 5).
- *"Some values are below 0.5. Is that a bug?"* No.
  - Below 0.5 means a slight negative relationship. For example, Pearson AUROC is 0.39 on Myeloid cells, so co-expressed TF pairs are slightly *less* often linked there.
  - Averaged over the 12 cell types, they come out at about 0.5.
- *"Is 0.53 significantly above chance?"* It may be statistically detectable, but it is practically useless. Every one of these scores is far below any model that learns from the training links.
- *"Why does GENIE3 not beat correlation?"* GENIE3 is better at finding *which TF predicts which gene's expression*. That is still a co-expression question, so it inherits the same mismatch with the motif answer key. The scMultiomeGRN paper also reports GENIE3 last among its baselines (AUROC 0.54).

---

## Slide 4. Simple supervised models on scMultiomeGRN's inputs

### Step 1: one vector per TF

For each input, we start from a TF × feature matrix:

| Input | Size | What a row holds |
|---|---|---|
| MAESTRO | TF × cells | the TF's accessibility "regulatory potential" in each cell |
| GRNBoost2 | TF × genes | the TF's importance for each gene |
| Geneformer | TF × 768 | the TF's pretrained embedding |
| Promoter openness | TF × 4 | summaries of how open the TF's promoter is |

These matrices are too wide to use directly, so we compress each one with a **truncated SVD**. This uses no labels.

$$M - \bar M = U \Sigma V^\top, \qquad x_u = (U_{u,1:32}\,\Sigma_{1:32}) \ \text{(standardised)}$$

### Step 2: pair features that ignore order

The task is undirected, so the pair feature must not change when u and v swap:

$$\varphi(u,v) = \big[\, x_u \odot x_v ,\ \ |x_u - x_v| \,\big]$$

- $\odot$ is element-wise multiplication. It is large when both TFs are high on the same dimension, like a similarity.
- $|x_u - x_v|$ measures how different the two TFs are.
- Both terms stay the same when u and v are swapped.

### Step 3: logistic regression

$$P(\text{link} \mid u,v) = \sigma\big(w^\top \varphi(u,v) + b\big), \qquad \sigma(z) = \frac{1}{1+e^{-z}}$$

It is fitted by minimising the penalised log-loss over the training pairs (labels $y_i \in \{0,1\}$):

$$\min_{w,b}\ -\sum_i \big[y_i\log p_i + (1-y_i)\log(1-p_i)\big] + \frac{1}{2C}\lVert w\rVert^2$$

### The variants your advisor asked for

| Name in the deck | Feature used | Test AUPR |
|---|---|---|
| Side by side ("given independently") | $[\varphi_M(u,v), \varphi_G(u,v)]$, the MAESTRO and GRNBoost2 blocks concatenated | **0.788** |
| Multiplied per TF ("cross-modal product") | first $z_u = m_u \odot g_u$ per TF, then $\varphi$ on z | 0.721 |
| Within-pair product only | $[m_u\odot m_v,\ g_u\odot g_v]$, with no difference terms | 0.523 |
| Linear regression | ridge regression on the 0/1 label, used as a score | 0.788 |
| Random forest | 300 trees on the side-by-side features | 0.836 |

**Why side by side beats multiplied.**

- Concatenation keeps every piece of information, so the model can learn its own weights.
- Multiplying $m_u \odot g_u$ before learning throws information away. Where either input is near zero, the product is near zero whatever the other one says.
- Logistic regression can only add terms. It cannot undo a product that was taken before it saw the data.

**Why the within-pair product alone is near chance (0.523).**

- $x_u \odot x_v$ only measures **similarity** between two TFs.
- Similarity is close to what the slide-3 scores measure, and that carries no signal here.
- The difference term $|x_u - x_v|$ and the per-TF scale are what let the model learn which individual TFs are "link-prone". In effect, they let it learn degree.

**Why a random forest helps (0.836).** Trees can learn interactions between features and thresholds, for example "u's accessibility is above some value AND v's is below some value". A linear model cannot.

**Likely questions**

- *"So does the cell data have signal or not?"* It has signal once a model can learn **which TFs tend to have links**. The features act as a proxy for a TF's propensity to be linked. They do not reveal the specific pairs. That is why they are useful here (0.70–0.84) but useless as label-free similarity (slide 3).
- *"Why use SVD to 32 dimensions?"*
  - With thousands of feature columns and only a few hundred training links, the model would overfit.
  - The SVD uses no labels, so it adds no leakage.
- *"Why use logistic regression and not linear regression on 0/1?"*
  - Logistic regression models probabilities properly.
  - Linear regression on the label (a "linear probability model") gives an almost identical ranking, 0.788 against 0.788.
  - AUROC and AUPR only depend on the ranking, so the two score the same.

---

## Slide 5. Counting the known links (link-prediction heuristics)

All of these are computed from $A$ (training links only). These are the standard heuristics from the link-prediction literature.

### Degree product and degree sum

$$s_{\text{prod}}(u,v) = (d_u+1)(d_v+1), \qquad s_{\text{sum}}(u,v) = d_u + d_v$$

The degree product is called **preferential attachment**: well-connected nodes keep getting more links. It works because of the hub structure in 0.2.

### Logistic regression on the link counts

With $\ell_u = \log(1+d_u)$:

$$\varphi(u,v) = \big[\ell_u + \ell_v,\ |\ell_u-\ell_v|,\ \ell_u\ell_v\big], \qquad P(\text{link}) = \sigma(w^\top\varphi + b)$$

This is a learned version of "product or sum": it decides how much each combination should matter. It reaches 0.931 AUPR, beating scMultiomeGRN (0.918).

**An important technical detail: computing training-pair features with the pair's own link held out.**

- A **training** link $u - v$ adds 1 to both $d_u$ and $d_v$. A **test** link is never in $A$.
- If we computed the training features on the full $A$, the model would learn "a link exists when the counts include this edge", a pattern that never occurs at test time.
- So we split the training links into 5 chunks. For pairs in chunk k, we compute features on $A$ with chunk k's links removed.
- This matches how test pairs look.
- Without this fix, the spectral-embedding model scored an AUROC of 0.06, almost perfectly backwards. That is the tell-tale sign of this mismatch.

### Neighbourhood rules

**Common neighbours:**
$$s_{\text{CN}}(u,v) = |N(u)\cap N(v)| = (A^2)_{uv}$$

**Jaccard index:**
$$s_{\text{Jac}}(u,v) = \frac{|N(u)\cap N(v)|}{|N(u)\cup N(v)|}$$

**Adamic–Adar and resource allocation:**
$$s_{\text{AA}}(u,v) = \sum_{w \in N(u)\cap N(v)} \frac{1}{\log d_w}, \qquad s_{\text{RA}}(u,v) = \sum_{w \in N(u)\cap N(v)} \frac{1}{d_w}$$

Both give more weight to a shared neighbour w that is *not* a hub. Sharing a rare friend says more than sharing a popular one.

**Katz index:** count all paths between u and v, with longer paths down-weighted. $(A^l)_{uv}$ is the number of walks of length l from u to v. We use $\beta = 0.05$ and stop at length 4:

$$s_{\text{Katz}}(u,v) = \sum_{l=1}^{4} \beta^l (A^l)_{uv}$$

**Random walk with restart (personalised PageRank):** a walker starts at u. At each step it moves to a random neighbour, or with probability $\alpha = 0.15$ jumps back to u. $R_{uv}$ is the long-run fraction of time it spends at v. With the transition matrix $P = D^{-1}A$:

$$R = \alpha\,\big(I - (1-\alpha)P\big)^{-1}, \qquad s(u,v) = R_{uv} + R_{vu}$$

**Low-rank reconstructions:**

- **SVD:** $A \approx U_k \Sigma_k V_k^\top$ with $k = 16$. The score is the reconstructed entry.
- **Spectral:** the same, using the top 16 eigenvectors of the normalised matrix $D^{-1/2} A D^{-1/2}$.

### How to explain the results

- **Degree rules win** (0.924–0.931) because the answer key is hub-structured.
- **Random walk with restart scores about the same as degree** (0.924). A random walker spends more time on high-degree nodes, so these scores are mostly degree in disguise.
- **Katz** (0.887) is weaker than degree. Test pairs are never training links, so the length-1 term is zero for all of them, and only paths of length 2–4 count. Those path counts grow with degree, which is why it still does well.
- **Common neighbours and Adamic–Adar are weaker** (0.78–0.81). They only see shared neighbours, and many test pairs have few.
- **The Jaccard index is backwards** (0.410). It divides by the size of the **union** of neighbourhoods. Test links mostly involve hubs, whose union is large, so Jaccard pushes exactly the true links down.
- **Low-rank SVD and spectral reconstruction are weak** (0.72 and 0.66 AUPR). They are above chance, but far below degree. Both test links and test non-links are non-edges of $A$. A low-rank fit of $A$ is built to reproduce the edges it has seen, not to rank pairs it has never seen, and it is dominated by a few large components.

**Likely questions**

- *"Why does adding features to the link counts not help?"*
  - The features mainly tell the model which TFs are link-prone, and the counts already say that directly.
  - Adding about 130 more dimensions just adds variance.
  - Result: between 0.926 and 0.931, no better than the counts alone (0.931).
- *"Then what is the benchmark really measuring?"* Mostly whether you can tell which TFs are hubs. That is why we added the degree-matched test (slide 7).
- *"Is the 'Hidden sampling biases' paper peer-reviewed?"* No, it is a bioRxiv preprint from December 2025. Our own data show the same effect here.

---

## Slide 6. Taking scMultiomeGRN apart

### How scMultiomeGRN works, in equations

**Message-passing graph.** The model passes messages over the training links, plus a self-loop on every node:

$$\hat A = D^{-1/2}(A+I)D^{-1/2}$$

**Edge features.** For every edge of that graph, the model computes a 16×16 joint histogram of the two TFs' $\log_{10}$ expression across cells:

$$E_{uv} = \frac{\log_{10}(H_{uv}/C + 10^{-4}) + 4}{4}$$

**GraFrank encoder.** For each modality k (scRNA features, ATAC features), there are two attention-based message-passing layers. Each layer:

- computes a message from neighbour j and the edge features: $m_{ij} = W[h_j \,\|\, \text{conv}(E_{ij})]$;
- weights the messages by attention, $\alpha_{ij} = \text{softmax}_j(\text{LeakyReLU}(a^\top m_{ij} + c_i))$;
- updates node i by combining the weighted sum with its own projection.

**Fusion and output.**

- **Cross-modality attention** merges the two modality embeddings into one vector per TF, through a softmax weighting across modalities.
- Then comes one graph convolution, $Z = \hat A\, H\, W$.
- The pair decoder scores $(z_u, z_v)$ through a small network on the two TF vectors.

**Loss.** A weighted binary cross-entropy, where $w_+$ (`pos_weight`) is the ratio of non-links to links:

$$\mathcal L = -\sum_{(u,v)} \big[w_+\, y\log\hat y + (1-y)\log(1-\hat y)\big]$$

### What each ablation switches off

| Ablation | What changes | Test AUPR |
|---|---|---|
| Unchanged | nothing (checked to reproduce the original code exactly) | 0.918 |
| Without the TF-TF graph | $A = 0$, so $\hat A = I$: self-loops only, no message passing between TFs | 0.826 (lung 0.906, blood 0.585) |
| STRING instead of the known links | $A$ = STRING v12 protein interactions, combined score ≥ 700 | 0.781 |
| Without pair histograms | the edge-feature channel removed completely | 0.920 |
| ATAC only, RNA only | one modality, GraFrank kept | 0.916, 0.916 |
| With Geneformer | Geneformer added as a third modality | 0.917 |

In every ablation, the **training labels** stay the same. Only the inputs or the graph change.

### How to explain the results

- **Without the graph, the model loses its degree information.** Message passing over training links tells each TF how many links it has. Without it, the model must guess link propensity from features alone.
- **Lung keeps 0.906 without the graph, but blood falls to 0.585.**
  - In the lung dataset the "RNA" part is gene-body **accessibility**, so both feature types come from ATAC and are closely tied to promoter openness, a core part of the answer key.
  - In blood the RNA comes from a separate experiment, so the features are weaker.
  - Present this as a likely explanation, not a proven one.
- **STRING is worse than no graph.**
  - STRING contains only 2–6% of the answer-key links.
  - Passing messages over unrelated edges mixes the representations of TFs that are not linked, which blurs exactly the degree signal the model needs.
- **Geneformer does not help scMultiomeGRN**, even though it is MeVD-GRN's strongest single part. A possible reason is that scMultiomeGRN's score is already set by the graph, so one more node feature adds nothing. This is a hypothesis, not tested.
- **We checked that Geneformer was really used.** The per-fold scores change slightly (for example Myeloid fold 1: 0.9616 → 0.9642). The input is used, it just doesn't help.

**Likely questions**

- *"Why not use their own variant class for single-modality runs?"* The authors' `scMultiomeGRNVariant` drops the whole GraFrank encoder when only one modality is used. That would remove more than one thing at once, so we kept GraFrank and only removed the modality.
- *"How do you know your patch didn't change the model?"* With every switch at its default, the patched code gives identical metrics to the original on the same seed and epochs, and the full rerun gives 0.918, the same as before.
- *"Why only one seed (666)?"*
  - We agreed on one seed beforehand.
  - Earlier runs showed the official code is nearly deterministic across seeds (Myeloid, 6 seeds, all within 0.001).

---

## Slide 7. The degree-matched test

### The idea

On the official test, non-links are **random pairs**. Most TFs have few links, so a random pair usually involves two low-degree TFs. True links usually involve hubs. A degree rule separates the two easily.

The fix, following the "Hidden sampling biases" preprint (bioRxiv 2025), is to choose test non-links **with the same degrees as the test links**. Then degree alone can no longer separate them.

### The exact rule

1. Take the fold's test links. These are the positives.
2. Describe each pair by its two TFs' log training degrees, $(\ell_{\min}, \ell_{\max}) = (\min(\ell_u,\ell_v), \max(\ell_u,\ell_v))$.
3. For each positive, find the **3 official test non-links** closest to it in L1 distance:
$$\text{dist} = |\ell_{\min} - \ell'_{\min}| + |\ell_{\max} - \ell'_{\max}|$$
4. Pick one of those 3 at random, with replacement and a fixed seed per fold.
5. Score with the same AUROC and AUPR. The test is still balanced, one non-link per link.

**Why only the official test non-links?** The split files cover **every** TF pair (0.3). The test non-links are the only non-links that no model trained on. Using them guarantees that no model has seen the pairs. So the new test is a **re-weighting of the official test set** toward the hard non-links.

### How to explain the results

| Model | Official AUPR | Degree-matched AUPR |
|---|---|---|
| Product of the link counts | 0.924 | 0.724 |
| Best simple model (counts and neighbour scores) | 0.932 | 0.759 |
| scMultiomeGRN | 0.918 | 0.705 |
| MeVD-GRN | 0.935 | **0.878** |

- **The counting rule does not fall all the way to 0.5.** In these dense graphs, 14–37% of all pairs are links, so a non-link between two hubs is rare. The nearest available non-link often has lower degrees than the link it is matched to, and the degree rule can still use that gap.
- So read 0.72 as the **floor set by degree** on this test. A model shows signal beyond hubness only by being clearly above it.
- **scMultiomeGRN and every model that counts links end up near that floor** (0.70–0.76).
- **MeVD-GRN stays far above it** (0.878): it ranks hub–hub non-links low.

**Likely questions**

- *"Could MeVD-GRN have seen these non-links in training?"* No.
  - MeVD-GRN draws its training negatives only from the fold's training non-links.
  - The degree-matched non-links are official test non-links, never in any training file.
  - Every MeVD-GRN result file passed the zero-leak check when the tables were compiled.
- *"Why would MeVD-GRN do better?"*
  - It sees per-TF information that the counting rules don't have: the Geneformer embedding and accessibility-based candidate graphs.
  - Possibly it learns which specific *pairs* go together, not only which TFs are hubs.
  - Which part carries the gap is not tested yet. That is the next step (slide 8). Be honest that this is open.
- *"Isn't this test small?"* Yes, and that is the honest limit.
  - The non-links come only from the official test set, whose size equals the number of test links (27 per fold in Myeloid cells, about 500 in Vascular endothelial cells).
  - Sampling with replacement means a few hard non-links can appear several times.
  - That is why we propose a held-out-TF split next.
- *"Why match on degree and not on something else?"* Degree is the specific shortcut we showed (slide 5). Matching on it tests whether anything is left once that shortcut is removed.

---

## Slide 8. What we take from this and what comes next

**The four conclusions, each tied to its evidence:**

1. **The benchmark mostly rewards knowing hub TFs.** Two numbers per pair, the link counts, beat the published model (0.931 vs 0.918, p = 2e-13).
2. **Label-free network inference does not fit a motif-based answer key.** All 12 scores are at 0.50–0.53.
3. **scMultiomeGRN's extra inputs add nothing on top of its graph.** All feature ablations change AUPR by 0.002 or less, with p ≥ 0.29.
4. **MeVD-GRN is the only model that holds up once hubs are controlled for.** It scores 0.878, against 0.705 for scMultiomeGRN.

**Next steps, and why:**

- **A held-out-TF split.** Hide all links of some TFs. No model can count a hidden TF's links, so the degree shortcut disappears completely, and the gap from slide 7 can be confirmed or ruled out.
- **Degree-matched scores for the MeVD-GRN ablations.** This tells us which part carries the gap: Geneformer, ATAC, the graph network, or the candidate graphs.

**Likely questions**

- *"If the benchmark is flawed, why use it?"*
  - It is the published head-to-head setting, so we must report it.
  - Showing *where* it is flawed, and that MeVD-GRN still holds up on the harder version, is itself a contribution.
- *"What would convince you MeVD-GRN learns real biology?"*
  - A held-out-TF split where it still beats the degree floor.
  - An answer key that is not built from motifs, for example ChIP-seq co-binding.

---

## Cheat sheet: numbers to remember

| What | Number |
|---|---|
| Random AUROC or AUPR on this test | 0.5 |
| Best label-free statistic | 0.529 (Geneformer cosine / MAESTRO correlation) |
| GRNBoost2 / GENIE3 | 0.508 / 0.504 |
| Logistic regression, MAESTRO + GRNBoost2 side by side | 0.788 |
| Random forest, same features | 0.836 |
| Degree product / logistic regression on counts | 0.924 / 0.931 |
| scMultiomeGRN | 0.918 (degree-matched 0.705) |
| MeVD-GRN | 0.935 (degree-matched 0.878) |
| scMultiomeGRN without the graph / with STRING | 0.826 / 0.781 |
| STRING links that are in the answer key | 2–6% of the key |
| Link density of the graphs | 14–37% of all TF pairs |
| Split | test 10% of links plus equal non-links, validation 18%, training 72%, 10 folds, seed 666 |
