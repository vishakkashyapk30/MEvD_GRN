# MEvD-GRN: Multi-Evidence Distillation for Gene Regulatory Network Inference
## Complete Theory & Implementation Plan

---
**Target agent**: Claude in VS Code (or equivalent)
**Format**: Theory first, then numbered implementation steps per module. No code in this document.
**Goal**: Implement a novel, publishable GRN inference method using SC-MO-GRN-DB as the data foundation.
---

## Part 0 — Scientific Framing

### 0.1 What This Method Is

MEvD-GRN (Multi-Evidence Distillation for Gene Regulatory Network Inference) is a compact graph neural network trained to predict transcription factor (TF) to target gene regulatory edges from paired scRNA-seq and scATAC-seq data. Its distinguishing feature is a three-stage curriculum training scheme that exploits SC-MO-GRN-DB's tiered reference networks — localization, perturbation, and dual-evidence — as progressively stricter supervision signals.

### 0.2 Why This is Novel

SC-MO-GRN-DB (Valensi et al., iScience 2026, DOI: 10.1016/j.isci.2026.115323) explicitly constructs three reference network types representing a biological confidence hierarchy:

- **Localization networks**: Built from ChIP-seq and ChIP-chip data. They record where a TF physically binds in the genome (proximity to target TSS). These networks have high recall but many false positives because TF binding does not always mean transcriptional regulation.
- **Perturbation networks**: Built from knockdown, knockout, and CRISPR-based TF perturbation experiments. If removing TF_i causes differential expression of gene_j, that is functional causal evidence of regulation. These networks have lower recall but much higher precision.
- **Dual-evidence networks**: The intersection of localization and perturbation networks for the same cell type. A TF-target pair must have both physical binding evidence AND functional perturbation evidence. These are the highest-confidence regulatory interactions in the database.

The key insight is that this hierarchy is a natural curriculum: localization is easy to learn (many examples, noisy signal), perturbation is harder (fewer examples, stronger biological signal), and dual-evidence is the hardest and most meaningful supervision (fewest examples, highest confidence). No existing GRN method uses this three-tier evidence structure as a training curriculum. That is the publishable contribution.

### 0.3 The Scientific Claim to Be Demonstrated

After training, the model trained via evidence-aware curriculum (localization → perturbation → dual-evidence) should achieve strictly higher AUPR and early precision on the dual-evidence test set than the same architecture trained on any single tier or on all tiers simultaneously without curriculum ordering. Secondary claims: it is computationally efficient (trains in hours, not days), and it generalizes across cell types.

---

## Part 1 — Database & Data Reference

### 1.1 SC-MO-GRN-DB

The database is publicly available at: https://scmogrndb.psu.edu  
Source repository: https://github.com/UzunLab/SC-MO-GRN-DB

It contains:
- 30 reference networks spanning 16 tissue types (human and mouse)
- 1,616 transcription factors, 63,512 target genes, 22,389,021 regulatory edges total
- 31 single-cell datasets, 9 tissues, 6 molecular modalities, 2,322,285 cells total
- Modalities: scRNA-seq, scATAC-seq, scChIP-seq, scDNA methylation, scHi-C, scCRISPR-seq
- For each cell type, datasets are matched to reference networks of the same cell type

### 1.2 Cell Types to Use

For Phase 1 (required, enough data for all three evidence tiers):
- **K562** (human erythroleukemia cell line): 3 reference networks, 8 datasets — this is the primary training cell type because it has all three evidence tiers available
- **ESC** (embryonic stem cells, human + mouse): 6 reference networks, 8 datasets — second training cell type and also used for cross-cell-type transfer evaluation

For Phase 2 (if time permits):
- HSC (hematopoietic stem cells): 2 networks, 4 datasets
- GM12878 (B-lymphoblastoid cell line): 1 network, 3 datasets

### 1.3 Data to Download

From the SC-MO-GRN-DB web interface, download the following:

**Reference networks** (two-column TSV files, format: TF_name TAB Target_name):
- K562 Localization network
- K562 Perturbation network
- K562 Dual-evidence network
- ESC Localization network
- ESC Perturbation network
- ESC Dual-evidence network

**Single-cell datasets** (h5ad files):
- K562 scRNA-seq dataset (prefer the one with most cells, around 8 available)
- K562 scATAC-seq dataset
- ESC scRNA-seq dataset
- ESC scATAC-seq dataset

If jointly-assayed (RNA+ATAC paired from the same cells) datasets are available for K562 or ESC in the database, use those in preference to separate modality files. Jointly-assayed data eliminates the need for modality alignment and provides paired observations per cell.

### 1.4 File Organization

Organize downloaded data in this directory structure before starting implementation:
```
data/raw/networks/   — all downloaded TSV reference network files
data/raw/single_cell/ — all downloaded h5ad single-cell files
```

---

## Part 2 — Biological & Mathematical Problem Formulation

### 2.1 What Gene Regulatory Network Inference Is

A GRN is a directed graph where nodes are genes and directed edges represent regulatory relationships: a directed edge from gene A to gene B means that gene A (acting as a transcription factor) regulates the expression of gene B (as a target). The key computational task is to infer which edges exist given observations of gene expression (and optionally, chromatin accessibility data) across many individual cells.

This is a supervised learning problem: we have positive edges from the reference networks, and we treat all other possible TF-to-gene pairs as potential negatives (the exact set of negatives depends on the prior graph construction, described in Part 4).

### 2.2 Why This is a Bipartite Directed Link Prediction Problem

The TF-to-Target relationship is inherently bipartite:
- **Node type 1 (TF nodes)**: A constrained set of known transcription factors. These are genes that encode proteins capable of binding DNA and regulating transcription. The SC-MO-GRN-DB paper identifies 1,616 such TFs.
- **Node type 2 (Target Gene nodes)**: All genes in the expressed gene universe, including TF genes themselves (TFs can regulate each other).
- **Directed edges**: Regulatory edges go FROM TF nodes TO target gene nodes. The direction matters biologically: TF_A may regulate TF_B, but TF_B may not regulate TF_A.

The task is therefore: given the node feature representations of all genes (derived from scRNA and scATAC data), predict which directed edges (TF → Target) exist in the GRN. This is formally equivalent to link prediction in a directed bipartite graph.

Bipartite structure matters for the GNN design because standard GNNs assume homogeneous nodes. In a bipartite graph, messages should flow FROM targets TO TFs and FROM TFs TO targets via separate update rules, because the two node types carry different biological semantics.

### 2.3 Mathematical Notation

Let:
- G = (V, E_prior) be the prior candidate graph (defined in Part 4)
- V = V_TF ∪ V_Gene with |V_TF| = N_TF, |V_Gene| = N_G (TF nodes are a subset of all gene nodes)
- X_RNA ∈ ℝ^(N_G × d_rna) = gene-level RNA features
- X_ATAC ∈ ℝ^(N_G × d_atac) = gene-level ATAC features
- E_loc, E_pert, E_dual ⊂ V_TF × V_Gene = positive edge sets from the three evidence tiers
- The goal is to learn a function f(X_RNA, X_ATAC, G) → S ∈ ℝ^(N_TF × N_G) where S_ij is the predicted probability of a regulatory edge from TF_i to Target_j

### 2.4 Why scRNA-seq and scATAC-seq are Complementary

scRNA-seq measures mRNA abundance across cells. It captures which genes are being expressed, at what level, and how variable that expression is across cells. From the GRN perspective, the co-expression pattern of a TF and its targets provides a correlational signal for regulation.

scATAC-seq measures chromatin accessibility (open vs. closed chromatin) at specific genomic loci across cells. Biologically, transcription can only occur when the chromatin around a gene's promoter and enhancers is "open" (accessible to the transcriptional machinery). For TF regulation specifically:
- If a TF's binding motif overlaps an ATAC-seq peak near a target gene's TSS, there is physical evidence that the TF can access and bind that region
- The accessibility level of a target gene's promoter region predicts how amenable the gene is to regulation by TFs that can bind nearby

The two modalities therefore provide orthogonal but complementary signals: RNA tells us about current expression states, and ATAC tells us about the regulatory landscape (what is structurally possible to regulate). Methods that integrate both modalities outperform single-modality methods, as shown by scMultiomeGRN (Xu et al., Nucleic Acids Research 2025).

---

## Part 3 — Data Preprocessing Theory

### 3.1 scRNA-seq Preprocessing Theory

Raw scRNA-seq count data is count-valued, highly sparse (due to biological and technical dropout), and varies in total counts across cells (library size variation). Standard preprocessing applies the following transformations before any model input:

**Quality control**: Filter out cells with very few detected genes (likely empty droplets or dead cells) and genes detected in very few cells (likely noise). SC-MO-GRN-DB data has already been filtered to wild-type/control cells only, so minimal additional QC is needed.

**Library size normalization**: Divide each cell's counts by the total counts in that cell, then multiply by a scaling constant (conventionally 10,000 = "counts per 10k"). This removes the effect of different total sequencing depths between cells, making cells comparable.

**Log transformation**: Apply log(x + 1) to the normalized counts. This compresses the dynamic range (gene expression spans several orders of magnitude), makes the data more normally distributed, and stabilizes variance for downstream modeling.

**Feature aggregation to gene level**: The MEvD-GRN model operates at the gene level, not the cell level. The GNN nodes are genes, not cells. Therefore, the preprocessed cell-by-gene matrix must be reduced to gene-level statistics. Compute, for each gene: (1) the mean log-normalized expression across all cells (captures average activity), (2) the variance of log-normalized expression across all cells (captures regulatory variability / cell-to-cell heterogeneity), and (3) the detection rate — the fraction of cells in which the gene is detected (non-zero) — which distinguishes broadly-expressed housekeeping genes from sparsely-expressed cell-state genes. These three numbers form the per-gene RNA feature vector of dimension 3.

**Why mean + variance + detection**: The mean captures the baseline expression level of a gene and its average regulatory activity. The variance captures how dynamically regulated the gene is across cells — highly variable genes are often cell-identity or state-determining genes, which are more likely to be key regulatory nodes. The detection rate separates ubiquitously-on genes from bursty/rare ones. Together they give the model richer per-node information than mean/variance alone.

**Co-expression signatures (preserving the pairwise signal)**: Per-gene summary statistics discard the single most important correlational signal in GRN inference — whether a TF and a candidate target co-vary across cells. Classical methods (GENIE3, PIDC) are built almost entirely on this signal. To retain it without carrying the full cell matrix into the model, compute a per-gene co-expression signature vector s_g ∈ ℝ^d (d ≈ 50) via a truncated SVD of the column-standardized, mean-centered expression matrix W (where W = (X − column_mean) / (column_std · √C), so that WᵀW is exactly the gene-gene Pearson correlation matrix). With W = U S Vᵀ, set s_g = (V S)_g. Then the dot product s_i · s_j ≈ (V S² Vᵀ)_ij ≈ (WᵀW)_ij = corr(i, j). Centering is applied implicitly through a `LinearOperator` so the sparse matrix is never densified. These signatures are used two ways downstream: (a) to build the co-expression message-passing graph (Part 4), and (b) as an explicit edge-level co-expression term in the decoder (Part 5.6).

**Highly variable gene selection**: If the gene universe is very large (>10,000 genes), restrict to the top 3,000 highly variable genes to reduce computation. Highly variable genes (HVGs) are those with high normalized dispersion — they contribute the most biological signal.

### 3.2 scATAC-seq Preprocessing Theory

scATAC-seq produces a peaks-by-cells (or cells-by-peaks) count matrix where each peak represents a genomic interval of open chromatin. This data is even sparser than scRNA-seq (typically >95% zeros).

**Quality filtering**: Remove cells with very low TSS enrichment scores (indicating poor chromatin accessibility signal) and very few detected peaks.

**Normalization**: Apply library-size normalization (counts per 10k) followed by log(x+1), same as RNA. Optionally apply TF-IDF normalization (term frequency-inverse document frequency), which re-weights peaks by how specifically they occur in individual cells vs. across all cells.

**Gene Activity Score construction**: This is the most critical preprocessing step for integrating ATAC with the gene-level GNN. The idea is to convert the peaks-by-cells matrix into a genes-by-cells accessibility matrix by mapping each ATAC peak to the gene(s) it is likely to regulate.

The mapping rule: for each gene, identify all ATAC peaks whose midpoint falls within a genomic window of ±100 kb around the gene's Transcription Start Site (TSS). The 100 kb window is biologically motivated because most cis-regulatory elements (promoters, proximal and distal enhancers) act within this range. For each gene, sum the accessibility values of all peaks within this window. This sum is the Gene Activity Score for that gene in a given cell.

After computing Gene Activity Scores for all genes across all cells, aggregate to gene-level statistics: compute mean and variance of the Gene Activity Score across cells. This gives a 2-dimensional ATAC feature per gene.

**Locus openness (a motif-free accessibility proxy)**: In addition to the activity-score mean/variance, compute a per-gene openness scalar = log1p(number of accessible peaks within ±window of the gene's TSS), normalized by the maximum across genes to lie in [0, 1]. Whereas the Gene Activity Score summarizes how much accessibility signal a locus carries, openness summarizes how many independent cis-regulatory elements are open near the locus — a direct, motif-independent proxy for how regulatable the locus is. Openness is consumed as an explicit edge-level term in the decoder (Part 5.6) and as a filter when constructing the TF-candidate graph (Part 4). Genes with no mapped peaks get openness 0.

**Why this mapping works**: Chromatin accessibility at a gene's locus reflects the likelihood that TF binding sites near the gene are accessible and that the gene can be regulated. High accessibility near a target gene's TSS increases the probability that TFs with binding motifs in those accessible regions can bind and regulate that gene. Crucially, accessibility is a *precondition* for regulation, not an interchangeable substitute for expression — this asymmetry is what motivates the role-aware integration in Part 5.4.

**Genome coordinates**: Use GENCODE v38 annotations for human (hg38) and GENCODE vM25 or equivalent for mouse (mm39). These define TSS positions for each gene. Use the pyranges library for fast interval operations.

### 3.3 Gene Universe Alignment

The RNA and ATAC preprocessing produce two feature matrices indexed by genes. These must be aligned to the same gene list. The gene universe is determined by the intersection of:
1. Genes detected in the scRNA-seq data (after HVG selection)
2. Genes covered by ATAC peaks (after gene activity score computation)
3. Genes appearing in any reference network (TFs or target genes)

Compute the union of TF and Target gene names appearing in any of the three evidence tier networks for the cell type. Intersect with scRNA-expressed genes. This defines the final gene universe of size N_G. All feature matrices are indexed by this gene list. TF nodes are the subset of the gene universe that appear as the source (column 0) in any reference network file, intersected with the known TF list from AnimalTFDB 4.0.

---

## Part 4 — Graph Construction Theory

### 4.1 The Message-Passing Graphs

The GNN performs message passing over graphs that must carry real biological structure. Note a separation of concerns clarified in this revision: these graphs are used **only** for representation learning (message passing). They do **not** define the candidate edge space — the decoder can score any (TF, gene) pair, positives come from the reference networks, and negatives are sampled from the full bipartite space (Part 4.3). Decoupling the message-passing graph from the candidate/label space avoids the accessibility-shortcut artifact documented during implementation.

Two complementary graphs are built:

**(a) Co-expression kNN graph (gene–gene, undirected).** For each gene, connect it to its top-k most co-expressed genes, where co-expression is the signature dot product s_i · s_j ≈ corr(i, j) from Part 3.1. k ≈ 20 by default. Genes that co-vary across cells belong to shared regulatory programs; message passing over this graph lets each node aggregate its co-regulation module. This is the primary source of biological topology and is entirely expression-derived.

**(b) TF-candidate graph (TF → target, directed).** For each TF, connect it to the top-K target genes that are **both** (i) accessible (positive openness, i.e. the locus has open cis-regulatory elements near its TSS) **and** (ii) most co-expressed with that specific TF (highest s_TF · s_target). K = 500 default. This yields a **different neighborhood per TF**, which is the key correction over the earlier design: the previous version ranked genes by accessibility once *globally* and connected every TF to the same top-K genes, so the graph carried no TF-specific information and the GNN had almost nothing discriminative to learn. Intersecting accessibility (can this locus be regulated at all?) with TF-specific co-expression (is this TF plausibly linked to this target?) produces a biologically meaningful, TF-specific candidate neighborhood.

Both graphs **exclude known positive edges** so that validation/test labels cannot leak into node embeddings through message passing. When co-expression signatures are unavailable (degenerate data), the TF-candidate graph falls back to the legacy global-accessibility ranking; this legacy mode is retained only for the graph ablation.

An alternative, more biologically precise construction of the TF-candidate graph is motif-based: for each TF, use its known DNA-binding motif (from JASPAR2024 or AnimalTFDB 4.0) to scan ATAC peaks for predicted binding sites, then link TF_i to gene_j if a peak near gene_j's TSS contains a motif site for TF_i. This is more selective but requires motif-scanning infrastructure; the accessibility × co-expression proxy above is the default.

### 4.2 Evidence Edge Graphs

For each evidence tier, load the positive edges from the SC-MO-GRN-DB reference network files and restrict to edges where both TF and target gene appear in the gene universe. The three resulting positive edge sets (E_loc, E_pert, E_dual) are nested by definition:

- E_dual ⊆ E_pert (dual is the intersection, so every dual edge is also a perturbation edge)
- E_dual ⊆ E_loc (same logic)
- |E_dual| < |E_pert| < |E_loc| (dual is sparsest, localization is densest)

Verify these relationships after loading. This nesting structure is the biological basis for the curriculum: the easiest supervision (E_loc, most edges) comes first, and the hardest (E_dual, fewest edges) comes last.

**Candidate/label space is independent of the message-passing graphs**: Because the decoder can score any (TF, gene) pair, positives are taken directly from the reference networks (restricted to the gene universe) and are **not** filtered to the message-passing graphs. This removes the artificial recall ceiling the earlier "only E_prior edges are scorable" design imposed, and eliminates the accessibility shortcut (if negatives were drawn only from the accessibility-ranked prior, the model could separate positives from negatives by an accessibility artifact rather than by learning regulation). Report the fraction of database edges within the gene universe as the coverage statistic.

### 4.3 Negative Edge Construction

The GRN positive/negative class imbalance is extreme: the vast majority of TF×gene pairs are not regulatory.

**Random negatives**: Sample random TF×gene pairs from the **full bipartite candidate space** that are positive in NO evidence tier (a universal negative pool, consistent across all curriculum stages). Sampling from the full space — rather than restricting to the accessibility-ranked message-passing graph — is essential: prior-restricted negatives let the model separate positives from negatives via an accessibility artifact, badly inflating AUPR. Random negatives are the field-standard choice.

**Hard negatives**: For Stage 2 (perturbation supervision), edges that are in E_loc but NOT in E_pert are biologically meaningful hard negatives: TF physically binds near the target but does not functionally regulate it. The model must learn to distinguish physical binding (easier) from functional regulation (harder). Similarly for Stage 3: edges in E_pert but not in E_dual are hard negatives. Using hard negatives (from the lower-confidence tier that are NOT in the current tier) alongside random negatives improves the model's discrimination of the confidence hierarchy.

**Class weighting**: During training, upweight positive edges by the ratio of negatives to positives (e.g., if 5 negatives per positive are sampled, use positive weight = 5 in the loss function). This corrects for the still-significant imbalance in the sampled mini-batch.

### 4.4 Train/Val/Test Split Strategy

Split must be performed at the **edge level**, not the node level. Node-level splitting would leak information (a gene that appears in training as a target also appears in validation as a TF, giving the model indirect access to validation interactions). Edge-level splitting is the correct approach for link prediction tasks.

Split the positive edges of each evidence tier independently:
- 70% training, 15% validation, 15% test
- Use a fixed random seed (42) for reproducibility
- Negative edges are sampled fresh for each split: training negatives, validation negatives, and test negatives are disjoint samples from the prior graph minus all positives

Ensure test edges from the dual-evidence tier do NOT appear in any earlier tier's training set (they won't, by construction, since dual-evidence positives are not in training negatives). But verify this.

---

## Part 5 — Model Architecture Theory

### 5.1 Design Philosophy

The architectural choices are governed by three constraints: computational efficiency (fit and train on 4× RTX 2080 Ti in <24 hours), biological interpretability (each component should correspond to a meaningful biological signal), and performance (competitive with state-of-the-art lightweight methods like GMFGRN and RegDiffusion).

The total model parameter count is intentionally small: targeting 100K–500K parameters. This is feasible because the database provides the biological knowledge through its reference networks — the model's job is to learn how to interpolate from the evidence tiers, not to memorize a massive regulatory atlas.

### 5.2 Module 1 — RNA Encoder

**What it does**: Maps the 3-dimensional RNA gene feature vector (mean expression, expression variance, detection rate) to a higher-dimensional latent embedding.

**Theory**: A two-layer MLP (Multi-Layer Perceptron) with LayerNorm and GELU activation. LayerNorm is preferred over BatchNorm because it normalizes over the feature dimension rather than the batch dimension, which is more stable for variable-size gene universes and small batches. GELU (Gaussian Error Linear Unit) is preferred over ReLU because it provides smooth gradients at zero, which is important for sparse gene expression data where many values are exactly or near zero.

**Input**: (N_G × 3) matrix — one row per gene in the gene universe
**Output**: (N_G × D) matrix — D = 128 hidden dimension
**Architecture depth**: 2 layers (input → 64 → 128) with LayerNorm after each linear layer and dropout (p=0.2) for regularization

**Why not deeper**: With only low-dimensional input features and a dataset of thousands of genes (not millions), a deeper encoder would overfit. The GNN backbone provides additional capacity through graph-structured feature propagation, and the co-expression signal is injected directly at the decoder rather than through this encoder.

### 5.3 Module 2 — ATAC Encoder

**What it does**: Maps the 2-dimensional ATAC gene activity feature vector (mean accessibility, accessibility variance) to the same D-dimensional latent space as the RNA encoder.

**Theory**: A single linear projection followed by LayerNorm and GELU. This is intentionally simpler than the RNA encoder because ATAC gene activity scores are already a derived/aggregated signal (not raw counts), and because ATAC data is much sparser and noisier. A single projection is sufficient to embed it into the shared latent space.

**Input**: (N_G × 2) matrix
**Output**: (N_G × D) matrix, same dimension as RNA encoder output

**Why separate encoder for ATAC**: The RNA and ATAC features have different statistical distributions and encode different biological information. Treating them as a single 4-dimensional input vector would implicitly assume that their features are directly comparable, which they are not. Separate encoders allow the model to learn independent representations before fusion.

### 5.4 Module 3 — Role-Aware Modality Integration (replaces early fusion)

**Motivation — why not fuse the two modalities into one node vector.** The earlier design collapsed RNA and ATAC into a single per-gene embedding via a convex gate, h_fused = g ⊙ h_RNA + (1 − g) ⊙ h_ATAC. Because the weights sum to one per dimension, this is a *soft either/or*: it treats accessibility and expression as **interchangeable substitutes** — competing readouts of the same quantity — and closing the gate toward one modality suppresses the other. That is biologically backwards. scATAC measures whether a locus is *accessible* (whether regulation is physically possible); scRNA measures the *expression state* and TF–target *co-expression*. Accessibility is a **precondition** for expression-based regulation, not an alternative measurement of it. The two are complementary and asymmetric, so they must be kept as separate channels and combined by their **biological role**, not averaged.

Two further problems with early fusion: (1) it happens at the node level, erasing the asymmetry of a regulatory edge — on the TF side what matters is the regulator's activity (an RNA property), while on the target side what matters is whether the locus is accessible (an ATAC property); (2) the ablations confirmed the substitutive gate was not helping (plain concatenation, which at least preserves both channels, matched or beat it).

**Design.** There is no early fusion module. Instead:
- The RNA encoder output and ATAC encoder output are each carried through the GNN as **separate channels** (h_RNA and h_ATAC), producing two context-aware embeddings per gene.
- The modalities are integrated **at the edge**, inside the decoder (Part 5.6), with explicit roles: the target's RNA embedding is gated multiplicatively by an accessibility gate derived from the target's ATAC embedding (accessibility as a regulatability precondition), the TF's RNA embedding supplies regulator activity, and explicit co-expression and locus-openness terms are added.

**Input**: two (N_G × D) matrices from the RNA and ATAC encoders
**Output**: two (N_G × D) channel embeddings (no collapse)

**Ablations retained**: the legacy convex gated fusion (`integration='gated'`) and plain concatenation-then-projection (`integration='concat'`) are kept as ablations. Each collapses the two channels into one fused embedding and uses the plain bilinear decoder, so the comparison `role_aware` vs `gated` vs `concat` directly isolates the value of role-aware integration.

### 5.5 Module 4 — GNN Backbone

**What it does**: Performs graph-based message passing over the prior candidate graph, allowing each gene node to aggregate regulatory context from its neighbors (co-regulated genes and co-regulating TFs).

**Theory — Message Passing Neural Networks (MPNNs)**: In a standard MPNN, the update rule at layer l for node v is:

h_v^(l) = UPDATE(h_v^(l-1), AGG({h_u^(l-1) : u ∈ N(v)}))

Where N(v) is the neighborhood of v in the graph, AGG is a permutation-invariant aggregation function (typically mean, max, or sum), and UPDATE combines the node's own previous state with the aggregated neighbor information.

**Why GraphSAGE for this task**: GraphSAGE (Hamilton et al., 2017) uses mean aggregation and concatenates the node's own embedding with the aggregated neighbor embedding before a linear transformation. This design has two advantages for GRN inference: (1) it preserves the node's own information independently of neighborhood aggregation (important for isolated TF nodes with few prior graph connections), and (2) mean aggregation is robust to degree variation (TFs regulate different numbers of targets).

**Two relations, summed per layer**: The GNN aggregates over the two message-passing graphs from Part 4.1 — the co-expression kNN graph and the TF-candidate graph. Each layer runs one SAGEConv per relation and sums their outputs before normalization and activation, so a node's update mixes its co-expression module (gene–gene) with its regulatory context (TF↔target). The TF-candidate relation is made bidirectional for representation learning (TF nodes aggregate from candidate targets and vice versa); directionality is reintroduced only at the decoder. In the `role_aware` integration the RNA and ATAC channels are each passed through their own GNN backbone (two backbones); in the fused ablations a single backbone processes the fused embedding.

**Two-layer stack**: Two message-passing layers allow each node to aggregate information from 2-hop neighborhoods. After Layer 1, a gene knows its direct co-expression neighbors and (for TFs) its candidate targets; after Layer 2, it also incorporates the neighbors of those neighbors (e.g., co-regulated modules and shared-target TFs). Two layers strike the right balance between expressive power and over-smoothing for our graph sizes.

**Residual connections**: After the first GNN layer, add the input to the output of the second GNN layer (skip connection). This prevents over-smoothing (the tendency of deep GNNs to make all node embeddings identical due to repeated averaging) and helps gradient flow during backpropagation.

**Dropout in GNN layers**: Apply dropout (p=0.2) after each GNN layer's activation. This regularizes the graph representations and prevents the GNN from over-relying on any specific neighbor.

**Input**: (N_G × D) channel gene representations + the two message-passing edge_index tensors
**Output**: (N_G × D) context-aware gene embeddings, per channel

### 5.6 Module 5 — Role-Aware Decoder

**What it does**: Scores each candidate TF→Target pair by integrating the two modality channels according to their biological roles, rather than scoring a single fused embedding.

**Theory**: For a pair (TF_i, Target_j), with RNA-channel embeddings h_RNA and ATAC-channel embeddings h_ATAC:

- a_j = sigmoid(gate(h_ATAC[j])) ∈ [0,1]^D — a per-dimension **accessibility gate** on the target locus. This encodes the biology that a target must be accessible to be regulated: accessibility multiplicatively modulates (rather than substitutes for) the target's expression representation.
- z_j = h_RNA[j] ⊙ a_j — the target's **regulatable-expression state**.
- raw = (h_RNA[i] · W_dec) · z_j — an **asymmetric** TF→target compatibility, with the TF's RNA embedding acting as regulator activity and W_dec a D×D learnable interaction matrix.
- score(TF_i, Target_j) = sigmoid( raw + w_coexp · (s_i · s_j) + w_open · openness_j + b )

where s_i · s_j is the co-expression signature dot product (≈ Pearson correlation, Part 3.1), openness_j is the target-locus openness scalar (Part 3.2), and w_coexp, w_open, b are learnable scalars.

**Why this is the right structure**: Each modality plays a distinct, non-interchangeable role — ATAC gates the *target's regulatability*, RNA supplies *TF activity* and *TF–target co-expression*, and openness contributes a *prior on how regulatable the target locus is*. This directly encodes the complementarity the earlier convex fusion violated. The bilinear core keeps the asymmetry needed for directed TF→target prediction (Distmult, Yang et al. 2014; Decagon, Zitnik et al. 2018) while adding the two biologically explicit terms that re-introduce the co-expression and accessibility evidence.

**Input**: two (N_G × D) channel embeddings + a batch of (TF_index, Target_index) pairs + signatures + openness
**Output**: (batch_size,) regulatory logits (sigmoid at loss/eval time)

**Ablation head**: the fused integrations (`gated`/`concat`) instead use a plain bilinear decoder on the single fused embedding, matching the original design.

### 5.7 Why Not Transformers, Mamba, or VAEs

These alternatives are explicitly ruled out for this problem and hardware budget:

- **Transformers**: Full self-attention over all gene pairs is O(N_G²) in time and memory. For N_G = 3000 genes, this is 9 million attention entries per layer. On 2080 Ti with 11GB VRAM, fitting a full transformer over the gene graph is infeasible without extensive engineering. Moreover, transformers do not natively encode the bipartite TF-target graph structure that is central to GRN biology.

- **Mamba / State-Space Models**: Designed for long sequential inputs. GRN inference is fundamentally a graph-structured problem, not a sequential one. Mapping genes to a sequence order is arbitrary and biologically unmotivated.

- **Variational Autoencoders (VAEs)**: Useful for learning cell embeddings (e.g., scVI), but not naturally suited to the edge-prediction formulation of GRN inference. They add reconstruction loss complexity and a latent sampling step that increases training variance without clear benefit for link prediction.

The compact GNN with bilinear decoder is the state-of-the-art choice for efficient, biologically principled GRN inference, as supported by GMFGRN (Li et al., 2024, PubMed 38261340) and GNNLink (Briefings in Bioinformatics 2023).

---

## Part 6 — Multi-Evidence Curriculum Learning Theory

### 6.1 What Curriculum Learning Is

Curriculum learning (Bengio et al., 2009) is the training strategy of presenting samples to a model in a meaningful order, typically from easier to harder, rather than in random order. The biological motivation is that learning broad, approximate patterns first gives the model a good initialization for learning the harder, more selective patterns. This mirrors how humans learn: general concepts before fine distinctions.

In graph neural network settings, curriculum learning has been shown to improve generalization and final performance. CurGraph (Wang et al., WWW 2021) demonstrated curriculum graph classification; CLNode (IEEE TKDE 2022) demonstrated curriculum node classification. The key principle is the same: ordering matters, and training on progressively more demanding objectives improves final model quality.

### 6.2 How the Three Evidence Tiers Map to Curriculum

**Stage 1: Localization Supervision (E_loc)**

Biological meaning: Use ChIP-seq / ChIP-chip TF binding data as supervision. The model is trained to distinguish TF-target pairs where the TF is known to physically bind near the target's TSS from pairs where there is no such binding evidence. This is the easiest signal to learn because there are many positive examples and the signal (TF binding = likely nearby target) is relatively direct.

Training setup: 30 epochs, learning rate 1e-3, standard BCE with weighted positive class. Negative ratio 5:1. Train all model components from random initialization.

What the model learns at this stage: The general spatial relationship between TF expression/accessibility and target gene accessibility. The GNN learns co-regulatory context (TFs that bind similar targets cluster together in the embedding space). The bilinear decoder learns a rough directional compatibility between TF and target embeddings.

**Stage 2: Perturbation Supervision Fine-Tuning (E_pert)**

Biological meaning: The model now receives supervision from perturbation-based networks: TF-target pairs where experimentally knocking out or knocking down TF_i causes differential expression of Target_j. This is functional/causal evidence. The positive set (E_pert) is a strict subset of E_loc typically, but not always — some perturbation edges do not have corresponding ChIP-seq data.

Training setup: 15 epochs of fine-tuning, reduced learning rate (3e-4), all model parameters trainable. Introduce hard negatives from E_loc \ E_pert (localization edges NOT in perturbation — TF binds near target but does not regulate it).

What fine-tuning accomplishes: The model refines its embeddings to distinguish physical binding (easy, learned in Stage 1) from functional regulation (harder). The gate in the fusion module may shift to emphasize ATAC accessibility more (functional regulation requires accessible chromatin) or RNA expression more (regulatory edges tend to be between co-expressed genes). The GNN refines its co-regulatory context toward functionally related genes.

**Stage 3: Dual-Evidence Supervision Fine-Tuning (E_dual)**

Biological meaning: The model now receives supervision from only the highest-confidence interactions: TF-target pairs with BOTH binding AND perturbation evidence. This is the most selective, highest-precision positive set. E_dual is a strict subset of both E_loc and E_pert.

Training setup: 10 epochs of fine-tuning, very low learning rate (1e-4), encoder weights frozen (only GNN and decoder are updated — the encoder representations are already well-learned). Negative ratio increased to 10:1 to enforce tighter precision. Hard negatives from E_pert \ E_dual.

Why freeze encoders in Stage 3: The encoder representations of gene expression and chromatin accessibility are already well-calibrated from Stage 1 and 2. Continuing to update them at this stage would be prone to catastrophic forgetting (overwriting useful representations with signal from very few positive examples). The GNN backbone and bilinear decoder have more capacity to absorb the Stage 3 signal without needing to update the raw feature embeddings.

What final fine-tuning accomplishes: The decoder bilinear matrix W_dec is refined to give highest scores to interactions that satisfy both binding AND functional evidence. This sharpens the precision at the top of the ranking, which is what matters biologically — the model's top-ranked predictions should align with the highest-confidence known regulatory interactions.

### 6.3 Memory Replay (Optional)

A risk of sequential fine-tuning is catastrophic forgetting: Stage 2 fine-tuning might reduce performance on Stage 1's localization test set. To mitigate this, memory replay can be enabled: during Stage 2 training, a small fraction (10%) of each mini-batch is drawn from the Stage 1 positive set at a reduced loss weight (0.1). Similarly during Stage 3, replay from both Stage 1 and Stage 2 edges. This is analogous to experience replay in reinforcement learning and elastic weight consolidation (EWC).

Implement memory replay as an optional flag (use_memory_replay: True/False in config). Run ablations both with and without it to determine if it helps for this dataset scale.

### 6.4 Curriculum vs. Multi-Task (the ablation)

An alternative to curriculum learning is multi-task learning: train simultaneously on all three evidence tier losses, with separate loss weights for each tier. This is the "all-at-once" ablation. The curriculum hypothesis predicts that the ordered curriculum (localization → perturbation → dual) will outperform multi-task training on the dual-evidence test set, because the ordered training provides a better initialization landscape for learning the hardest signal.

---

## Part 7 — Training Theory

### 7.1 Loss Function

**Primary loss: Weighted Binary Cross-Entropy (BCE)**

For each mini-batch of edge pairs (positive and negative), compute:

L = −(1/N_batch) × Σ [ w_pos × y_i × log(p_i) + (1−y_i) × log(1−p_i) ]

Where:
- y_i ∈ {0, 1} is the label (1 for positive edges, 0 for negative edges)
- p_i ∈ [0,1] is the model's predicted probability
- w_pos is the positive class weight, set to the negative sampling ratio (e.g., 5 if 5 negatives per positive are sampled in this mini-batch)

Weighted BCE corrects for the class imbalance within mini-batches. Without weighting, the loss would be dominated by the numerically larger negative class, and the model would learn to predict everything as negative (the trivial solution with low loss but no discriminative power).

**Alternative: Focal Loss (for Stage 3)**

Focal loss (Lin et al., 2017) adds a modulating factor (1 − p_t)^γ to the standard BCE:

L_focal = −(1/N) × Σ α × (1 − p_t)^γ × log(p_t)

Where p_t is the predicted probability of the true class. This down-weights easy examples (edges the model already predicts correctly with high confidence) and focuses gradient updates on hard examples (edges near the decision boundary). For Stage 3's small, high-confidence positive set, focal loss (with α=0.75, γ=2.0) may improve the model's ability to correctly rank the small number of dual-evidence positives above the large number of negatives. Test both BCE and focal loss for Stage 3 and report in ablations.

### 7.2 Mini-Batch Graph Training

**Problem**: The full prior graph may have hundreds of thousands of edges. Computing GNN message passing over the full graph for every training step is computationally expensive and may not fit in VRAM.

**Solution: Edge mini-batching with full-graph GNN encoding**

For our scale (N_G ≈ 3,000 genes, N_prior ≈ 1.5 million candidate edges), the gene node feature matrix fits trivially in VRAM (3000 × 128 × 4 bytes ≈ 1.5 MB). Only the GNN message passing over E_prior is the bottleneck.

For K562 with ~3,000 genes and 500 prior targets per TF (~200 TFs in universe → 100,000 prior edges), full-graph GNN encoding at each step is entirely feasible on the 2080 Ti. Perform ONE full-graph GNN encode per epoch (or per gradient step, whichever is cheaper) to compute all node embeddings h_v for all genes. Then sample random mini-batches of edges for the bilinear decoder and loss computation.

**When to use ClusterGCN / GraphSAINT**: If the prior graph grows large (>500k edges, e.g., when using HSC or GM12878 with more TFs), switch to mini-batch graph training:
- **ClusterGCN** (Chiang et al., 2019): Partition the gene nodes into K clusters using METIS or random partitioning. Each mini-batch processes one or more clusters and only computes message passing within the cluster's subgraph. This drastically reduces memory and allows parallel processing of multiple clusters. Cluster size controls the trade-off between variance (small clusters) and memory (large clusters).
- **GraphSAINT** (Zeng et al., 2020): Instead of clustering nodes, sample random subgraphs by sampling edges with importance weights (based on node degree, random walk probabilities, or edge degree). Each mini-batch processes the sampled subgraph. GraphSAINT corrects for sampling bias via per-node and per-edge normalization. It is generally better for heterogeneous degree distributions.

For the default configuration (K562, N_G ≈ 3000), full-graph encoding is recommended. Implement ClusterGCN as a fallback option that is activated by a config flag.

### 7.3 Optimizer and Learning Rate Schedule

Use AdamW optimizer (Adam with decoupled weight decay, Loshchilov & Hutter 2019). AdamW is the standard choice for modern deep learning because it correctly separates the adaptive learning rate from L2 regularization. Weight decay λ = 1e-4.

Learning rate schedule: CosineAnnealingLR with T_max equal to the number of epochs in the current stage. Cosine annealing starts at the stage's peak LR and smoothly decays to near zero, which allows the model to explore early in each stage and converge tightly to a local minimum by the end. This is particularly important for fine-tuning stages where overstepping would destroy the prior stage's learned representations.

Gradient clipping: Apply max-norm gradient clipping (clip_norm = 1.0) before each optimizer step. This prevents exploding gradients, which are common when fine-tuning across stages with very different loss scales (localization has many positive edges, dual-evidence has very few).

### 7.4 Early Stopping

Within each curriculum stage, apply early stopping based on validation AUPR. If validation AUPR does not improve for patience=10 consecutive validation checks (evaluated every 5 epochs), terminate the current stage early and advance to the next. Save the best checkpoint (highest validation AUPR within the stage) rather than the final epoch's checkpoint.

---

## Part 8 — Evaluation Theory

### 8.1 Why AUPR is the Primary Metric

AUPR (Area Under the Precision-Recall Curve) is the primary evaluation metric for GRN inference for the following reasons:

1. **Class imbalance**: GRN networks are extremely sparse. Even with a prior graph of 100,000 candidate edges and 1,000 positive edges (a generous estimate), the positive fraction is 1%. In such heavily imbalanced settings, the Precision-Recall curve is far more informative than the ROC curve. A method can achieve AUROC > 0.9 by slightly doing better than random, which is misleading. AUPR directly measures the model's ability to place true positives at the top of the ranked list without false positives.

2. **Biological relevance**: In practice, researchers use GRN inference to generate a ranked list of regulatory edges to validate experimentally. They only validate the top few dozen or hundred predictions. The precision of those top-ranked predictions (the "early precision") directly determines the utility of the method.

3. **Field standard**: The BEELINE benchmarking framework (Pratapa et al., Nature Methods 2020) and subsequent GRN papers report AUPR as the primary metric. Using AUPR ensures comparability with the literature.

### 8.2 AUROC as Secondary Metric

AUROC (Area Under the ROC Curve) measures the probability that a randomly chosen positive edge is ranked above a randomly chosen negative edge. It is less sensitive to class imbalance than AUPR but is still reported for completeness and cross-study comparison.

### 8.3 Early Precision and EPR

**Early Precision (EP)**: Precision@K where K = number of true positives in the test set. This measures: among the model's top-K predictions, what fraction are true regulatory edges? This is the most practically relevant metric because it evaluates the quality of the model's top-ranked predictions.

**Enrichment over Random (EPR)**: EP divided by the random baseline precision (which is: number of positives / total number of candidates). EPR > 1 means the model is better than random. EPR > 2 means predictions are at least twice as enriched for true positives as random. Report EPR alongside absolute EP.

### 8.4 Per-Tier Evaluation Protocol

Evaluate the final model (after all three curriculum stages) against the test split of each evidence tier separately:

- Evaluation on E_loc test set: measures ability to recover TF binding sites
- Evaluation on E_pert test set: measures ability to recover functional regulatory edges
- Evaluation on E_dual test set: measures ability to recover highest-confidence causal edges (the most biologically meaningful test)

This per-tier evaluation directly substantiates the curriculum's benefit: the model should improve on E_dual by going through Stage 2 and Stage 3, while maintaining reasonable performance on E_loc.

### 8.5 Cross-Cell-Type Transfer Evaluation

After training on K562, evaluate (without any retraining) on:
1. ESC scRNA + scATAC data with ESC reference networks

This tests generalizability: a model that has learned general principles of TF-target regulatory logic (rather than memorizing K562-specific patterns) should transfer to a different cell type. Report AUPR on ESC E_loc, E_pert, E_dual test sets using the model trained on K562.

Cross-cell-type transfer is a particularly compelling result for journal reviewers because it demonstrates that the model has learned something biologically meaningful, not just fit noise.

### 8.6 Runtime Measurement

Measure and report:
- Preprocessing time (scRNA + scATAC, from raw h5ad to processed features)
- Training time per stage
- Total training time (all three stages)
- Inference time (scoring all prior graph edges once)
- Peak GPU memory during training

This supports the efficiency claim in the paper abstract.

---

## Part 9 — Baseline Methods Theory

### 9.1 GRNBoost2 / GENIE3

GRNBoost2 (part of the Arboreto library) is a tree-ensemble method based on gradient boosting. It treats GRN inference as a regression problem: for each target gene j, train a regression model to predict the expression of gene j from the expression of all TF genes. The importance of each TF_i in predicting the expression of Target_j is interpreted as the regulatory edge weight. This is done in parallel for all target genes.

Why include: It is the most widely used classical baseline for GRN inference and is included in every major benchmarking study. It requires only scRNA-seq input (no graph structure, no ATAC). It is fast (~5-30 minutes on commodity hardware). Include it to demonstrate that the graph-based and multimodal approach provides improvements beyond classical co-expression-based methods.

Installation: pip install arboreto

### 9.2 RegDiffusion

RegDiffusion is a diffusion-model-based GRN inference method (Luyuan et al., PubMed 39387266). It learns to model the noise process over gene expression and uses the diffusion trajectory to infer regulatory relationships. It is reported to run in under 5 minutes on large single-cell datasets, making it the fastest deep learning GRN method currently available.

Why include: RegDiffusion establishes the lower bound on runtime for deep learning GRN methods. If MEvD-GRN achieves significantly better AUPR than RegDiffusion on the dual-evidence test set, it justifies the additional training time of the curriculum approach.

Installation: pip install regdiffusion

### 9.3 GMFGRN

GMFGRN (Li et al. 2024, GitHub: https://github.com/Lishuoyy/GMFGRN) is a graph neural network method that uses GNN-based matrix factorization for GRN inference. It learns gene embeddings via message passing and uses them to score TF-target pairs, similar to MEvD-GRN's architecture. It is reported as dramatically more memory-efficient than competing GNN methods and achieves competitive AUPR.

Why include: GMFGRN is the most direct architectural competitor to MEvD-GRN. Both use GNN + bilinear-style decoding. The difference is that GMFGRN does not use ATAC data, does not use multi-evidence curriculum training, and does not use SC-MO-GRN-DB. Showing MEvD-GRN outperforms GMFGRN on the SC-MO-GRN-DB evaluation directly establishes the value of the curriculum approach.

---

## Part 10 — Ablation Study Design Theory

Each ablation tests one specific hypothesis about the model design.

| Ablation Name | What is Changed | Hypothesis Tested |
|---|---|---|
| full_curriculum | Full model (role-aware integration), all three stages | N/A — this is the full method |
| loc_only | Train only on localization, skip Stages 2–3 | Does curriculum help beyond Stage 1 alone? |
| pert_only | Train only on perturbation (Stage 2 only) | Is Stage 1 initialization necessary? |
| dual_only | Train only on dual-evidence (Stage 3 only) | Can the model learn from only the hardest supervision? |
| all_at_once | Train simultaneously on union of all tiers, equal weight | Is curriculum ordering better than multi-task? |
| rna_only | Zero all ATAC-derived signal (features, openness, gate) | What does chromatin accessibility contribute? |
| gated_fusion | Replace role-aware integration with legacy convex gated fusion + bilinear decoder | Does the role-aware integration beat substitutive fusion? |
| concat_fusion | Replace role-aware integration with concatenation-then-projection + bilinear decoder | Does role-aware integration beat naive concatenation? |
| no_gnn | Remove GNN backbone; decode from encoder features directly | What do the co-expression / TF-candidate graphs contribute? |
| with_replay | Full curriculum + memory replay from earlier stages | Does replay prevent catastrophic forgetting? |

Each ablation uses identical hyperparameters to the full model except the one changed component. Run all ablations on K562, evaluate on K562 dual-evidence test set. The `gated_fusion` and `concat_fusion` ablations together isolate the contribution of the role-aware integration (accessibility gate + co-expression + openness) versus collapsing the two modalities into one fused embedding.

---

## Part 11 — Repository Structure

Create the following directory and file structure. Do not deviate from this structure, as subsequent parts reference specific paths.

```
MEvD-GRN/
├── README.md
├── plan.md
├── configs/
│   ├── default.yaml
│   ├── k562.yaml
│   └── esc.yaml
├── data/
│   ├── raw/
│   │   ├── networks/
│   │   └── single_cell/
│   ├── processed/
│   │   ├── K562/
│   │   └── ESC/
│   └── splits/
├── src/
│   ├── __init__.py
│   ├── data/
│   │   ├── __init__.py
│   │   ├── preprocessing.py
│   │   ├── feature_extraction.py
│   │   ├── graph_builder.py
│   │   └── dataset.py
│   ├── models/
│   │   ├── __init__.py
│   │   ├── encoders.py
│   │   ├── fusion.py
│   │   ├── gnn.py
│   │   ├── decoder.py
│   │   └── mevd_grn.py
│   ├── training/
│   │   ├── __init__.py
│   │   ├── curriculum.py
│   │   ├── sampler.py
│   │   ├── losses.py
│   │   └── trainer.py
│   ├── evaluation/
│   │   ├── __init__.py
│   │   ├── metrics.py
│   │   └── benchmarker.py
│   ├── baselines/
│   │   ├── __init__.py
│   │   ├── grnboost2_wrapper.py
│   │   ├── regdiffusion_wrapper.py
│   │   └── gmfgrn_wrapper.py
│   └── utils/
│       ├── __init__.py
│       ├── io.py
│       └── visualization.py
├── scripts/
│   ├── 01_download_guide.md
│   ├── 02_preprocess.py
│   ├── 03_train.py
│   ├── 04_evaluate.py
│   ├── 05_run_baselines.py
│   └── 06_ablation.py
├── slurm/
│   ├── train.sh
│   ├── baselines.sh
│   └── ablation.sh
├── notebooks/
│   ├── 01_data_exploration.ipynb
│   └── 02_results_analysis.ipynb
├── results/
│   ├── checkpoints/
│   └── figures/
└── requirements.txt
```

---

## Part 12 — Detailed Implementation Steps

Work through these steps in order. Each step references the relevant source file and the theory section that justifies the design.

---

### Step 1 — Environment Setup

1.1. Create a new conda environment named `mevd-grn` with Python 3.10.

1.2. Install PyTorch 2.1.0 with CUDA 11.8 support using the official PyTorch wheels.

1.3. Install PyTorch Geometric 2.4.0 and its required sparse dependencies (torch_scatter, torch_sparse) from the PyG wheel index for torch-2.1.0+cu118.

1.4. Install the following packages via pip: scanpy, anndata, muon, numpy, scipy, pandas, scikit-learn, matplotlib, seaborn, tqdm, pyyaml, pyranges, pybiomart, arboreto, tqdm.

1.5. Optionally install wandb for experiment tracking if desired.

1.6. Verify CUDA setup: confirm that `torch.cuda.device_count()` returns 4 and `torch.cuda.get_device_name(0)` returns the expected GPU name.

1.7. Create `requirements.txt` with pinned versions of all installed packages.

1.8. Create `setup.py` with the `mevd_grn` package declaration pointing to the `src/` directory.

---

### Step 2 — Configuration File

2.1. Create `configs/default.yaml` with the following top-level sections and values. All training scripts must read from this config file rather than using hardcoded values.

**data section**: cell_types list (K562, ESC), evidence_tiers list (localization, perturbation, dual_evidence), n_hvg=3000, atac_window_bp=100000, top_k_prior_targets=500, train_ratio=0.70, val_ratio=0.15, test_ratio=0.15, neg_train_ratio=5, seed=42.

**model section**: rna_in_dim=2, atac_in_dim=2, hidden_dim=128, n_gnn_layers=2, dropout=0.2, gnn_type=sage.

**curriculum section**: Three stage entries (Stage1_Localization with n_epochs=30 and lr=1e-3; Stage2_Perturbation with n_epochs=15 and lr=3e-4; Stage3_DualEvidence with n_epochs=10 and lr=1e-4). freeze_encoder_stage3=true. use_memory_replay=false. replay_weight=0.1.

**training section**: optimizer=adamw, weight_decay=1e-4, batch_size=2048, clip_grad_norm=1.0, early_stopping_patience=10, loss_type=bce_weighted, neg_stage3_ratio=10.

**evaluation section**: metrics list (auroc, aupr, early_precision, epr), test_on_all_tiers=true, cross_cell_type=true.

**hardware section**: device=cuda:0, n_dataloader_workers=4.

2.2. Create `configs/k562.yaml` and `configs/esc.yaml` that inherit from default.yaml and override only cell-type-specific settings (e.g., file paths, any cell-type-specific HVG thresholds).

---

### Step 3 — Data Download Guide

3.1. Create `scripts/01_download_guide.md` as a plain text guide (not a script) documenting the exact steps to download from SC-MO-GRN-DB:

- Navigate to https://scmogrndb.psu.edu
- Go to the Reference Networks section
- Use the filter dropdowns to select Cell Type = K562 and Evidence Type = Localization. Download the TSV file and save as `data/raw/networks/K562_localization.tsv`
- Repeat for Perturbation → `data/raw/networks/K562_perturbation.tsv`
- Repeat for Dual evidence → `data/raw/networks/K562_dual_evidence.tsv`
- Repeat all three for ESC
- Go to the Single-Cell Datasets section
- Filter by Cell Type = K562, Modality = scRNA, download and save as `data/raw/single_cell/K562_scRNA.h5ad`
- Filter by Modality = scATAC, save as `data/raw/single_cell/K562_scATAC.h5ad`
- Repeat for ESC
- Note: if any cell type has a joint assay (RNA + ATAC from the same cells, e.g., 10x Multiome), prefer that file over separate modality files

3.2. Document the expected file sizes and column formats:
- Networks: two-column TSV, no header, column 0 = TF name (gene symbol), column 1 = target gene name (gene symbol)
- Single-cell: .h5ad format (AnnData), rows = cells, columns = genes or peaks depending on modality
- Check that the organism label in the AnnData metadata matches expectations (human for K562, human or mouse for ESC)

---

### Step 4 — Preprocessing Module (`src/data/preprocessing.py`)

This module implements the theory from Part 3.

4.1. Implement a `preprocess_scrna(adata_path, n_hvg, output_dir)` function with the following internal steps:

- Load the h5ad file using `scanpy.read_h5ad`
- Run QC filtering: remove cells with fewer than 200 detected genes; remove genes detected in fewer than 3 cells
- Apply `sc.pp.normalize_total` with target_sum=10000 to library-size normalize
- Apply `sc.pp.log1p` for log(x+1) transformation
- Run `sc.pp.highly_variable_genes` with n_top_genes=n_hvg and flavor='seurat_v3' to identify the top HVGs
- Extract the normalized log-count matrix (cells × HVGs)
- Compute per-gene mean and variance across all cells; stack into a (n_hvg × 2) numpy array
- Save: the (n_hvg × 2) feature array as `gene_rna_features.npy`; the list of gene names as `gene_names.json`; the HVG boolean mask as `hvg_mask.npy`
- Return the gene names list (this defines the gene universe for this cell type)

4.2. Implement a `preprocess_scatac(adata_path, gene_names, genome, atac_window_bp, output_dir)` function with these internal steps:

- Load the scATAC h5ad file using `scanpy.read_h5ad`. The feature names (adata.var_names) should be genomic peak coordinates in the format "chrN:start-end".
- Run QC filtering: remove cells with fewer than 100 detected peaks; remove peaks detected in fewer than 3 cells
- Apply `sc.pp.normalize_total` with target_sum=10000 and `sc.pp.log1p`
- Parse peak coordinates from adata.var_names: extract chromosome, start position, end position for each peak. Compute peak midpoint = (start + end) / 2
- Load TSS annotations for all genes in `gene_names` using the `pybiomart` library. Query the Ensembl BioMart for the genome assembly (hsapiens_gene_ensembl for hg38, mmusculus_gene_ensembl for mm39). Retrieve: external_gene_name, chromosome_name, transcription_start_site. Prepend 'chr' to chromosome names if not present. If BioMart is unavailable, warn and proceed with zero ATAC features (handled downstream).
- For each gene in `gene_names`, identify which peak indices have their midpoint within atac_window_bp of the gene's TSS and on the same chromosome. This is the peak-to-gene mapping.
- Compute the Gene Activity Score: for each cell and each gene, sum the log-normalized accessibility values of all peaks mapped to that gene
- Compute per-gene mean and variance of the Gene Activity Score across all cells; stack into (n_genes × 2) numpy array
- If a gene has no mapped peaks, set its ATAC features to (0.0, 0.0)
- Save the (n_genes × 2) ATAC feature array as `gene_atac_features.npy`
- Report: fraction of genes with non-zero ATAC features, mean number of peaks per gene

4.3. Implement a `align_features(rna_features, atac_features, gene_names, network_dir, cell_type)` function:

- Load all reference network files for the cell type from `network_dir`, collect all unique TF and target gene names
- Compute the final gene universe: intersection of gene_names (from scRNA), genes with non-zero ATAC features, and genes appearing in any reference network
- Filter rna_features and atac_features to only the rows corresponding to genes in the final gene universe
- Create and save `gene_index.json`: a dictionary mapping each gene name to its integer row index in the feature matrices
- Create and save `tf_indices.json`: the list of integer indices (in gene_index) of genes that appear as TFs in any reference network
- Save aligned `rna_features_aligned.npy` and `atac_features_aligned.npy`

4.4. Create `scripts/02_preprocess.py` as the command-line entry point that:
- Parses arguments: --cell_types (list), --config (path to config yaml)
- Calls preprocess_scrna and preprocess_scatac for each cell type
- Calls align_features to finalize the gene universe
- Saves all processed files to `data/processed/{cell_type}/`
- Prints summary statistics: number of genes, number of TFs, coverage of reference network genes

---

### Step 5 — Graph Construction Module (`src/data/graph_builder.py`)

5.1. Implement `build_prior_graph(gene_index, tf_indices, atac_features, top_k_targets)` that:

- For each TF index in tf_indices, rank all other genes by their ATAC mean accessibility value (atac_features[:, 0])
- Exclude the TF's own index from the target list
- Select the top_k_targets (default 500) highest-accessibility genes as candidate targets for that TF
- Build two parallel lists: source TF indices and destination target indices
- Return a (2, N_prior_edges) integer tensor (edge_index in PyTorch Geometric convention) and the total number of prior edges
- This prior graph structure is used BOTH as the GNN message-passing graph AND as the candidate edge space for scoring

5.2. Implement `load_network_edges(network_path, gene_index)` that:

- Reads a two-column TSV network file
- Maps each TF name and target gene name to their integer index using gene_index
- Skips rows where either name is not in gene_index (out-of-universe genes)
- Returns a (2, N_positive) integer tensor of positive edges and a count of skipped rows
- Report the coverage rate: what fraction of the database edges are within the gene universe

5.3. Implement `filter_to_prior(positive_edges, prior_edge_set)` that:

- Converts prior_edges (2 × N_prior) to a Python set of (src, dst) tuples
- Filters positive_edges to only those that also appear in prior_edge_set
- Returns filtered positive edges and the count of edges lost due to prior filtering

5.4. Implement `create_negative_pool(prior_edges, all_positive_sets, gene_index, tf_indices)`:

- Creates the pool of candidate negative edges: prior graph edges NOT in ANY positive edge set (not in E_loc, E_pert, or E_dual)
- This ensures negatives are consistent across all training stages — an edge that is positive in any tier is NEVER used as a negative
- Returns a list or tensor of candidate negative edge pairs

5.5. Implement `create_edge_splits(positive_edges, negative_pool, train_ratio, val_ratio, seed)`:

- Shuffles positive edges with the fixed seed
- Splits into train/val/test by the given ratios
- Samples negative edges from negative_pool for each split at the appropriate ratio (neg_train_ratio for training, 5:1 for val and test)
- Returns a dict with keys 'train', 'val', 'test', each containing a dict with 'pos' and 'neg' edge tensors

5.6. Implement `save_all_evidence_splits(data_dir, cell_type, evidence_graphs, prior_edges, config)` as a coordinator that:

- Calls create_negative_pool once to build the universal negative pool
- For each evidence tier, calls create_edge_splits and saves the result to `data/splits/{cell_type}_{tier}_splits.pt`
- Saves prior_edges to `data/processed/{cell_type}/prior_edges.pt`
- Saves negative_pool to `data/processed/{cell_type}/negative_pool.pt`

---

### Step 6 — Model Implementation

**6.1. RNA Encoder (`src/models/encoders.py`)**

Implement a class `RNAEncoder` (inherits torch.nn.Module) with:
- `__init__`: Creates a sequential network: Linear(rna_in_dim → hidden_dim//2), LayerNorm(hidden_dim//2), GELU(), Dropout(dropout), Linear(hidden_dim//2 → hidden_dim), LayerNorm(hidden_dim). All dimensions are constructor arguments.
- `forward(x)`: Applies the sequential net to input x of shape (N_genes, rna_in_dim), returns (N_genes, hidden_dim)

**6.2. ATAC Encoder (`src/models/encoders.py`)**

Implement a class `ATACEncoder` (inherits torch.nn.Module) with:
- `__init__`: Creates: Linear(atac_in_dim → hidden_dim), LayerNorm(hidden_dim)
- `forward(x)`: Applies linear, LayerNorm, then GELU. Input (N_genes, atac_in_dim), output (N_genes, hidden_dim)

**6.3. Gated Fusion Module (`src/models/fusion.py`)**

Implement a class `GatedFusion` (inherits torch.nn.Module) with:
- `__init__`: Creates a gate network: Linear(hidden_dim × 2 → hidden_dim), followed by Sigmoid. The gate input is the concatenation of both modal embeddings.
- `forward(h_rna, h_atac)`: Concatenate h_rna and h_atac along the feature dimension to get (N_genes, 2×hidden_dim). Pass through gate network to get g of shape (N_genes, hidden_dim) with values in [0,1]. Return g * h_rna + (1 − g) * h_atac, which is the fused representation of shape (N_genes, hidden_dim).

**6.4. GNN Backbone (`src/models/gnn.py`)**

Implement a class `GNNBackbone` (inherits torch.nn.Module) with:
- `__init__`: Creates a list of `n_gnn_layers` SAGEConv layers (from torch_geometric.nn). For bipartite message passing, pass a tuple `(hidden_dim, hidden_dim)` as `in_channels` to SAGEConv to allow different source and destination node dimensions. Creates matching LayerNorm layers. Creates a Dropout layer.
- `forward(x, edge_index)`: Applies each layer sequentially: SAGEConv → LayerNorm → GELU → Dropout. After the first layer, add a residual connection (x + h_new) if the dimensions match. The edge_index should encode BOTH directions (TF→Target AND Target→TF) so both node types aggregate from each other. Construct the bidirectional edge_index by concatenating [src, dst] with [dst, src] if not already bidirectional.
- Note on bipartite: In the prior graph, most edges are TF→Target. To allow target genes to aggregate from TFs and TFs to aggregate from targets, ensure the message passing is bidirectional. This is not the same as treating the edge as undirected for the DECODER (where directionality matters for scoring). The GNN uses bidirectional propagation only for representation learning, not for the final prediction.

**6.5. Bilinear Decoder (`src/models/decoder.py`)**

Implement a class `BilinearDecoder` (inherits torch.nn.Module) with:
- `__init__`: Creates a learnable parameter matrix W of shape (hidden_dim, hidden_dim) initialized with Xavier uniform initialization, and a scalar bias parameter b initialized to 0.
- `forward(h_tf, h_target)`: Both inputs are (batch_size, hidden_dim). Compute h_transformed = h_tf @ W of shape (batch_size, hidden_dim). Compute the dot product between h_transformed and h_target element-wise and sum along the feature dimension to get a (batch_size,) score vector. Add the bias b. Apply sigmoid to return predicted probabilities in [0,1].

**6.6. Main Model Class (`src/models/mevd_grn.py`)**

Implement a class `MEvDGRN` (inherits torch.nn.Module) with:
- `__init__`: Instantiates RNAEncoder, ATACEncoder, GatedFusion, GNNBackbone, BilinearDecoder as submodules. All hyperparameters (hidden_dim, n_gnn_layers, dropout) are constructor arguments.
- `encode(rna_features, atac_features, prior_edge_index)`: Sequential call: rna_enc = RNAEncoder(rna_features); atac_enc = ATACEncoder(atac_features); h_fused = GatedFusion(rna_enc, atac_enc); h_out = GNNBackbone(h_fused, prior_edge_index). Returns h_out of shape (N_genes, hidden_dim).
- `decode(node_embeddings, tf_indices, target_indices)`: Gathers h_tf = node_embeddings[tf_indices] and h_target = node_embeddings[target_indices]. Calls BilinearDecoder(h_tf, h_target). Returns scores of shape (batch_size,).
- `forward(rna_features, atac_features, prior_edge_index, tf_indices, target_indices)`: Calls encode then decode.
- Implement a `count_parameters()` method that returns the total number of trainable parameters.

---

### Step 7 — Curriculum and Sampler Modules

**7.1. Curriculum Manager (`src/training/curriculum.py`)**

Implement a `CurriculumStage` dataclass with fields: name (str), evidence_tier (str), n_epochs (int), learning_rate (float), neg_sampling_ratio (int), freeze_encoder (bool).

Define a `CURRICULUM_STAGES` list of three CurriculumStage instances matching the values in the config: Stage1_Localization (30 epochs, lr=1e-3, freeze=False), Stage2_Perturbation (15 epochs, lr=3e-4, freeze=False), Stage3_DualEvidence (10 epochs, lr=1e-4, freeze=True).

Implement a `CurriculumManager` class with:
- `__init__`: Takes the list of stages as input, initializes current_stage_idx=0, initializes a history dict for logging.
- `current_stage` property: Returns the stage at current_stage_idx.
- `advance()`: Increments current_stage_idx if not at the last stage.
- `is_done()`: Returns True when all stages are complete.
- `get_encoder_freeze_state()`: Returns the freeze_encoder value of the current stage.

**7.2. Negative Sampler (`src/training/sampler.py`)**

Implement `sample_negatives(negative_pool, n_to_sample, seed)`:
- Takes the pre-built negative pool (tensor of candidate negative edges), a requested count, and a seed
- Randomly samples n_to_sample edges from the pool without replacement
- Returns a (2, n_to_sample) tensor

Implement `get_hard_negatives(evidence_graphs, current_tier, tier_hierarchy, n_hard)`:
- The tier_hierarchy defines the order: ['localization', 'perturbation', 'dual_evidence']
- For the current tier, identify the immediately preceding tier in the hierarchy
- Hard negatives are edges in the preceding tier's positive set that are NOT in the current tier's positive set
- Convert both edge sets to Python sets of tuples for set difference operation
- Sample n_hard edges from the hard negative candidates
- Return a (2, n_hard) tensor (or None if no preceding tier or no candidates)

**7.3. Loss Functions (`src/training/losses.py`)**

Implement `bce_weighted_loss(scores, labels, pos_weight)`:
- Scores and labels are (batch_size,) tensors
- Create a per-sample weight tensor: elements where labels==1 get pos_weight, others get 1.0
- Apply `torch.nn.functional.binary_cross_entropy` with the per-sample weight tensor
- Return scalar loss

Implement `focal_loss(scores, labels, alpha, gamma)`:
- Compute standard BCE per sample (reduction='none')
- Compute p_t: predicted probability of the true class (scores where label=1, 1-scores where label=0)
- Compute focal weight: alpha × (1 − p_t)^gamma
- Return mean of (focal_weight × bce_per_sample)

---

### Step 8 — Trainer Module (`src/training/trainer.py`)

Implement a class `MEvDTrainer` with:

**`__init__`**: Takes model, data dict (containing rna_features, atac_features, prior_edges), evidence_graphs dict (tier_name → positive edge tensor), config dict, device string. Moves all data to device. Initializes best_val_aupr = -1.0 and patience counter.

**`train_stage(stage, splits)`**: Implements one complete curriculum stage:

- Create AdamW optimizer with stage.learning_rate and weight_decay from config
- Create CosineAnnealingLR scheduler with T_max = stage.n_epochs
- If stage.freeze_encoder is True: set `requires_grad = False` for all parameters in model.rna_encoder and model.atac_encoder
- Otherwise: ensure all parameters have `requires_grad = True`
- For each epoch:
  - Set model to train mode
  - Concatenate splits['train']['pos'] and splits['train']['neg'] along dim=1; create matching label tensor (ones for pos, zeros for neg)
  - Permute the concatenated edges and labels with a random permutation
  - Iterate over edge mini-batches of size config['batch_size']:
    - Call model.encode (full graph, not mini-batched) to get node embeddings h
    - Extract the batch's TF and target indices from the batch edge tensor
    - Call model.decode(h, tf_batch, target_batch) to get scores
    - Compute loss using bce_weighted_loss (or focal_loss if stage==Stage3 and config['loss_type']=='focal')
    - Call loss.backward()
    - Clip gradients (torch.nn.utils.clip_grad_norm_ with max_norm = config['clip_grad_norm'])
    - Call optimizer.step() and optimizer.zero_grad()
  - Call scheduler.step() once per epoch
  - Every 5 epochs, call evaluate_split(splits['val']) to get val AUPR; log to console; check for early stopping; save checkpoint if best
- Return history dict of {epoch: {train_loss, val_auroc, val_aupr}}

**`run_full_curriculum(splits_per_tier, stages, use_replay)`**: Orchestrates all stages:

- Initialize epoch_offset = 0
- For each stage in stages:
  - Print stage name
  - Call train_stage with the tier's splits
  - Advance the curriculum manager
  - Increment epoch_offset by stage.n_epochs
  - If use_replay and not the first stage: add replay edges from earlier tiers into the next stage's training splits at reduced weight (this requires modifying the splits dict before calling train_stage for the next stage)

**`evaluate_split(split)`**: Evaluates on a val or test split:

- Set model to eval mode, disable gradients
- Concatenate positive and negative edges from split; create labels
- Call model.encode then model.decode on ALL edges (not just a batch) — for evaluation, we want scores for all edges
- Compute AUROC and AUPR using sklearn.metrics.roc_auc_score and average_precision_score
- Compute Early Precision: sort edges by score descending; compute precision in the top K predictions where K = number of positives in split
- Compute EPR: EP / (n_positives / n_total_candidates)
- Return dict: {'auroc': float, 'aupr': float, 'early_precision': float, 'epr': float}

**`save_checkpoint(path)`**: Save model state dict and best_val_aupr to the given path using torch.save.

**`load_checkpoint(path)`**: Load a previously saved checkpoint.

---

### Step 9 — Evaluation Module (`src/evaluation/metrics.py` and `benchmarker.py`)

**metrics.py**:

Implement `compute_auroc(labels, scores)` using sklearn.metrics.roc_auc_score.

Implement `compute_aupr(labels, scores)` using sklearn.metrics.average_precision_score.

Implement `compute_early_precision(labels, scores, k=None)`: Sort scores descending; if k is None, set k = number of positive labels; return fraction of top-k predictions that are true positives.

Implement `compute_epr(early_precision, n_positives, n_total)`: epr = early_precision / (n_positives / n_total). Return EPR value.

Implement `compute_all_metrics(labels, scores, n_total_candidates)` that calls all four metric functions and returns a dict.

**benchmarker.py**:

Implement `run_all_baselines(cell_type, data, config, output_dir)` that calls each baseline wrapper sequentially, evaluates each on the same test splits as MEvD-GRN, and saves results.

Implement `compile_results_table(results_dir)` that loads all JSON result files from results_dir and compiles them into a single pandas DataFrame suitable for export to CSV for the paper table.

---

### Step 10 — Baseline Wrappers (`src/baselines/`)

**grnboost2_wrapper.py**:

Implement `run_grnboost2(scrna_path, gene_names, tf_names, n_workers, seed)`:
- Load the scRNA h5ad, apply the same normalization and log transform as in preprocessing
- Filter to only the genes in gene_names (the training gene universe)
- Call `arboreto.algo.grnboost2` with the expression DataFrame and TF name list
- Return a DataFrame with columns ['TF', 'target', 'importance'] sorted descending by importance

Implement `grnboost2_to_scored_edges(grnboost_df, gene_index, prior_edge_set)`:
- Convert the GRNBoost2 ranked edge list to a score array aligned with the prior edge candidate set
- For each candidate edge in the prior edge set, assign the GRNBoost2 importance score if found, else 0
- Return a (N_prior_edges,) score array for evaluation

**regdiffusion_wrapper.py**:

Implement `run_regdiffusion(scrna_path, gene_names, tf_mask, device)`:
- Install regdiffusion via pip if not available (catch ImportError and print instructions)
- Load and preprocess the scRNA data (same normalization as main pipeline)
- Filter to gene_names
- Create a boolean TF mask array of shape (n_genes,) — True for TFs
- Instantiate and run the RegDiffusion model to obtain a TF × Target score matrix
- Return the score matrix

**gmfgrn_wrapper.py**:

Implement `run_gmfgrn(rna_features, prior_edges, n_genes, device)` that:
- Instantiates GMFGRN from the upstream repository (https://github.com/Lishuoyy/GMFGRN)
- Clones/imports the GMFGRN code into `src/baselines/gmfgrn/` from the GitHub repository
- Trains GMFGRN using the same scRNA features and the same prior graph as MEvD-GRN for fair comparison
- Note: GMFGRN does NOT use ATAC data; set ATAC=zeros for this baseline OR train on RNA only to match its input specification
- Returns scored edge index tensor

---

### Step 11 — Main Training Script (`scripts/03_train.py`)

Implement as a command-line script with the following flow:

11.1. Parse arguments: --config (yaml path), --cell_type (K562 or ESC), --device.

11.2. Load config yaml.

11.3. Load preprocessed features from `data/processed/{cell_type}/`: rna_features_aligned.npy, atac_features_aligned.npy, gene_index.json, tf_indices.json, prior_edges.pt.

11.4. Load evidence graphs for the cell type: for each evidence tier, load the positive edge tensor from `data/processed/{cell_type}/`.

11.5. Load train/val/test splits from `data/splits/`.

11.6. Instantiate MEvDGRN with dimensions from config. Print total parameter count.

11.7. Instantiate MEvDTrainer with model, data, evidence_graphs, config, device.

11.8. Call trainer.run_full_curriculum with the three curriculum stages defined in CURRICULUM_STAGES.

11.9. After training completes, call trainer.evaluate_split on the TEST split of each evidence tier and print final results.

11.10. Save final evaluation results to `results/{cell_type}_results.json`.

---

### Step 12 — Evaluation and Baseline Scripts

**scripts/04_evaluate.py**:

Command-line script that loads a saved checkpoint and runs comprehensive evaluation:
- Evaluate on all three evidence tier test sets for the training cell type
- Evaluate on the other cell type's test sets (cross-cell-type transfer) without retraining
- Save all metrics to JSON

**scripts/05_run_baselines.py**:

Command-line script that runs all three baselines and evaluates them on the same test splits:
- Run GRNBoost2, RegDiffusion, GMFGRN
- Score each on the test edge splits using compute_all_metrics
- Save to `results/baselines/{baseline_name}_{cell_type}_results.json`

**scripts/06_ablation.py**:

Command-line script accepting a --ablation argument from the list: full_curriculum, loc_only, pert_only, dual_only, all_at_once, rna_only, concat_fusion, no_gnn, with_replay.

For each ablation, modify the config or model accordingly:
- loc_only: Run only Stage 1 (limit curriculum to localization stage)
- pert_only: Run only Stage 2 from random initialization (skip Stage 1)
- dual_only: Run only Stage 3 from random initialization
- all_at_once: Modify trainer to train simultaneously on the union of all three positive sets, equal weight; no stage-based structure
- rna_only: Before training, replace atac_features with a zero tensor of the same shape
- concat_fusion: Replace GatedFusion with a simple Linear(hidden_dim × 2 → hidden_dim) followed by GELU
- no_gnn: Remove GNNBackbone, use fused representation directly as input to BilinearDecoder
- with_replay: Run full curriculum with use_memory_replay=True

Save results to `results/ablations/{ablation_name}_results.json`.

---

### Step 13 — SLURM Scripts

**slurm/train.sh**:

SBATCH directives: job-name=MEvD-GRN, time=20:00:00, nodes=1, ntasks=1, cpus-per-task=8, mem=64G, gres=gpu:4, partition=gpu.

Body: activate conda environment, load CUDA module, call preprocessing script if `data/processed/K562/` does not exist, then call main training script.

**slurm/baselines.sh**:

SBATCH directives: time=04:00:00, gres=gpu:1, array=0-2.

Array job with three baselines indexed by SLURM_ARRAY_TASK_ID. Each task runs one baseline.

**slurm/ablation.sh**:

SBATCH directives: time=06:00:00, gres=gpu:1, array=0-8.

Array job with nine ablation variants indexed by SLURM_ARRAY_TASK_ID mapped to ablation names.

---

### Step 14 — Results and Visualization

**src/utils/visualization.py**:

Implement `plot_pr_curves(method_results_dict, tier, output_path)`:
- Takes a dict of {method_name: (labels, scores)} pairs
- Plots Precision-Recall curves for all methods on the same axes
- Adds AUPR in the legend
- Saves to output_path

Implement `plot_ablation_bars(ablation_results_dict, metric='aupr', output_path)`:
- Plots a horizontal bar chart comparing ablation variants
- Highlights the full_curriculum bar in a distinct color
- Saves to output_path

Implement `plot_training_curves(history_dict, output_path)`:
- Plots training loss and validation AUPR across all epochs, with vertical dashed lines marking stage transitions

**notebooks/02_results_analysis.ipynb**:

Create this notebook with cells that:
- Load all result JSON files
- Compile into pandas DataFrames
- Call the visualization functions to generate Figures 2, 3, 4, 5, 6 (as described in Part 15)
- Print the main results table in a format suitable for the paper

---

### Step 15 — Sanity Checks and Validation

Before running full training, verify the following:

15.1. **Data shapes**: Assert that rna_features has shape (N_G, 2) and atac_features has shape (N_G, 2) with the same N_G. Assert that prior_edges has shape (2, K) for some K.

15.2. **Evidence nesting**: Assert that |E_dual| ≤ |E_pert| and |E_dual| ≤ |E_loc|. Print the three sizes and coverage statistics.

15.3. **No data leakage**: Assert that no edge in the test split of any tier appears in the training split of the same tier. Also verify that the universal negative pool contains no edges that are positive in any tier.

15.4. **Forward pass**: Instantiate MEvDGRN with default config. Run encode on small subsets of rna_features and atac_features (first 100 rows). Assert output shape is (100, 128). Run decode with dummy index tensors. Assert scores are in [0,1].

15.5. **Loss computation**: Compute weighted BCE loss with a batch of scores and labels. Assert it is a scalar with finite, positive value.

15.6. **One training step**: Run one mini-batch forward pass and backward pass. Assert gradients are non-None and non-NaN for all model parameters.

15.7. **Checkpoint round-trip**: Save and reload a model checkpoint. Assert that the reloaded model produces identical outputs to the original on the same input.

15.8. **GPU memory check**: After moving the model and a full batch to GPU, print peak GPU memory usage. Verify it is under 8 GB per GPU to leave headroom for multi-GPU scenarios.

---

## Part 13 — Expected Output Files

After successful training, the following files should exist in `results/`:

- `results/checkpoints/best_model_Stage1_Localization.pt` — model checkpoint after Stage 1
- `results/checkpoints/best_model_Stage2_Perturbation.pt` — model checkpoint after Stage 2
- `results/checkpoints/best_model_Stage3_DualEvidence.pt` — final best model
- `results/K562_results.json` — metrics dict: {tier: {auroc, aupr, early_precision, epr}}
- `results/ESC_transfer_results.json` — cross-cell-type transfer metrics
- `results/baselines/grnboost2_K562_results.json`
- `results/baselines/regdiffusion_K562_results.json`
- `results/baselines/gmfgrn_K562_results.json`
- `results/ablations/{ablation_name}_results.json` (9 files)
- `results/figures/figure2_training_curves.pdf`
- `results/figures/figure3_pr_curves.pdf`
- `results/figures/figure4_ablation_bars.pdf`
- `results/figures/figure5_transfer_bars.pdf`
- `results/figures/figure6_runtime_scatter.pdf`

---

## Part 14 — Key Configuration Reference

This section defines exact hyperparameter values. These were chosen based on the literature (GMFGRN, RegGAIN, GNNLink) and the hardware constraints of RTX 2080 Ti. Do not change these values in the main experiment; only the ablations are permitted to modify them.

| Parameter | Value | Justification |
|---|---|---|
| n_hvg | 3000 | Standard for single-cell analysis; captures most biologically relevant variation |
| atac_window_bp | 100,000 | Covers promoters + proximal and distal enhancers; standard in ATAC-RNA integration |
| top_k_prior_targets | 500 | Balances prior graph coverage vs. computation |
| hidden_dim | 128 | Sufficient expressive power; keeps model under 500K parameters |
| n_gnn_layers | 2 | Captures 2-hop neighborhood; avoids over-smoothing |
| dropout | 0.2 | Mild regularization for small gene universe |
| Stage 1 n_epochs | 30 | Enough for convergence on large E_loc set |
| Stage 2 n_epochs | 15 | Fine-tuning; shorter since initialization is from Stage 1 |
| Stage 3 n_epochs | 10 | Minimal fine-tuning on very small E_dual set |
| Stage 1 lr | 1e-3 | Standard initial LR for AdamW on small models |
| Stage 2 lr | 3e-4 | Reduced 3× for fine-tuning |
| Stage 3 lr | 1e-4 | Minimal perturbation to avoid catastrophic forgetting |
| weight_decay | 1e-4 | Standard L2 regularization |
| batch_size (edges) | 2048 | Fits on 2080 Ti with full-graph GNN encoding |
| clip_grad_norm | 1.0 | Prevents gradient explosion across stage boundaries |
| neg_ratio (Stages 1-2) | 5:1 | Moderate class balance correction |
| neg_ratio (Stage 3) | 10:1 | Stricter balance for the very small dual-evidence positive set |
| early_stopping_patience | 10 | In units of validation checks (every 5 epochs = 50 epoch patience) |

---

## Part 15 — Timeline for 24-Hour Budget

| Clock Time | Action | Duration |
|---|---|---|
| T+0:00 | Manual data download from SC-MO-GRN-DB (follow Step 3) | 30 min |
| T+0:30 | Run `scripts/02_preprocess.py` for K562 and ESC | 45 min |
| T+1:15 | Graph construction (runs as part of preprocessing) | 30 min |
| T+1:45 | Submit main training job via SLURM: Stage 1 begins | — |
| T+4:15 | Stage 1 completes; Stage 2 begins automatically | — |
| T+5:45 | Stage 2 completes; Stage 3 begins automatically | — |
| T+6:30 | Stage 3 completes; evaluation runs automatically | — |
| T+7:00 | Submit baseline array job + ablation array job in parallel | — |
| T+9:00 | All baselines complete | — |
| T+10:00 | All ablations complete | — |
| T+10:00 | Run notebook to compile results and generate figures | 1 h |
| T+11:00 | All results locked | — |
| T+11:00–T+24:00 | Buffer: reruns, debugging, additional cell types | 13 h |

---

## Part 16 — Key References

All papers that directly inform implementation decisions. Implementer should read the abstracts of papers marked [READ].

| Reference | Relevance | Priority |
|---|---|---|
| Valensi et al., iScience 2026 (DOI: 10.1016/j.isci.2026.115323) | SC-MO-GRN-DB: all data comes from here | **[READ]** |
| Li et al. 2024, PubMed 38261340 (GMFGRN) | Direct architectural competitor + baseline | **[READ]** |
| Luyuan et al. 2024, PubMed 39387266 (RegDiffusion) | Fastest baseline; speed benchmark | **[READ]** |
| Xu et al. 2025, NAR 2025 gkaf138 (scMultiomeGRN) | RNA+ATAC fusion justification | [READ] |
| Pratapa et al. 2020, Nature Methods 17:147 (BEELINE) | Benchmarking methodology | [READ] |
| Hamilton et al. 2017 (GraphSAGE) | GNN backbone algorithm | [READ] |
| Guan et al. 2025, Advanced Science (RegGAIN) | Contrastive learning for GRN | Skim |
| Wang et al. 2021, WWW (CurGraph) | Curriculum learning for graphs | Skim |
| Chiang et al. 2019 (ClusterGCN) | Mini-batch graph training | Skim if large graph needed |
| Zeng et al. 2020 (GraphSAINT) | Mini-batch graph training (alternative) | Skim if large graph needed |
| PyTorch Geometric Heterogeneous Graph Tutorial | Bipartite SAGEConv implementation | **[READ]** |
| Loshchilov & Hutter 2019 (AdamW) | Optimizer theory | Reference only |

---

*End of MEvD-GRN Implementation Plan*

*This document contains all theory and steps needed to implement the method from scratch. Implement modules strictly in the order of Steps 1–14. Each step references the theory section justifying its design. No deviations from the architecture or hyperparameters described here should be made in the main experiment; deviations belong in the ablation study.*
