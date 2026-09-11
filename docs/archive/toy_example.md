# Toy walkthrough: what happens to the data inside MEvD-GRN

This note follows one cell type (K562) through the pipeline. The goal is
simple: see how raw RNA and ATAC files become a score for
"does TF *i* regulate gene *j*?"

We use real K562 sizes where it helps, and tiny made-up matrices where it
helps you *see* the format.

**Symbols used throughout**

| Symbol | Meaning | K562 value |
|---|---|---|
| N | genes in the final universe | 22,943 |
| T | transcription factors among them | 225 |
| C_rna | RNA cells | 953 |
| C_atac | ATAC cells (post-QC) | 1,126 |
| d | hidden dimension | 128 |
| d_sig | co-expression signature dim | 50 |
| B | edge batch size | 8,192 |

Candidate space: T × N = 225 × 22,943 ≈ **5.2 million** TF→gene pairs.
Training never scores all of them at once — only batches of edges. Gene
embeddings are still built for all N genes in one shot.

---

## 1. What the raw files look like

### 1.1 RNA expression matrix (genes × cells)

K562 path: `data/raw/single_cell/K562_scRNA/**/*.out`

SC-MO-GRN-DB stores RNA as **features × cells** (gene rows, cell columns),
tab- or comma-delimited text. A tiny toy file would look like:

```text
Symbol,cell_001,cell_002,cell_003,cell_004
GATA1,2.31,0.00,1.85,0.42
MYC,0.00,3.10,2.74,1.05
ACTB,5.20,4.80,5.01,4.95
HBB,4.10,3.90,0.00,4.50
```

So each **row** is a gene, each **column** is a cell, and the entry is
expression. Real files have ~20k+ gene rows and hundreds of cell columns
(953 for K562 DS019).

Gene symbols in RNA are mixed-case (`Gata1`); network files use UPPERCASE
(`GATA1`). Everything is uppercased when matching.

### 1.2 ATAC peak matrix (peaks × cells)

K562 uses a 10x multiome `.h5` (`data/raw/single_cell/K562_multiome_h5/*.h5`).
Other datasets use text peak matrices. Same idea either way: rows are
**peaks** (genomic intervals), columns are cells, values are accessibility
counts.

Toy text form:

```text
peak_coord	cell_001	cell_002	cell_003
chr1:1000-1500	0	3	1
chr1:2000-2500	5	0	2
chr11:5248000-5249000	0	1	4
chr8:127735000-127736000	2	0	1
```

A peak like `chr11:5248000-5249000` is just an interval of open chromatin.
It is **not** yet attached to a gene name.

### 1.3 Reference networks (labels)

```text
Source	Target	Relationship
GATA1	HBB	Activation
MYC	CDK4	Activation
GATA1	ALAS2	Activation
```

These are the **supervision** edges (localization / perturbation /
dual-evidence). They are *not* the message-passing graphs.

---

## 2. RNA: genes × cells → one vector per gene

The model treats **genes as nodes**, not cells. So the big cell matrix is
compressed.

**Steps**

1. Load genes × cells, transpose → cells × genes  
2. QC (drop sparse cells/genes), CP10K normalize, log1p  
3. For each gene, keep three numbers across cells  
4. Also build a co-expression signature (length 50)

**Toy conversion for `GATA1`** (after log-norm, 4 cells):

```text
values:           [0.8,  0.0,  1.2,  0.5]
mean:             0.625
variance:         0.187
detection rate:   0.75     (3 of 4 cells non-zero)

→ rna_features row:  [0.625, 0.187, 0.75]
```

Do that for every gene:

| Stage | Shape (general) | Shape (K562) |
|---|---|---|
| Raw file on disk | genes × C_rna | ~20k+ × 953 |
| After load + transpose | (C_rna, N_raw) | (953, N_raw) |
| After QC + normalize | (C_rna, N') | (953, N') |
| `rna_features` | (N, 3) | **(22,943, 3)** |
| `signatures` | (N, d_sig) | **(22,943, 50)** |

The 953 cells are gone after this. Only per-gene stats and signatures remain.
The signature is set up so that `signatures[i] · signatures[j]` ≈ how
correlated genes i and j are across cells.

---

## 3. ATAC: peaks → gene activity + openness

Peaks are intervals. Genes have a transcription start site (TSS). We ask:
which peaks fall within **±100 kb** of this gene's TSS?

That builds a sparse peak→gene **incidence** matrix (0/1). Then:

```text
gene_activity = peak_matrix @ incidence
                (C_atac, P)  ×  (P, N)  →  (C_atac, N)
```

From gene activity, keep per gene:

- mean and variance across ATAC cells → **ATAC features**
- how many peaks sit near the TSS → **openness** in [0, 1]

**Toy for `HBB`**

```text
peaks near HBB TSS:           3
gene activity (3 cells):      [2.1,  0.0,  1.4]
mean / var:                   → atac_features row [mean, var]
openness:                     log1p(3) / max_over_all_genes   ∈ [0, 1]
```

| Stage | Shape (general) | Shape (K562) |
|---|---|---|
| Peak matrix | (C_atac, P) | (1,126, many peaks) |
| Incidence | (P, N) | sparse 0/1 |
| Gene activity | (C_atac, N) | (1,126, 22,943) |
| `atac_features` | (N, 2) | **(22,943, 2)** |
| `openness` | (N,) | **(22,943,)** |

Again the cells disappear. Each gene now has a tiny ATAC summary plus one
number saying "how open is this locus?"

---

## 4. Gene universe + labels

Intersect RNA genes with genes that appear in any reference network:

```text
universe = RNA genes ∩ network genes     →  N = 22,943
TFs      = sources in networks ∩ universe →  T = 225
```

Evidence edges that fall inside the universe become labeled positives, e.g.
`evidence_localization.pt` with shape `(2, E)`.

---

## 5. Two graphs for message passing

These are **not** the labels. They only tell the GNN who should talk to whom.
Both graphs sit on the same N gene nodes. Message passing runs on **both**.

### Co-expression graph (undirected)

For each gene, take its **20** most co-expressed neighbors (signature dot
product). Add edges both ways.

```text
coexpr_edges ≈ (2, N × 20) ≈ (2, ~460,000)
```

Point: genes that move together share regulatory context.

### TF-candidate graph (directed TF → gene)

For each of the 225 TFs, keep the top **500** genes that are (a) accessible
(`openness > 0`) and (b) most co-expressed with that TF. Known positive
edges are excluded so labels cannot leak into the graph.

```text
tf_candidate_edges ≈ (2, T × 500) ≈ (2, ~112,000)
```

Point: a sparse "plausible regulation" skeleton, without scoring all 5.2M pairs.

---

## 6. What preprocess writes to disk

After `scripts/02_preprocess.py`, training loads from
`data/processed/K562/`:

| File | What it is | Shape |
|---|---|---|
| `rna_features_aligned.npy` | mean, var, detection | (22,943, 3) |
| `atac_features_aligned.npy` | mean, var | (22,943, 2) |
| `rna_signature.npy` | co-expression signatures | (22,943, 50) |
| `openness.npy` | locus openness | (22,943,) |
| `coexpr_edges.pt` | gene–gene kNN | (2, ~460k) |
| `tf_candidate_edges.pt` | TF→candidate | (2, ~112k) |
| `evidence_*.pt` | labeled positives | (2, E_tier) |
| `gene_index.json` | `"GATA1" → 1234` | — |

Those matrices are what enter the neural net.

---

## 7. Encoders (still separate RNA and ATAC tracks)

RNA and ATAC are **not** fused early. Each gene gets its own RNA embedding
and ATAC embedding.

### RNA encoder (2-layer MLP)

```text
(N, 3)
  → Linear(3 → 64) → LayerNorm → GELU → Dropout
  → Linear(64 → 128) → LayerNorm
  → (N, 128)

K562:  (22,943, 3)  →  (22,943, 128)
```

### ATAC encoder (one projection)

```text
(N, 2)
  → Linear(2 → 128) → GELU → LayerNorm
  → (N, 128)

K562:  (22,943, 2)  →  (22,943, 128)
```

After encoding:

```text
h_rna:   (22,943, 128)
h_atac:  (22,943, 128)
```

Each row is still one gene — just in a richer space.

---

## 8. GNN: message passing on both graphs

There are two GraphSAGE backbones with separate weights:

- `gnn_rna` updates `h_rna`
- `gnn_atac` updates `h_atac`

Each has **2 layers**. In every layer, each gene aggregates from **both**
graphs; the two messages are added:

```text
for each GNN layer:
    msg_co = SAGEConv(h, coexpr_edges)           # (N, 128)
    msg_tf = SAGEConv(h, tf_candidate_edges)     # (N, 128)
    h      = Dropout(GELU(LN(msg_co + msg_tf)))  # (N, 128)
```

Shape never changes:

```text
(N, 128) → layer 1 → (N, 128) → layer 2 → (N, 128)
```

After both GNNs:

```text
h_rna:   (22,943, 128)    expression context + neighbors
h_atac:  (22,943, 128)    accessibility context + neighbors
```

A gene's vector is no longer only its own mean/var. It also carries
information from co-expressed genes and TF-candidate neighbors.

---

## 9. Decoder: one score per TF → gene pair

Training samples a batch of B = 8,192 edges `(tf_idx, target_idx)`.

Lookup:

```text
h_tf_rna   = h_rna[tf_idx]           # (B, 128)
h_tg_rna   = h_rna[target_idx]       # (B, 128)
h_tg_atac  = h_atac[target_idx]      # (B, 128)
coexpr     = signatures[tf] · signatures[tg]   # (B,)
open_j     = openness[target_idx]              # (B,)
```

Role-aware decoder:

```text
1. a_j   = sigmoid(Linear(h_tg_atac))     # (B, 128)  accessibility gate
2. z_tg  = h_tg_rna * a_j                 # (B, 128)  gated target RNA
3. raw   = sum( (h_tf_rna @ W) * z_tg )   # (B,)
4. logit = raw + w_coexp * coexpr + w_open * open_j + b
5. prob  = sigmoid(logit)                  # (B,) in (0, 1)
```

One number per candidate edge: how likely is this regulation?

---

## 10. Dimensional analysis (full pipe)

### General

| Stage | Shape |
|---|---|
| Raw RNA file | genes × C_rna |
| RNA after load | (C_rna, N_raw) |
| RNA features | (N, 3) |
| Signatures | (N, d_sig) |
| Peak matrix | (C_atac, P) |
| ATAC features | (N, 2) |
| Openness | (N,) |
| Co-expression edges | (2, ~ N × 20) |
| TF-candidate edges | (2, ~ T × 500) |
| After RNA / ATAC encoders | (N, d) and (N, d) |
| After GNN (2 layers) | still (N, d) and (N, d) |
| Edge batch | B pairs |
| Decoder logits | (B,) |

### K562 numbers

```text
RNA:   genes×953 cells  →  (22943, 3)  → encoder → (22943, 128)
ATAC:  peaks×1126 cells →  (22943, 2)  → encoder → (22943, 128)
                                              ↓ GNN × 2  (both graphs)
                                      (22943, 128) each channel
                                              ↓ lookup batch of edges
                        h_tf, h_tg_rna, h_tg_atac: (8192, 128)
                                              ↓ role-aware decode
                                      logits: (8192,)
```

Pattern: after preprocessing almost everything is **genes × 128**, until
the last step picks pairs and collapses to one score per pair.

---

## 11. One sentence version

Compress RNA and ATAC into a few numbers per gene, embed those into 128-d
vectors, let each gene talk to co-expression and TF-candidate neighbors
through a GNN on **both** graphs, then ask — for a specific TF and target —
whether the TF's RNA state matches a target whose RNA state has been gated
by how open its chromatin is.
