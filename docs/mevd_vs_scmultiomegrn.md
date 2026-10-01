# MeVD-GRN vs scMultiomeGRN: comparative study and novelty audit

Written 2026-10-01 for the paper's Related Work and Contributions sections.

**Scope.**
- MeVD-GRN as implemented on branch `leak-fix-and-benchmarks`: `src/models/`, `src/training/`, `src/data/`, `configs/`, and the result JSONs in `results/ablations/`.
- scMultiomeGRN (Xu *et al.*, *NAR* 53(5):gkaf138, 2025) from three sources: the paper (`docs/reference/scMultiomeGRN.pdf`), the official Zenodo code (`~/.cache/mevd_scmgrn/official/ScmultiomeGRN-main/`), and our adapter (`src/baselines/scmultiomegrn_wrapper.py`).
- Where the code and a document disagree, **the code wins**. Every such disagreement is listed in Section 8.

**Naming.** Code and paper spell the model "MEvD-GRN"; newer docs write "MeVD-GRN". This note uses MeVD-GRN.

**What "MeVD-GRN" means here.** It means the paper's headline configuration, recovered from the result files rather than the shipped default config. That configuration is `scripts/06_ablation.py --ablation all_at_once` on `configs/sweep/k562_fm_h384_l2.yaml`. The result files are `results/ablations/all_at_once_fm_h384l2_seed4{2..6}_K562.json` (3,257,479 params, empty `relation_weights`). Concretely:

| Setting | Headline value | Shipped default (`configs/default.yaml`) |
|---|---|---|
| Protocol | `all_at_once`: one stage on loc ∪ pert, 45 epochs, lr 1e-3 (hard-coded in `build_all_at_once_stage`) | `sequential` (loc 30 ep, then pert 15 ep) |
| Hidden / GNN layers | 384 / 2 | 128 / 2 |
| Geneformer FM embedding | on | off |
| Relations | 2 (co-expression kNN + TF-candidate), unweighted `sum` | same |
| Motif relation | **off** | off |
| Gated relation combiner | **off** (`combine_mode: sum`) | off |
| Hard negatives | **off** (forced off by `all_at_once`) | on |
| Memory replay | **not applicable** (single stage) | on |
| Decoder | role-aware (ATAC-gated bilinear + linear co-expression/openness terms) | same |

So several components that the docs present as MeVD-GRN's contributions are not in the model that produces the headline numbers. These are the curriculum, replay, hard negatives, the motif graph and the gated combiner. This shapes the novelty analysis below.

---

## 1. Summary

Both methods encode RNA-derived and ATAC-derived node features in separate modality-specific graph encoders and train with a class-weighted BCE link-prediction loss. The shared idea is GraFRank's modality-specific aggregation, which scMultiomeGRN took over directly (Sankar *et al.*, WWW 2021). They differ in almost everything else:

- **Task.** scMultiomeGRN predicts an *undirected TF–TF* graph. MeVD-GRN predicts *directed TF→any-gene* edges (~225 TFs × ~23k genes in K562).
- **Labels.** scMultiomeGRN's labels are motif hits (FIMO/HOCOMOCO) inside promoters, computed from the *same* scATAC data it takes as input. MeVD-GRN's labels are experimental: ChIP localization, knockdown/knockout perturbation, and their intersection, all from SC-MO-GRN-DB.
- **Input graph.** scMultiomeGRN passes messages over the *training-positive adjacency*, with attention over CNN-encoded joint-expression-histogram edge features, then fuses modalities per node with cross-modal attention. MeVD-GRN passes messages with mean-GraphSAGE over label-free structural graphs (co-expression kNN, accessible and co-expressed TF-candidates) and keeps the modalities apart until the decoder, where the target's ATAC embedding gates its RNA embedding.
- **Features.** scMultiomeGRN keeps per-cell feature vectors (MAESTRO RP per cell; GRNBoost2 importance vectors). MeVD-GRN collapses each modality to a few pseudobulk statistics per gene and relies mostly on a frozen Geneformer gene-token embedding.
- **Protocol.** MeVD-GRN's distinctive protocol is a hierarchy-consistent global edge split across nested evidence tiers, with the intersection tier held out as a zero-shot test. It also includes a small multi-cell-type joint-training and transfer study.

Most individual architectural pieces of MeVD-GRN are borrowed or adapted. What is new is mainly in the supervision and evaluation protocol and in the empirical findings, and some of those claims carry open validity risks (Section 6).

---

## 2. Side-by-side comparison

| Dimension | scMultiomeGRN (paper + official code) | MeVD-GRN (code, headline config) |
|---|---|---|
| **Task** | Semi-supervised link completion on an **undirected TF–TF** graph. Nodes are TFs only; non-TF targets are excluded, which the authors list as a limitation. Scores are symmetrised: code `(out + out.T)/2`. | Supervised, **directed TF→target** link prediction over the TF × all-genes space: K562 has 22,943 genes and ~225 TFs. Direction enters only at the decoder; GNN graphs are made bidirectional (`gnn.py: to_bidirectional`). |
| **Inputs and modality pairing** | scATAC + scRNA per cell type. Lung (main benchmark, GSM4508936) is from the sci-ATAC-seq3 atlas; per the GEO description (see `docs/experiments/scmultiomegrn_generalization.md` §2), its "scRNA" assay is a gene-body *accessibility* count matrix, so both modalities come from ATAC. Hematopoiesis and PBMC use separately profiled RNA and ATAC (unpaired). | Per-gene features for one cell type. **K562 RNA (DS019, 953 cells) and ATAC (DS025 multiome) come from different datasets**, so they are not cell-paired (`configs/k562.yaml`: "gene-level aggregation makes the RNA/ATAC pairing irrelevant"). Macrophage (DS026) and MCF7 (DS027) are paired multiome, but only their gene-level summaries are used. |
| **Gene/TF features** | ATAC: MAESTRO "Enhanced" RP score of each TF **in every cell** (dim = #cells). RNA: GRNBoost2 TF→gene importance vector (dim = #genes). Both modalities are fed to the encoder. | RNA: `[mean, var, detection]` (3-d) after CP10K + log1p. ATAC: `[mean, var, detection]` of RP-weighted gene activity (exponential TSS decay, 10 kb, ±100 kb window) (3-d). A separate 4-d "locus-shape" descriptor (mean/max RP of the top-10 peaks, mean signed distance, peak count) is used only at the decoder. **FM:** Geneformer-V2-104M static input-token embedding (768-d, frozen, 77% gene coverage in K562), passed through a 2-layer MLP and added to the RNA channel. A 50-d co-expression signature (truncated SVD, s_i·s_j ≈ corr) is used only for graph construction and the decoder's co-expression term. It is **not** an encoder input. |
| **Prior / input graphs** | A single graph: the **training-positive adjacency** plus self-loops, GCN-normalised, i.e. the known edges themselves. Edge attribute is a 16×16 log joint histogram of the two TFs' expression (CNNC-style, Eqs 1–2). | Two relations, rebuilt every layer by per-relation SAGEConv: (1) undirected co-expression kNN, k = 20, 871,450 edges; (2) directed TF→candidate graph: top-500 targets per TF that are proximally accessible (max RP > 0.1) and most co-expressed, 112,500 edges, **excluding every known positive of every tier and split** (`prior_exclude_positives: true`). Optional third relation (off in the headline): 28,204 JASPAR PWM-hit edges. Because the scan is restricted to TF-candidate edges, this relation is a *subset* of relation (2). |
| **Encoder** | Per modality: 2 GraFrankConv layers. Each message = Linear([x_j ‖ CNN(edge histogram)]); neighbour attention = LeakyReLU softmax; update = lin_l(agg) + lin_r(self). Width 512. Then one GCN layer (512→256) over the normalised adjacency. | Per modality: 2-layer MLP encoder (LayerNorm, GELU), then its **own 2-layer relational GraphSAGE tower** (one mean-SAGEConv per relation, summed, LayerNorm + GELU + dropout, residual from layer 2). Width 384. The RNA and ATAC towers share the graphs but not their weights. |
| **Cross-modal fusion** | **Per-node cross-modal attention** (β_k = softmax over modalities of a 2-layer tanh MLP; h = Σ β_k W z_k), applied *before* the decoder. | **No fusion in the encoder.** The only interaction is at the decoder: a_j = σ(W_g h^atac_j) gates the target's RNA embedding (z_j = h^rna_j ⊙ a_j). The FM embedding is added into the RNA channel. The legacy `gated`/`concat` fusion variants exist only as ablations. |
| **Decoder** | Paper Eq. 11 says an inner product σ(hhᵀ). **The code instead uses a concatenation MLP** (`InnerProductDecoder`: Linear(512→256) → Linear(256→1) → sigmoid, with no nonlinearity in between) over [h_u ‖ h_v], symmetrised. | Asymmetric full bilinear (h^rna_i)ᵀ W z_j (RESCAL-style, not DistMult), plus w_c · coexpr(i,j), plus Linear(open_j ∈ ℝ⁴), plus bias. An optional edge-MLP variant exists (off). |
| **Supervision / labels (ground truth)** | "Initial adjacency" from each cell type's own scATAC: (1) keep peaks in > 10% of cells; (2) FIMO over HOCOMOCO v11 (code: `--thresh 1e-4`, then p ≤ 1e-6; the paper says FIMO v5.1 with `-thresh 1e-6`); (3) edge TF_a–TF_b if a TFBS of b lies in a's promoter (TSS ± 2 kb). Motif- and ATAC-derived, and circular with the ATAC input. The ground-truth construction follows DeepTFni (Li *et al.*, *Nat Mach Intell* 2022). | SC-MO-GRN-DB (Valensi *et al.*, *iScience* 2026) K562 reference networks: localization (RN117, ChIP-seq), perturbation (RN118, KO) and dual evidence (RN119, from a different study). About 818k / 177k / ~8k in-universe edges; 100% of dual is inside perturbation and 77.5% inside localization. Macrophage and MCF7 have localization only. |
| **Training objective and negatives** | BCE-with-logits with a VGAE-style `pos_weight` and `norm` (Kipf & Welling). The code applies it to the *post-sigmoid* decoder output, so the sigmoid is applied twice (quirk in `model_helper.loss_rec`). Train negatives = **all** remaining non-edges; val/test negatives are drawn 1:1. Adam, lr 1e-5, ≤ 2000 full-batch epochs, early stopping on val loss (patience 100), ReduceLROnPlateau. | Weighted BCE-with-logits, pos_weight = n_neg/n_pos recomputed every epoch. 5 uniform-random negatives per positive, resampled each epoch from a 1M-pair pool. Since 2026-09-30 that pool excludes every val/test negative (`restrict_negative_pool`); **all reported numbers predate this fix**. AdamW (wd 1e-4), cosine LR, batch 8,192, early stopping on val AUPR (10 checks × 5 epochs). The GNN re-encodes the full graph for every mini-batch. |
| **Curriculum / multi-evidence handling** | None. One flat label set per cell type. Robustness comes from **ensembling**: 10 runs, each binarised at the test-median score, keeping edges present in ≥ 6/10 runs. | Tier-aware machinery is available: sequential loc→pert curriculum, replay of 10% of earlier-tier positives at loss weight 0.1, and "hard" negatives = lower-tier train positives absent from the current tier. **The headline uses none of these**: `all_at_once` pools loc ∪ pert and disables hard negatives. The part of tier handling the headline does keep is the *protocol*: the global split and the dual tier held out as a zero-shot test. |
| **Split protocol and leakage controls** | 10-fold KFold over positive TF pairs (seed 666). Test = one fold + an equal number of random non-edges. Training adjacency = the remaining positives, which are also the message-passing graph. The paper's reported AUROC/AUPR are the full-matrix `all` set, which includes training edges (see our generalization doc §3), i.e. transductive. There is no TF-held-out protocol. | A **global, hierarchy-consistent edge split**: each unique (TF, gene) pair is assigned to train/val/test once across all tiers (70/15/15), which blocks cross-tier nesting leaks. The code checks this (`verify_no_cross_tier_leakage`, plus an assertion in `build_all_at_once_stage`). Hard negatives use only train splits (fix for the earlier 0.33-AUROC bug). Val/test negatives are excluded from the training pool (2026-09-30 fix). **Open issue:** the TF-candidate and motif graphs exclude *all* positives, including test ones, so absence from a graph is label-correlated. The split is by random edge, so TFs are shared across splits. |
| **Evaluation metrics** | AUROC; AUPR computed as trapezoidal `auc(recall, precision)` (not AP); accuracy at the test-median threshold; KL divergence between positive and negative score distributions. Reported on the full N × N matrix, mean ± std over folds. | AUPR (sklearn average precision), AUROC, early precision at k = #positives, EPR relative to the eval set. Test sets use 1:5 positive:negative, so the random AUPR baseline is ≈ 0.167. The dual tier is scored on val + test pooled. Headline: 5 seeds (42–46), std ≤ 0.0017. |
| **Scale: datasets** | Fetal lung: 9 cell types, 544–43,289 cells each, 72,622 in total; TF graphs of 542–9,609 "initial" interactions. Hematopoiesis: HSC/CMP/GMP. PBMC (5 cell types) and AD microglia are used qualitatively only. | K562 is the only cell type with all three tiers. Macrophage (16,202 genes, 22 TFs, 112k loc edges) and MCF7 (16,730 genes, 250 TFs, 1.2M loc edges) are used for transfer and joint training, localization only. |
| **Scale: parameters** | ~4.08M in our TF→gene adaptation (`paper/main.tex`), per trained tier. In the official code, the first-layer width equals #cells (ATAC RP) and #genes (GRNBoost2), so the size grows with the dataset. By our arithmetic from the layer definitions (not measured), the ATAC branch alone exceeds 40M parameters for the 43k-cell stromal type. | 315,015 (base, h128, no FM); 3,257,479 (headline h384/l2 + FM); 4,438,663 with the motif relation. |
| **Compute** | ≤ 2000 epochs per model × 10 runs per cell type. Our K562 adaptation took 5h58m on 4 × 2080 Ti for 2 tiers. | 45 epochs in total; 16.6 min for the base model on an RTX 4060. Not compute-matched to the baseline. |
| **Code availability** | Zenodo 10.5281/zenodo.14848389. No license declared, and no processed lung/hematopoiesis data or splits are published. | This repository (branch `leak-fix-and-benchmarks`). No public release or license yet. |

---

## 3. Provenance of each MeVD-GRN component

Labels: **Borrowed** = used essentially as published. **Adapted** = a known idea with a specific modification. **Novel** = no precedent found in the local knowledge base (the ~130-paper survey in `docs/reference/`) or in targeted searches. "Novel" is always relative to that search, not a proof of priority.

| # | Component | What it does | Source idea (venue, year) | What we changed | Label | In headline? |
|---|---|---|---|---|---|---|
| C1 | Directed TF→gene link prediction over TF × all genes | Scores any (TF, gene) pair | VGAE link prediction (Kipf & Welling, arXiv 2016); GENELink / GNNLink family of supervised TF→target link prediction (*Bioinformatics* 2022; *Brief Bioinform* 2023) | Genome-scale candidate space (~5.2M pairs) instead of BEELINE's TFs+500/1000 | Borrowed | yes |
| C2 | RNA pseudobulk features `[mean, var, detection]` | Per-gene expression summary | Standard scRNA QC statistics | — | Borrowed | yes (adds ~nothing once FM is present: `fm_only` 0.963 = `with_fm` 0.963) |
| C3 | Co-expression signature (truncated SVD, s_i·s_j ≈ Pearson r) | Correlation without an N×N matrix | Standard low-rank correlation factorisation; co-expression as a regulatory cue (GENIE3, *PLoS ONE* 2010; GRNBoost2, *Bioinformatics* 2019) | Used for the kNN graph, the TF-candidate ranking and a decoder term | Borrowed | yes |
| C4 | RP-weighted ATAC gene activity (exponential TSS decay) | Peak→gene aggregation | MAESTRO / BETA regulatory potential (Wang *et al.*, *Genome Biol* 2020); Signac gene activity. **scMultiomeGRN uses MAESTRO RP too, so this is shared.** | Collapsed to per-gene mean/var/detection instead of per-cell RP; decay 10 kb instead of MAESTRO's 1 kb (`genedistance=1000` in scMultiomeGRN) | Borrowed | yes |
| C5 | 4-d locus-shape descriptor (mean/max RP, mean signed distance, peak count) | Target-openness term in the decoder; gate for the TF-candidate graph | Descriptive statistics over RP weights | Our own engineering; small | Adapted (minor) | yes |
| C6 | Geneformer gene-token embedding as an auxiliary node feature | Pretrained gene prior | Geneformer (Theodoris *et al.*, *Nature* 2023); **scRegNet** feeds frozen Geneformer/scBERT/scFoundation embeddings into a GNN for supervised TF→target prediction (Kommu *et al.*, *Bioinformatics* 2025, ISMB) | Uses the *static input-token* table of the V2-104M checkpoint (a context-free gene-identity vector), not cell-contextual embeddings mean-pooled over cells as in scRegNet; added into the RNA tower before message passing | Borrowed (closest precedent: scRegNet) | yes |
| C7 | Separate RNA and ATAC encoders + GNN towers | Keeps modalities apart | Modality-specific neighbour aggregation: GraFRank (Sankar *et al.*, WWW 2021), used for GRNs by **scMultiomeGRN** (*NAR* 2025) | Mean-SAGE towers over shared structural graphs instead of edge-feature attention over the training adjacency; no fusion before the decoder | Adapted from scMultiomeGRN / GraFRank | yes |
| C8 | Relational GraphSAGE (one SAGEConv per relation, summed) | Multi-relation message passing | GraphSAGE (Hamilton *et al.*, NeurIPS 2017); per-relation sum as in R-GCN (Schlichtkrull *et al.*, ESWC 2018) | Plain composition | Borrowed | yes |
| C9 | Co-expression kNN graph | Co-regulation modules as structure | Co-expression gene graphs (WGCNA tradition; GMFGRN, *Brief Bioinform* 2024) | — | Borrowed | yes |
| C10 | TF-candidate graph (proximally accessible × top co-expressed, per TF) | A label-free, TF-specific "plausible target" relation | Accessibility-constrained candidate sets, e.g. CellOracle's base GRN (*Nature* 2023) and PECA-style models (Duren *et al.*, *PNAS* 2017) | A per-TF top-k co-expression ranking inside accessible genes, used **as a message-passing relation rather than as labels or a hard mask**. Excludes all positives, which is a leakage risk (Section 6) | Adapted | yes |
| C11 | Structural (non-label) input graphs instead of the training adjacency | The encoder never sees training edges as structure | Contrast with scMultiomeGRN, DeepTFni, GENELink and scRegNet, which all message-pass over known edges; prior-graph GNNs such as GMFGRN use co-expression | A design choice, not a new mechanism; weakened by the positive exclusion in C10 | Adapted (differentiator vs scMultiomeGRN) | yes |
| C12 | Role-aware decoder: ATAC-gated target RNA state, asymmetric bilinear, plus co-expression and openness terms | ATAC acts as a multiplicative precondition | Biology: regulation ∝ TF expression × RE accessibility (PECA; Duren *et al.*, *PNAS* 2017); hard accessibility masks (CellOracle). Mechanism: sigmoid gating (GLU; Dauphin *et al.*, ICML 2017). Bilinear scoring: RESCAL (Nickel *et al.*, ICML 2011) | A learned soft per-dimension gate from the ATAC tower applied to the target only; explicit asymmetric TF/target roles | Adapted. The specific formulation is new to our knowledge, but ablations show **no accuracy gain** over gated/concat fusion (0.881 vs 0.883 / 0.890 dual AUPR, small model, n = 1) | yes |
| C13 | Optional edge-MLP on [coexpr, openness] | Nonlinear edge term | MPNN edge features (Gilmer *et al.*, ICML 2017); scMultiomeGRN's edge features | Five scalars, not histograms | Borrowed | no (neutral, 0.882) |
| C14 | Learned gated relation combiner (per-layer softmax over relations) | Interpretable relation weights | Semantic-level attention over relations/meta-paths (HAN; Wang *et al.*, WWW 2019); scMultiomeGRN's own softmax over modalities. `docs/citations.md` cites the Relational Graph Transformer (arXiv 2025), which is a looser match | A global scalar per (layer, relation), not node-dependent | Borrowed (simplified) | **no**. With 2 relations and FM at h384 it learns near-uniform weights (0.43–0.58) and gives the same accuracy (0.962 vs 0.963 dual AUPR). The 3-relation run was never done. |
| C15 | JASPAR PWM motif relation | Sequence-backed TF→gene edges | Motif-in-accessible-peak priors (scMultiomeGRN/DeepTFni ground truth; CellOracle; SCENIC+); FIMO (Grant *et al.*, *Bioinformatics* 2011); JASPAR 2024 | A NumPy log-odds scan (hit at ≥ 75% of max score, no p-values), limited to TF-candidate pairs, so it is a subset of C10; used as a relation, not as labels | Adapted | **no** (`use_motif: false` everywhere; `motif_graph_only` alone reaches 0.955 dual AUPR, sequential protocol) |
| C16 | Global hierarchy-consistent edge split | Prevents a test edge in one tier from being a training edge in another | Standard deduplication before splitting; OGB / link-prediction split hygiene | Applied to nested evidence tiers, with runtime leak assertions | Adapted (standard hygiene, carefully applied) | yes |
| C17 | Evidence-tier curriculum (loc → pert) | Easy/broad → specific supervision | Curriculum learning (Bengio *et al.*, ICML 2009); tier definitions from SC-MO-GRN-DB | Stages defined by *experimental evidence type* rather than a difficulty score | Novel as an application, to our search | **no** (headline is `all_at_once`) |
| C18 | Memory replay across tiers | Limits forgetting of the earlier tier | Experience replay (Rolnick *et al.*, NeurIPS 2019); catastrophic interference (McCloskey & Cohen 1989). **In GRN inference, LINGER already uses lifelong learning** (EWC from bulk atlas to single cell; Yuan & Duren, *Nat Biotechnol* 2024) | Replay of 10% of earlier-tier positives at weight 0.1 | Borrowed | **no** |
| C19 | Tier-difference hard negatives ("in the weaker tier, absent from the stronger") | "Binds but does not regulate" contrast | Hard-negative mining (FaceNet, CVPR 2015); GRN hard negatives from other TFs' targets (HGATLink, *BMC Bioinformatics* 2025) | Defined by evidence-tier set difference | Novel construction, to our search | **no**. There is no dedicated ablation, so its effect has never been isolated. |
| C20 | Zero-shot evaluation on the intersection (dual-evidence) tier | Tests transfer to the highest-confidence labels | Intersection gold standards (ChIP ∩ perturbation) are used *for evaluation* by scMTNI (Zhang *et al.*, *Nat Commun* 2023) and McCalla *et al.* (*G3* 2023). The dual tier itself is defined by SC-MO-GRN-DB | Train only on the constituent tiers and hold the intersection out, under the global split | Adapted (the protocol framing is new; intersection gold standards are not) | yes |
| C21 | Multi-cell-type joint training + zero-shot transfer to a third cell type | Cross-cell-type generalisation | Multi-task learning (Caruana 1997); cross-cell-type blind tests (GRNFormer, Hegde & Cheng, *Bioinformatics* 2026, RNA-only); unseen-dataset GRN benchmarks with scFMs (Qi *et al.*, arXiv 2605.08128, 2026); transfer from external atlases (LINGER; CellPolaris, *Adv Sci* 2026) | Supervised, ChIP-labelled, multiome gene-level features across K562/MCF7/Macrophage from SC-MO-GRN-DB | Adapted. It is new on SC-MO-GRN-DB as far as we can tell (BEAR-GRN, from the same lab, benchmarks unsupervised methods rather than training jointly) | separate study (h128, sequential/loc-only, pre-fix) |
| C22 | Weighted BCE, AdamW, cosine LR, early stopping on val AUPR | Optimisation | Standard | — | Borrowed | yes |
| C23 | Metrics: AUPR, AUROC, EP, EPR | Evaluation | Saito & Rehmsmeier (*PLoS ONE* 2015); BEELINE EPR (*Nat Methods* 2020) | EPR normalised to the eval set (bug fix) | Borrowed | yes |

---

## 4. MeVD-GRN's novel contributions, ranked

The ranking weighs how defensible the novelty is against how much the headline evidence depends on it. None of the items is an architectural invention on the scale of a new layer type. The independent review reached the same view (`docs/critical_review_independent.md` §9).

1. **An evidence-tier-aware supervision and evaluation protocol for multiome GRN inference: a hierarchy-consistent split, training on the constituent tiers, and zero-shot scoring on their intersection.**
   - *Evidence it is new.*
     - In the local 130-paper survey and targeted web searches, we found no method that *trains* on separate experimental evidence types of one reference resource and holds their intersection out.
     - The surveyed methods pool labels into one positive set, including scMultiomeGRN, scRegNet, GENELink, EpiAwareNet and LINGER.
     - SC-MO-GRN-DB provides the tier annotations but no training protocol.
   - *Closest prior work.*
     - scMTNI (*Nat Commun* 2023) and McCalla *et al.* (*G3* 2023) evaluate against ChIP, perturbation and ChIP ∩ perturbation gold standards, but only as separate *evaluation* sets.
     - DoRothEA's confidence levels are another tiered label resource.
     - Curriculum learning (Bengio 2009), multi-fidelity and weak-to-strong supervision are the generic ML precedents.
   - *Caveats.*
     - The headline configuration does not use the curriculum itself (`all_at_once`), so the claim must be about the protocol, not about "curriculum training improves accuracy".
     - The zero-shot dual set is a subset of held-out perturbation edges (Section 6, risk R4).

2. **An empirical finding: a frozen FM gene embedding is the largest in-domain gain but *hurts* single-source cross-cell-type transfer, and joint multi-cell-type training shrinks that penalty.**
   - *Numbers.*
     - In-domain K562: dual AUPR 0.881 → 0.951 (`fm_only`).
     - Single-source transfer: K562→Macrophage 0.607 → 0.447 and K562→MCF7 0.758 → 0.597 once FM is added.
     - Joint training (Macrophage held out): 0.743 without FM, 0.703 with FM.
     - With K562 held out the sign flips: 0.704 vs 0.744.
   - *Evidence it is new.*
     - scRegNet and the FM-benchmark papers report in-domain gains only.
     - arXiv 2605.08128 argues that scFM features *should* generalise, and our result is a counterpoint in a supervised multiome setting.
     - We found no report of FM embeddings degrading supervised cross-cell-type GRN transfer.
   - *Closest prior work.* scRegNet (*Bioinformatics* 2025); Qi *et al.* (arXiv 2605.08128); GRNFormer's cross-cell-type test (*Bioinformatics* 2026).
   - *Caveats.*
     - All n = 1, at h128.
     - Run before the eval-negative leak fix.
     - Macrophage and MCF7 have localization labels only.
     - The effect is holdout-dependent.
     - The obvious mechanism is gene-identity memorisation (Section 6, risk R6). That makes this a finding about a failure mode rather than a method contribution.

3. **Multi-cell-type joint training with zero-shot transfer to a third, unseen cell type on SC-MO-GRN-DB.**
   - *Numbers.* It beats single-source transfer in both holdouts that have a baseline: 0.743 vs 0.607 and 0.837 vs 0.758 AUPR.
   - *Evidence it is new.* Joint training plus third-cell-type transfer on SC-MO-GRN-DB is unreported as far as we can tell, and BEAR-GRN (*Nat Commun* 2026) benchmarks unsupervised methods on the same resource.
   - *Closest prior work.* Multi-task learning (Caruana 1997); scMTNI (multi-task across a lineage); GRNFormer's blind cell-type test; LINGER and CellPolaris (transfer from external atlases).
   - *Caveats.*
     - Joint multi-task training is a generic technique, so the novelty is the setting and the result, not the method.
     - The same caveats as item 2 apply, plus shared TF/gene identities across cell types.

4. **A role-aware, accessibility-gated decoder over separate modality towers with label-free structural graphs, as an explicit alternative to scMultiomeGRN's attention fusion over the training adjacency.**
   - *Evidence it is new.* We found no GRN GNN that applies a learned ATAC-derived gate only to the target's RNA state inside an asymmetric bilinear scorer.
   - *Closest prior work.*
     - scMultiomeGRN (cross-modal attention).
     - PECA (Duren *et al.*, *PNAS* 2017: TF expression × accessibility).
     - CellOracle (hard accessibility mask).
     - GLU gating.
     - scTFBridge (motif-masked decoder, *Nat Commun* 2025).
   - *Caveat.* There is no measurable accuracy gain: gated fusion, concat fusion, edge-MLP and role-aware all lie within 0.01 dual AUPR. Present it as an interpretability and biological-prior choice, not as a performance contribution. The "label-free graph" half is undermined by the positive exclusion (risk R3).

5. **Minor engineering contributions.**
   - The RP locus-shape descriptor, which fixed a 96%-redundant ATAC feature.
   - A dependency-free motif relation.
   - An interpretable relation-weight combiner.
   - *Caveat.* These are Adapted/Borrowed, and neither the motif relation nor the combiner is in the headline model. List them as implementation details, or as ablations with honest null results, not as contributions.

**Not claimable as novel:**
- Using FM embeddings (scRegNet did it).
- MAESTRO-style RP features (shared with scMultiomeGRN).
- Modality-specific towers (GraFRank / scMultiomeGRN).
- GraphSAGE / R-GCN message passing.
- Hard-negative mining in general.
- Experience replay or continual learning in GRN inference (LINGER).
- Evaluation against an intersection gold standard in itself (scMTNI, McCalla).

---

## 5. Where scMultiomeGRN is stronger, or does something we don't

1. **Cell-level signal.**
   - scMultiomeGRN keeps per-cell RP vectors and GRNBoost2 importance profiles. MeVD-GRN reduces each gene to 3 + 3 + 4 numbers.
   - Our own `fm_only` ≈ `with_fm` result shows that our hand-built RNA features carry almost no marginal information.
   - *Adopt:* metacell-level or per-cell summaries (e.g. SEACells or cluster-level pseudobulk vectors) as encoder inputs, or contextual FM embeddings pooled per cell type (scRegNet-style) instead of the static token table.
2. **Rich pairwise edge features.**
   - scMultiomeGRN encodes a 16 × 16 joint expression histogram (CNNC; Yuan & Bar-Joseph, *PNAS* 2019) with a small CNN inside every message.
   - Our pairwise information is one scalar (the signature dot product), plus target-only openness.
   - *Adopt:* a histogram or copula feature for decoded (TF, target) pairs. This is feasible because only mini-batched pairs are decoded.
3. **Attention with edge features, and node-wise modality weighting.**
   - scMultiomeGRN learns *per-node* modality weights, which are interpretable ("which TFs are ATAC-driven").
   - Our gate is per-dimension on the target only, and our relation weights are global scalars.
   - *Adopt:* report per-TF gate statistics, or add node-wise modality attention as an ablation arm.
4. **Using the observed network.**
   - Message passing over known edges is a strong signal that scMultiomeGRN exploits (transductively). We deliberately withhold it.
   - The scMultiomeGRN-benchmark pipeline already has an `--obs_graph` variant (`scripts/17_scmgrn_train_eval.py`).
   - *Adopt:* report an obs-graph variant on SC-MO-GRN-DB as well, built from train positives only, so reviewers see what the label-free design costs.
5. **Final network output and ensembling.**
   - scMultiomeGRN outputs a binarised GRN by a 10-run majority vote.
   - MeVD-GRN outputs only rankings and has never been used to produce a network.
   - *Adopt:* seed-ensemble plus vote to produce a final K562 network, and report its stability (Jaccard across seeds, BEELINE-style).
6. **Robustness analyses we lack.**
   - scMultiomeGRN has masked-positive recovery (Fig. 2D–E), cell-number downsampling (Fig. 2F, which plateaus above ~50 cells), and KL separation of score distributions (Fig. 3C).
   - *Adopt:* downsample K562 cells, and mask positives then measure recovery. Both are cheap.
7. **Biological case studies.**
   - scMultiomeGRN shows cell-type-specific TF subnetworks (SPI1 in monocytes; KLF12, GATA3 and EGR4 in lymphocytes; SPI1/RUNX1 in AD microglia).
   - MeVD-GRN has no biological validation beyond the metrics.
   - *Adopt:* K562 case studies for erythroid/megakaryocytic regulators such as GATA1 and TAL1, checked against held-out perturbation data.
8. **Breadth of cell types.** Nine lung cell types (from 544 to 43k cells) plus three hematopoietic ones, against our three cell lines, only one of which has all three tiers.
9. **Released code.** Their code is public (Zenodo) and includes a demo. Ours is not yet released.
10. **Training to convergence.**
    - The baseline trains ≤ 2000 epochs per model; we train 45.
    - This favours scMultiomeGRN on compute, and our adapter's perturbation run had not converged at 2000 epochs. The fairness question goes both ways (risk R9).

---

## 6. Weaknesses and risks in our novelty and performance claims

- **R1. The headline model does not contain most of the "contributions."**
  - `all_at_once` removes the curriculum, replay and hard negatives. The motif relation and the gated combiner are off.
  - The paper abstract and introduction still present the evidence-tier *curriculum* as the central contribution. Under the headline configuration the tier structure only drives the split and the held-out evaluation.
  - Reframe the contribution as a protocol, or show a tier-aware training variant that beats `all_at_once`. None does at h384 so far: sequential 0.751/0.834/0.963 and `all_at_once` 0.968/0.696/0.952 trade tiers against each other.
- **R2. Every reported MeVD-GRN number predates the 2026-09-30 eval-negative fix, and the comparison with scMultiomeGRN is asymmetric.**
  - Before the fix, 100% of val/test negatives were in MeVD-GRN's training pool (`docs/experiments/scmultiomegrn_generalization.md` §10).
  - Our scMultiomeGRN adapter trains only on `split["train"]["neg"]`.
  - So the "5 of 6 metrics" win compares a leaky MeVD-GRN against a non-leaky baseline. The headline must be rerun with `exclude_eval_negatives: true` before any claim.
- **R3. Label-dependent graph construction.**
  - The TF-candidate and motif graphs exclude every known positive of every split. A test positive is therefore systematically *absent* from its TF's top-500 candidate list, while a test negative with high co-expression can be present.
  - The effect is weak at 23k genes, but it is a leak by construction.
  - The fix is already used for the scMultiomeGRN benchmark (`prior_exclude_positives: false`): exclude train positives only, or build the graph label-free.
- **R4. "Zero-shot" dual evidence is not a distribution shift.**
  - Dual edges are 100% nested in perturbation. Under the global split, the dual val/test positives are exactly those perturbation val/test positives that also carry localization evidence. These are edges supported by two assay types, and plausibly the easiest held-out edges (dual AUPR 0.95 vs perturbation 0.70).
  - Dual *val* positives are also part of the merged val set used for early stopping.
  - Dual AUPR drove about six design decisions (plan.md §9).
  - Describe it as "held-out intersection-tier edges", not as generalisation to a new label type, and add an untouched final-check split.
- **R5. Degree / hub bias (InfoSEM, Stock *et al.*, Yılmaz *et al.*).**
  - Splits are random edge splits with every TF in train and test, and negatives are sampled uniformly over TFs and genes.
  - InfoSEM (ICML 2025) shows a gene-ID-only logistic regression matches GENELink/scGREAT under such splits, with supervised methods dropping 42–79% on unseen genes.
  - Stock *et al.* (bioRxiv 2025) show an out-degree sorter matches GNNs, and degree-matched negatives push them to near-random.
  - Yılmaz *et al.* (*PNAS* 2025) show degree-only predictors beat link predictors in network biology.
  - MeVD-GRN has no degree-only or gene-ID-only baseline and no TF-disjoint split. The 1:5 evaluation ratio (random AUPR ≈ 0.167) also inflates AUPR relative to the real candidate density.
  - (A TF-disjoint, degree-matched protocol is being built for BEAR-GRN in `src/benchmarks/bear_protocol.py`. It is in progress and not yet part of any reported result.)
- **R6. The FM embedding is a gene-identity lookup.**
  - The Geneformer input-token table is the same for a gene in every context. A model given it can memorise per-TF and per-target base rates, which is exactly the shortcut InfoSEM describes.
  - That would explain why it is the biggest in-domain win (`fm_only` ties `with_fm`) and hurts transfer to cell types with different hub structure.
  - Required control: a learned one-hot gene-ID embedding of the same width, and a TF-disjoint split with and without FM.
- **R7. ATAC may add little under realistic evaluation.**
  - geneRNIB, BEAR-GRN and PEREGGRN find motif/ATAC methods no better than, or worse than, expression-only ones (`docs/reference/multiome_grn_benchmark_consensus.md` §5).
  - Our own `rna_only` gap is small: dual AUPR 0.856 vs 0.881 in results.md (0.859 in the JSON). That was measured at the 315k, no-FM, sequential size, single seed, before the fix. It has **never been measured at the headline size**.
  - K562's RNA and ATAC are not even from the same cells.
  - The "dual-modality" framing needs a clean `rna_only` ablation at h384 + FM, with multiple seeds. The result may be null.
  - scMultiomeGRN's strong ATAC-only ablation (0.79 AUROC) is not evidence the other way: its labels are derived from the same ATAC.
- **R8. The curriculum idea has precedent under other names.**
  - Curriculum learning, multi-fidelity or noisy-to-clean supervision, weak-to-strong generalisation.
  - Continual learning is already used for GRNs in LINGER.
  - Intersection gold standards are used by scMTNI and McCalla.
  - Claim novelty only for the *specific use of experimental evidence types of one resource as training stages and held-out targets*, and only "to our knowledge".
- **R9. The baseline is an adaptation, not scMultiomeGRN as published.**
  - It uses our pseudobulk features instead of MAESTRO per-cell RP and GRNBoost2 vectors.
  - Its neighbours are capped at 500 per TF, with a 200k training-edge cap per epoch.
  - Its symmetric decoder is applied to a *directed* task, which is a structural handicap.
  - It is single-seed and was not converged at 2000 epochs.
  - The fair comparison on scMultiomeGRN's own benchmark (`configs/scmgrn/`) is pre-registered but **not yet run**.
  - Until then, claims should read "outperforms an adapted scMultiomeGRN on our benchmark".
- **R10. Null architectural ablations.** The role-aware decoder and the gated combiner show no accuracy gain. The motif relation was never tested in the headline setting or with 3 gated relations.
- **R11. Narrow evidence base.**
  - One cell line (K562) has all three tiers.
  - The transfer and joint-training numbers are n = 1, localization-only for the target cell types, and at h128.
- **R12. Compute mismatch.** 45 vs 2000 epochs, and 3.26M vs ≈4.08M parameters. We cannot claim architectural superiority from accuracy alone.

---

## 7. Suggested honest wording for the paper

- *Related work.* "scMultiomeGRN (Xu *et al.*, 2025) is the closest multi-omic graph model. It adopts GraFRank's modality-specific neighbour aggregation and cross-modal attention, and predicts undirected TF–TF links on motif-derived labels. scRegNet (Kommu *et al.*, 2025) previously combined frozen single-cell foundation-model gene embeddings with a GNN for supervised TF→target prediction. LINGER (Yuan & Duren, 2024) uses lifelong learning to transfer regulatory knowledge from external bulk data. Intersection gold standards (ChIP ∩ perturbation) have been used for evaluation (scMTNI; McCalla *et al.*)."
- *Contributions.* Lead with the protocol (R1, R4) and the cross-cell-type and FM findings (items 2–3), each with its caveats. Describe the decoder as an interpretable biological prior with no measured accuracy gain. Do not describe the curriculum, replay, motif relation or gated combiner as components of the headline model.

---

## 8. Code-vs-documentation discrepancies found

| # | Where | Document says | Code / data says |
|---|---|---|---|
| D1 | `paper/main.tex` §Model encoders | RNA input = 3 stats **concatenated with the 50-d co-expression signature** | `rna_in_dim: 3`. The signature is used only in graph construction and the decoder's `coexpr` term (`mevd_grn.py: decode`). |
| D2 | `paper/main.tex` §Model encoders | ATAC ∈ ℝ², single Linear + LN + GELU | ATAC is 3-d RP `[mean, var, detection]` through a 2-layer MLP (`encoders.py`). The paper's own parameter table says 3→64→128, so the paper is internally inconsistent. |
| D3 | `paper/main.tex` §Role-aware decoder | Scalar `open_j` with weight w_o | 4-d locus descriptor through `Linear(4,1)` (`decoder.py`). |
| D4 | `paper/main.tex` §GNN towers; `docs/citations.md` §8 | Graphs are "never built from labeled positive edges" / "so labels cannot leak" | Graphs are built by **excluding all** positives, which makes them label-dependent. The repo itself says so in `docs/experiments/scmultiomegrn_generalization.md` §10 and `configs/scmgrn/mevd_base.yaml`. |
| D5 | `paper/main.tex` abstract + introduction | The curriculum is the central contribution; hard negatives and replay are described as part of training | The headline is `all_at_once` (`build_all_at_once_stage` disables hard negatives; there is no replay). |
| D6 | `paper/main.tex` abstract | "no result in this paper yet has a repeated-seed variance estimate" | A 5-seed headline exists (results.md §7), and results.md/plan.md say the paper was updated with it. The abstract was not. |
| D7 | `paper/main.tex` §Architectural ablation hooks | "Three flags" | There are also `graph_mode`, `combine_mode`, `use_fm`, `use_motif` and `use_edge_mlp`. |
| D8 | `configs/default.yaml`, `configs/k562*.yaml`, `configs/sweep/k562_fm_h384_l2.yaml` | Default protocol `sequential`, h128, no FM; the `protocol` comment says sequential beats all_at_once | The headline comes from `06_ablation.py --ablation all_at_once` on the h384 sweep config. No shipped config reproduces it directly. plan.md §10 item 1 acknowledges this. |
| D9 | All reported numbers (results.md, paper) | Presented as current | They predate the `restrict_negative_pool` fix (2026-09-30), and `exclude_eval_negatives: true` is now the default. They need a rerun. |
| D10 | `docs/citations.md` §4 | Asymmetric bilinear scoring attributed to DistMult | DistMult is diagonal and **symmetric**. Our W is a full matrix, i.e. RESCAL (Nickel *et al.*, ICML 2011). |
| D11 | `docs/citations.md` §2, §6; `paper/references.bib` `duren2018integrative` | Duren *et al.*, "Modeling gene regulation from paired expression and chromatin accessibility data", PNAS 2018, doi 10.1073/pnas.1802973115 | **That DOI does not resolve in Crossref.** The title is PNAS 2017, 114(25) E4914, doi 10.1073/pnas.1704553114 (Duren, Chen, Jiang, Wang, Wong). The bib's author list and volume 115(30) belong to the 2018 coupled-NMF paper, doi 10.1073/pnas.1805681115. |
| D12 | `docs/citations.md` §17 | The gated combiner is motivated by the Relational Graph Transformer | The mechanism is a global per-layer softmax scalar per relation. HAN's semantic attention (WWW 2019) is the closer precedent. It is off in the headline, and at h384 + FM with 2 relations it learns near-uniform weights. |
| D13 | `docs/citations.md` §19; `plan.md` §2, §4 | The motif graph is "a genuinely different third relation" | The scan is restricted to TF-candidate edges (`candidate_edges=`), so the motif relation is a **subset** of the TF-candidate relation. `use_motif` is false in every config except `k562_fm_h384_l2_motif.yaml`, and that config has no result. |
| D14 | `docs/citations.md` §10, §12 | Replay "on by default"; half the negatives are hard | True for `default.yaml`, not for the headline. |
| D15 | `docs/reference/literature_knowledge_base.md` header | MeVD-GRN predicts edges "from paired scRNA-seq + scATAC-seq" | K562 RNA is DS019 and ATAC is DS025, from different datasets (`configs/k562.yaml`). Features are gene-level pseudobulk. |
| D16 | `paper/main.tex` §Related Work | scMultiomeGRN is "a supervised graph neural network" | The authors call it semi-supervised generative. It is supervised on its own motif-derived labels. |
| D17 | scMultiomeGRN paper vs its own code | Decoder σ(hhᵀ) (Eq. 11); FIMO v5.1 `-thresh 1e-6` | The code uses a concatenation-MLP decoder, symmetrised, and applies BCE-with-logits to an already-sigmoided output. The pipeline runs FIMO at 1e-4 and then filters to p ≤ 1e-6 (MEME 5.4.1 in our reproduction). |
| D18 | results.md §2 `rna_only` | 0.856 / 0.962 dual | `results/ablations/rna_only_K562.json`: 0.859 / 0.963. A minor mismatch, probably a rerun. |

---

## 9. References

- Xu J, Lu C, Jin S, *et al.* Deep learning-based cell-specific gene regulatory networks inferred from single-cell multiome data (scMultiomeGRN). *Nucleic Acids Res* 53:gkaf138 (2025). doi:10.1093/nar/gkaf138. Code: doi:10.5281/zenodo.14848389
- Valensi H, Karamveer K, Moeller E, *et al.* SC-MO-GRN-DB: a comprehensive repository for single-cell multiomic gene regulatory networks. *iScience* 29:115323 (2026). doi:10.1016/j.isci.2026.115323
- Karamveer K, Moeller E, Valensi H, *et al.* BEAR-GRN: systematic assessment of single-cell multi-omics-based GRN inference methods. *Nat Commun* (2026). doi:10.1038/s41467-026-77838-w
- Sankar A, Liu Y, Yu J, Shah N. Graph neural networks for friend ranking in large-scale social platforms (GraFRank). *WWW* 2021. doi:10.1145/3442381.3450120
- Li H, Sun Y, Hong H, *et al.* Inferring transcription factor regulatory networks from single-cell ATAC-seq data based on graph neural networks (DeepTFni). *Nat Mach Intell* (2022). doi:10.1038/s42256-022-00469-5
- Kommu S, Wang Y, Wang Y, Wang X. Prediction of gene regulatory connections with joint single-cell foundation models and graph-based learning (scRegNet). *Bioinformatics* 41 (ISMB 2025). doi:10.1093/bioinformatics/btaf217
- Theodoris CV, Xiao L, Chopra A, *et al.* Transfer learning enables predictions in network biology (Geneformer). *Nature* 618 (2023). doi:10.1038/s41586-023-06139-9
- Qi J, Li H, Cui Y, Zheng Y, Huang J. Towards universal gene regulatory network inference: unlocking generalizable regulatory knowledge in single-cell foundation models. arXiv:2605.08128 (2026)
- Yuan Q, Duren Z. Inferring gene regulatory networks from single-cell multiome data using atlas-scale external data (LINGER). *Nat Biotechnol* (2024). doi:10.1038/s41587-024-02182-7
- Bravo González-Blas C, De Winter S, Hulselmans G, *et al.* SCENIC+: single-cell multiomic inference of enhancers and gene regulatory networks. *Nat Methods* 20 (2023). doi:10.1038/s41592-023-01938-4
- Wang F-a, Yi C, Chen J, *et al.* scTFBridge. *Nat Commun* 16:9166 (2025). doi:10.1038/s41467-025-64227-y
- Xu T, *et al.* EpiAwareNet: prior-guided multi-omic transformers for single-cell GRN inference. KDD 2026, arXiv:2606.00685
- Chen G, Liu Z-P. Graph attention network for link prediction of gene regulations from single-cell RNA-sequencing data (GENELink). *Bioinformatics* 38(19) (2022). doi:10.1093/bioinformatics/btac559
- Mao G, Pang Z, Zuo K, *et al.* Predicting gene regulatory links from single-cell RNA-seq data using graph neural networks (GNNLink). *Brief Bioinform* 24(6) (2023). doi:10.1093/bib/bbad414
- Zhang S, Pyne S, Pietrzak S, *et al.* Inference of cell type-specific gene regulatory networks on cell lineages from single cell omic datasets (scMTNI). *Nat Commun* 14 (2023). doi:10.1038/s41467-023-38637-9
- McCalla SG, Fotuhi Siahpirani A, Li J, *et al.* Identifying strengths and weaknesses of methods for computational network inference from single-cell RNA-seq data. *G3* 13 (2023). doi:10.1093/g3journal/jkad004
- Hegde A, Cheng J. GRNFormer: accurate gene regulatory network inference using graph transformer. *Bioinformatics* 42(4) (2026). doi:10.1093/bioinformatics/btag144
- Sun, Gao, *et al.* HGATLink. *BMC Bioinformatics* (2025). doi:10.1186/s12859-025-06071-x
- Cui T, Xu S-J, Moskalev A, *et al.* InfoSEM: a deep generative model with informative priors for gene regulatory network inference. ICML 2025, arXiv:2503.04483
- Stock *et al.* Hidden sampling biases inflate performance in gene regulatory network inference. bioRxiv (2025). doi:10.64898/2025.12.19.695616
- Yılmaz S, Yorgancıoğlu K, Koyutürk M, *et al.* Bias-aware training and evaluation of link prediction algorithms in network biology. *PNAS* (2025). doi:10.1073/pnas.2416646122
- Kamimoto K, Stringa B, Hoffmann CM, *et al.* Dissecting cell identity via network inference and in silico gene perturbation (CellOracle). *Nature* 614 (2023). doi:10.1038/s41586-022-05688-9
- Duren Z, Chen X, Jiang R, Wang Y, Wong WH. Modeling gene regulation from paired expression and chromatin accessibility data (PECA). *PNAS* 114(25):E4914 (2017). doi:10.1073/pnas.1704553114
- Duren Z, Chen X, Zamanighomi M, *et al.* Integrative analysis of single-cell genomics data by coupled nonnegative matrix factorizations. *PNAS* 115(30) (2018). doi:10.1073/pnas.1805681115
- Wang C, Sun D, Huang X, *et al.* Integrative analyses of single-cell transcriptome and regulome using MAESTRO. *Genome Biol* 21 (2020). doi:10.1186/s13059-020-02116-x
- Yuan Y, Bar-Joseph Z. Deep learning for inferring gene relationships from single-cell expression data (CNNC). *PNAS* 116 (2019). doi:10.1073/pnas.1911536116
- Moerman T, Aibar Santos S, Bravo González-Blas C, *et al.* GRNBoost2 and Arboreto. *Bioinformatics* 35 (2019). doi:10.1093/bioinformatics/bty916
- Huynh-Thu VA, Irrthum A, Wehenkel L, Geurts P. Inferring regulatory networks from expression data using tree-based methods (GENIE3). *PLoS ONE* 5 (2010). doi:10.1371/journal.pone.0012776
- Li *et al.* GMFGRN. *Brief Bioinform* 25 (2024). doi:10.1093/bib/bbad529
- Zhu H, Slonim D. From noise to knowledge: diffusion probabilistic model-based neural inference of GRNs (RegDiffusion). *J Comput Biol* (2024). doi:10.1089/cmb.2024.0607
- Rauluseviciute I, *et al.* JASPAR 2024. *Nucleic Acids Res* 52 (2024). doi:10.1093/nar/gkad1059
- Grant CE, Bailey TL, Noble WS. FIMO: scanning for occurrences of a given motif. *Bioinformatics* 27 (2011). doi:10.1093/bioinformatics/btr064
- Hamilton W, Ying R, Leskovec J. Inductive representation learning on large graphs (GraphSAGE). NeurIPS 2017, arXiv:1706.02216
- Schlichtkrull M, *et al.* Modeling relational data with graph convolutional networks (R-GCN). ESWC 2018, arXiv:1703.06103
- Wang X, Ji H, Shi C, *et al.* Heterogeneous graph attention network (HAN). *WWW* 2019. doi:10.1145/3308558.3313562
- Dwivedi VP, *et al.* Relational Graph Transformer. arXiv:2505.10960 (2025)
- Gilmer J, *et al.* Neural message passing for quantum chemistry (MPNN). ICML 2017, arXiv:1704.01212
- Kipf TN, Welling M. Variational graph auto-encoders. arXiv:1611.07308 (2016)
- Nickel M, Tresp V, Kriegel H-P. A three-way model for collective learning on multi-relational data (RESCAL). ICML 2011 (no DOI)
- Yang B, *et al.* Embedding entities and relations for learning and inference in knowledge bases (DistMult). ICLR 2015, arXiv:1412.6575
- Dauphin YN, *et al.* Language modeling with gated convolutional networks (GLU). ICML 2017, arXiv:1612.08083
- Bengio Y, Louradour J, Collobert R, Weston J. Curriculum learning. ICML 2009. doi:10.1145/1553374.1553380
- McCloskey M, Cohen NJ. Catastrophic interference in connectionist networks. *Psychol Learn Motiv* 24 (1989). doi:10.1016/S0079-7421(08)60536-8
- Rolnick D, *et al.* Experience replay for continual learning. NeurIPS 2019, arXiv:1811.11682
- Kirkpatrick J, *et al.* Overcoming catastrophic forgetting in neural networks (EWC). *PNAS* 114 (2017). doi:10.1073/pnas.1611835114
- Caruana R. Multitask learning. *Mach Learn* 28 (1997). doi:10.1023/A:1007379606734
- Schroff F, Kalenichenko D, Philbin J. FaceNet. CVPR 2015. doi:10.1109/CVPR.2015.7298682
- Saito T, Rehmsmeier M. The precision-recall plot is more informative than the ROC plot… *PLoS ONE* 10 (2015). doi:10.1371/journal.pone.0118432
- Pratapa A, Jalihal AP, Law JN, *et al.* Benchmarking algorithms for GRN inference from single-cell transcriptomic data (BEELINE). *Nat Methods* 17 (2020). doi:10.1038/s41592-019-0690-6

DOIs above were checked against the Crossref API and arXiv IDs against the arXiv API on 2026-10-01. Exceptions: the Zenodo DOI (DataCite, not Crossref), RESCAL (no DOI), and the EpiAwareNet, HGATLink and GMFGRN author lists, which were taken from the local knowledge base.
