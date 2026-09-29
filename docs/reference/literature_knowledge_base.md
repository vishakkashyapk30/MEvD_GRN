# Literature knowledge base for MeVD-GRN

Purpose: this maps which datasets each competing paper trained and evaluated on, and which numbers MeVD-GRN has to beat. MeVD-GRN is a dual-modality GNN that predicts TF→target edges from paired scRNA-seq + scATAC-seq, trained on SC-MO-GRN-DB.

Conventions:
- Every number is copied from the paper, with the table or figure reference.
- "figure only, approx. X" means the value was read off a plot and is not printed in the paper.
- "not in paper" means the paper does not state it.
- Links marked "(found online, not in paper)" were found by web search and are not in the source PDF.
- "(reviewer note)" marks our own critical observations, not the authors' claims.
- Source PDFs are in `/home/vishak/research/` (top level, `papers/`, `ref-theses/`). The project's own paper (`mevd_grn/paper/main.pdf`) is excluded.

Status: complete (2026-09-30). All 17 source PDFs were read. Main-results figures without printed tables were checked on high-resolution page renders: ATFGRN Fig. 2, SMOGT Fig. 2, scTFBridge Fig. 5.

---

## 1. Overview table

| # | Paper (file) | Year / venue | Task | Input modalities | Model family | Code |
|---|---|---|---|---|---|---|
| 2.1 | scUniGP (`8402_Unifying_Graph_Based_and_.pdf`) | 2025, under review at ICLR 2026 | supervised TF→target link prediction | scRNA + prior GRN | GAT/GCN + pairwise Transformer, multi-scale fusion | not public |
| 2.2 | ATFGRN (`ATFGRN.pdf`) | 2026, Brief Bioinform 27(1) bbaf733 | supervised TF→target link prediction | scRNA + prior GRN | SEAL-style subgraph GCN + TransformerConv + KNN/Node2Vec/GAT, attention fusion | https://github.com/Evert-zrt/ATFGRN |
| 2.3 | CaHoT-GRN (`CaHoT-GRN.pdf`) | 2026, Brief Bioinform 27(2) bbag202 | supervised TF→target link prediction | scRNA + DNABERT + ESM-2 embeddings + HuRI PPI + prior GRN | multi-view GAT with similarity co-attention | https://github.com/ydkvictory/CaHoT-GRN |
| 2.4 | scRegNet (`scregnet.pdf`) | 2025, Bioinformatics 41 (ISMB) | supervised TF→target link prediction | scRNA (+ frozen scFM) + prior GRN | scBERT / Geneformer / scFoundation embeddings + GCN | https://github.com/sindhura-cs/scRegNet |
| 2.5 | Kommu MS thesis (`ref-theses/thesis_found_model_scRNA.pdf`) | 2025, Virginia Tech | = scRegNet + GRN-guided Geneformer pretraining | scRNA | scFM + GNN | (scRegNet repo) |
| 2.6 | SMOGT (`SCMOGT.pdf`) | 2025, Brief Bioinform 26(6) bbaf664 | TF→CRE and CRE–CRE edges, hierarchical TF–CRE–TG network | **paired scRNA + scATAC** + priors | heterogeneous graph transformer, semi-supervised | https://github.com/YuHongHuang-lab/SMOGT |
| 2.7 | scMultiomeGRN (`scMultiomeGRN.pdf`) | 2025, NAR 53 gkaf138 | semi-supervised TF–TF link prediction (motif-derived labels) | scRNA + scATAC (mostly unpaired) | modality-specific GAT aggregators + cross-modal attention, graph autoencoder | https://doi.org/10.5281/zenodo.14848389 |
| 2.8 | scTFBridge (`scTFBridge.pdf`) | 2025, Nat Commun 16:9166 | unsupervised TF→TG and RE→TG regulatory scores | **paired scRNA + scATAC** | disentangled multimodal VAE, motif-masked decoder, SHAP | https://github.com/FengAoWang/scTFBridge |
| 2.9 | EpiAwareNet (`epiawarenet.pdf`) | 2026, KDD (arXiv 2606.00685) | weakly supervised regulator→target ranking | **paired scRNA + scATAC** + curated prior | multi-omic Transformer (gene–peak cross-attention) + MLP head | https://github.com/tianyang-x/EpiAwareNet_pub |
| 2.10 | SC-MO-GRN-DB (`SC_MO_GRN_Db.pdf`) | 2026, iScience 29 115323 | resource: reference GRNs + single-cell multiomic datasets | 6 modalities | – | https://github.com/UzunLab/SC-MO-GRN-DB ; https://scmogrndb.psu.edu |
| 2.11 | EpiXFormer (`EpiXFormer.pdf`) | 2026, Brief Bioinform 27(1) bbaf721 | TF binding-site (region) classification | DNA sequence + bulk epigenome | cross-attention Transformer | not stated |
| 2.12a | DECODE (`decode.pdf`) | 2026, Nat Methods | cell-type deconvolution (**out of scope**) | bulk + sc references (RNA, protein, metabolite) | deep deconvolution | Zenodo |
| 2.12b | Maizels & Briscoe (`Gene regulatory networks .pdf`) | 2026, Nat Rev Genet | perspective on GRN causality (**out of scope**, no experiments) | – | – | – |
| 2.12c | HyperCLSA (`papers/HyperCLSA.pdf`) | IIIT-H preprint / proceedings | breast-cancer subtyping (**out of scope**) | bulk mRNA + methylation + miRNA | hypergraph contrastive learning | https://github.com/Gaurav2543/HyperCLSA |
| 2.12d | Shetty et al. pRCC (`papers/s00438-023-02022-4.pdf`) | 2023, Mol Genet Genomics 298 | survival subtyping (**out of scope**) | bulk CNV / RNA | network diffusion + clustering | – |
| 2.12e | Essential genes (`papers/2024.10.09.616990v1.full.pdf`) | 2024, bioRxiv | essential-gene classification (**out of scope**) | PPI networks | node embeddings + classifiers | not in paper |
| 2.12f | Aswin Jose MS thesis (`ref-theses/Thesis_Aswin.pdf`) | 2025, IIIT-H | cancer subtyping / prognosis (**out of scope**) | bulk TCGA mutations / CNV | pyNBS; graph autoencoder | – |

**How close each paper is to MeVD-GRN.** Closest, since they use paired RNA+ATAC and gene-level targets: EpiAwareNet, then scTFBridge, then scMultiomeGRN. Next are SMOGT (peak-level targets) and the BEELINE RNA-only family: scUniGP > scRegNet > CaHoT-GRN ≈ ATFGRN.

---

## 2. Per-paper sections

### 2.0 Shared protocol: the BEELINE / GENELink scRNA-seq benchmark

Five of the papers use the same benchmark: scUniGP, ATFGRN, CaHoT-GRN, scRegNet and the Kommu thesis. It is described once here.

- **Datasets.** Seven BEELINE scRNA-seq datasets (Pratapa et al., Nat Methods 2020), all **scRNA-seq only, with no ATAC**:

| Dataset | GEO | Species | Cells | Genes (TFs+500 / TFs+1000)† |
|---|---|---|---|---|
| hESC | GSE75748 | human | 758 | 910 / 1410 |
| hHEP | GSE81252 | human | 425 | 948 / 1448 |
| mDC | GSE48968 | mouse | 383 | 821 / 1321 |
| mESC | GSE98664 | mouse | 421 | 1120 / 1620 |
| mHSC-E | GSE81682 | mouse | 1071 | 704 / 1204 |
| mHSC-GM | GSE81682 | mouse | 889 | 632 / 1132 |
| mHSC-L | GSE81682 | mouse | 847 | 560 / 692 |

  †Cells and gene counts are from scUniGP Tables 4/5.

- **Gene sets.** "TFs+500" and "TFs+1000" are the significantly varying TFs plus the 500 or 1000 most-variable genes (BEELINE preprocessing: variance p < 0.01, Bonferroni).
- **Ground-truth networks (4 types).**
  - STRING (v11).
  - Non-specific ChIP-seq (DoRothEA / RegNetwork / TRRUST).
  - Cell-type-specific ChIP-seq (ChIP-Atlas / ENCODE / ESCAPE).
  - LOF/GOF perturbation (ESCAPE), mESC only.

  That gives 7+7+7+1 = **22 dataset×GT combinations per gene scale, 44 in total**.
- **Network density** (scUniGP Table 4, TFs+500):
  - STRING / non-specific: 0.013–0.048 (sparse).
  - Cell-type-specific: 0.082–0.578 (dense).
  - LOF/GOF: 0.158.

  AUPRC therefore depends strongly on the GT type, because random AUPRC ≈ density.
- **Split (GENELink, `Code/Train_Test_Split.py` in https://github.com/zpliulab/GENELink; found online, not in the papers).**
  - Per TF, 2/3 of positives go to train (1/5 of that to validation) and 1/3 to test.
  - Train/val negatives are 1:1 from the same TF's non-targets.
  - Test negatives are sampled so the positive fraction ≈ network density.
  - A "hard-negative" mode (0.67/0.10/0.23 of all pairs per TF) is used for the cell-type-specific networks.
  - This is a **transductive edge split: every TF and gene appears in both train and test**. There is no held-out-TF protocol.
- **Data sources.** Processed files are in the GENELink repo under `Dataset/Benchmark Dataset/{STRING, Non-Specific, Specific, Lofgof} Dataset/<cell>/TFs+{500,1000}/`, and raw BEELINE data is on Zenodo record 3701939 (both found online, not in papers). Download effort: very low (tens of MB).
- **What this means for MeVD-GRN.** No dataset has paired (or any) scATAC. MeVD-GRN can compete only in an RNA-only mode (ATAC branch ablated), or with an unpaired chromatin prior such as ENCODE/ChIP-Atlas bulk ATAC/DNase for H1, mESC, etc. Such a prior must be declared as external information.
- **Leakage warning (reviewer note).** An ATFGRN code audit (§2.2) found that trivial degree-based scorers reach AUROC 0.87–0.99 on some of these splits, and that STRING test positives often include reversed training edges. Any comparison should also report a leakage-controlled variant.

---

### 2.1 scUniGP: "Unifying Graph-Based and Pairwise-Based Representations for GRN Inference from scRNA-seq Data" (`8402_Unifying_Graph_Based_and_.pdf`)

**Citation.** Anonymous; under review at ICLR 2026 (submission 8402, PDF dated 2025-09-25). OpenReview PDF: https://openreview.net/pdf/bfe5ab916595dd5453bc226437aaabb52ffab65b.pdf (found online, not in paper). **Code:** not public. The paper says it is in the supplemental material; no GitHub repo was found.

**Task / modalities / family.** Supervised TF→target link prediction; scRNA-seq only; hybrid GNN plus pairwise Transformer.

**Problem.** Graph-based methods (GENELink, GNNLink) capture topology but overfit and miss local pairwise signal. Pairwise methods (GNE, CNNC, scGREAT) ignore global topology. scUniGP fuses the two.

**Method (enough to reimplement).**
- **Node features.** Each gene's raw feature is its expression vector across cells, e = X[:,k] ∈ R^c, standardized.
- **Global branch.** A GAT (the headline model) or GCN is pre-trained on the undirected training-positive graph with BCE.
  - Main text: 2 layers, 128→64-d.
  - Appendix C: 3 hidden layers 128/64/32 → 16-d output; GAT with 3 heads; dropout 0.01; Adam lr 3e-3; batch 256; early stopping with patience 5; max 30 epochs.
  - "Expert score": s_ij = h_i^(L)·h_j^(L).
- **Pairwise branch.** A 2-token sequence [ℓ(x_TF)+p0 ; ℓ(x_target)+p1] goes through a 4-layer Transformer encoder (8 heads, embedding 1024). p0/p1 are learnable role embeddings. Up to 5 neighbouring pairs are sampled for context; the mechanism is not detailed.
- **Multi-scale fusion.** At layer l, z^(l) = m_f(z^(l-1), h_i^(l), h_j^(l), 1{l=L}·s_ij), combined by concatenation. The classifier is 3 residual layers (BatchNorm + PReLU) with a sigmoid output, and the GNN score is late-fused.
- **Training.**
  - Stage 1: GNN pre-training. Stage 2: joint fine-tuning.
  - Adam, lr 5e-6, weight decay 1e-5, StepLR (γ = 0.999, step 10), dropout 0.2, L2 0.01, batch 512, ≤200 epochs, early stopping on validation AUROC (patience 8).
- **Loss.** BCE only.

**Datasets.** The 7 BEELINE datasets × {TFs+500, TFs+1000} × 4 GTs (§2.0). Accessions are given in the paper (Sec. 7). Density statistics: Tables 4/5.

**Protocol.**
- GENELink per-TF split: 2/3 train+val (9:1), 1/3 test, same-TF 1:1 hard negatives. The text also says "67%/6.7%/23.3%", which does not sum to 100 (reviewer note).
- Mean of 5 seeds (Appendix D says 3). **No std reported.**
- Model selection on validation AUROC; "both the best and average" test performance is reported.

**Main results. Table 1: AUROC, TFs+500, mean of 5 seeds.** Best baseline in brackets.

| GT | hESC | hHEP | mDC | mESC | mHSC-E | mHSC-GM | mHSC-L |
|---|---|---|---|---|---|---|---|
| STRING | 0.948 (GNNLink 0.921) | 0.941 (GNNLink 0.929) | 0.956 (GENELink 0.941) | 0.951 (scGREAT 0.934) | 0.942 (scGREAT 0.924) | 0.937 (scGREAT 0.920) | 0.882 (GNNLink 0.851) |
| Non-specific ChIP | 0.896 (scGREAT 0.882) | 0.906 (scGREAT 0.886) | 0.926 (scGREAT 0.907) | 0.928 (GENELink 0.887) | 0.893 (scGREAT 0.874) | 0.882 (scGREAT 0.880) | 0.842 (scGREAT 0.802) |
| Cell-type-specific ChIP | 0.895 (scGREAT 0.890) | 0.910 (scGREAT 0.908) | 0.813 (scGREAT 0.808) | 0.941 (scGREAT 0.930) | 0.930 (scGREAT 0.927) | 0.937 (scGREAT 0.928) | 0.885 (scGREAT 0.876) |
| LOF/GOF | – | – | – | 0.891 (scGREAT 0.888) | – | – | – |

- **Averages over the 22 combinations (Table 1):** scUniGP **0.911**, scGREAT 0.892, GNNLink 0.869, GENELink 0.865, GNE 0.715, GENIE3 0.596, DeepSEM 0.591, GRNBoost2 0.586, PCC 0.584, MI 0.583.
- **TFs+1000 AUROC average:** **0.921** (text, Sec. 4.3; Fig. 14 violin). scGREAT 0.908 is figure only (Fig. 14). Per-dataset TFs+1000 AUROC is **not reported numerically**: Figs 5/6, labelled as AUROC heatmaps, are duplicates of the AUPRC heatmaps (Figs 7/8) (reviewer note).
- **AUPRC averages:**
  - TFs+500: **0.522** (text; Fig. 15). scGREAT 0.463, GENELink 0.445.
  - TFs+1000: **0.543**, figure only (Fig. 16). scGREAT 0.538. Paired t-test vs scGREAT: p = 0.484, not significant (Fig. 12).
- **Per-dataset AUPRC, TFs+500, figure only (numbers printed in Fig. 7 heatmap cells):**
  - STRING: hESC 0.51, hHEP 0.52, mDC 0.66, mESC 0.39, mHSC-E 0.52, mHSC-GM 0.51, mHSC-L 0.35.
  - Non-specific: hESC 0.24, hHEP 0.23, mDC 0.38, mESC 0.32, mHSC-E 0.29, mHSC-GM 0.36, mHSC-L 0.34.
  - Specific: hESC 0.62, hHEP 0.85, mDC 0.18, mESC 0.88, mHSC-E 0.94, mHSC-GM 0.94, mHSC-L 0.89.
  - LOF/GOF mESC: 0.55.
- **Per-dataset AUPRC, TFs+1000, figure only (Fig. 8 cells):**
  - STRING: hESC 0.53, hHEP 0.54, mDC 0.67, mESC 0.60, mHSC-E 0.54, mHSC-GM 0.58, mHSC-L 0.30.
  - Non-specific: hESC 0.28, hHEP 0.24, mDC 0.42, mESC 0.26, mHSC-E 0.37, mHSC-GM 0.48, mHSC-L 0.21.
  - Specific: hESC 0.65, hHEP 0.85, mDC 0.18, mESC 0.89, mHSC-E 0.96, mHSC-GM 0.95, mHSC-L 0.88.
  - LOF/GOF mESC: 0.56.
  - scGREAT beats scUniGP on 6/7 STRING AUPRC cells at TFs+1000 (reviewer note).
- **Ablation (Table 2, AUROC TFs+500, 7-dataset means).**
  - No GNN (pairwise Transformer only): STRING 0.918 / Non-specific 0.867.
  - With GCN: 0.929 / 0.879.
  - With GAT: 0.937 / 0.896.
- **Fusion ablation (Table 3).** Early fusion 0.916/0.863; late fusion 0.931/0.886; the proposed fusion 0.937/0.896.
- **Runtime (Table 6).** TF+500: 18m21s vs scGREAT 15m10s. TF+1000: 41m08s.

**Baselines.** MI, PCC, GRNBoost2 and GENIE3 (reproduced following the scGREAT protocol); DeepSEM, GNE, GENELink, GNNLink and scGREAT (reimplemented). CNNC appears in the runtime table only.

**Limitations.**
- Author-stated: no reliable negatives; RNA only (chromatin accessibility is listed as future work).
- (Reviewer notes)
  - The figures are duplicated (see above).
  - The "40/42" and "33/42" win counts are wrong: there are 44 combinations, and by the printed values it is best or tied on 30/44 AUPRC.
  - No variances are reported.
  - It reports the "best" test run.
  - Appendix C says the GNN is pre-trained "on the ground-truth TF–target interaction graph", which is a possible leak.
  - The †-baseline values (0.650, 0.610, …) look copied from earlier papers rather than re-run.

**Relevance to MeVD-GRN.** This is the strongest RNA-only number set on BEELINE: AUROC 0.911 / 0.921 and AUPRC 0.522 / 0.543. It is easy to reproduce with the GENELink split files, but only as an RNA-only comparison. There is no code, so the numbers must be taken from the paper. A useful idea is that a pairwise TF–target Transformer alone already reaches AUROC 0.918 on STRING (Table 2, TFs+500 mean), so a pairwise head could be a strong addition to MeVD-GRN.

---

### 2.2 ATFGRN: "Revealing hidden regulatory dependencies: multi-perspective graph learning for single-cell GRN inference" (`ATFGRN.pdf`)

**Citation.** He W., Zhang R., Zhu Y., Zhou H., Zuo Y., Bai Y., Yang L., Guo F. *Briefings in Bioinformatics* 27(1), bbaf733, 2026. doi:10.1093/bib/bbaf733.
- **Code and data (in paper):** https://github.com/Evert-zrt/ATFGRN. The repo ships the expression matrices and fixed Train/Validation/Test CSVs (checked online in the earlier reading pass).

**Task / modalities / family.** Supervised TF→target link prediction; scRNA-seq plus a prior GRN; **no ATAC**. It is a three-branch graph model with attention fusion.

**Method.**
1. **Local-subgraph branch.** For each candidate edge, a 2nd-order enclosing subgraph is extracted from the prior adjacency A, labelled with DRNL (SEAL-style), encoded by a GCN, SortPooled to the top-k nodes and flattened.
2. **Global branch.** A TransformerConv (self-attention GNN) runs over the prior GRN with normalized expression features. Best settings: 8 heads, 2 layers (Fig. 5).
3. **Similarity branch.** A KNN graph is built on Euclidean distances between expression profiles. Node2Vec gives 32-d initial embeddings, followed by a 2-layer GAT and an LSTM across layers.
4. **Fusion.** Shared-query attention weights the three edge representations, and an MLP gives the score.
- **Loss:** L_total = L_fusion + L1 + L2 + L3, each BCE (Eq. 15–16).
- **Learning rate:** 0.001 (Fig. 5c).

**Datasets.** The 7 BEELINE sets × {TFs+500, TFs+1000} (§2.0). The paper gives no GEO accessions, only a reference to BEELINE. **Table 1** gives TF / gene / edge / density statistics identical to scUniGP's Tables 4/5.

**GT.**
- STRING.
- Non-specific ChIP-seq (DoRothEA / RegNetwork / TRRUST).
- Cell-type-specific ChIP-seq (ChIP-Atlas / ENCODE).
- LOF/GOF (ESCAPE), mESC only.

**Protocol.**
- The paper defers grouping and negative sampling to the Supplementary file, which was not available.
- The earlier code audit of the repo split files found:
  - Specific networks: per-TF split ≈67/10/23 with natural class imbalance.
  - STRING / Non-specific / LOF/GOF: ≈54/10/35 of positives, 1:1 negatives in train/val, and test negatives at (1 − density)/density drawn from all TFs.
- **(Reviewer note: serious problems found in that audit)**
  1. `train1.py` keeps the checkpoint with the best **test** AP, i.e. selection on the test set.
  2. On STRING, 25–49% of test positives are reversed training edges, and the model treats edges as undirected.
  3. Hundreds of training positives reappear as test negatives (377 in STRING hESC-500).
  4. 95% of LOF/GOF test negatives come from TFs with zero positives.
  5. A degree-only scorer with no expression gets AUROC 0.92 vs ATFGRN 0.81 on Non-specific hESC-500, and 0.99 vs 0.92 on LOF/GOF. **Treat ATFGRN numbers as inflated.**
  - These audit figures come from an earlier reading pass whose scratch notes were lost, and have not been re-verified here.

**Results. Fig. 2 heatmaps: figure only, but numbers are printed in each cell (2 decimals, no std).** Read from a 300-dpi render of PDF p. 7. Format: ATFGRN / best of GATCL, scMGATGRN, GENELink, DeepSEM (TFs+500 panel) or "GNNLink" (TFs+1000 panel), PCC, SCODE, GENIE3.

| GT | Dataset | AUROC 500 | AUROC 1000 | AUPRC 500 | AUPRC 1000 |
|---|---|---|---|---|---|
| STRING | hESC | 0.95 / scMGATGRN 0.93 | 0.95 / scMGATGRN 0.92 | 0.53 / GATCL 0.29 | 0.55 / GATCL 0.30 |
| STRING | hHEP | 0.94 / 0.91 | 0.95 / 0.92 | 0.54 / GATCL 0.29 | 0.56 / GATCL 0.29 |
| STRING | mDC | 0.95 / 0.92 | 0.95 / 0.93 | 0.62 / 0.32 | 0.64 / scMGATGRN 0.35 |
| STRING | mESC | 0.97 / 0.92 | 0.96 / 0.92 | 0.63 / GATCL 0.22 | 0.62 / GATCL 0.25 |
| STRING | mHSC-E | 0.94 / 0.91 | 0.92 / 0.91 | 0.54 / GATCL 0.35 | 0.53 / GATCL 0.38 |
| STRING | mHSC-GM | 0.92 / 0.91 | 0.90 / scMGATGRN 0.90 | 0.50 / GATCL 0.43 | 0.49 / GATCL 0.40 |
| STRING | mHSC-L | 0.87 / scMGATGRN 0.85 | 0.91 / 0.85 | 0.36 / scMGATGRN 0.35 | 0.37 / scMGATGRN 0.34 |
| Non-spec | hESC | 0.81 / 0.74 | 0.79 / 0.72 | 0.10 / 0.06 | 0.07 / 0.06 |
| Non-spec | hHEP | 0.82 / 0.74 | 0.80 / 0.73 | 0.11 / 0.06 | 0.09 / 0.06 |
| Non-spec | mDC | 0.87 / 0.82 | 0.86 / 0.82 | 0.29 / 0.17 | 0.26 / 0.16 |
| Non-spec | mESC | 0.88 / 0.83 | 0.88 / 0.82 | 0.23 / 0.09 | 0.23 / 0.09 |
| Non-spec | mHSC-E | 0.82 / 0.80 | 0.83 / 0.78 | 0.22 / GATCL 0.21 | 0.26 / GATCL 0.22 |
| Non-spec | mHSC-GM | 0.81 / 0.79 | 0.84 / 0.80 | 0.27 / GATCL 0.27 | 0.35 / GATCL 0.31 |
| Non-spec | mHSC-L | 0.75 / 0.71 | 0.76 / 0.71 | 0.13 / GATCL,PCC 0.14 | 0.12 / GATCL 0.16 |
| LOF/GOF | mESC | 0.92 / GATCL 0.83 | 0.90 / scMGATGRN 0.81 | 0.68 / GATCL 0.63 | 0.67 / scMGATGRN 0.59 |
| Specific | hESC | 0.87 / 0.87 | 0.88 / 0.88 | 0.61 / scMGATGRN 0.60 | 0.60 / scMGATGRN 0.61 |
| Specific | hHEP | 0.90 / 0.89 | 0.90 / 0.89 | 0.84 / 0.81 | 0.84 / 0.81 |
| Specific | mDC | 0.75 / 0.76 | 0.78 / scMGATGRN 0.79 | 0.13 / 0.13 | 0.15 / 0.14 |
| Specific | mESC | 0.93 / 0.91 | 0.93 / scMGATGRN 0.94 | 0.86 / 0.81 | 0.86 / 0.84 |
| Specific | mHSC-E | 0.92 / 0.91 | 0.93 / 0.92 | 0.94 / 0.93 | 0.95 / 0.91 |
| Specific | mHSC-GM | 0.92 / 0.92 | 0.94 / 0.92 | 0.93 / 0.92 | 0.95 / 0.93 |
| Specific | mHSC-L | 0.87 / 0.85 | 0.88 / scMGATGRN 0.88 | 0.87 / 0.85 | 0.87 / 0.85 |

- **Text claims.**
  - +5.09% average AUROC over the second-best method (scMGATGRN) on the 15 STRING/Non-specific/LOF-GOF TFs+500 datasets.
  - TFs+1000 average AUROC "exceeding 0.88".
  - STRING average AUPRC 0.53 vs GATCL 0.32.
  - Highest AUROC on "all 30 datasets"; best AUPRC on 27/30.
  - Specific networks: AUROC ≥ scMGATGRN on 11/14; average AUPRC +2.67%.
  - The paper publishes no overall averages.
- **(Reviewer note) Column label.** The fifth column is labelled "DeepSEM" in the TFs+500 panels but "GNNLink" in the TFs+1000 panels, and the values differ. Which method the TFs+1000 column is remains unclear.
- **Ablation (Fig. 4, STRING, text only).** Average AUROC drop: without the subgraph branch −4.27%, without TransformerConv −1.61%, without KNN −1.74%.
- **Imputation and training-size experiments (Figs 6–7):** figure only.

**Baselines.** PCC, GENIE3, SCODE, DeepSEM, GENELink, GATCL, scMGATGRN, and "GNNLink" in the TFs+1000 panels.

**Limitations.**
- Author-stated: the negative-sampling scheme may include false negatives.
- (Reviewer notes)
  - Leakage and test-set selection, as above.
  - No std is reported.
  - The results are figure-only heatmaps.

**Relevance to MeVD-GRN.** ATFGRN's STRING TFs+500 AUPRC averages 0.53 (text; the 7 cells average 0.531). That is above scUniGP (0.494 heatmap mean) and CaHoT-GRN (0.363 recomputed mean). Yet GATCL scores about the same average in both ATFGRN (0.32) and CaHoT-GRN (0.320). So the split files look comparable, and the ATFGRN margin may partly come from the leakage issues above. If MeVD-GRN is compared on ATFGRN's released split files, **also report a leakage-controlled re-run of ATFGRN**: validation-based selection, reverse-edge dedup, negatives only from GT TFs. The subgraph (SEAL/DRNL) branch is the component its ablation shows matters most.

---

### 2.3 CaHoT-GRN: context-aware high-order topology learning (`CaHoT-GRN.pdf`)

**Citation.** Yao D., Zhang B., Zhan X., Wang W., Liang N. *Briefings in Bioinformatics* 27(2), bbag202, 2026 (published 30 Apr 2026). doi:10.1093/bib/bbag202.
- **Code (in paper):** https://github.com/ydkvictory/CaHoT-GRN. The repo holds only a README and a .rar archive (found online, not in paper). The README example is `--epochs 20 --batch_size 256 --net Specific --num 500 --data hHEP`, which is GENELink's CLI.
- **Data link:** the "Zenodo" hyperlink in the paper is broken (its target is the literal string "Zenodo").

**Task / modalities / family.** Supervised TF→target link prediction.
- Inputs: scRNA-seq expression, frozen DNABERT embeddings of gene DNA sequences, frozen ESM-2 embeddings of protein sequences, the HuRI human PPI (as a filter) and the prior GRN topology. **No ATAC.**
- Model: multi-view GAT with similarity co-attention; GENELink family.

**Problem.** Expression-plus-topology models ignore sequence context and cooperative, long-range (TF-complex) regulation, which leads to false positives.

**Method.**
1. **Sequence features.** DNA sequences (NCBI) → DNABERT; protein sequences (UniProt) → ESM-2. Both are average-pooled and frozen. The DNA region used and the embedding sizes are not specified.
2. **High-order graph (HIN).**
   - A2 = A·A. A2n = binary(A2) − A, where binary() is printed as a threshold > 1.
   - HIN = A_PPI ∘ A2n: an intersection with a dataset-specific HuRI PPI adjacency.
   - 2-hop is best (Table 4).
3. **GAT per modality** (expression, DNA, protein), with separate masked attention over GRN neighbours and HIN neighbours.
4. **Similarity co-attention.** co_att = softmax(H_G H_Gᵀ + H_H H_Hᵀ + H_G H_Hᵀ). This is followed by HAN-style view weights β1/β2 and fusion β1·H_GRN + β2·H_HIN.
5. **Output.** K = 3 heads are concatenated (Table 5), the modalities are concatenated (Table 6), then an MLP and a dot-product decoder.
- **Loss:** BCE.
- **Not stated in paper:** lr, optimizer, hidden dims, dropout.

**Datasets.** The 7 BEELINE sets × {TFs+500, TFs+1000} × 4 GTs = 44 scenarios (§2.0). Accessions are in the paper: GSE81252, GSE75748, GSE98664, GSE48968, GSE81682. There is **no dataset-statistics table**. A bladder-cancer spatial transcriptomics set (GSE171351) is used as a qualitative case study only, with hTFtarget priors and no metrics.

**GT.**
- STRING.
- Non-specific ChIP-seq (TRRUST v2, RegNetwork, DoRothEA).
- Cell-type-specific ChIP-seq (RegNetwork, ChIP-Atlas, ENCODE).
- LOF/GOF (ESCAPE), mESC only.

**Protocol.**
- "10 independent tests", means only, **no std**, paired t-tests.
- **Split ratios, negative sampling and prior-graph masking are not stated.** The GENELink split is probably used (reviewer inference).

**Results.** Tables 1 (AUROC) and 2 (AUPRC); CaHoT-GRN value / best of 7 baselines.

| GT | Dataset | AUROC TFs+500 | AUROC TFs+1000 | AUPRC TFs+500 | AUPRC TFs+1000 |
|---|---|---|---|---|---|
| STRING | hESC | 0.91 / GATCL 0.91 | 0.92 / GCLink,GATCL 0.91 | 0.33 / GATCL 0.29 | 0.34 / GATCL 0.30 |
| STRING | hHEP | 0.93 / GCLink 0.91 | 0.92 / GCLink 0.92 | 0.38 / 0.30 | 0.33 / GCLink 0.33 |
| STRING | mDC | 0.92 / 0.90 | 0.93 / GCLink 0.91 | 0.46 / GCLink 0.41 | 0.50 / GCLink 0.43 |
| STRING | mESC | 0.93 / GCLink 0.91 | 0.93 / 0.92 | 0.34 / GCLink 0.27 | 0.34 / GATCL 0.25 |
| STRING | mHSC-E | 0.91 / 0.90 | 0.90 / 0.90 | 0.28 / GATCL 0.35 | 0.28 / GATCL 0.38 |
| STRING | mHSC-GM | 0.89 / 0.89 | 0.86 / 0.88 | 0.49 / GATCL 0.43 | 0.37 / GATCL 0.40 |
| STRING | mHSC-L | 0.81 / GENELink 0.82 | 0.81 / GeneLink+ 0.84 | 0.26 / GeneLink+ 0.42 | 0.24 / GeneLink+ 0.42 |
| Non-spec ChIP | hESC | 0.72 / GCLink 0.73 | 0.68 / 0.72 | 0.05 / 0.06 | 0.03 / GATCL 0.06 |
| Non-spec ChIP | hHEP | 0.74 / 0.73 | 0.72 / GATCL 0.73 | 0.05 / GeneLink+ 0.07 | 0.04 / 0.06 |
| Non-spec ChIP | mDC | 0.81 / 0.81 | 0.79 / GATCL 0.81 | 0.21 / scMGATGRN 0.19 | 0.17 / GATCL 0.16 |
| Non-spec ChIP | mESC | 0.79 / 0.80 | 0.82 / scMGATGRN 0.82 | 0.06 / GATCL 0.09 | 0.07 / scMGATGRN 0.09 |
| Non-spec ChIP | mHSC-E | 0.79 / 0.76 | 0.76 / GATCL 0.77 | 0.15 / GATCL 0.21 | 0.19 / GATCL 0.22 |
| Non-spec ChIP | mHSC-GM | 0.76 / GENELink 0.76 | 0.74 / GENELink 0.77 | 0.27 / GATCL 0.27 | 0.27 / GATCL 0.31 |
| Non-spec ChIP | mHSC-L | 0.73 / GENELink 0.67 | 0.73 / 0.67 | 0.11 / GeneLink+ 0.17 | 0.10 / GeneLink+ 0.43 |
| LOF/GOF | mESC | 0.84 / GATCL 0.83 | 0.79 / 0.80 | 0.63 / GATCL 0.63 | 0.51 / GeneLink+ 0.53 |
| Cell-type-spec ChIP | hESC | 0.88 / 0.87 | 0.88 / 0.88 | 0.61 / scMGATGRN 0.60 | 0.62 / scMGATGRN 0.61 |
| Cell-type-spec ChIP | hHEP | 0.90 / 0.89 | 0.90 / 0.89 | 0.83 / 0.82 | 0.83 / GeneLink+ 0.82 |
| Cell-type-spec ChIP | mDC | 0.75 / GCLink 0.77 | 0.76 / GCLink 0.80 | 0.18 / GATCL 0.13 | 0.16 / GCLink 0.16 |
| Cell-type-spec ChIP | mESC | 0.92 / GeneLink+ 0.92 | 0.93 / GeneLink+ 0.93 | 0.83 / GeneLink+ 0.85 | 0.85 / GeneLink+ 0.85 |
| Cell-type-spec ChIP | mHSC-E | 0.92 / GeneLink+ 0.92 | 0.93 / GeneLink+ 0.93 | 0.94 / GeneLink+ 0.94 | 0.95 / GeneLink+ 0.95 |
| Cell-type-spec ChIP | mHSC-GM | 0.93 / GeneLink+ 0.93 | 0.93 / GeneLink+ 0.93 | 0.94 / GeneLink+ 0.94 | 0.94 / 0.93 |
| Cell-type-spec ChIP | mHSC-L | 0.87 / GeneLink+ 0.87 | 0.88 / 0.86 | 0.86 / GeneLink+ 0.87 | 0.86 / GeneLink+ 0.88 |

**Averages over 22 scenarios (Table 1/2 "Average" rows; the Table 2 row is mislabelled "Average AUROC").**

| Method | AUROC 500 | AUROC 1000 | AUPRC 500 | AUPRC 1000 |
|---|---|---|---|---|
| **CaHoT-GRN** | **0.846** | **0.842** | **0.420** | **0.409** |
| GATCL | 0.827 | 0.831 | 0.393 | 0.391 |
| GCLink | 0.818 | 0.818 | 0.385 | 0.375 |
| scMGATGRN | 0.815 | 0.82 | 0.353 | 0.368 |
| GENELink | 0.815 | 0.817 | 0.342 | 0.334 |
| GeneLink+ | 0.811 | 0.805 | 0.375 | 0.380 |
| GNE | 0.696 | 0.676 | 0.223 | 0.222 |
| DeepRIG | 0.564 | 0.624 | 0.239 | 0.237 |

- **Paired t-test p-values:** all < 0.05 except AUPRC TFs+1000 vs GeneLink+ (4.37E-01) and vs GATCL (3.33E-01).
- **Ablations.**
  - Table 4 (HIN order), TFs+500 AUROC/AUPRC: S-GRN 0.832/0.406, PPI 0.834/0.396, 1-hop 0.837/0.407, 2-hop 0.846/0.420, 3-hop 0.834/0.407.
  - Table 5 (heads): K = 1 0.824/0.385 … K = 3 0.846/0.420.
  - Table 6 (fusion): max 0.827/0.408, mean 0.820/0.378, concat 0.846/0.420.
  - The Fig. 2 component ablation is figure only. Approximate STRING AUROC: full ≈0.901, without sequence features ≈0.878.

**Baselines.** GeneLink+, scMGATGRN, GCLink, GATCL, GENELink, GNE, DeepRIG, all at default parameters.

**Limitations.**
- Author-stated: rigid hard-intersection masking causes artefacts on sparse networks; similarity attention can under-fuse.
- (Reviewer notes)
  - The split, negative sampling and whether test edges are removed from the prior graph A (which feeds both the GAT and A²) are unspecified, so leakage is possible.
  - No std is reported.
  - Wins over the best baseline, counted cell by cell: TFs+500 AUROC 10 win / 8 tie / 4 loss; TFs+1000 AUROC 6/7/9; TFs+500 AUPRC 9/4/9; TFs+1000 AUPRC 7/4/11.
  - The "% improvements" in the text are absolute point differences.
  - The GNE TFs+1000 average printed as 0.676 recomputes to 0.705.
  - HuRI is human-only, yet 5/7 datasets are mouse; ortholog mapping is unexplained.

**Relevance to MeVD-GRN.** Same benchmark as scUniGP but lower numbers. **scUniGP's averages dominate CaHoT's**, so beating scUniGP on BEELINE also beats CaHoT-GRN, provided the splits are equivalent. There is no ATAC data. A transferable idea is adding frozen DNA and protein language-model embeddings as node features.

---

### 2.4 scRegNet: "Prediction of gene regulatory connections with joint single-cell foundation models and graph-based learning" (`scregnet.pdf`)

**Citation.** Kommu S., Wang Y., Wang Y., Wang X. *Bioinformatics* 41 (ISMB/ECCB 2025 Supplement), i619–i627, 2025. doi:10.1093/bioinformatics/btaf217. **Code (in paper):** https://github.com/sindhura-cs/scRegNet.

**Task / modalities / family.** Supervised TF→target link prediction; scRNA-seq only; frozen single-cell foundation-model (scFM) gene embeddings plus a GNN over the prior graph.

**Method.**
- **Gene embeddings.** Three scFM backbones are used, all **frozen**: scBERT (200-d), Geneformer (20 layers, 896-d) and scFoundation (512-d). Gene embeddings are mean-pooled over all cells of the dataset.
- **Graph branch.** A GNN (GCN in the headline; GraphSAGE and GAT compared) runs over the training-edge prior graph with expression-derived initial features.
- **Fusion and decoder.** Z_joint = scFM ⊕ GNN embedding. Output P̂ = Softmax(FCN(MLP(Z_joint[i]) ⊕ MLP(Z_joint[j]))) (Eq. 6).
- **Loss and tuning.** BCE; hyperparameters tuned with Optuna. Settings are in Supp. Table S2, which is not in the PDF.

**Datasets.** The 7 BEELINE sets × {TFs+500, TFs+1000} (§2.0). Accessions in paper: GSE81252, GSE75748, GSE98664, GSE48968, GSE81682. **Only the cell-type-specific ChIP-seq GT is used.** No STRING, non-specific or LOF/GOF GTs.
- **Table 1 (verbatim).** Cells are listed as hESC 759, hHEP 426, mDC 384, mESC 422, mHSC-E 1072, mHSC-GM 890, mHSC-L 848. These are one more than other papers report, probably counting a header row.
- **TFs:** 34 (34), 30 (31), 20 (21), 88 (89), 29 (33), 22 (23), 16 (16).
- **Train / test sizes, TFs+500:** 20 677 / 7142 (hESC), 19 002 / 6563, 10 969 / 3792, 65 895 / 22 736, 13 632 / 4718, 9280 / 3216, 5976 / 2076.

**GT.** Cell-type-specific ChIP-seq (ENCODE, ChIP-Atlas, ESCAPE), following GNNLink (Mao et al. 2023).

**Protocol.**
- GENELink split: per TF, positives and negatives go 67% train / 33% test; 10% of train is used as validation for early stopping.
- The negatives are all non-target genes of each TF (hard-negative mode).
- Every TF appears in both train and test.
- Results are averaged over **50 independent evaluations**; std is printed as ±0.00–0.01.

**Results. Table 2 (TFs+500) and Table 3 (TFs+1000), cell-type-specific GT.** Best scRegNet variant vs best baseline (GENELink or GNNLink).

| Dataset | AUROC 500: scRegNet (Geneformer) / best baseline | AUPRC 500 | AUROC 1000 | AUPRC 1000 |
|---|---|---|---|---|
| hESC | 0.89±0.00 / GNNLink 0.85 | 0.62±0.00 / GNNLink 0.52 | 0.88±0.00 / GENELink 0.83 | 0.62±0.00 / GNNLink 0.51 |
| hHEP | 0.90±0.00 / GENELink 0.84 | 0.84±0.00 / GNNLink 0.75 | 0.90±0.00 (scFoundation 0.91) / GENELink 0.85 | 0.84±0.00 (scFoundation 0.85) / GNNLink 0.78 |
| mDC | 0.81±0.00 / GENELink 0.71 | 0.17±0.00 / **GNNLink 0.25** | 0.84±0.00 / GNNLink 0.78 | 0.17±0.01 / **GNNLink 0.21** |
| mESC | 0.93±0.00 / GENELink 0.88 | 0.86±0.00 / GNNLink 0.76 | 0.93±0.00 / GENELink 0.89 | 0.87±0.00 / GNNLink 0.78 |
| mHSC-E | 0.92±0.00 / GENELink 0.87 | 0.94±0.00 / GENELink 0.89 | 0.94±0.00 / GENELink 0.90 | 0.95±0.00 / GNNLink 0.93 |
| mHSC-GM | 0.93±0.00 / 0.89 | 0.94±0.00 / 0.89 | 0.94±0.00 / GNNLink 0.92 | 0.95±0.00 / GNNLink 0.93 |
| mHSC-L | 0.88±0.00 / GNNLink 0.84 | 0.88±0.00 / GNNLink 0.85 | 0.88±0.00 / GNNLink 0.86 | 0.87±0.00 / GNNLink 0.86 |

- **scBERT variant, TFs+500:** AUROC 0.88/0.90/0.75/0.92/0.92/0.92/0.85; AUPRC 0.61/0.83/0.12/0.84/0.94/0.93/0.85.
- **scFoundation variant, TFs+500:** AUROC 0.89/0.90/0.81/0.93/0.92/0.93/0.88; AUPRC 0.62/0.83/0.15/0.86/0.94/0.94/0.88. The full table is in the paper.
- **Text claims.**
  - Geneformer variant, TFs+500: +7.4% / +6.9% AUROC and +18.6% / +4.1% AUPRC over GNNLink / GENELink.
  - TFs+1000: +6.2% / +7.5% AUROC and +16.5% / +3.9% AUPRC over GENELink / GNNLink.
  - (Reviewer note) The baseline labels are swapped between the two sentences, and the claim of beating the baselines "across all" datasets fails for mDC AUPRC, where GNNLink wins: 0.25 vs 0.17 at TFs+500 and 0.21 vs 0.17 at TFs+1000.
- **Other results.** The GNN-architecture comparison (GCN / SAGE / GAT; supplementary, and thesis Table 3.6) and the label-noise robustness experiment (Fig. 3) are not reproduced here. The latter is figure only.

**Baselines.** GRNBoost2, GENIE3, PCC, GRN-Transformer, DeepDRIM, CNNC, GNE, GENELink, GNNLink.

**Limitations.**
- (Reviewer notes)
  - Only the dense cell-type-specific GT is used (density 0.08–0.58), which favours high AUPRC.
  - The ±0.00 std over 50 runs is implausibly small.
  - The split is transductive, with all TFs seen in training.
  - The scFM embeddings are mean-pooled, so cell-level signal is lost.

**Relevance to MeVD-GRN.** On the cell-type-specific GT, scRegNet (Geneformer) comes second only to scUniGP. scUniGP is equal or higher on every TFs+500 Specific cell: AUROC 0.813–0.941 vs 0.81–0.93, and AUPRC equal or higher. scRegNet is the stronger reference for the TFs+1000 Specific AUROC cells, which scUniGP does not report numerically. The comparison is RNA-only. Frozen scFM gene embeddings as node features would be a cheap addition to MeVD-GRN.

---

### 2.5 Kommu MS thesis: "Towards Network-Guided Large-Scale Foundation Models on Single-Cell Transcriptomics" (`ref-theses/thesis_found_model_scRNA.pdf`), brief

**Citation.** Sindhura Kommu, MS in Computer Science and Applications, Virginia Tech, May 2025. Advisor: Xuan Wang.

**Contents.**
- **Ch. 3 is scRegNet (§2.4).** Tables 3.2/3.4/3.5 are identical to the paper's Tables 1–3.
  - The thesis adds, in the main body, the GNN-architecture comparison (Table 3.6: GCN/SAGE/GAT, 3 decimals) and the noise-robustness experiment (figures only).
  - (Reviewer note) Table 3.6 disagrees with Table 3.5 in two cells. At 1% label noise, the robustness figure shows hESC AUROC ≈0.59 (figure only), far below the 0.89 noise-free value. The robustness runs may have used a different setup.
- **Ch. 4, scNetFormer (new, not in any paper).** Continued pretraining of a 6-layer Geneformer, injecting GRAND tissue GRNs with early / intermediate / late fusion.
  - Evaluated on Geneformer-style gene-classification tasks, **not edge prediction**.
  - Best preliminary AUC 0.87 vs 0.81 for the baseline Geneformer (Tables 4.1/4.2).
  - Downstream data sources are not given. This includes about 30k heart endothelial cells with no accession.

**Datasets.** The same 7 BEELINE sets with the cell-type-specific GT (§2.0), all RNA only. GRAND networks are used for pretraining (the list of 18 is not given).

**Relevance to MeVD-GRN.** There are no new benchmark numbers beyond scRegNet. It supports the idea of injecting GRN priors into a foundation model.

---

### 2.6 SMOGT (file `SCMOGT.pdf`): "Deciphering hierarchical regulatory network of cell fate via an epigenetics-informed heterogeneous graph transformer on single-cell multi-omics data"

**Citation.** Huang Y., Liu C. (equal contribution), Yang Z., Liu B., Zhai X., Zheng J., Xiao J., Song T. *Briefings in Bioinformatics* 26(6), bbaf664, 2025. doi:10.1093/bib/bbaf664. **Code (in paper):** https://github.com/YuHongHuang-lab/SMOGT.

**Task / modalities / family.**
- **Task.** Build a hierarchical TF-TF → TF-CRE → CRE-CRE → CRE-TG network ("HRNet"). **The benchmarked edges are TF→CRE (peak) binding and CRE-CRE (long-range contacts), not TF→target-gene edges.**
- **Modalities.** Paired scRNA-seq + scATAC-seq, plus priors. The priors include ENCODE epigenomic resources and STRING for TG-TG; seven resource types in total (Fig. 1A).
- **Model family.** Heterogeneous graph transformer (HGT) with meta-path-guided message passing, trained semi-supervised.

**Method.**
- **Node types:** TF, CRE and TG. Initial features are the paired expression and accessibility profiles on SEACells metacells, after a Seurat pipeline.
- **Encoders:**
  - TF-TF: GCN, which was empirically better than GAT or HGT.
  - TF-CRE, CRE-CRE, CRE-TG: HGT layers with gated residual H = θ·ReLU(H̃) + (1−θ)·H^(l−1).
  - Loops: GraphSAGE.
  - Information flows along the meta-path TF-TF → TF-CRE → CRE-CRE → CRE-TG.
- **Dual-task decoder:**
  - Edge reconstruction by dot product plus sigmoid.
  - A 4-layer FCNN that regresses TG expression (embedding_dim 16, hidden_dim 32).
- **Loss:**
  - L = λ·MSE(expression) + (1−λ)·BLCE(edges), where BLCE is class-weighted BCE with weight α (Eqs 5–7).
  - Semi-supervised with adversarial training. Pseudo-positive TF-CRE and CRE-CRE edges are highly correlated pairs above a per-dataset PCC threshold (Supp. Note S8). Pseudo-negatives are added at controlled ratios.
- **Downstream modules:** BioStreamNet expression prediction and in-silico perturbation, MRWR_Pert and MRWR_CTS random walks for driver regulators, and Louvain co-CRE modules.

**Datasets.**
- Paired scRNA+scATAC from GEO, 10x Genomics and ENCODE:
  - (i) bone marrow (BM)
  - (ii) healthy PBMC
  - (iii) human cerebral cortex
  - (iv) leukemia stem cells
  - (v) K562 and HCT116 cell lines
  - (vi) A549 and GM12878 cell lines
  - (vii) K562 Perturb-seq
- **Accessions, cell counts and preprocessing are only in Supplementary Table S1 / Note S1, which are not in the PDF.** Accessions: not in paper.
- Auxiliary data: GWAS SNPs (OpenGWAS, GWAS Catalog), PBMC eQTL (eQTLGen, FDR < 1e-150), cortex eQTL (GTEx), super-enhancers (SEanalysis), TCGA AML (prognosis case study).

**GT.** ChIP-seq for TF-CRE (BM, K562, HCT116, A549, GM12878; Supp. Table S2) and Hi-C for CRE-CRE (Supp. Table S3). Sources and thresholds are in Supp. Note S1, which is not available.

**Protocol.**
- Per-TF AUPR and max-F1 of TF-CRE predictions against ChIP-seq.
- Precision and recall at top-N pairs per TF vs STREAM.
- CRE-CRE AUC/AUPR against Hi-C.
- Thresholds: TF-CRE 0.8, CRE-CRE 2% (for the expression-prediction comparison).
- It is not a held-out-edge supervised split; pseudo-labels come from correlation.

**Results (all figure only).**
- **TF-CRE AUPR / max-F1 per TF (Fig. 2A, boxplots)** on BM, K562, HCT116, A549 and GM12878 vs regX, LINGER, TRIPOD, Pando, REUNION and Pearson. Not reported numerically. The text claims "significant superiority"; visually, SMOGT's box sits highest in each panel, with y-axes spanning 0–0.8.
- **Top-N precision vs STREAM (Fig. 2B).** SMOGT is better at top-200/300/40 in BM, A549 and GM12878, matches STREAM at top 60–200 in K562, and is worse in HCT116 (text).
- **CRE-CRE vs Hi-C (Fig. 2C).** The legends print AUCs; read from a 500-dpi render, so figure only, but the numbers are printed:

| Cell line | ROC-AUC SMOGT / Cicero / ArchR | PR-AUC SMOGT / Cicero / ArchR |
|---|---|---|
| K562 | 0.787 / 0.576 / 0.585 | 0.698 / 0.340 / 0.231 |
| HCT116 | 0.861 / 0.538 / 0.64 | 0.219 / 0.042 / 0.035 |
| A549 | 0.776 / 0.512 / 0.580 | 0.567 / 0.294 / 0.241 |
| GM12878 | 0.681 / 0.521 / (illegible) | 0.615 / 0.437 / 0.378 |

- **Other results.**
  - TG-expression PCC vs regX, SCARlink and TSS-proximal CREs on BM, PBMC and A549 (Fig. 3D–E): figure only.
  - Driver-TF identification vs SCENIC+, CEFCON and STREAM on BM: counts of known lineage TFs in the top 5/10 (text), e.g. HSC top 5 has 4 specific TFs for SMOGT vs 3 for SCENIC+.

**Baselines.** regX, LINGER, TRIPOD, Pando, REUNION, Pearson (TF-CRE); STREAM (top-N); Cicero, ArchR (CRE-CRE); SCARlink, regX (expression); SCENIC+, CEFCON, STREAM (driver TFs).

**Limitations.**
- (Reviewer notes)
  - No TF→gene benchmark.
  - All quantitative results are figure-only.
  - Datasets and accessions are hidden in unavailable supplements.
  - Pseudo-labels come from correlation, and the model is scored partly against correlated structure.
  - Heavy downstream biology with little ablation. The homogeneous-graph ablation is in Supp. Fig. S3.

**Relevance to MeVD-GRN.**
- Its datasets are paired RNA+ATAC, which suits MeVD-GRN. PBMC and BM overlap with scTFBridge, scMultiomeGRN and EpiAwareNet. K562 and GM12878 overlap with SC-MO-GRN-DB (DS025 K562 scRNA/scATAC; RN117/RN201 ChIP-seq networks).
- Its benchmark targets peaks, not genes, so a direct numeric comparison needs a TF→peak head.
- If MeVD-GRN adds peak-level edges, the ChIP-seq TF-CRE AUPR (per TF) on K562/GM12878 would be the comparable protocol. SC-MO-GRN-DB provides K562 and GM12878 ChIP-seq networks.

---

### 2.7 scMultiomeGRN: "Deep learning-based cell-specific gene regulatory networks inferred from single-cell multiome data" (`scMultiomeGRN.pdf`)

**Citation.** Xu J., Lu C. (equal contribution), Jin S., Meng Y., Fu X., Zeng X., Nussinov R., Cheng F. *Nucleic Acids Research* 53, gkaf138, 2025. doi:10.1093/nar/gkaf138. **Code (in paper):** https://doi.org/10.5281/zenodo.14848389.

**Task / modalities / family.**
- **Task.** Semi-supervised link prediction on a **TF–TF network** (nodes are TFs only; non-TF target genes are excluded, which the authors list as a limitation).
- **Modalities.** scRNA-seq + scATAC-seq.
- **Model family.** GCN/GAT-style graph autoencoder with modality-specific neighbour aggregators and cross-modal attention.

**Method.**
- **"Ground-truth" / initial adjacency.**
  1. Remove peaks detected in < 10% of cells.
  2. Scan the remaining peaks with FIMO v5.1 (`--thresh 1e-6`).
  3. Draw edge TF_u → TF_v if a binding site of u lies in the promoter of v (TSS ± 2 kb). Redundant motifs are merged.
  4. This gives a binary |N_TF| × |N_TF| matrix.
  - **The labels are therefore motif- and ATAC-derived, not ChIP-seq.** This is despite the Methods sentence "we used the following three standard ChIP-seq data".
- **Node features.**
  - ATAC: MAESTRO regulatory-potential scores.
  - RNA: GRNBoost2 TF–gene importance vectors.
- **Edge features.** A 16×16 joint histogram of log expression for each TF pair (CNNC-style), log-normalised (Eqs 1–2).
- **Encoder.**
  - For each modality k: messages m = W1·x_v + W_e·e_uv + b, with attention α = LeakyReLU(aᵀ[W2 x_v ‖ m]) and softmax over neighbours.
  - Self plus neighbourhood representations are concatenated and passed through a dense layer with ELU. L layers are stacked.
  - Cross-modal attention: β_k = softmax(a_mᵀ W_m z_k + b_m); h = Σ β_k z_k.
- **Decoder.** Â = σ(h hᵀ). **Loss:** BCE.
- **Training.** PyTorch 2.4, Adam, lr 1e-5, ≤ 2000 iterations, early stopping.
- **Output.** 10 runs; each is binarised at the median test weight; edges kept if present in ≥ 6/10 runs.

**Datasets.**

| Dataset | Source / accession (as in paper) | Species / tissue | Modalities | Cells | Use |
|---|---|---|---|---|---|
| Human fetal lung | GEO **GSM4508936** | human fetal lung, 59 samples, 89–125 days post-conception, 9 cell types | scATAC + scRNA ("lung single-cell multiome dataset") | 72 622 | main benchmark (Fig. 2) |
| Human adult hematopoietic differentiation (Buenrostro et al. 2018) | scRNA: Data S2 of Buenrostro 2018 (Cell); scATAC: github.com/pinellolab/scATAC-benchmarking/…/Buenrostro_2018 | human HSC, CMP, GMP | scATAC + separately profiled RNA (**unpaired**) | not stated | second benchmark (Supp. Fig. S1) |
| 10x PBMC | scATAC: 10x "10k PBMCs from a healthy donor" (ATAC v1); scRNA: 10x `pbmc_10k_v3` | human PBMC; CD14+ Mono 3447, CD4.Memory 874, CD4.Naive 589, CD8.effector 454, CD8.Naive 361 (scATAC clusters) | **two separate 10x datasets, not paired multiome** | 7189 ATAC / 9129 RNA | qualitative cell-type-specific GRNs (SPI1, KLF12, GATA3, EGR4) |
| AD brain snMultiome | Synapse **syn52293417** | human prefrontal cortex, 48 non-AD / 29 early / 15 late AD individuals | snRNA + snATAC | 414 000 snRNA / 437 000 snATAC nuclei; microglia 8600 | qualitative (SPI1, RUNX1) |

- **Preprocessing.** Peaks and genes detected in < 10% of cells are removed. PBMC annotation uses Seurat v3.
- (Reviewer note) GSM4508936 is a single GEO sample and is probably from the fetal sci-ATAC-seq3 atlas (Domcke et al. 2020). If so, the lung RNA and ATAC come from separate atlases and are **not cell-paired**. This has not been verified against GEO.

**Protocol.**
- "Tenfold" scheme: positives of the motif adjacency are split into 10 subsets. One subset plus an **equal number of random negatives (1:1)** forms the test set; everything else is training.
- The model is run per cell type. Metrics: accuracy, AUROC, AUPR, "average of multiple runs".
- There is no held-out-TF protocol, and the labels come from the same ATAC data used as features (reviewer note: circular).

**Results.**
- **Human fetal lung, averaged over 9 cell types** (text; per-cell-type values are figure only in Fig. 2B / Supp. Table S1, which is not in the PDF):

| Method | AUROC | AUPR |
|---|---|---|
| **scMultiomeGRN** | **0.924** | **0.790** |
| GENELink | 0.900 (P = 8.7e-03) | 0.758 (P = 2.3e-02) |
| DeepTFni | 0.871 (P = 1.1e-03) | 0.721 (P = 2.4e-02) |
| SCRIP | 0.638 | 0.644 |
| DeepWalk | 0.614 | 0.292 |
| scMTNI | 0.553 | 0.290 |
| GENIE3 | 0.540 | 0.242 |

- Test accuracy is "exceeding 0.83" on every lung cell type (Supp. Table S1).
- **Hematopoietic (HSC/CMP/GMP).** It beats the second-ranked DeepTFni "by at least 2.5% and 4.6% in AUROC and AUPR". Absolute values are figure only (Supp. Fig. S1), not reported numerically in the main text.
- **Omics ablation on lung lymphoid cells (Fig. 3).** The fused model beats RNA-only and ATAC-only; values are figure only. KL divergence between positive and negative score distributions: 0.7410 and 0.8572 (single-omics) vs 1.7182 (fused).
- **Robustness to masked positives and to cell number (Figs 2D–F).** Figure only. Performance plateaus above about 50 cells.

**Baselines.** DeepTFni, DeepWalk, SCRIP (ATAC-based); GENELink, GENIE3 (RNA-based); scMTNI.

**Limitations.**
- Author-stated:
  - Relies on data quality.
  - Filtering peaks in < 10% of cells sharply reduces the number of TFs.
  - Non-TF genes and causal TF→target direction are not modelled.
  - No spatial data.
- (Reviewer notes)
  - The labels are motif scans on the same scATAC, which is circular.
  - Balanced 1:1 test sets inflate AUPR.
  - The "multiome" datasets are mostly unpaired.
  - No per-cell-type numbers are in the main text.

**Relevance to MeVD-GRN.**
- This is the nearest architectural relative (dual-modality graph model).
- Its headline 0.924 / 0.790 is on a **motif-derived TF–TF network with 1:1 test negatives**. To compare, MeVD-GRN would have to reproduce that label construction (FIMO 1e-6, TSS ± 2 kb, TF–TF only), which is easy to script.
- (Reviewer note) Beating it on its own protocol is feasible but scientifically weak. Report it alongside a ChIP-seq-labelled evaluation such as SC-MO-GRN-DB, or scTFBridge's Cistrome PBMC set.
- Data effort is low for PBMC (two 10x downloads, a few hundred MB) and lung (one GEO sample). The AD Synapse set is large and requires registration.

---

### 2.8 scTFBridge: "a disentangled deep generative model informed by TF-motif binding for gene regulation inference in single-cell multi-omics" (`scTFBridge.pdf`)

**Citation.** Wang F.-a., Yi C., Chen J., He R., Liu J., Li Y. *Nature Communications* 16:9166, 2025. doi:10.1038/s41467-025-64227-y.
- **Code (in paper):** https://github.com/FengAoWang/scTFBridge (package) and https://github.com/FengAoWang/scTFBridge_reproduce.
- **TF-motif prior:** https://doi.org/10.5281/zenodo.15954450.

**Task / modalities / family.**
- **Task.** Unsupervised / self-supervised inference of cis (RE→TG) and trans (TF→TG) regulatory scores, plus cross-modal generation.
- **Modalities.** Paired scRNA-seq + scATAC-seq (10x Multiome-type).
- **Model family.** Disentangled multimodal VAE.

**Method.**
- **Latent spaces.** A modality-shared latent z_s plus modality-private latents z_pr (RNA) and z_pa (ATAC).
- **Shared latent = TFs.** The shared latent has **128 dimensions, one per TF**: the 128 most highly expressed TFs in scRNA.
- **Motif-masked ATAC decoder.** The ATAC decoder is a single linear layer masked by a TF×RE motif binding-score matrix B, so W'_ij = W_ij·B_ij (Eq. 18).
  - Motifs: 713 PWMs from the PECA2 repository (JASPAR / TRANSFAC / UniPROBE).
  - Peaks scanned with HOMER.
- **Disentanglement.** Mutual information between shared and private latents, and between the two private latents, is minimised with the CLUB upper bound (Eqs 10–17). Cross-modal generation uses a Product-of-Experts.
- **Loss.** L = L_VAE + β·L_MI + γ·L_cross-modal (Eq. 21). A contrastive CMI term is described in Results but has no equation (reviewer note).
- **Regulatory scores.**
  - Trans (TF→TG): SHAP attributions from the TF latents to imputed TG expression.
  - Cis: REscore_ij = ρ_ij · O_j · E_i · e^(−d_ij/d0) (Eq. 22), where ρ is SHAP, O mean accessibility, E mean expression and d the RE–TSS distance.
- **Features.** 3000 target genes; peaks open in ≥ 1% of cells; cells present in both modalities.
- **Training.** 5-fold CV (64/16/20 per the earlier reading pass); SHAP computed on test cells.

**Datasets (Data availability).**

| Dataset | Accession (as in paper) | Tissue | Modalities | Cells |
|---|---|---|---|---|
| Human BMMC (NeurIPS 2021 multimodal benchmark) | GEO **GSE194122** | bone marrow mononuclear cells | paired 10x Multiome RNA+ATAC | not in paper |
| Human PBMC 10k | https://scglue.readthedocs.io/zh-cn/latest/data.html (the 10x `pbmc_granulocyte_sorted_10k` Multiome demo, identified online, not in paper) | PBMC | **paired** 10x Multiome | not in paper |
| Rheumatoid-arthritis synovium | GEO **GSE243917** | synovial tissue | paired multiome | not in paper |
| Lung, ever- and never-smokers | printed as "**GSE2414468**" (reviewer note: invalid length; probably GSE241468) | tumour-distant normal lung | paired multiome | not in paper |

- **GWAS:** RA summary statistics GCST90132222 (for S-LDSC heritability).

**GT (PBMC only).**
- **TF→TG:** Cistrome ChIP-seq regulatory scores for **17 TFs across 4 blood cell types**. The IDs are in Supplementary Data 6, which is not in the PDF.
- **Cis RE→TG:** single-cell eQTL from scQTLbase for 8 immune cell types (Supp. Data 5), and promoter-capture Hi-C (Javierre et al. 2016) from 3 primary blood cell types.

**Protocol.**
- Top-k predicted targets per TF are scored against the ChIP-derived targets, over the intersection of genes kept by all methods.
- Metrics: AUC, "AUPR ratio" (AUPR divided by the positive fraction), and F1 (vs SCENIC+, which gives unranked pairs).
- Cis evaluation is stratified by TSS distance bins.

**Results.**
- **TF→TG, STAT1 in CD14 monocytes (Fig. 5b/c; numbers printed in the legends, verified from a 300-dpi render):**

| Method | AUC | AUPR ratio |
|---|---|---|
| scTFBridge | 0.693 | **2.221** |
| LINGER | **0.719** | 2.105 |
| scMTNI | 0.561 | 1.476 |
| scMultiomeGRN | 0.556 | 1.226 |
| SCENIC+ | 0.509 | 1.041 |

- **TF→TG across all 17 TFs (Fig. 5d/e boxplots): figure only.**
  - Approximate medians: AUC scTFBridge ≈0.63 vs **LINGER ≈0.67**; AUPR ratio ≈1.6 vs **LINGER ≈1.95**.
  - The significance brackets (** between Ours and LINGER) mark a difference, and **LINGER is the higher method**. SCENIC+, scMTNI and scMultiomeGRN medians are ≈0.51–0.53 AUC.
  - (Reviewer note) The text's claim that scTFBridge "consistently surpassed the other baseline methods" holds only against those three, not LINGER. The authors note LINGER uses external bulk pretraining.
- **Cis RE→TG (Fig. 4): figure only.**
  - vs eQTL in naive CD4 T cells: highest AUC and AUPR ratio across all TSS-distance groups; F1 > SCENIC+ (paired t-test p = 0.00051, n = 6 distance groups).
  - Other 7 cell types: highest overall ranking (dot plots, not reported numerically); F1 vs SCENIC+ p = 0.0016.
  - vs pcHi-C: highest AUC across distance groups.
- **Integration, generation and ablation results:** figure-only or in the supplement, which was not available.

**Baselines.**
- TF→TG: LINGER, scMTNI, scMultiomeGRN, SCENIC+.
- Cis: PCC, Distance, Random, SCENIC+.
- Integration: scGLUE, scJoint, MultiVI, etc. (named in the intro).

**Limitations.**
- (Reviewer notes)
  - Only 17 TFs have ChIP ground truth.
  - The shared space is limited to 128 TFs.
  - Nearly all numbers are figure-only.
  - Cell counts are not stated.
  - One accession is malformed.
  - It is unsupervised, so it is not directly comparable with supervised edge prediction.

**Relevance to MeVD-GRN.**
- This is **the most relevant ChIP-labelled paired-multiome TF→TG benchmark in the corpus**: 10x PBMC 10k Multiome with Cistrome ChIP for 17 TFs in 4 blood cell types. MeVD-GRN needs both RNA and ATAC, which this dataset has, and the download is small (about 100s of MB).
- Because MeVD-GRN is supervised, the fair protocol is **held-out TFs**: train on TFs other than the 17, then score the 17 TFs' target rankings with AUC / AUPR ratio.
- The bar to beat is **LINGER** (STAT1 CD14-Mono AUC 0.719), not scTFBridge (0.693).
- BMMC GSE194122 is also paired, but scTFBridge reports no GRN ground truth for it. SMOGT also uses a bone-marrow multiome set, but its accession is only in SMOGT's unavailable Supp. Table S1, so it may not be the same data.

---

### 2.9 EpiAwareNet: "Prior-Guided Multi-Omic Transformers for Single-Cell Gene Regulatory Network Inference" (`epiawarenet.pdf`)

**Citation.** Xu T., Liu T., Rayamajhi N., Patrick R., Varala K., Li Y., Gao J. (Purdue). KDD 2026 (Jeju), arXiv:2606.00685 (30 May 2026). **Code (in paper):** https://github.com/tianyang-x/EpiAwareNet_pub.

**Task / modalities / family.**
- **Task.** Weakly-supervised regulator→target ranking over the full candidate space E = T × G.
- **Modalities.** **Paired** scRNA + scATAC (same cells).
- **Model family.** Two-stage multi-omic Transformer.

**Method.**
- **Stage 1 (backbone f_θ, trained per dataset from scratch).**
  - Gene–gene attention plus sparse, candidate-constrained gene–peak cross-attention with Top-K routing (K = 8). Candidate peaks are within **±5 kb** of the gene, ranked by distance, capped at 50 per gene.
  - Pretrained by masked-gene reconstruction (15% masked) with a negative-binomial loss.
  - AdamW, lr 1e-4, weight decay 1e-2, batch 64, seed 42, early stopping on reconstruction loss.
- **Stage 2 (backbone frozen).**
  - A one-hidden-layer MLP head (256, ReLU, dropout 0.1) on concatenated regulator–target embeddings.
  - BCE with prior train edges as positives and **10:1 uniformly sampled unlabeled pairs as negatives**.
  - AdamW, lr 3e-4, weight decay 1e-2, 20 epochs.
- **nnPU variant.** A non-negative positive–unlabeled loss.

**Datasets. Table 1 (verbatim); all paired multiome.**

| Dataset | Species | Cells N | Genes G | Peaks P | Prior edges \|E\| | Regulators \|T\| | Source |
|---|---|---|---|---|---|---|---|
| mN (tomato root, low nitrogen) | tomato | 8,748 | 34,074 | 53,570 | 170,300 | 2,381 | Patrick et al. 2026 bioRxiv 2026.02.06.704465; accession not in paper |
| pN (tomato root, normal N) | tomato | 12,448 | 34,074 | 55,930 | 170,300 | 2,381 | same |
| PBMC | human | 12,012 | 36,600 | 111,856 | 5,092 | 118 | 10x "10k Human PBMCs Multiome v1.0, Chromium Controller" |
| Mouse brain | mouse | 4,881 | 32,284 | 42,755 | 5,064 | 819 | 10x "Fresh Embryonic E18 Mouse Brain 5k" Multiome |

- **Preprocessing.** RNA: library-size normalisation plus log1p. ATAC: dataset-specific normalisation plus log.

**GT / priors.**
- PBMC: **DoRothEA**.
- Mouse brain: **TRRUST**.
- Tomato: a "bulk-derived" regulator–target prior shipped with the tomato benchmark (construction not documented).

**Protocol.**
- A **random 80/20 split of prior edges** (train / "val"). Metrics are on the held-out 20%, ranked against **all T × G pairs**.
- No checkpoint selection on val edges; the final checkpoint is reported.
- Seed 42. Table 7 adds 3 more seeds.
- This is **not a held-out-TF split**: regulators are shared.

**Results. Table 2 (±5 kb default).**

| Dataset | EpiAwareNet AUPRC | AUROC | P@100 | AUPRC ratio | Best baseline (AUPRC / AUROC) |
|---|---|---|---|---|---|
| PBMC | **0.00244** | **0.8218** | 0.03 | 8.0970 | Pando 0.00064 / GRNBoost2 AUROC 0.5713; RNA-only ablation 0.00165 / 0.7782 |
| Mouse brain | **0.00057** | **0.7263** | 0.01 | 6.4543 | scGPT+head 0.00015 / 0.5936 |
| mN tomato | **0.08458** | **0.7415** | 0.13 | 2.6939 | GRNBoost2 0.04505 / 0.6534 |
| pN tomato | 0.07451 | **0.7340** | 0.16 | 2.3756 | GRNBoost2 0.04435 / 0.6435; its own RNA-only ablation has higher AUPRC, 0.08173 |

- **Table 6 (gene–peak window).**
  - PBMC AUPRC / AUROC: ±5 kb 0.002440 / 0.8218; ±50 kb 0.002864 / 0.8673; ±100 kb **0.002933 / 0.8708**.
  - Mouse: ±5 kb 0.000570 / 0.7263; ±50 kb **0.001018** / 0.7662; ±100 kb 0.000779 / **0.8002**.
- **Table 7 (seeds 43–45).** PBMC AUPRC 0.002475 ± 4.55e-5, AUROC 0.8232 ± 7.57e-4. Mouse brain 0.000563 ± 1.0e-6, 0.7317 ± 1.63e-3.
- **Baselines (Table 2):** WGCNA, GRNBoost2, scGPT + matched head, scGLUE + matched head, SCENIC+, Pando. SCENIC+ and Pando score ≈0.50 AUROC because unreturned edges get score 0.

**Limitations.**
- (Reviewer notes)
  - The labels are curated priors (DoRothEA/TRRUST) with a random edge split, so the model learns the prior's hub structure.
  - AUPRC is tiny because of the huge candidate space.
  - The implied scored-pair count (AUPRC / AUPRC-ratio) is smaller than T × G in Table 1, so the exact evaluation space should be taken from the repo.
  - The scGPT baseline gives identical numbers on mN and pN (0.03139 / 0.4998 / 0.9999), which looks like a failed run.

**Relevance to MeVD-GRN.** This is **directly compatible**: paired RNA+ATAC, a public 10x PBMC 10k Multiome and E18 mouse brain 5k (small downloads), and simple priors. MeVD-GRN can be trained on the same 80% DoRothEA / TRRUST edges and scored the same way. The targets to beat are PBMC AUROC 0.8218 (±5 kb) or 0.8708 (±100 kb), and mouse brain AUROC 0.7263 / 0.8002. PBMC 10k Multiome is the same 10x demo family as scTFBridge's PBMC (not necessarily the identical file).

---

### 2.10 SC-MO-GRN-DB: "A comprehensive repository for single-cell multiomic gene regulatory networks" (`SC_MO_GRN_Db.pdf`), the project's training database

**Citation.** Valensi H., Karamveer K., Moeller E., Ozdogan S.E., Edwards R.M., Uzun Y. *iScience* 29, 115323, 17 Apr 2026. doi:10.1016/j.isci.2026.115323.
- **Portal:** https://scmogrndb.psu.edu.
- **Build scripts:** https://github.com/UzunLab/SC-MO-GRN-DB.

**Type.** Resource paper. **It reports no GRN-method benchmark numbers** (no AUROC/AUPRC anywhere in the text), so there is no number to beat. Its value is data plus reference networks.

**Contents.**
- **Reference networks.**
  - More than 22 million TF→TG edges.
  - Evidence types: localization (ChIP-seq/ChIP-chip), perturbation (knockdown/knockout/overexpression), dual evidence (localization ∩ perturbation in the same cell type), literature-curated, text-mined and PPI.
  - Human and mouse; non-specific plus cell-type-specific.
  - Cell types include ESC, HSC, DC, HepG2, K562, H1, BJ, GM12878, MCF7, macrophage, and literature networks for HSC/VSC/GSD/PSC/CAD/iPSC.
- **Single-cell datasets.**
  - More than 2 million cells across six modalities: scRNA, scATAC, scChIP, scDNAme, scHi-C, scCRISPR.
  - Seven joint assays.
  - Cell-type counts: ESC 8, K562 8, HSC 4, GM12878 3, BJ 2, MCF7 2, and one each of DC, H1, HepG2, macrophage (Fig. 4D).
- **Download sizes (from `docs/reference/dataset_info.md`):**
  - REFERENCE_NETWORKS_Human.zip 56.02 MB; REFERENCE_NETWORKS_Mouse.zip 14.27 MB.
  - SINGLE_CELL_DATASETS_Human.zip 17 481.61 MB; SINGLE_CELL_DATASETS_Mouse.zip 7188.55 MB.

**How the reference networks were built (STAR Methods).**
- **Localization.** Sources: BEELINE, ChIP-Atlas, Cistrome, ESCAPE and individual papers (Table S3). Each peak is assigned to the **nearest TSS** (edgeR `nearestTSS`, GENCODE v38; hg38 / mm39).
- **Perturbation.** From KnockTF 2.0, ESCAPE and primary papers.
- **Other networks.** PPI from STRING (via BEELINE); text-mined from TRRUST; literature Boolean models from BEELINE.
- **Filtering.** All networks are filtered to known TFs (AnimalTFDB-type list, ref. 50).
- **Datasets.** Only wild-type or control cells, from BEELINE, GEO, ArrayExpress, scPerturb and papers (Tables S4/S5).
- **Pairing caveat.** "These resources are not intended to be paired in a one-to-one" alignment; networks and datasets are matched by cell type.

**Key cross-links to the competitor benchmarks (from `docs/reference/dataset_info.md`).**
- **DS001–DS007 are the BEELINE scRNA sets:**
  - DS001 mESC, 422 cells
  - DS002 mDC, 384
  - DS003–005 mHSC, 1072 / 890 / 848
  - DS006 hESC, 759
  - DS007 hHEP, 426
  - These cell counts match scRegNet's Table 1.
- **BEELINE-derived reference networks:**
  - RN101 hESC ChIP-seq; RN102 HepG2 ChIP-seq; RN110 mouse DC ChIP-seq.
  - RN111 mESC ChIP-seq; RN112 mESC LOF/GOF; RN113 mouse HSC ChIP-seq.
  - RN005/RN006 human / mouse non-specific ChIP-seq; RN001/RN002 STRING-type PPI.
  - So the competitors' BEELINE benchmark is **already inside SC-MO-GRN-DB**, though the GENELink TFs+500/1000 split files must still be taken from the GENELink repo.
- **Paired RNA+ATAC datasets in SC-MO-GRN-DB:**
  - DS010 mESC scRNA/scATAC/scHiC, 9021 cells
  - DS011 mESC scRNA/scATAC/HiC, 36 535
  - DS012 mESC, 930
  - DS014 mESC, 68 794
  - DS025 K562, 434
  - DS026 macrophage, 3114
  - DS027 MCF7, 2995
  - With matching ChIP networks: RN111 / RN114–116 (mESC), RN117–119 (K562), RN204 (macrophage), RN205 (MCF7).
  - **None of these paired datasets is used by any competitor paper in this corpus.** That makes them the natural "bonus" results.

**Relevance to MeVD-GRN.** This is the training source. For head-to-head comparisons, the BEELINE subsets (DS001–007 with RN101/102/110–113) connect SC-MO-GRN-DB to scUniGP, CaHoT-GRN, ATFGRN and scRegNet, but only as RNA-only data. The K562/GM12878 networks (RN117, RN201) are the same cell lines SMOGT benchmarks on.

---

### 2.11 EpiXFormer: cross-attention network for cell-type-specific TF binding sites (`EpiXFormer.pdf`), adjacent, not GRN edge inference

**Citation.** Peng Y., Liu X., Wu J., Lin S., Zhan S., Li H., Wang J., et al. *Briefings in Bioinformatics* 27(1), bbaf721, 2026. doi:10.1093/bib/bbaf721.

**Task.** Binary classification of whether a ~200 bp genomic region is bound by a given DNA-binding protein in a given cell line.

**Inputs.** DNA sequence cross-attended with **bulk** epigenomic tracks (chromatin accessibility / histone). No single-cell data and no TF→gene edges.

**Data.**
- **199 protein × cell-line combinations** (32 TFs plus 11 other DNA-binding proteins, 43 DBPs in total) across 7 human ENCODE cell lines: GM12878, HCT116, HepG2, IMR-90, K562, MCF-7, SK-N-SH (Supp. Table S2).
- Plus bulk GEO sets **GSE167979, GSE147716, GSE90895, GSE132440, GSE216090** for external tests.
- ENCODE accessions are in Supp. Table S6, which is not in the PDF.

**Labels.** ChIP-seq peaks are positives. Negatives are random genomic regions plus unbound motif hits.

**Protocol.** A random 80/20 split of regions, plus ENCODE-DREAM-style cross-cell-line tests.

**Headline.** Mean **AUROC 0.9905, AUPRC 0.9608** over the 199 combinations (text). Per-pair values are in Supp. Table S3, which is not available. Baseline comparisons (DeepBind, FIMO, other methods; Fig. 2C–F) are figure only.

**Code.** The earlier reading pass found none stated in the paper.

**Relevance.** Out of scope for the benchmark targets. The idea of combining sequence with accessibility could serve MeVD-GRN as a TF→peak binding prior, which scMultiomeGRN and scTFBridge approximate with motif scans.

---

### 2.12 Out-of-scope sources (not GRN inference): short summaries

These sources were read in full. **None provides a dataset or number that MeVD-GRN can benchmark against.**

**(a) DECODE: "deep learning-based common deconvolution framework for various omics data"** (`decode.pdf`). Zhao T., Liu R., Sun Y., … Wang Y. *Nature Methods* 2026, doi:10.1038/s41592-026-03007-y.
- **Task.** Estimate cell-type and cell-state proportions in bulk samples across transcriptomics, proteomics and metabolomics, with scRNA/CITE-seq references.
- **Data (examples).**
  - Breast cancer scRNA GSE176078; mouse islet GSE211799; PBMC CITE-seq GSE253721 (RNA + protein, not ATAC).
  - Proteomics MSV000086809 / MSV000084110; single-cell metabolomics including PR001858.
  - Further GEO sets for application cohorts: GSE161865, GSE184869, GSE196941, GSE200356, GSE222550, GSE243906, GSE245467, GSE253217, GSE256501, GSE267916, GSE269058.
- **Ground truth.** Simulated pseudo-tissue proportions, or measured blood compositions.
- **Metrics.** CCC / RMSE / Pearson, all **figure only**.
- **Code.** Zenodo (the ID differs between the text and the reporting summary, per the earlier reading pass).
- **No scATAC and no GRN.** Out of scope.

**(b) "Gene regulatory networks: from correlative models to causal explanations"** (`Gene regulatory networks .pdf`). Maizels R.J., Briscoe J. *Nature Reviews Genetics* 2026, doi:10.1038/s41576-026-00939-1.
- A perspective with no experiments. Arguments relevant to MeVD-GRN's evaluation:
  1. Expression-based GRN inference often does no better than random or simple baselines (it cites the BEELINE study of Pratapa et al. 2020, Chen & Mar 2018, and Kernfeld et al. 2024).
  2. Correlation cannot establish causality. Perturbation data are needed to validate edges.
  3. Multimodal (RNA+ATAC) methods show limited robustness, are sensitive to parameters, and predict perturbation effects poorly (Badia-i-Mompel et al. 2025).
  4. Synthetic GRNs with known causal structure are proposed as future ground truth.
- **Implications for MeVD-GRN:**
  - Always report trivial baselines (degree, correlation, motif-only).
  - Validate on perturbation-derived networks; SC-MO-GRN-DB has KO/LOF networks, e.g. RN118 K562 KO.
  - Guard against leakage between motif/ATAC-derived features and motif-derived labels.
- Out of scope for numbers.

**(c) HyperCLSA: "Breast Cancer Subtyping with HyperCLSA: A Hypergraph Contrastive Learning Pipeline for Multi-Omics Data Integration"** (`papers/HyperCLSA.pdf`). Bhole G., HC P., JR M., Vinod P.K., Bhimalapuram P. (IIIT Hyderabad).
- **Code:** https://github.com/Gaurav2543/HyperCLSA.
- **Task.** PAM50 subtyping of **TCGA-BRCA, 850 bulk tumours** with mRNA, 450K methylation and miRNA (Broad GDAC Firehose; no accession given).
- **Method.** Hypergraph contrastive learning.
- **Results (Table 1, 5-fold CV).** Accuracy **0.901 ± 0.007**, F1-macro **0.866 ± 0.019**, F1-weighted **0.901 ± 0.007**, vs HyperTMO at 0.858 / 0.821 / 0.863.
- Out of scope.

**(d) "Network diffusion-based approach for survival prediction and identification of biomarkers using multi-omics data of papillary renal cell carcinoma"** (`papers/s00438-023-02022-4.pdf`). Shetty K.S., Jose A., Bani M., Vinod P.K. *Mol Genet Genomics* 298:871–882, 2023. doi:10.1007/s00438-023-02022-4.
- **Data.** Bulk TCGA-KIRP CNV (288 samples × 5913 genes), TCGA RNA-seq (GDC), microarray **GSE2748** for validation.
- **Method.** Network propagation (pyNBS-style) followed by consensus clustering.
- **Results (Table 1).** Best network GeneNet: log-rank p = 0.002, silhouette 0.95, CCC 0.98, PAC 15%.
- Out of scope.

**(e) "Deep Prediction of Human Essential Genes using Weighted Protein-Protein Interaction Networks"** (`papers/2024.10.09.616990v1.full.pdf`). Mehrpou S., Mansoori E.G. bioRxiv 2024, doi:10.1101/2024.10.09.616990.
- **Task.** Essential vs non-essential gene classification from node embeddings on two human PPI networks: Reactome FIs (12,277 genes, 230,243 interactions) and InWeb-IM (17,428 genes, 625,641 interactions). Labels from Guo et al. 2017.
- **Protocol.** 5-fold CV.
- **Results (Table 2 onward, ANN classifier).**
  - FIs: ACC 0.871 (under-sampling), 0.909 (GAN over-sampling).
  - InWeb: 0.868 and 0.923.
- **Code:** not in paper.
- Out of scope.

**(f) Aswin Jose MS thesis: "Network-Based Approaches for Cancer Subtype Identification and Prognosis"** (`ref-theses/Thesis_Aswin.pdf`). IIIT Hyderabad, advisor P.K. Vinod.
- **Ch. 2 (pyNBS on TCGA-KIRC/KIRP CNV and somatic mutations).** Uses PCNet and CRN networks. This is the thesis counterpart of paper (d).
- **Ch. 3 (DeepGraphMut).**
  - A graph autoencoder: TransformerConv + GraphNorm with focal loss, on an NCG cancer-gene PPI subnetwork (3217 nodes) with binary mutation features.
  - Covers **16 TCGA cancers, 7352 samples**.
  - Cox C-index is reported in Table 3.2, an image in the PDF. The text gives the range "0.5 to 0.8"; the earlier visual reading gave GBM 0.56 and LGG 0.80.
- **Data.** All bulk TCGA. No single-cell data and no GRN.
- **Ideas that could carry over:** focal loss for sparse positives; strong sensitivity to the choice of prior network.
- Out of scope.

---

## 3. Cross-paper dataset catalogue

The "Paired RNA+ATAC?" column is what matters for MeVD-GRN, which needs both modalities from the same cells. Cell counts are as printed in the papers; "n.s." means not stated.

### 3a. Datasets used by 2 or more sources (highest-value targets)

| Dataset | Accession / URL | Species | Cell type | Modalities | Paired RNA+ATAC? | #cells | Ground truth used | Used by |
|---|---|---|---|---|---|---|---|---|
| BEELINE hESC | GEO GSE75748 | human | embryonic stem cells | scRNA | **No (RNA only)** | 758 (759 in scRegNet / SC-MO-GRN-DB DS006) | STRING, non-spec ChIP, cell-type-specific ChIP | scUniGP, ATFGRN, CaHoT-GRN, scRegNet, Kommu thesis; in SC-MO-GRN-DB (DS006, RN101) |
| BEELINE hHEP | GEO GSE81252 | human | hepatocytes | scRNA | **No** | 425 (426) | same 3 GTs | same 5 + SC-MO-GRN-DB (DS007, RN102 HepG2) |
| BEELINE mDC | GEO GSE48968 | mouse | dendritic cells | scRNA | **No** | 383 (384) | same 3 GTs | same 5 + SC-MO-GRN-DB (DS002, RN110) |
| BEELINE mESC | GEO GSE98664 | mouse | embryonic stem cells | scRNA | **No** | 421 (422) | STRING, non-spec, specific, **LOF/GOF** | same 5 + SC-MO-GRN-DB (DS001, RN111, RN112) |
| BEELINE mHSC-E / -GM / -L | GEO GSE81682 | mouse | HSC erythroid / granulocyte-monocyte / lymphoid lineages | scRNA | **No** | 1071 / 889 / 847 (1072 / 890 / 848) | same 3 GTs | same 5 + SC-MO-GRN-DB (DS003–005, RN113) |
| 10x human PBMC (several different 10x demo files) | EpiAwareNet: 10x "10k Human PBMCs Multiome v1.0 Chromium Controller"; scTFBridge: scGLUE data page (10x `pbmc_granulocyte_sorted_10k` Multiome, identified online); scMultiomeGRN: 10x "10k PBMCs healthy donor" scATAC v1 + `pbmc_10k_v3` scRNA; SMOGT: "healthy PBMC" (Supp. Table S1, not available) | human | PBMC | RNA + ATAC | **Yes** for EpiAwareNet and scTFBridge; **No** for scMultiomeGRN (two separate datasets); SMOGT n.s. | 12,012 (EpiAwareNet); 7189 ATAC / 9129 RNA (scMultiomeGRN); n.s. (scTFBridge, SMOGT) | DoRothEA (EpiAwareNet); Cistrome ChIP, 17 TFs × 4 cell types + sc-eQTL + pcHi-C (scTFBridge); none quantitative (scMultiomeGRN); eQTL for expression tasks (SMOGT) | EpiAwareNet, scTFBridge, scMultiomeGRN, SMOGT (4) |
| K562 | SMOGT: Supp. Table S1 (n.s.); EpiXFormer: ENCODE bulk; SC-MO-GRN-DB DS025 (PMID 35501385) | human | CML cell line | SMOGT: paired sc; EpiXFormer: bulk; DS025: scRNA/scATAC | **Yes** (SMOGT, DS025) | DS025: 434 | ChIP-seq (TF-CRE) and Hi-C (CRE-CRE) in SMOGT; ENCODE ChIP peaks in EpiXFormer; RN117 ChIP / RN118 KO / RN119 dual evidence in SC-MO-GRN-DB | SMOGT, EpiXFormer, SC-MO-GRN-DB |
| GM12878 | SMOGT: Supp. Table S1; EpiXFormer: ENCODE bulk; SC-MO-GRN-DB DS024/DS028/DS033 (none RNA+ATAC paired) | human | lymphoblastoid cell line | SMOGT: paired sc | Yes in SMOGT; **no paired set in SC-MO-GRN-DB** | n.s. | ChIP-seq + Hi-C (SMOGT); ENCODE ChIP (EpiXFormer); RN201 ChIP-seq (SC-MO-GRN-DB) | SMOGT, EpiXFormer, SC-MO-GRN-DB |
| HCT116 | SMOGT: Supp. Table S1; EpiXFormer: ENCODE | human | colorectal cancer cell line | SMOGT: paired sc; EpiXFormer: bulk | Yes (SMOGT) | n.s. | ChIP-seq + Hi-C (SMOGT) | SMOGT, EpiXFormer |
| Human bone marrow | scTFBridge: BMMC GEO **GSE194122** (NeurIPS 2021); SMOGT: "BM" (Supp. Table S1) | human | BM mononuclear / HSC niche | paired 10x Multiome | **Yes** | n.s. in paper | none for GSE194122 in scTFBridge; ChIP-seq TF-CRE in SMOGT | scTFBridge, SMOGT (possibly different data) |

### 3b. Datasets used by one source

| Dataset | Accession / URL | Species | Cell type | Modalities | Paired RNA+ATAC? | #cells | Ground truth used | Used by |
|---|---|---|---|---|---|---|---|---|
| 10x E18 mouse brain 5k Multiome | 10x "Fresh Embryonic E18 Mouse Brain 5k" | mouse | embryonic brain | RNA + ATAC | **Yes** | 4,881 | TRRUST (5,064 edges, 819 regulators) | EpiAwareNet |
| Tomato root pN / mN | Patrick et al. 2026 bioRxiv (accession not in paper) | tomato | root, normal / low nitrogen | RNA + ATAC | **Yes** | 12,448 / 8,748 | "bulk-derived" prior (170,300 edges) | EpiAwareNet |
| Human fetal lung | GEO GSM4508936 | human | fetal lung, 9 cell types | scATAC + scRNA | **Unclear / likely No** (reviewer note: probably separate atlases) | 72,622 | motif-derived TF–TF network (FIMO 1e-6, TSS ± 2 kb) | scMultiomeGRN |
| Human hematopoiesis (Buenrostro 2018) | scRNA: Cell 2018 Data S2; scATAC: github.com/pinellolab/scATAC-benchmarking (Buenrostro_2018) | human | HSC / CMP / GMP | scATAC + separate RNA | **No** | n.s. | motif-derived TF–TF | scMultiomeGRN |
| AD brain snMultiome | Synapse syn52293417 | human | prefrontal cortex; microglia 8,600 | snRNA + snATAC | Yes (per paper) | 414,000 snRNA / 437,000 snATAC | none (qualitative) | scMultiomeGRN |
| RA synovium | GEO GSE243917 | human | synovium | paired multiome | **Yes** | n.s. | none (GWAS heritability only) | scTFBridge |
| Lung, ever/never smokers | printed "GSE2414468" (probably GSE241468) | human | normal lung | paired multiome | **Yes** | n.s. | none (qualitative) | scTFBridge |
| A549 | SMOGT Supp. Table S1 | human | lung cancer cell line | paired sc | Yes | n.s. | ChIP-seq + Hi-C | SMOGT |
| Human cerebral cortex; leukemia stem cells; K562 Perturb-seq | SMOGT Supp. Table S1 | human | – | paired sc (Perturb-seq: RNA) | Yes / Yes / No | n.s. | eQTL/GWAS/SE enrichment; perturbation DEGs | SMOGT |
| Bladder cancer spatial | GEO GSE171351 | human | bladder tumour | spatial transcriptomics | No | n.s. | hTFtarget (qualitative) | CaHoT-GRN |
| SC-MO-GRN-DB paired sets | https://scmogrndb.psu.edu: DS010 (PMID 34579774), DS011 (35977485), DS012 (29987051), DS014 (bioRxiv 2022.06.15.496239), DS025 (35501385), DS026 (38747151), DS027 (36318267) | mouse ESC (DS010–014); human K562, macrophage, MCF7 | – | scRNA + scATAC (+ Hi-C for DS010/011) | **Yes** | 9,021; 36,535; 930; 68,794; 434; 3,114; 2,995 | RN111/RN114–116 (mESC), RN117–119 (K562), RN204 (macrophage), RN205 (MCF7) | SC-MO-GRN-DB only; **no competitor uses these** (bonus results) |
| ENCODE 199 DBP × cell-line ChIP sets (+ GSE167979, GSE147716, GSE90895, GSE132440, GSE216090) | ENCODE (Supp. Table S6); GEO | human | 7 cell lines | bulk ChIP / ATAC / histone | n/a (bulk) | – | ChIP peaks | EpiXFormer |
| Out-of-scope bulk cohorts | TCGA-BRCA (HyperCLSA); TCGA-KIRP + GSE2748 (pRCC paper); 16 TCGA cancers (Aswin thesis); Reactome FIs / InWeb-IM PPI (essential genes); DECODE's GEO/MassIVE/PDC sets | human / mouse | – | bulk | No | – | – | out-of-scope sources |

---

## 4. Benchmark targets (numbers MeVD-GRN must beat)

**Read this first.**
- AUPRC is **not comparable across papers** unless the test prevalence (negatives per positive) is identical. On the same BEELINE Non-specific hESC TFs+500 cell, scUniGP reports 0.24, ATFGRN 0.10 and CaHoT-GRN 0.05 AUPRC, with AUROC 0.896 / 0.81 / 0.72. That gap reflects different split files and prevalence more than model quality.
- Always state the split source (GENELink `Train_Test_Split.py` or the ATFGRN repo split CSVs), the negative ratio and the seeds.

### 4a. Paired RNA+ATAC benchmarks (MeVD-GRN can use its full model)

| Paper | Dataset | Metric | Reported value | Protocol notes |
|---|---|---|---|---|
| EpiAwareNet | 10x PBMC 10k Multiome, DoRothEA prior (118 TFs, 5,092 edges) | AUROC / AUPRC | **0.8218 / 0.00244** (±5 kb); **0.8708 / 0.002933** (±100 kb); 3-seed 0.8232 ± 7.57e-4 / 0.002475 ± 4.55e-5 | random 80/20 edge split; scored over all T×G; 10:1 negatives in training; Tables 2, 6, 7 |
| EpiAwareNet | 10x E18 mouse brain 5k Multiome, TRRUST (819 TFs, 5,064 edges) | AUROC / AUPRC | **0.7263 / 0.00057** (±5 kb); 0.8002 AUROC (±100 kb); 0.001018 AUPRC (±50 kb) | same; Tables 2, 6 |
| EpiAwareNet | tomato root mN / pN | AUROC / AUPRC | 0.7415 / 0.08458 ; 0.7340 / 0.07451 (pN RNA-only ablation AUPRC 0.08173) | same; Table 2; no accession |
| scTFBridge (vs LINGER) | 10x PBMC Multiome, Cistrome ChIP, STAT1 in CD14 Mono | AUC / AUPR ratio | scTFBridge 0.693 / **2.221**; **LINGER 0.719** / 2.105 | legend values in Fig. 5b/c; unsupervised; top-k targets vs ChIP |
| scTFBridge (vs LINGER) | same, all 17 TFs × 4 cell types | AUC / AUPR ratio | figure only, approx. median AUC 0.63 (scTFBridge) vs **0.67 (LINGER)**; AUPR ratio approx. 1.6 vs **1.95** | Fig. 5d/e boxplots |
| scMultiomeGRN | human fetal lung (GSM4508936), 9 cell types, average | AUROC / AUPR | **0.924 / 0.790** (GENELink 0.900 / 0.758; DeepTFni 0.871 / 0.721) | motif-derived TF–TF labels; one positive fold + 1:1 random negatives; text (Results) |
| scMultiomeGRN | human hematopoiesis HSC/CMP/GMP | AUROC / AUPR | "≥ 2.5% / 4.6% above DeepTFni"; absolute values figure only (Supp. Fig. S1) | same |
| SMOGT | K562 / HCT116 / A549 / GM12878, CRE–CRE vs Hi-C | ROC-AUC / PR-AUC | 0.787 / 0.698 ; 0.861 / 0.219 ; 0.776 / 0.567 ; 0.681 / 0.615 | legend values in Fig. 2C; **peak–peak task** |
| SMOGT | BM, K562, HCT116, A549, GM12878, TF–CRE vs ChIP-seq | per-TF AUPR, max-F1 | not reported numerically (Fig. 2A boxplots) | **TF–peak task** |

### 4b. BEELINE scRNA-only benchmark (MeVD-GRN must run RNA-only or declare an external chromatin prior)

**Averages over the 22 dataset×GT combinations.**

| Paper | TFs+500 AUROC | TFs+1000 AUROC | TFs+500 AUPRC | TFs+1000 AUPRC | Notes |
|---|---|---|---|---|---|
| **scUniGP** | **0.911** (Table 1) | **0.921** (text / Fig. 14) | **0.522** (text) | 0.543 (figure only, Fig. 16) | 5 seeds (or 3), no std; best baseline scGREAT 0.892 / 0.908 / 0.463 / 0.538 |
| CaHoT-GRN | 0.846 | 0.842 | 0.420 | 0.409 | 10 runs, no std; Tables 1–2 |
| ATFGRN | not printed (heatmap only) | "exceeding 0.88" (text) | STRING-only average 0.53 (text) | – | Fig. 2 heatmaps; leakage concerns (§2.2) |
| scRegNet | – (specific GT only) | – | – | – | see per-cell table below |

**Best reported value per cell, TFs+500 (AUROC | AUPRC), with source.** S = scUniGP (AUROC from Table 1; AUPRC from Fig. 7 heatmap cells), A = ATFGRN (Fig. 2 heatmap cells), C = CaHoT-GRN (Tables 1–2), R = scRegNet (Table 2, Geneformer).

| GT | hESC | hHEP | mDC | mESC | mHSC-E | mHSC-GM | mHSC-L |
|---|---|---|---|---|---|---|---|
| STRING | 0.95 A (S 0.948) \| 0.53 A | 0.941 S \| 0.54 A | 0.956 S \| 0.66 S | 0.97 A \| 0.63 A | 0.942 S \| 0.54 A | 0.937 S \| 0.51 S | 0.882 S \| 0.36 A |
| Non-specific ChIP | 0.896 S \| 0.24 S | 0.906 S \| 0.23 S | 0.926 S \| 0.38 S | 0.928 S \| 0.32 S | 0.893 S \| 0.29 S | 0.882 S \| 0.36 S | 0.842 S \| 0.34 S |
| Cell-type-specific ChIP | 0.895 S \| 0.62 S/R | 0.910 S \| 0.85 S | 0.813 S \| 0.18 S/C | 0.941 S \| 0.88 S | 0.930 S \| 0.94 S/C/R | 0.937 S \| 0.94 S/C/R | 0.885 S \| 0.89 S |
| LOF/GOF | – | – | – | 0.92 A \| 0.68 A | – | – | – |

- The table only lists values from each paper's *proposed* method. On Specific mDC AUPRC, the baseline **GNNLink scores 0.25** (TFs+500, in both scUniGP Fig. 7 and scRegNet Table 2), which beats every proposed method. Use 0.25 as the bar there.
- The ATFGRN cells are heatmap readings (2 decimals); the scUniGP AUPRC cells are figure-printed values.

**Cell-type-specific GT, TFs+1000 (scRegNet Table 3, Geneformer; AUROC / AUPRC):** hESC 0.88/0.62, hHEP 0.90/0.84, mDC 0.84/0.17, mESC 0.93/0.87, mHSC-E 0.94/0.95, mHSC-GM 0.94/0.95, mHSC-L 0.88/0.87. scUniGP's TFs+1000 Specific AUPRC (figure only) is higher or equal on hESC 0.65, hHEP 0.85, mESC 0.89, mHSC-E 0.96 and mHSC-L 0.88, and lower on mDC 0.18 vs GNNLink 0.25.

### 4c. Numbers that are *not* targets
- **SC-MO-GRN-DB** publishes no method benchmark.
- **Out-of-scope papers:** HyperCLSA acc 0.901; pRCC log-rank p = 0.002; essential genes ACC 0.871–0.923; DECODE and EpiXFormer (AUROC 0.9905 on region classification). These cannot be compared with TF→target edge prediction.

---

## 5. Recommendations

Ranking criteria:
- **(a)** how many competitor papers report numbers on the dataset;
- **(b)** whether paired RNA+ATAC is available (MeVD-GRN needs both);
- **(c)** download size and effort.

Sizes marked "est." are our estimates and were not verified in the papers.

| Rank | Dataset / benchmark | (a) Competitors with numbers | (b) Paired RNA+ATAC | (c) Effort / size | What to run |
|---|---|---|---|---|---|
| 1 | **10x PBMC 10k Multiome** (EpiAwareNet file; also the scTFBridge `pbmc_granulocyte_sorted_10k` file) | 2 with numbers: EpiAwareNet (DoRothEA), scTFBridge / LINGER (Cistrome, 17 TFs). A different, unpaired PBMC is used by scMultiomeGRN, and SMOGT's is unknown. | **Yes** | Low: public 10x demo downloads (est. ≲ 1 GB each); DoRothEA and Cistrome are public | (i) EpiAwareNet protocol: 80/20 DoRothEA edge split, rank all T×G, report AUROC, AUPRC, AUPRC ratio, P@100. Targets: AUROC 0.8218 (±5 kb) and 0.8708 (±100 kb). (ii) scTFBridge protocol: hold out the 17 Cistrome TFs and score top-k targets per cell type. Target: beat LINGER (STAT1 CD14-Mono AUC 0.719; 17-TF median approx. 0.67, figure only). |
| 2 | **BEELINE 7 datasets** × TFs+500/1000 × 4 GTs | **5**: scUniGP, ATFGRN, CaHoT-GRN, scRegNet, Kommu thesis. Also in SC-MO-GRN-DB (DS001–007). | **No** (RNA only) | Very low: GENELink repo split files (est. tens of MB) | RNA-only MeVD-GRN (ATAC branch off) on the GENELink splits. Targets: average AUROC 0.911 / 0.921 and AUPRC 0.522 / 0.543 (scUniGP), plus the per-cell table in §4b. Also report a leakage-controlled variant: held-out TFs, reverse-edge dedup, degree-baseline sanity check. This gives the most competitor comparisons but cannot show the benefit of ATAC. Optionally add an external bulk chromatin prior (ENCODE H1 / mESC ATAC/DNase) as a clearly labelled extra variant. |
| 3 | **10x E18 mouse brain 5k Multiome** + TRRUST | 1 (EpiAwareNet) | **Yes** | Low: 10x demo download (est. < 1 GB) | EpiAwareNet protocol. Targets: AUROC 0.7263 (±5 kb) / 0.8002 (±100 kb); AUPRC 0.00057 / 0.001018 (±50 kb). |
| 4 | **Human fetal lung GSM4508936** (+ Buenrostro hematopoiesis) | 1 (scMultiomeGRN) | Probably **not** cell-paired (verify) | Medium: GEO sample, 72,622 cells; motif-label construction (FIMO) must be scripted | Reproduce the motif-derived TF–TF labels with 1:1 test negatives. Target: AUROC 0.924 / AUPR 0.790. Low scientific value, because the labels come from the same ATAC; report it alongside a ChIP-labelled benchmark. |
| 5 | **K562 / GM12878 / HCT116 / A549** (SMOGT) | 1 with numbers (SMOGT; peak-level) + EpiXFormer (bulk) | Yes (SMOGT's own data; accessions unavailable) + SC-MO-GRN-DB DS025 K562 (434 cells) | Medium–high: SMOGT accessions only in its supplement; Hi-C/ChIP processing needed | Only if MeVD-GRN gets a TF→peak head. Use SC-MO-GRN-DB RN117 (K562 ChIP) and RN201 (GM12878 ChIP) as gene-level labels instead. |
| 6 | **BMMC GSE194122** | 0 with GRN numbers (scTFBridge uses it for integration only) | **Yes** | Medium: NeurIPS 2021 multiome, large (est. several GB) | Optional extra paired set, labelled with ChIP-Atlas / SC-MO-GRN-DB HSC networks (RN113 is mouse; human networks would be needed). |
| 7 | **SC-MO-GRN-DB paired sets** (DS010/011/012/014 mESC; DS025 K562; DS026 macrophage; DS027 MCF7) | 0 competitors | **Yes** | Human zip 17.5 GB / mouse 7.2 GB (from `dataset_info.md`); could be fetched per dataset | The planned "bonus" results. mESC (RN111 ChIP, RN112/RN115 LOF/GOF, RN116 dual evidence) links to the BEELINE mESC cell type. |
| – | Tomato pN/mN; AD syn52293417; RA GSE243917; lung "GSE2414468" | ≤ 1, qualitative only (or no accession) | Yes | High, or accession missing | Skip. |

**Protocol hygiene for any comparison (drawn from the reviewer notes above and the NRG perspective):**
1. Select models on validation only. ATFGRN's code selects on test.
2. Remove reverse duplicates of training edges from STRING test sets.
3. Sample negatives only from TFs that have GT edges.
4. Report prevalence-aware metrics (AUPRC ratio = AUPRC / positive fraction, EPR).
5. Report trivial baselines: node degree, PCC, motif-only.
6. Add a held-out-TF split, since all BEELINE-family papers use transductive per-TF edge splits.
7. For ATAC/motif-labelled benchmarks (scMultiomeGRN), flag the circularity between features and labels.
8. Report mean ± std over ≥ 3 seeds. Most competitors report no std.
