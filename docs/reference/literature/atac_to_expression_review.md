# From open chromatin to gene expression across cell types: a short review

*Literature review, October 2026. Written for the mevd_grn project.*

## Abstract

A gene can only be read if its DNA is open, so open chromatin (measured by ATAC-seq or DNase-seq) should tell us something about expression. This review asks how much of the expression difference between cell types is explained by the difference in open chromatin, and how well models transfer to cell types they have never seen. We cover early regression models, gene activity scores, models that combine DNA sequence with accessibility, single-cell scATAC-to-scRNA translation, and studies comparing changes in accessibility with changes in expression. The short answer: accessibility predicts *which genes are high and which are low* well, and this transfers to new cell types. It predicts *how one gene changes between cell types* much less well. Simple baselines, such as the mean expression over training cell types, already score very high on the most common metric. Few papers report the harder metric, and fewer still the change between two cell types. That gap is where new work can help.

## 1. Why this question matters

Every cell has the same genome, but cells use different genes. Part of this choice is made by chromatin, the way DNA is packed around proteins. Where the packing is loose ("open" or "accessible"), transcription factors (proteins that switch genes on or off) can bind. ATAC-seq measures these open regions genome-wide. It is cheap and works in single cells.

If open chromatin could predict expression in a new cell line, we could learn the rules that link regulatory DNA to genes. In our setting, changing the cell line acts like a perturbation. We want f with f(ATAC of cell line c) close to the expression of c, for cell lines not used in training. There are two targets: (a) absolute expression in the new cell line, and (b) the log fold change between two cell lines, predicted from the change in ATAC.

## 2. Early models: open chromatin as a predictor of expression

The first large studies used DNase-seq, an older cousin of ATAC-seq. Thurman et al. (2012) mapped about 2.9 million open sites in 125 cell and tissue types and linked about 580,000 distal sites to promoters by checking which opened and closed together. Sheffield et al. (2013) correlated openness with expression across 112 samples from 72 cell types to predict target genes for more than 500,000 open sites.

Natarajan et al. (2012) asked directly whether open chromatin can predict cell-type-specific expression. Across 19 cell types, they classified genes as up-regulated, down-regulated or always on. Using transcription factor motifs inside open regions worked much better than using the promoter sequence alone. They also found that the CG content of the promoter matters a lot, and that any baseline should include it.

Dong et al. (2012) built regression models on ENCODE data from seven cell lines, using eleven histone marks plus DNase. Within one cell line, predicted and measured expression correlated at about r = 0.8 to 0.9 across genes. A model trained on K562 and applied to four other cell lines still reached r = 0.82 to 0.87. This shows a pattern we will meet again: ranking genes from high to low is easy, and the rule transfers. (This model used histone marks too, not accessibility alone.)

Later work focused on linking distal elements (enhancers) to their genes. The activity-by-contact (ABC) model (Fulco et al., 2019) scores an enhancer by its activity times its 3D contact with the promoter. It beat earlier methods on thousands of CRISPR tests and was later used to build enhancer–gene maps in 131 cell types (Nasser et al., 2021). The newer ENCODE-rE2G model (Gschwind et al., 2026) is trained on CRISPR data. Its main version needs only DNase-seq plus an average Hi-C map, so it can be run on new cell types. Among 523 ENCODE assays, DNase-seq was the best single proxy for enhancer activity, and DNase-seq models beat ATAC-seq models. PECA (Duren et al., 2017) instead modelled transcription factors, regulatory elements and target genes jointly from paired expression and accessibility data.

## 3. Gene activity scores: the simplest bridge

Single-cell ATAC tools turn peaks into a per-gene number called a "gene activity score". Signac counts reads over the gene body plus 2 kb upstream (Stuart et al., 2021; Signac docs). ArchR tested 56 rules and chose one that sums nearby signal with weights that fall with distance and stop at neighbouring genes; adding distance-weighted distal elements helped in every test (Granja et al., 2021). One caveat: that benchmark matched scATAC and scRNA from *different* cells by computational integration, not true paired data. MAESTRO uses a similar "regulatory potential" model (Wang et al., 2020). Cicero links peaks that are open together in the same cells and uses those links to score genes (Pliner et al., 2018).

How well do these scores track RNA? Not very well at the single-cell level. In BABEL's held-out GM12878 cells, ArchR and MAESTRO scores reached Pearson r = 0.15 and 0.19 against measured expression, while BABEL reached 0.59 (Wu et al., 2021). A 2026 benchmark on nine paired RNA + ATAC datasets found that the per-gene correlation between gene activity and expression across cells was below 0.2 for most datasets and methods. Yet the same scores still grouped similar cells together well (Naqing et al., 2026). So gene activity scores are good for labelling cell types, but they are weak predictors of how much of each gene a cell makes.

## 4. Adding DNA sequence: models that transfer to new cell types

Sequence models such as Enformer (Avsec et al., 2021) and Borzoi (Linder et al., 2025) predict expression tracks straight from DNA. They learn one output per training cell type, so they cannot predict a new cell type by themselves. Karollus et al. (2023) also showed that Enformer captures promoter effects well but mostly misses distal enhancers, even though it sees far enough along the DNA. ChromBPNet (Pampari et al., 2025, preprint) predicts accessibility, not expression, from sequence.

A newer group of models feeds accessibility in as the cell-type signal, next to the sequence. Because accessibility is an input, the same model can run on a new cell type.

- **EPCOT** (Zhang et al., 2023) takes sequence plus accessibility and predicts many outputs, including CAGE and RNA-seq. A model trained on four cell lines (GM12878, K562, HUVEC, IMR-90) predicted CAGE for H1 and HeLa, and the authors report cross-cell-type accuracy often as good as other models' cross-chromosome accuracy. We did not extract exact numbers.
- **EpiBERT** (Javed et al., 2025) first learns to fill in masked ATAC signal, then is fine-tuned to predict expression (RAMPAGE). Helpfully, it reports two metrics. "Across genes" asks whether high genes are predicted high. "Across cell types" asks whether each gene's change from one cell type to another is right. In training cell types (on held-out genes), EpiBERT scored mean r = 0.88 across genes and 0.72 across cell types. A simple linear baseline on gene accessibility scored 0.42 and 0.20, and fine-tuned Enformer scored 0.87 and 0.57. On 8 cell types never seen in training, EpiBERT scored 0.85 across genes but only 0.60 across cell types.
- **GET** (Fu et al., 2025) is pretrained on accessibility and motifs from 213 human cell types, then fine-tuned on paired scATAC and scRNA. On left-out fetal astrocytes, its correlation with measured expression was 0.94 (R² = 0.88), for genes with accessible promoters. Two simple baselines were also strong: the most similar training cell type reached r = 0.83, and the mean expression over training cell types reached r = 0.78. Without pretraining, GET fell to about r = 0.6. The authors also say GET predicts log fold change between two cell types, even when one is unseen, but the numbers are in an extended figure that we did not check.

Taken together, these papers show real transfer to new cell types. They also show that the usual "across genes" correlation is inflated by gene-level averages. Simply knowing which genes are usually high gives r near 0.8.

## 5. Single-cell translation: scATAC to scRNA

Paired assays such as SHARE-seq (Ma et al., 2020) and 10x Multiome measure ATAC and RNA in the same cell, and have led to many "translation" models. BABEL (Wu et al., 2021) uses two linked autoencoders. On a held-out cell cluster it reached r = 0.62 from ATAC to RNA, against 0.27 for a nearest-neighbour baseline, and 0.59 on never-seen GM12878 cells. This correlation pools every gene in every cell, so it rewards getting each gene's average level right. Polarbear (Zhang et al., 2022) adds unpaired single-assay data to make up for small paired datasets. CMOT (Alatkar & Wang, 2023) uses optimal transport (a way of matching two clouds of points) to map cells across modalities. scButterfly (Cao et al., 2024) was tested on held-out batches and held-out cell types, but scored mainly by whether predicted profiles still cluster by cell type, not by gene-level accuracy. scMoGNN (Wen et al., 2022), a graph neural network, was an official overall winner of modality prediction at NeurIPS 2021.

That challenge is the closest thing to a shared benchmark (Luecken et al., 2021; Lance et al., 2022). It used bone marrow from 12 donors processed at 4 sites, with held-out test data. The write-up is frank. Most entries scored about the same, and rankings were unstable: rerunning one winner four times moved it by up to 8 places. For ATAC to RNA, some entries clearly beat the baseline; the winner was a plain 4-hidden-layer network. Expression scales differed by site, which limited transfer to an unseen site, and RMSE (root mean squared error) did not reward getting the pattern right. The organisers suggested a correlation metric for future rounds.

## 6. Predicting changes: does a change in ATAC predict a change in RNA?

This is target (b), and the evidence is thinner and more mixed. Starks et al. (2019) compared promoter accessibility and expression in placenta and ten other cell types. Genes with high openness and high expression were mostly housekeeping genes. Many tissue-specific genes were highly expressed with only medium or low promoter openness. Some genes had open promoters but low expression, and these looked actively repressed (they carried a repressive histone mark). So an open promoter is often needed but is not enough.

For short-term stimuli, the link between fold changes can be very weak. After insulin treatment of fly S2 cells, the correlation between ATAC and RNA log fold changes across all peaks was only R = 0.05. It was a little higher for promoter peaks (Merrill et al., 2022). Between stable cell types the link is likely stronger, but we found few studies that report a single genome-wide number for it.

Timing is one reason. In SHARE-seq data, accessibility at dense regulatory regions (called DORCs) opened *before* the gene was expressed during lineage choice; the gap, called "chromatin potential", predicted future cell fate (Ma et al., 2020). MultiVelo (Li et al., 2023) fitted time-lag models and found genes whose chromatin closes before transcription stops, genes where it closes after, and cell states where chromatin and RNA are out of step. So open chromatin may mark genes that are "primed" (ready but not yet on), which limits what a single ATAC snapshot can say about RNA at the same moment.

Direct models of differential expression exist mostly for histone marks. DeepDiff, for example, predicts expression differences between pairs of Roadmap cell types (Sekhon et al., 2018). GET's log fold change result, noted above, is the main accessibility-based example we found.

## 7. Known limits, and how strong are simple baselines?

The evidence points to five limits.

1. **Most of the "explained variance" is between genes, not between cell types.** Across-gene correlations of 0.8 to 0.94 are common (Dong et al., 2012; Fu et al., 2025; Javed et al., 2025), but the training-mean baseline already gives r = 0.78 (GET). Across-cell-type correlations are much lower: 0.60 for EpiBERT on unseen cell types, and 0.20 for a linear accessibility baseline even on training cell types.
2. **Promoters are easier than enhancers.** Sequence models mostly capture promoter effects (Karollus et al., 2023). Adding distal peaks helps gene scores (Granja et al., 2021), but linking enhancers to genes still needs CRISPR-trained models (Gschwind et al., 2026).
3. **Open is not the same as active.** Open promoters can be repressed or primed (Starks et al., 2019; Ma et al., 2020; Li et al., 2023). Transcription factor levels and repressive marks are not visible in ATAC alone.
4. **Post-transcriptional effects.** RNA levels also depend on processing and mRNA stability. Dong et al. (2012) found CAGE (which measures transcription start) easier to predict than RNA-seq, and cytoplasmic and nuclear RNA differed in how predictable they were. Borzoi models some of these layers from sequence (Linder et al., 2025).
5. **Batch and scale shifts.** In NeurIPS 2021, site-specific expression scales hurt transfer to an unseen site (Lance et al., 2022).

On baselines, the record is mixed. Learned models clearly beat gene activity scores and kNN (Wu et al., 2021; Naqing et al., 2026). GET beat the nearest-cell-type and mean baselines (0.94 vs 0.83 and 0.78), a smaller gap than the headline suggests. In the NeurIPS challenge, most entries were close to the baseline on several subtasks. We found no paper that tests a nearest-cell-type baseline on log fold change.

## 8. Common datasets

- **ENCODE** cell lines (K562, GM12878, HepG2, and others): bulk DNase/ATAC with RNA-seq, CAGE and RAMPAGE; registry of candidate regulatory elements (ENCODE Project Consortium, 2020).
- **10x Multiome PBMC** (for example the 10k sorted PBMC set): paired scATAC and scRNA from one donor.
- **SHARE-seq** mouse skin and other tissues (Ma et al., 2020).
- **NeurIPS 2021 BMMC**: 12 donors across 4 sites, with a nested donor and site batch structure (Luecken et al., 2021).
- **TCGA ATAC**: 410 tumour samples across 23 cancer types, with matched TCGA expression (Corces et al., 2018).
- **CCLE**: has RNA-seq, methylation and other data for 1,072 cell lines (Ghandi et al., 2019), but the 2019 release lists no ATAC-seq. We did not find a CCLE-scale ATAC panel.

## 9. Open problems

Several gaps stand out. Most papers report only across-gene correlation; few report per-gene correlation across held-out cell types, and almost none test log fold change against a nearest-cell-type baseline. "Held-out cell type" means different things: in GET only the expression is hidden, and the cell type's ATAC may be seen in pretraining, while EpiBERT holds out both genes and cell types. Nobody has a clear estimate of the ceiling, meaning how much of the between-cell-type change any ATAC-based model could explain, given priming, repression and RNA decay. Single-cell translation methods are mostly scored by clustering, not gene-level accuracy. And there is no shared bulk benchmark of many cell lines with paired ATAC and RNA and a fixed leave-cell-line-out split.

## 10. What this means for our problem

For absolute expression in an unseen cell line, our across-gene correlation will look high almost whatever we do. So we should always report the training-mean and nearest-cell-line baselines next to the model, plus the per-gene correlation across held-out cell lines (EpiBERT's "across cell types" metric). For log fold change, the mean baseline predicts zero change, so this target is a cleaner test of whether ATAC carries cell-type-specific signal. But the evidence that ATAC changes track RNA changes is mixed, and promoter peaks carry most of the signal. A good model will likely need promoter plus distance-weighted distal accessibility, and ideally motif or sequence features, which drove the gains in GET and EpiBERT. Splits must hold out whole cell lines, with no test-line expression used anywhere in training or pretraining, and we should check for batch scale shifts between cell lines. Some errors will be expected and informative: genes with open but silent promoters may be primed or repressed, which ATAC alone cannot see.

## References

- Alatkar SA, Wang D (2023). CMOT: Cross-Modality Optimal Transport for multimodal inference. *Genome Biology* 24:163. https://doi.org/10.1186/s13059-023-02989-8
- Avsec Ž, et al. (2021). Effective gene expression prediction from sequence by integrating long-range interactions. *Nature Methods* 18:1196–1203. https://doi.org/10.1038/s41592-021-01252-x
- Cao Y, et al. (2024). scButterfly: a versatile single-cell cross-modality translation method via dual-aligned variational autoencoders. *Nature Communications* 15:2973. https://doi.org/10.1038/s41467-024-47418-x
- Corces MR, et al. (2018). The chromatin accessibility landscape of primary human cancers. *Science* 362. https://doi.org/10.1126/science.aav1898
- Dong X, et al. (2012). Modeling gene expression using chromatin features in various cellular contexts. *Genome Biology* 13:R53. https://doi.org/10.1186/gb-2012-13-9-r53
- Duren Z, et al. (2017). Modeling gene regulation from paired expression and chromatin accessibility data. *PNAS* 114:E4914–E4923. https://doi.org/10.1073/pnas.1704553114
- ENCODE Project Consortium, et al. (2020). Expanded encyclopaedias of DNA elements in the human and mouse genomes. *Nature* 583:699. https://doi.org/10.1038/s41586-020-2493-4
- Fu X, et al. (2025). A foundation model of transcription across human cell types (GET). *Nature*. https://doi.org/10.1038/s41586-024-08391-z (preprint: https://doi.org/10.1101/2023.09.24.559168)
- Fulco CP, et al. (2019). Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. *Nature Genetics* 51:1664–1669. https://doi.org/10.1038/s41588-019-0538-0
- Ghandi M, et al. (2019). Next-generation characterization of the Cancer Cell Line Encyclopedia. *Nature* 569:503. https://doi.org/10.1038/s41586-019-1186-3
- Granja JM, et al. (2021). ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. *Nature Genetics* 53:403–411. https://doi.org/10.1038/s41588-021-00790-6
- Gschwind AR, et al. (2026). An encyclopedia of human enhancer–gene regulatory interactions (ENCODE-rE2G). *Nature*. https://doi.org/10.1038/s41586-026-10781-4 (preprint: https://doi.org/10.1101/2023.11.09.563812)
- Javed N, et al. (2025). A multi-modal transformer for cell type-agnostic regulatory predictions (EpiBERT). *Cell Genomics* 5:100762. https://doi.org/10.1016/j.xgen.2025.100762
- Karollus A, Mauermeier T, Gagneur J (2023). Current sequence-based models capture gene expression determinants in promoters but mostly ignore distal enhancers. *Genome Biology* 24:56. https://doi.org/10.1186/s13059-023-02899-9
- Lance C, et al. (2022). Multimodal single cell data integration challenge: results and lessons learned. *Proceedings of the NeurIPS 2021 Competitions and Demonstrations Track*, PMLR vol. 176. https://proceedings.mlr.press/v176/lance22a.html (preprint: https://doi.org/10.1101/2022.04.11.487796)
- Li C, Virgilio MC, Collins KL, Welch JD (2023). Multi-omic single-cell velocity models epigenome–transcriptome interactions and improves cell fate prediction (MultiVelo). *Nature Biotechnology* 41:387–398. https://doi.org/10.1038/s41587-022-01476-y
- Linder J, et al. (2025). Predicting RNA-seq coverage from DNA sequence as a unifying model of gene regulation (Borzoi). *Nature Genetics* 57:949–961. https://doi.org/10.1038/s41588-024-02053-6
- Luecken MD, et al. (2021). A sandbox for prediction and integration of DNA, RNA, and proteins in single cells. *NeurIPS Datasets and Benchmarks Track*. https://datasets-benchmarks-proceedings.neurips.cc/paper/2021/hash/158f3069a435b314a80bdcb024f8e422-Abstract-round2.html
- Ma S, et al. (2020). Chromatin potential identified by shared single-cell profiling of RNA and chromatin (SHARE-seq). *Cell* 183:1103–1116. https://doi.org/10.1016/j.cell.2020.09.056
- Merrill CB, et al. (2022). Harnessing changes in open chromatin determined by ATAC-seq to generate insulin-responsive reporter constructs. *BMC Genomics*. https://doi.org/10.1186/s12864-022-08637-y
- Naqing F, Yuan Q, Duren Z (2026). Benchmarking component choices for unpaired single cell RNA and epigenomic integration. *Genome Biology* 27. https://doi.org/10.1186/s13059-026-04071-5
- Nasser J, et al. (2021). Genome-wide enhancer maps link risk variants to disease genes. *Nature*. https://doi.org/10.1038/s41586-021-03446-x
- Natarajan A, Yardimci GG, Sheffield NC, Crawford GE, Ohler U (2012). Predicting cell-type-specific gene expression from regions of open chromatin. *Genome Research* 22:1711–1722. https://doi.org/10.1101/gr.135129.111
- Pampari A, et al. (2025). ChromBPNet: bias factorized, base-resolution deep learning models of chromatin accessibility reveal cis-regulatory sequence syntax, transcription factor footprints and regulatory variants. *bioRxiv* (preprint). https://doi.org/10.1101/2024.12.25.630221
- Pliner HA, et al. (2018). Cicero predicts cis-regulatory DNA interactions from single-cell chromatin accessibility data. *Molecular Cell* 71:858–871. https://doi.org/10.1016/j.molcel.2018.06.044
- Sekhon A, Singh R, Qi Y (2018). DeepDiff: DEEP-learning for predicting DIFFerential gene expression from histone modifications. *Bioinformatics*. https://doi.org/10.1093/bioinformatics/bty612
- Sheffield NC, et al. (2013). Patterns of regulatory activity across diverse human cell types predict tissue identity, transcription factor binding, and long-range interactions. *Genome Research* 23:777. https://doi.org/10.1101/gr.152140.112
- Starks RR, Biswas A, Jain A, Tuteja G (2019). Combined analysis of dissimilar promoter accessibility and gene expression profiles identifies tissue-specific genes and actively repressed networks. *Epigenetics & Chromatin* 12:16. https://doi.org/10.1186/s13072-019-0260-2
- Stuart T, Srivastava A, Madad S, Lareau CA, Satija R (2021). Single-cell chromatin state analysis with Signac. *Nature Methods* 18:1333–1341. https://doi.org/10.1038/s41592-021-01282-5 (GeneActivity docs: https://stuartlab.org/signac/reference/geneactivity)
- Thurman RE, et al. (2012). The accessible chromatin landscape of the human genome. *Nature* 489:75. https://doi.org/10.1038/nature11232
- Wang C, et al. (2020). Integrative analyses of single-cell transcriptome and regulome using MAESTRO. *Genome Biology* 21:198. https://doi.org/10.1186/s13059-020-02116-x
- Wen H, et al. (2022). Graph neural networks for multimodal single-cell data integration (scMoGNN). *KDD 2022, Applied Data Science Track*. https://arxiv.org/abs/2203.01884
- Wu KE, Yost KE, Chang HY, Zou J (2021). BABEL enables cross-modality translation between multiomic profiles at single-cell resolution. *PNAS* 118:e2023070118. https://doi.org/10.1073/pnas.2023070118
- Zhang R, Meng-Papaxanthos L, Vert J-P, Noble WS (2022). Multimodal single-cell translation and alignment with semi-supervised learning (Polarbear). *Journal of Computational Biology* 29:1198–1212. https://doi.org/10.1089/cmb.2022.0264
- Zhang Z, Feng F, Qiu Y, Liu J (2023). A generalizable framework to comprehensively predict epigenome, chromatin organization, and transcriptome (EPCOT). *Nucleic Acids Research* 51:5931–5947. https://doi.org/10.1093/nar/gkad436
- 10x Genomics. PBMC from a healthy donor, granulocytes removed through cell sorting (10k), Single Cell Multiome ATAC + Gene Expression dataset. https://www.10xgenomics.com/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0
