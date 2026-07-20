iScience
Article
SC-MO-GRN-DB: A comprehensive repository for
single-cell multiomic gene regulatory networks
Graphical abstract Authors
HannahValensi,KaramveerKaramveer,
EricMoeller,SamiErenOzdogan,
RachaelM.Edwards,YasinUzun
Correspondence
yuzun@pennstatehealth.psu.edu
In brief
Medicine;Biologicalsciences;Genetics
Highlights
•
Asingle-cellmultiomicsgeneregulatorynetworkdatabase
foranalysisandunbiasedbenchmark
•
Evidence-basedreferencenetworksbuiltwithTFlocalization
andperturbationdata
•
Standardizedsingle-cellmultiomicdatasetspairedwith
referencenetworks
•
Intuitiveinterfacesupportingeasyaccessandexplorationof
data
Valensietal.,2026,iScience29,115323
April17,2026©2026TheAuthor(s).PublishedbyElsevierInc.
ll
https://doi.org/10.1016/j.isci.2026.115323

iScience
ll
OPEN ACCESS
Article
SC-MO-GRN-DB: A comprehensive repository
for single-cell multiomic gene regulatory networks
Hannah Valensi,1Karamveer Karamveer,1Eric Moeller,1Sami Eren Ozdogan,2Rachael M. Edwards,3and Yasin Uzun1,4,5,*
1Department of Pediatrics, Pennsylvania State University College of Medicine, Hershey, PA, USA
2Department of Computer Science and Engineering, Pennsylvania State University, State College, PA, USA
3Department of Neuroscience and Experimental Therapeutics, Pennsylvania State University College of Medicine, Hershey, PA, USA
4Department of Molecular and Precision Medicine, Pennsylvania State University College of Medicine, Hershey, PA, USA
5Lead contact
*Correspondence: yuzun@pennstatehealth.psu.edu
https://doi.org/10.1016/j.isci.2026.115323
SUMMARY
Gene regulatory networks (GRNs) are essential models for understanding gene expression, cell differentia-
tion, and cellular function. Existing GRN resources primarily focus on transcriptomic data and often overlook
the epigenetic mechanisms of gene regulation, limiting evaluation and improvement of GRN inference
methods. To address this challenge, we developed SC-MO-GRN-DB, a publicly accessible database of exper-
imentally validated GRNs alongside tissue-matched single-cell multiomic datasets across human and mouse
tissues. This repository includes ground-truth GRNs comprising over 22 million regulatory edges, curated
from high-confidence experimental datasets. In addition, it hosts multiomic single-cell datasets totaling
more than two million cells across six molecular modalities: single-cell RNA sequencing (scRNA-seq), chro-
matin accessibility (scATAC-seq), chromatin immunoprecipitation (scChIP-seq), DNA methylation (scDNA-
Met), chromatin conformation (scHi-C), and gene perturbation screens (scCRISPR-seq). By offering standard-
ized input datasets paired with experimentally supported ground-truth networks, SC-MO-GRN-DB provides a
platform for the development, benchmarking, and validation of GRNs with single-cell multiomic data.
INTRODUCTION A variety of GRN-focused databases have since been devel-
oped to aggregate TF-gene interactions, providing benchmarks
The regulation of gene expression is a fundamental process that and reference points for the field.14,16–19While these resources
determines how cells acquire their identity, respond to environ- have been instrumental in advancing systems-level studies,
mental cues, and maintain their functions. Central to this regula- most rely on bulk data or focus narrowly on transcriptomic-
tion are transcription factors (TFs), which recognize and bind to only data. Although bulk profiling captures genome-wide activity
specific DNA sequences to activate or repress their target genes patterns, it reflects population averages across millions of cells.
(TGs).1Rather than acting in isolation, TFs operate within inter- This masks cell type-specific regulation and obscures heteroge-
connected gene regulatory networks (GRNs) that capture the de- neity in the cell population.20,21Consequently, bulk-based net-
pendencies among regulators and their downstream targets.2 works are limited in representing the full spectrum of regulatory
These networks provide a framework for understanding how ge- mechanisms that shape gene expression, including epigenetic
netic programs are established,3 how cellular states transition modifications, chromatin accessibility, higher-order chromatin
during development or disease,2,4and how dysregulation of TF organization, and direct perturbation effects.
activity contributes to pathological processes.1 The introduction of single-cell technologies has provided
Over the past two decades, a diverse range of strategies for unprecedented resolution for characterizing TF activity and
inferring the structure and function of GRNs have been tested.5–10 network structure at the level of individual cells.22Both single-
Prior to the high-throughput sequencing (HTS) era, attempts to cell transcriptomics and epigenomics have been instrumental
reconstruct GRNs relied on experimentally validated TF-TG inter- in enabling the construction of GRNs that better capture
actions compiled from the literature, ChIP-chip, or small-scale dynamic regulatory programs and the molecular diversity of
molecular assays with limited size and scope.11–14The advent tissues.23–25 These advances have driven the development of
of HTS technologies enabled large-scale inference of GRNs us- computational approaches that leverage single-cell omics to
ing data types such as RNA-seq, ChIP-seq, and related assays, infer more precise and context-dependent regulatory net-
providing genome-wide views of transcriptional regulation.15 works.26 Although recent efforts have incorporated single-cell
Single-cell RNA-seq and multiomic assays further facilitated transcriptomic data into GRN resources, a collection of
the development of algorithms and resources that capture regu- multiple single-cell modalities with experimentally supported
latory programs with elevated power and resolution.7–10 ground-truth networks is lacking. This gap constrains both our
iScience 29, 115323, April 17, 2026 © 2026 The Author(s). Published by Elsevier Inc. 1
This is an open access article under the CC BY-NC-ND license (http://creativecommons.org/licenses/by-nc-nd/4.0/).

iScience
ll
OPEN ACCESS Article
Figure 1. Reference network curation pro-
A
cess and summary statistics
(A) TF localization network construction flow chart.
(B) TF perturbation network construction flow
chart.
B
(C) Dual evidence network construction flow chart.
(D) Network overall statistics table.
(E) Number of networks per organism pie chart.
(F) Number of networks per cell type pie chart.
C
methods. By bridging transcriptomic and
epigenomic layers of regulation, SC-MO-
GRN-DB establishes a unique foundation
for studying GRNs in diverse cellular con-
texts and for accelerating discoveries in
D E F systems biology, disease mechanisms,
and therapeutic target identification.
RESULTS
Current GRN research is built upon two
fundamental components: ground truth
(reference) networks and single-cell da-
tasets. Ground truth networks curated
through various types of genetic data
serve as a reference for evaluating in-
mechanistic understanding of transcriptional regulation and the ferred GRNs. Depending on the type of evidence, these net-
evaluation of new computational methods for GRN inference. works can either be tissue- or organism-specific. Single-cell
Evaluating and comparing GRN inference methods remains a sequencing datasets are used as input data for methods that
major challenge. Each computational study typically constructs computationally infer GRNs. We curated both types of data to
its own benchmarking datasets from disparate sources, provide a foundation for GRN research, especially in the context
requiring extensive preprocessing to match biological context of single-cell genomics.
or experimental modality. As each individual tool employs a
different set of criteria for selecting or filtering reference net- Reference networks
works, results are often incomparable across studies. This frag- To establish reliable benchmarks for GRN inference, we curated
mentation slows down methodological progress, as each group high-confidence, experimentally validated TF-TG interaction
must rebuild their evaluation framework from scratch rather networks. For curating tissue and cell type-specific networks,
than relying on standardized benchmarking. Although previous we employed complementary strategies. One approach that
efforts8,27have provided useful reference datasets, they require we utilized is to employ TF localization data from ChIP-Seq
substantial manual curation and are limited in their coverage of and ChIP-chip experiments to identify physical binding events
single-cell and multiomic modalities. The field has therefore (Figure 1A). Widely used for evaluating inferred networks,28–30
lacked a unified, comprehensive, and ready-to-use resource this approach assumes that TFs bind to the DNA regions that
that pairs high-confidence ground-truth networks with the corre- are in proximity to their TGs to initiate transcription.
sponding single-cell data necessary for rigorous benchmarking. Although ChIP-Seq and ChIP-chip data provide TF localization
To overcome these challenges, we developed SC-MO-GRN-DB information, it does not guarantee a regulatory interaction, lead-
(https://scmogrndb.psu.edu), a comprehensive database that ing to potential false positives. As an alternative, TF perturbation
unifies experimentally validated reference GRNs with harmonized assays coupled with gene expression profiling provide functional
single-cell multiomic datasets. This resource catalogs more than evidence of gene regulation. These assays reveal how altering TF
22 million high-confidence regulatory edges derived from experi- activity via genetic manipulation influences downstream gene
mental evidence and contains over 2 million cells across six molec- expression, and therefore, provide causal evidence of regula-
ular single-cell modalities: scRNA-seq, scATAC-seq, scChIP-seq, tion. More specifically, perturbation-based data from genetic
scDNA methylation, scHi-C, and scCRISPR-based perturbation knockdowns, knockouts, and overexpression assays capture
screens. All datasets have been systematically processed and the functional consequences of changes in TF activity on down-
standardized to facilitate cross-study and cross-modality compar- stream TG expression. For the genetic perturbation experiment
isons. In addition to ground truth GRNs, the database provides a of a specific TF, the genes showing significant differential
curated collection of single-cell multiomic data that supports expression in the readout as a response to the perturbation
reproducible analyses and benchmarking of computational event are considered as its targets. We also curated reference
2 iScience 29, 115323, April 17, 2026

iScience
ll
Article OPEN ACCESS
networks using TF perturbation data obtained from available re- cord (VSC), and cell lines BJ, GM12878, HepG2, K562, H1,
sources (Figure 1B). Although generating such data at scale re- and MCF7 (Figure 1F).
mains technically challenging, CRISPR-based technologies The sizes of curated reference networks vary markedly, re-
streamline the process, expanding the feasibility of systemati- flecting differences in experimental scope, evidence depth,
cally capturing TF-gene relationships across many conditions. and biological context (Figures 2and S2). The literature-based
To reduce potential false positives, we also built dual-evi- networks include the lowest numbers of TFs as regulators,
dence networks by intersecting TF localization- and perturba- limiting their usability for comprehensive studies whereas gen-
tion-based networks for the same cell type (Figure 1C). These eral, non-specific networks have greatest TFs, on the opposite
networks require both types of evidence (localization and end of the spectrum (Figure 2A). Tissue-specific networks that
perturbation) to hold at the same time for each individual TF- are based on experimental evidence span the space between
TG interaction, increasing the confidence. However, this these two groups, where the number of TFs tends to be higher
approach decreases network size, as for a TF to be included for heavily studied cell types and cell lines such as ESC, HSC,
for a specific cell type, experimental data need to be available and K562.
for both assays. It should be noted that this reduced coverage, In terms of overall regulon size (overall TG set size; Table S1),
while enhancing confidence, results in sparser TF regulons, literature-based networks are also limited (Figure 2B). Experi-
which may influence benchmarking and downstream analyses mentally derived networks are larger in size when compared to
(Figure S1). This design reflects the broader tradeoff between others (Figure 2B). As these networks are based on genome-
confidence and comprehensiveness. Dual-evidence networks wide experiments, they cover a large set of targets and have
prioritize high-confidence interactions, whereas single-evidence larger sizes than PPI and text mining-based networks. The total
networks offer more extensive coverage, and the choice be- number of edges (TF-to-TG interactions) increase linearly with
tween them may depend on the goals of a given analysis. target set sizes (Figure 2C). A similar trend is seen for the regu-
In addition to experimentally derived networks, SC-MO-GRN- lons of individual TFs. Literature-based and text mining-based
DB includes smaller sets of literature-curated, text-mined, and networks are limited in size compared to experimentally derived
protein-protein interaction (PPI)-based networks to broaden bio- reference networks with large regulons (Figure 2D).
logical and contextual coverage. Literature-curated networks SC-MO-GRN-DB enables researchers to explore reference
are assembled manually by expert curators who extract TF-TG networks directly through an intuitive web interface. The
relationships from published studies. These networks often cap- browsing functionality allows researchers to navigate by organ-
ture high-quality interactions with direct experimental evidence ism, tissue, and evidence type, presenting summaries of the size
but can be limited by human bias and incomplete coverage, as of the networks via TF, TG, and edge counts (Figure 3A). This hi-
not all regulatory relationships are reported or curated consis- erarchical view facilitates comparison across datasets and quick
tently. PPI-based networks leverage protein-protein interaction assessment of network scope and quality. Complementing this,
data31,32 to infer potential co-regulatory relationships among the search function enables targeted exploration of specific TFs,
TFs and other transcriptional regulators. Although PPI-based TGs, or TF-TG pairs across all networks (Figure 3B). Search re-
networks often lack tissue specificity, they expand SC-MO- sults display the supporting evidence types dynamically, allow-
GRN-DB’s coverage, allowing both detailed and exploratory an- ing the investigators to trace how each regulatory relationship
alyses across diverse biological conditions. is supported by localization, perturbation, dual evidence, or liter-
Automated text mining has also been used as an approach to ature-based sources. These interactive features streamline
derive gene networks from the published literature.33 These discovery and cross-validation of transcriptional regulatory infor-
methods extract TF-to-TG associations by processing vast mation within and across biological systems.
amounts of scientific literature and can be useful for hypothesis
generation and for supplementing areas with limited experi- Single-cell datasets
mental data. Although text mining-based networks can include As the complementary component to reference networks, SC-
indirect associations, they can be useful for some cell types MO-GRN-DB hosts curated single-cell datasets. Each single-
where tissue-specific experimental data are not available and cell dataset is paired with a corresponding tissue-type-specific
are used for studying genetic interactions and evaluating infer- reference network, enabling direct benchmarking of inferred
ence models.8 Therefore, they are included in SC-MO-GRN- GRNs against experimentally supported ground-truth interac-
DB, in addition to the cell type-specific reference networks sup- tions for the same biological context. The quality and composi-
ported with experimental data.31 tion of the single-cell data have a strong effect on the accuracy
Curated with the described approaches, SC-MO-GRN-DB of inferred GRNs. Large-scale single-cell datasets provide
hosts 30 curated reference networks spanning 16 tissues across elevated power to detect rare cell states and subtle regulatory
human and mouse in total. These networks contain 1,616 TFs, relationships; however, excessive cell numbers without appro-
63,512 TGs, and 22,389,021 regulatory edges (Figure 1D). Net- priate curation can introduce technical noise and increase
works are distributed across 18 human and 12 mouse tissues computational burden. Datasets that are smaller in size can fail
(1E) encompassing varied biological contexts including embry- to detect rare cell states and miss subtle regulatory relationships,
onic stem cells, hematopoietic stem cells (HSCs), cortical area but can still capture core regulatory programs and high-fre-
development (CAD), dendritic cells (DCs), gonadal sex determi- quency interactions. The size of the datasets in SC-MO-GRN-
nation (GSD), induced pluripotent stem cells (IPSCs), blood mac- DB is curated to ensure both sufficient cellular coverage and bio-
rophages (Macr.), pluripotent stem cells (PSCs), ventral spinal logical relevance.
iScience 29, 115323, April 17, 2026 3

1,000
100
10
1
In this context, we curated single-cell multiomic datasets span- using sequencing), scChIP-seq (single-cell chromatin immunopre-
ning six modalities: scRNA-seq (single-cell RNA sequencing), scA- cipitation sequencing), scDNA methylation (single-cell DNA
TAC-seq (single-cell assay for transposase-accessible chromatin methylation profiling), scHi-C (single-cell high-throughput
sFT
#
DAC-701 CSV-401 CSH-301 DSG-501 CSE-611 JB-002 rcaM-402 265K-911 1H-202 CSPI-801 CD-011 CSP-601 CSE-411 CSE-211 CSE-511 PEH-201 265K-811 1H-302 CSE-101 CSH-311 265K-711 MG-102 CSE-111 FCM-502 SN-300 SN-400 SN-200 SN-600 SN-100 SN-500
10,000
1,000
100
10
1
sGT
#
DAC-701 CSV-401 CSH-301 DSG-501 CSPI-801 CSP-601 SN-400 SN-300 CSE-611 SN-200 SN-100 265K-911 1H-202 CD-011 265K-811 CSE-211 PEH-201 SN-600 CSE-101 CSE-511 CSH-311 CSE-411 SN-500 CSE-111 MG-102 JB-002 rcaM-402 1H-302 265K-711 FCM-502
1,000,000
100,000
10,000
1,000
100
10
1
segde
#
DAC-701 CSV-401 CSH-301 DSG-501 CSPI-801 CSP-601 SN-400 CSE-611 SN-300 265K-911 CD-011 1H-202 SN-200 SN-100 SN-600 CSE-211 CSE-511 CSE-411 rcaM-402 JB-002 265K-811 PEH-201 SN-500 CSE-101 1H-302 CSE-111 CSH-311 MG-102 265K-711 FCM-502
Evidence: Localization Perturbation Dual evidence
Literature Text mining Protein−protein interaction
RN ID:
CSV-401 CSH-301 DAC-701 SN-300 SN-400 CSP-601 CSPI-801 DSG-501 SN-500 SN-200 SN-600 SN-100 CSE-611 CD-011 265K-911 1H-202 CSE-101 CSE-211 CSE-511 CSE-411 265K-811 CSE-111 PEH-201 rcaM-402 FCM-502 1H-302 MG-102 265K-711 CSH-311 JB-002
10,000
1,000
100
10
1
RN ID:
sezis
nolugeR
iScience
ll
OPEN ACCESS Article
A
B
C
D
Figure 2. Reference network sizes
Bars and boxes are colored with respect to the evidence type. The y axes are displayed in logarithmic scale for visibility purposes. Solid border lines indicate non-
specific networks. Networks derived from TF binding localization data and protein-protein interaction evidence tend to exhibit higher TF, TG, and edge counts,
reflecting broader regulatory coverage, whereas perturbation-based and literature-curated networks are generally sparser and more selective. Labels for mouse
datasets are colored with blue and human datasets are black. Middle lines show median values for each group. Box edges represent 25% and 75% quartiles.
Whiskers extend 1.5 times the inter-quartile range from the individual box edges within the limits of data points.
(A) Total number of TFs for each reference network.
(B) Overall TG set sizes for each reference network (each TG counted only once, regardless of the regulator).
(C) Total number of TF-to-TG edges for each reference network.
(D) Distribution of individual regulon sizes (TG counts) for individual TFs for each reference network.
4 iScience 29, 115323, April 17, 2026

| iScience |     |     |     |     |     |     |     | ll          |
| -------- | --- | --- | --- | --- | --- | --- | --- | ----------- |
| Article  |     |     |     |     |     |     |     | OPEN ACCESS |
A
|        |           |          | EvidenceType: | Organism:   | Cell Type:    |        |               |          |
| ------ | --------- | -------- | ------------- | ----------- | ------------- | ------ | ------------- | -------- |
|        |           | All      |               | All         | All           |        |               |          |
|        | Reference |          |               |             |               | Number |               |          |
|        |           | Evidence |               | Cell        | Number Number |        |               |          |
| Select | Network   |          | Organism      |             |               |        | of PMID       | Download |
|        | ID        | Type     |               | Type        | of TFs of TGs | Edges  |               |          |
|        | RN001     | PPI      | Human         | Nonspecific | 1489          | 8712   | 9914931907445 | Download |
|        | RN002     | PPI      | Mouse         | Nonspecific | 1350          | 7636   | 7856731907445 | Download |
RN003 TextMining Human Nonspecific 795 2492 938429087512 Download
B
|     | Transcription Factor (TF): |     |     | Target Gene (TG): |     |     |     |     |
| --- | -------------------------- | --- | --- | ----------------- | --- | --- | --- | --- |
Search
|        | RUNX1 |          |          | e.g., BRCA1 |       |        |          |          |
| ------ | ----- | -------- | -------- | ----------- | ----- | ------ | -------- | -------- |
|        |       | Evidence |          |             | TF    | TGs    |          |          |
| Select | RN ID |          | Organism | Tissue      |       |        | PMID     | Download |
|        |       | Type     |          |             | Name  | per TF |          |          |
|        | RN001 | PPI      | Human    | Nonspecific | RUNX1 | 161    | 31907445 | Download |
|        | RN002 | PPI      | Mouse    | Nonspecific | RUNX1 | 131    | 31907445 | Download |
RN003 TextMining Human Nonspecific RUNX1 40 29087512 Download
RN004 TextMining Mouse Nonspecific RUNX1 17 29087512 Download
Figure 3. Reference network investigation functionalities
(A) Reference network browsing functionality.
(B) Reference network search functionality.
chromosome conformation capture), and scCRISPR-seq (single-  depth, and biological system complexity (Figure 4E). Dataset
cell CRISPR-based perturbation sequencing) (Figure 4A). Partic- sizes range from as few as 63 cells in smaller targeted studies
ular emphasis was placed on joint assays that simultaneously  to nearly two million cells in larger-scale studies. This extensive
measure gene expression and chromatin accessibility, especially  range ensures that the repository includes both high-resolution,
scRNA-seq with scATAC-seq, providing complementary layers  small-scale datasets ideal for mechanistic validation and large-
of information critical for accurate GRN reconstruction. For each  scale datasets suited for benchmarking scalability and statistical
dataset, only wild-type or control cells were retained to ensure  robustness. Variation in dataset size also enables researchers to
that the data reflected baseline cellular states. When datasets con- assess how GRN inference methods perform across datasets
tained multiple annotated cell types, cells were separated into  with differing sparsity, depth, and biological diversity.
individual subsets so that each dataset corresponded to a single,  Similar to the case of reference networks, SC-MO-GRN-DB
untreated, unmodified cell type. This curation strategy ensures that  provides multiple ways to explore these curated datasets. The
single-cell data are optimally matched to the corresponding refer- individual datasets can be browsed by molecular modality, or-
ence networks, which are likewise cell type-specific and derived  ganism, and cell type, and or multiple datasets can be selected
from unperturbed contexts. as a batch through the interface (Figure 4F). Each dataset is sum-
Our repository contains 31 harmonized single-cell datasets  marized with its total number of cells, cell type, and source study.
covering nine tissues for human and mouse, totaling 2,322,285  This functionality enables researchers to quickly identify data-
cells,  and  spans  six  modalities  with  seven  joint  assays  sets relevant to their study design, whether they seek large-scale
(Figure 4B). Datasets cover both human and mouse tissues  multiomic assays or more narrowly defined, cell type-specific
| (Figure 4C) and cell type representation includes ESC (n = 8),  |     |     |     | datasets. |     |     |     |     |
| --------------------------------------------------------------- | --- | --- | --- | --------- | --- | --- | --- | --- |
K562 (n = 8), HSC (n = 4), GM12878 (n = 3), BJ (n = 2), MCF7
| (n = 2), and single datasets from DC, H1, HepG2, and macro- |     |     |     | DISCUSSION |     |     |     |     |
| ----------------------------------------------------------- | --- | --- | --- | ---------- | --- | --- | --- | --- |
phages (Figure 4D).
The number of cells per dataset in SC-MO-GRN-DB varies  Benchmarking GRN inference methods is challenging because
widely, reflecting differences in experimental design, sequencing  computational methods often train or validate models using
iScience 29, 115323, April 17, 2026  5

A
Single-cell dataset Wild-type/ Cell type filtering
Filtered single-
control selection
cell dataset
+
+
+ Cell annotation
B Overall Statistics C # datasets by organism D # single-cell datasets by cell type
Datasets: 31 H1 (1) HEP (1) Macr. (1)
Tissues: 9 DC (1)
Mouse MCF7 (2) ESC (8)
Modalities: 6 (11) Human BJ (2)
Joint Assays: 7 (20)
Cells: 2,322,285 GM (3)
HSC (4) K562 (8)
1,000,000
100,000
10,000
1,000
100
10
1
ground-truth networks derived from heterogeneous evidence
types.26 While this variability complicates direct comparisons,
SC-MO-GRN-DB mitigates this issue by providing paired refer-
ence networks and curated input datasets for each cell type
and modality, which facilitates fair and standardized bench-
marking. Furthermore, most current GRN inference tools rely pri-
marily on ChIP-seq or other localization-based ground-truth net-
works.34 By including TF knock-down/knockout perturbation
networks in addition to localization and perturbation networks,
SC-MO-GRN-DB presents alternative approaches for ground
truth selection.
Ground-truth networks differ not only in scope but also in the
type of evidence supporting each TF-TG interaction. Traditional
GRN resources rely primarily on TF localization assays such as
ChIP-seq, which identify physical DNA-binding events but
cannot confirm their regulatory effect. SC-MO-GRN-DB ex-
pands beyond this paradigm by including perturbation-based
networks. These networks are derived from knockdown,
knockout, and CRISPR interference studies, which reveal the
220SD 320SD 130SD 230SD 030SD 420SD 200SD 330SD 100SD 700SD 520SD 120SD 020SD 310SD 600SD 500SD 400SD 210SD 910SD 300SD 720SD 620SD 920SD 820SD 010SD 510SD 110SD 610SD 410SD 710SD 810SD
E
DS ID:
sllec#
iScience
ll
OPEN ACCESS Article
Figure 4. Single-cell dataset curation pro-
cess and summary statistics
(A) Single-cell datasets with different modalities
(scRNA-Seq, scATAC-Seq, scChIP-Seq, scHiC,
and scBS-Seq) matching with cell types of the
reference networks are curated from the literature.
When there are multiple cell types in the data, the
relevant cell types are filtered using cell type an-
notations associated with the dataset. Cleaned,
filtered single-cell datasets are deposited in the
database.
(B) Overall statistics for the single-cell datasets in
the repository.
(C) Number of single-cell datasets by organism.
(D) Number of single-cell datasets for each indi-
vidual cell type. ESC, embryonic stem cells; Hep,
HepG2 human hepatocellular carcinoma cell line;
HSC, hematopoietic stem cells; DC, dendritic
cells; K562, human erythroleukemia cell line; GM,
GM12878 human B-lymphoblastoid cell line; BJ,
human fibroblast cell line; H1, human embryonic
stem cell line; Macr, macrophage; MCF7, human
breast cancer cell line.
(E) Number of cells per single-cell dataset.
(F) Browsing interface for single-cell datasets us-
ing data modality (such as scRNA-Seq and scA-
TAC-Seq), organism (human/mouse) and cell type
(such as embryonic stem cells, HSCs, and
macrophage).
functional consequences of TF perturba-
tion. We also provide dual-evidence net-
F works, which intersect localization and
perturbation data to provide high-confi-
Data Modality:Organism: Cell Type:
All All All dence interactions. We include text-
mined, literature-curated, and PPI-based
SelectDataset ID Data Modality OrganismCell TypeNumber of Cells PMID or DOI Download
networks to contribute breadth and
DS001 scRNA Mouse ESC 422 31907445 Download
contextual depth. Literature-curated net-
DS002 scRNA Mouse DC 384 31907445 Download
works deliver precision through expert
DS003 scRNA Mouse HSC 1072 31907445 Download
annotation but remain limited by scope.
Text-mined networks offer a large span-
ning, robust network size, but require careful interpretation to
avoid inaccurate associations. PPI-based networks capture reg-
ulatory cofactors and signaling interactions that shape TF activ-
ity but may not directly indicate transcriptional regulation.
Collectively, these complementary evidence types make SC-
MO-GRN-DB a flexible benchmarking resource adaptable to
both high-confidence and exploratory GRN studies.
Together with these curated reference networks, the accom-
panying single-cell datasets provide the essential foundation
for benchmarking GRN inference methods across modalities
and cell types. The variety of tissues across both organisms al-
lows systematic evaluation of model generalizability, while the
inclusion of multiomic and joint assays enables integrative ap-
proaches that more faithfully capture the multilayered nature of
transcriptional regulation. The scale of over two million cells sup-
ports robust statistical inference and motivates the use of
advanced computational methods that can fully leverage the
richness of the data. By providing complementary, cell type-spe-
cific datasets aligned with curated reference networks, SC-MO-
6 iScience 29, 115323, April 17, 2026

iScience
ll
Article OPEN ACCESS
GRN-DB supports the construction and benchmarking of GRN ances these trade-offs by providing both diverse and cell type-
inference tools. specific networks. By providing a range of network sizes,
SC-MO-GRN-DB is therefore designed to support flexible SC-MO-GRN-DB enables users to tailor analyses for either
mapping between single-cell datasets and reference networks comprehensive exploration or mechanistic precision.
derived from independent sources. For a given tissue or cell The curated pairing of networks and single-cell datasets is cen-
type, the database may include multiple single-cell datasets tral to SC-MO-GRN-DB’s design. Each dataset has been
and multiple reference networks, reflecting differences in exper- matched to an experimentally validated reference network repre-
imental conditions, molecular modalities, and evidence types. senting the same tissue, cell type, and biological condition. This
These resources are not intended to be paired in a one-to-one alignment ensures that benchmarking remains interpretable, as
manner. Instead, compatible datasets and reference networks both the reference network and the corresponding dataset
can be combined to support GRN construction and evaluation. describe the same regulatory context. Such systematic curation
Our recommendation is to leverage all available datasets rather reduces ambiguity in evaluation and allows for fair comparison
than relying on a single pairing, as this enables more robust between methods, even across distinct molecular modalities.
and reproducible analyses. For benchmarking and method eval- Unlike repositories that rely heavily on large-scale text-mining
uation, an all-to-all strategy, testing each single-cell dataset or automated inference,8,14,19 SC-MO-GRN-DB emphasizes
against each reference network, allows for the assessment of curated, experimentally supported content (Table S2). Each
consistency and performance across data sources. For biolog- network and dataset has been34 standardized and annotated
ical discovery, results may be combined across datasets and to facilitate reproducibility and meaningful comparison. This
reference networks to construct consensus GRNs that capture design principle ensures that the database delivers resources
shared regulatory signals. This design provides multiple, inde- with a strong focus on accuracy and biological relevance,
pendent lines of evidence to support GRN construction, valida- directly supporting the construction and evaluation of GRN infer-
tion, and comparative analysis. ence tools.
In parallel with the rapid expansion of single-cell datasets, a A defining strength of SC-MO-GRN-DB is its coverage of
wide range of machine learning and deep learning approaches molecular modalities encompassed within the single-cell data-
have been developed to model regulatory and disease-associ- sets. While many current GRN inference tools have focused pri-
ated processes in genomic data.35–37 These methods span marily on transcriptomics or chromatin accessibility as input for
diverse applications, including regulatory interaction prediction, GRN inference,41 SC-MO-GRN-DB incorporates scRNA-seq,
GRN inference, cell clustering and annotation, and disease as- scATAC-seq, scChIP-seq, scDNA methylation, scHi-C, and
sociation modeling, often leveraging attention mechanisms, scCRISPR-seq datasets. This diversity enables researchers to
graph-based representations, or deep autoencoder frame- investigate regulatory mechanisms that extend beyond gene
works.38–40While such approaches have demonstrated strong expression alone, such as TF binding, epigenetic modifications,
performance within specific tasks and datasets, they are typi- three-dimensional chromatin organization, and causal perturba-
cally evaluated using custom data processing pipelines tion effects.42Joint assays, such as scRNA-seq with scATAC-
and heterogeneous reference networks, limiting reproducibility seq, provide complementary perspectives on gene regulation
and cross-method comparison. By providing standardized, by simultaneously measuring transcriptional output and chro-
experimentally supported reference networks paired with matin accessibility. Similarly, scHi-C captures long-range regu-
tissue-matched single-cell multiomic datasets, SC-MO-GRN- latory interactions, while scCRISPR-seq offers direct perturba-
DB offers a unifying resource to support the systematic tion-based evidence of causality. Together, these multiomic
evaluation and future development of these computational layers create new opportunities for tool development that more
approaches. faithfully capture the complexity of transcriptional regulation.
Beyond evidence type, the scale of both networks and data- By providing such a comprehensive resource, SC-MO-GRN-
sets strongly influences their interpretive and benchmarking DB positions itself as a foundation for next-generation GRN infer-
value. By encompassing over 1,000 TFs and tens of thousands ence methods that integrate across molecular modalities.
of TGs, SC-MO-GRN-DB captures broad regulatory landscapes By pairing high-confidence reference networks with curated
suitable for large-scale method evaluation. The database in- multiomic datasets, SC-MO-GRN-DB serves not only as a repos-
cludes networks ranging from small, focused systems with itory of curated data but also as a practical platform for advancing
only a few regulators and targets to extensive, high-density net- both the benchmarking and construction of GRNs in single-cell
works containing nearly 1,000 TFs, nearly 30,000 TGs, and more multiomic contexts. This resource reduces the burden on re-
than a million edges. This variability reflects the underlying searchers to locate and standardize datasets, which will help to
experimental scope and biological diversity of the source data. accelerate the development of more accurate multimodal GRN
Network properties such as the number of TFs, TGs, and edges inference tools. This comprehensive framework supports the
determine how much regulatory complexity is represented, while exploration of transcriptional regulation across multiple molecular
dataset properties, particularly cell count and gene coverage, layers, from gene expression and chromatin accessibility to epige-
define the resolution at which those networks can be inferred. netic modifications, three-dimensional genome architecture, and
Larger networks and datasets offer statistical robustness and perturbation-driven causality. As a result, SC-MO-GRN-DB en-
broader coverage but may introduce noise and redundancy. ables the systematic construction and evaluation of context-spe-
Smaller or more focused systems yield clearer signals but risk cific GRNs, empowering researchers to generate insights into
underrepresenting biological diversity. SC-MO-GRN-DB bal- cellular regulatory programs, disease mechanisms, and potential
iScience 29, 115323, April 17, 2026 7

iScience
ll
OPEN ACCESS Article
therapeutic targets in a reproducible and biologically meaningful ACKNOWLEDGMENTS
manner.
SC-MO-GRN-DB is curated from currently available public The research reported in this publication was supported by the National Insti-
tute of General Medical Sciences of the National Institutes of Health under
data sources, including reference networks and single-cell data-
award number R35GM150616.
sets, with a specific emphasis on multiomics. The single-cell ge-
nomics field is growing rapidly, and new datasets are being AUTHOR CONTRIBUTIONS
shared with the research community continuously. As a living
resource, SC-MO-GRN-DB is designed for continued expansion Conceptualization was performed by Y.U.; data curation was carried out by
as new tissues, modalities, and evidence types become avail- H.V., Y.U., and R.E.; formal analysis, project administration, visualization,
and the majority of the investigation were conducted by H.V.; methodology
able. In this sense, it serves as a seed repository, a scalable foun-
was developed by H.V. with input from Y.U.; software development was led
dation for GRN inference, capable of growing into a central hub
by H.V. with technical support from S.E.O.; validation and quality control
for the development, comparison, and validation of computa- were performed by H.V., K.K., and E.M.; resources were provided by Y.U.; su-
tional models of transcriptional regulation, similar to the UCI re- pervision was provided by Y.U.; writing of the original draft was performed by
pository for machine learning.43 H.V.; and writing-review and editing were contributed by H.V., K.K., E.M., and
Y.U.
Limitations of the study DECLARATION OF INTERESTS
The present release includes data from 16 tissue types and is
Y.U. is the owner of Systems Biology Consulting & Analytics LLC (Systems
biased toward well-characterized cell lines and abundant cell
Bio). However, Systems Bio did not provide any support for this manuscript
populations, while rare cell types and disease-associated states
and has no financial or other interests related to its content.
remain underrepresented. This limitation reflects the availability
of high-confidence reference networks and matched single- STAR★METHODS
cell multiomic datasets rather than constraints of the database
framework. In particular, the construction of ground truth refer- Detailed methods are provided in the online version of this paper and include
the following:
ence networks requires experimentally validated regulatory in-
teractions, and for many TFs, such evidence is only available • KEY RESOURCES TABLE
through dedicated perturbation or ChIP-based experiments. • EXPERIMENTAL MODEL AND STUDY PARTICIPANT DETAILS
As a result, ground truth coverage remains incomplete and un- • METHOD DETAILS
even across TFs and cell types. To address this limitation, we ○ Construction of reference networks
plan to continuously expand the database by incorporating ○ Curation of single-cell datasets
○ Interactive browsing
new ground truth networks and single-cell datasets. A key focus
○ Interactive search
of future development will be the inclusion of more truly paired ○ Analytical functions
multi-omic datasets. • QUANTIFICATION AND STATISTICAL ANALYSIS
While integration methods exist for unpaired modalities, such
approaches can introduce uncertainty and batch effects, poten- SUPPLEMENTAL INFORMATION
tially reducing the reliability of downstream GRN inference.
Therefore, prioritizing true multiomics (joint assay) datasets will Supplemental information can be found online at https://doi.org/10.1016/j.isci.
2026.115323.
provide better data alignment, reduce preprocessing burden
for users, and enhance the robustness of multi-modal analyses.
Received: November 13, 2025
Overall, the database is designed as a living resource, with Revised: January 9, 2026
ongoing updates aimed at improving data coverage, methodo- Accepted: March 6, 2026
logical reliability, and support for increasingly complex GRN Published: March 11, 2026
inference workflows.
REFERENCES
1.Lee, T.I., and Young, R.A. (2013). Transcriptional regulation and its misre-
RESOURCE AVAILABILITY
gulation in disease. Cell 152, 1237–1251.
Lead contact 2.Macneil, L.T., and Walhout, A.J.M. (2011). Gene regulatory networks and
Requests for further information should be directed to and will be fulfilled by the role of robustness and stochasticity in the control of gene expression.
the lead contact, Yasin Uzun (yuzun@pennstatehealth.psu.edu). Genome Res. 21, 645–657.
3.Feigin, C., Li, S., Moreno, J., and Mallarino, R. (2023). The GRN concept as
a guide for evolutionary developmental biology. J. Exp. Zool. B Mol. Dev.
Materials availability Evol. 340, 92–104.
This study did not generate new reagents.
4.Davidson, E.H., and Erwin, D.H. (2006). Gene regulatory networks and the
evolution of animal body plans. Science 311, 796–800.
Data and code availability 5.Ideker, T., Galitski, T., and Hood, L. (2001). A new approach to decoding
• SC-MO-GRN-DB is publicly accessible at https://scmogrndb.psu.edu. life: systems biology. Annu. Rev. Genomics Hum. Genet. 2, 343–372.
• The scripts that are used to build SC-MO-GRN-DB are deposited into 6.Davidson, E.H., Rast, J.P., Oliveri, P., Ransick, A., Calestani, C., Yuh, C.H.,
and publicly accessible at: https://github.com/UzunLab/SC-MO- Minokawa, T., Amore, G., Hinman, V., Arenas-Mena, C., et al. (2002). A
GRN-DB. genomic regulatory network for development. Science 295, 1669–1678.
8 iScience 29, 115323, April 17, 2026

iScience
ll
Article OPEN ACCESS
7.Huynh-Thu, V.A., Irrthum, A., Wehenkel, L., and Geurts, P. (2010). Inferring by Shared Single-Cell Profiling of RNA and Chromatin. Cell 183, 1103–
regulatory networks from expression data using tree-based methods. 1116.e20.
PLoS One 5, e12776.
26.Karamveer & Uzun, Y., and Uzun, Y. (2024). Approaches for Benchmarking
8.Pratapa, A., Jalihal, A.P., Law, J.N., Bharadwaj, A., and Murali, T.M. Single-Cell Gene Regulatory Network Methods. Bioinform. Biol. Insights
(2020). Benchmarking algorithms for gene regulatory network inference 18, 11779322241287120.
from single-cell transcriptomic data. Nat. Methods 17, 147–154.
27.Xu, H., Baroukh, C., Dannenfelser, R., Chen, E.Y., Tan, C.M., Kou, Y., Kim,
9.Bravo Gonza´lez-Blas, C., De Winter, S., Hulselmans, G., Hecker, N., Ma- Y.E., Lemischka, I.R., and Ma’ayan, A. (2013). ESCAPE: database for inte-
tetovici, I., Christiaens, V., Poovathingal, S., Wouters, J., Aibar, S., and grating high-content published data collected from human and mouse em-
Aerts, S. (2023). SCENIC+: single-cell multiomic inference of enhancers bryonic stem cells. Database 2013, bat045.
and gene regulatory networks. Nat. Methods 20, 1355–1367.
28.Marbach, D., Costello, J.C., Ku¨ffner, R., Vega, N.M., Prill, R.J., Camacho,
10.Yuan, Q., and Duren, Z. (2025). Inferring gene regulatory networks from D.M., Allison, K.R., DREAM5 Consortium; Kellis, M., Collins, J.J., and Sto-
single-cell multiome data using atlas-scale external data. Nat. Biotechnol. lovitzky, G. (2012). Wisdom of crowds for robust gene network inference.
43, 247–257. Nat. Methods 9, 796–804.
11.Keenan, A.B., Torre, D., Lachmann, A., Leong, A.K., Wojciechowicz, M.L., 29.Sikora-Wohlfeld, W., Ackermann, M., Christodoulou, E.G., Singaravelu,
Utti, V., Jagodnik, K.M., Kropiwnicki, E., Wang, Z., and Ma’ayan, A. (2019). K., and Beyer, A. (2013). Assessing computational methods for transcrip-
ChEA3: transcription factor enrichment analysis by orthogonal omics inte- tion factor target gene identification based on ChIP-seq data. PLoS Com-
gration. Nucleic Acids Res. 47, W212–W224. put. Biol. 9, e1003342.
12.Garcia-Alonso, L., Holland, C.H., Ibrahim, M.M., Turei, D., and Saez-Ro- 30.Aibar, S., Gonza´lez-Blas, C.B., Moerman, T., Huynh-Thu, V.A., Imrichova,
driguez, J. (2019). Benchmark and integration of resources for the estima- H., Hulselmans, G., Rambow, F., Marine, J.C., Geurts, P., Aerts, J., et al.
tion of human transcription factor activities. Genome Res. 29, 1363–1375. (2017). SCENIC: single-cell regulatory network inference and clustering.
13.Liu, Z.-P., Wu, C., Miao, H., and Wu, H. (2015). RegNetwork: an integrated Nat. Methods 14, 1083–1086.
database of transcriptional and post-transcriptional regulatory networks in 31.Szklarczyk, D., Kirsch, R., Koutrouli, M., Nastou, K., Mehryary, F., Hachilif,
human and mouse. Database 2015, bav095. R., Gable, A.L., Fang, T., Doncheva, N.T., Pyysalo, S., et al. (2023). The
14.Han, H., Cho, J.W., Lee, S., Yun, A., Kim, H., Bae, D., Yang, S., Kim, C.Y., STRING database in 2023: protein-protein association networks and func-
Lee, M., Kim, E., et al. (2018). TRRUST v2: an expanded reference data- tional enrichment analyses for any sequenced genome of interest. Nucleic
base of human and mouse transcriptional regulatory interactions. Nucleic Acids Res. 51, D638–D646.
Acids Res. 46, D380–D386. 32.Oughtred, R., Rust, J., Chang, C., Breitkreutz, B.J., Stark, C., Willems, A.,
15.Angelini, C., and Costa, V. (2014). Understanding gene regulatory mecha- Boucher, L., Leung, G., Kolas, N., Zhang, F., et al. (2021). The BioGRID
nisms by integrating ChIP-seq and RNA-seq data: statistical solutions to database: A comprehensive biomedical resource of curated protein, ge-
biological problems. Front. Cell Dev. Biol. 2, 51. netic, and chemical interactions. Protein Sci. 30, 187–200.
16.Zhang, Q., Liu, W., Zhang, H.M., Xie, G.Y., Miao, Y.R., Xia, M., and Guo, 33.Altman, R.B., Bergman, C.M., Blake, J., Blaschke, C., Cohen, A., Gannon,
A.Y. (2020). hTFtarget: A Comprehensive Database for Regulations of Hu- F., Grivell, L., Hahn, U., Hersh, W., Hirschman, L., et al. (2008). Text mining
man Transcription Factors and Their Targets. Genom. Proteom. Bio- for biology–the way forward: opinions from leading scientists. Genome
inform. 18, 120–128. Biol. 9, S7.
17.Ben Guebila, M., Lopes-Ramos, C.M., Weighill, D., Sonawane, A., Bur- 34.Loers, J.U., and Vermeirssen, V. (2024). A single-cell multimodal view on
kholz, R., Shamsaei, B., Platig, J., Glass, K., Kuijjer, M., and Quackenbush, gene regulatory network inference from transcriptomics and chromatin
J. (2022). GRAND: a database of gene regulatory network models across accessibility data. Brief. Bioinform. 25, bbae382.
human conditions. Nucleic Acids Res. 50, D610–D621.
35.Yuan, L., Zhao, L., Lai, J., Jiang, Y., Zhang, Q., Shen, Z., Zheng, C.H., and
18.Fang, L., Li, Y., Ma, L., Xu, Q., Tan, F., and Chen, G. (2021). GRNdb: de- Huang, D.S. (2024). iCRBP-LKHA: Large convolutional kernel and hybrid
coding the gene regulatory networks in diverse human and mouse condi- channel-spatial attention for identifying circRNA-RBP interaction sites.
tions. Nucleic Acids Res. 49, D97–D103. PLoS Comput. Biol. 20, e1012399.
19.Huang, X., Song, C., Zhang, G., Li, Y., Zhao, Y., Zhang, Q., Zhang, Y., Fan, 36.Yuan, L., Zhao, J., Shen, Z., Zhang, Q., Geng, Y., Zheng, C.H., and Huang,
S., Zhao, J., Xie, L., and Li, C. (2024). scGRN: a comprehensive single-cell D.S. (2023). iCircDA-NEAE: Accelerated attribute network embedding and
gene regulatory network platform of human and mouse. Nucleic Acids dynamic convolutional autoencoder for circRNA-disease associations
Res. 52, D293–D303. prediction. PLoS Comput. Biol. 19, e1011344.
20.Fiers, M.W.E.J., Minnoye, L., Aibar, S., Bravo Gonza´lez-Blas, C., Kalender 37.Yuan, L., Sun, S., Zhang, Q., Li, H.T., Shen, Z., Hu, C., Zhao, X., Ye, L.,
Atak, Z., and Aerts, S. (2018). Mapping gene regulatory networks from sin- Zheng, C.H., and Huang, D.S. (2024). Identification of ferroptosis-related
gle-cell omics data. Brief. Funct. Genomics 17, 246–254. lncRNAs for predicting prognosis and immunotherapy response in non-
21.Cha, J., and Lee, I. (2025). Single-cell network biology enabling cell-type- small cell lung cancer. Future Gener. Comput. Syst. 159, 204–220.
resolved disease genetics. Genomics Inform. 23, 10. 38.Yuan, L., Zhao, L., Jiang, Y., Shen, Z., Zhang, Q., Zhang, M., Zheng, C.H.,
22.Tabula Sapiens Consortium, Jones, R.C., Karkanias, J., Krasnow, M.A., and Huang, D.S. (2024). scMGATGRN: a multiview graph attention
Pisco, A.O., Quake, S.R., Salzman, J., Yosef, N., Bulthaup, B., Brown, network-based method for inferring gene regulatory networks from sin-
P., et al. (2022). The Tabula Sapiens: A multiple-organ, single-cell tran- gle-cell transcriptomic data. Brief. Bioinform. 25, bbae526.
scriptomic atlas of humans. Science 376, eabl4896. 39.Yuan, L., Xu, Z., Meng, B., and Ye, L. (2025). scAMZI: attention-based
23.Chen, S., Lake, B.B., and Zhang, K. (2019). High-throughput sequencing deep autoencoder with zero-inflated layer for clustering scRNA-seq
of the transcriptome and chromatin accessibility in the same cell. Nat. Bio- data. BMC Genom. 26, 350.
technol. 37, 1452–1457. 40.Yuan, L., Sun, S., Jiang, Y., Zhang, Q., Ye, L., Zheng, C.H., and Huang,
24.Liu, L., Liu, C., Quintero, A., Wu, L., Yuan, Y., Wang, M., Cheng, M., Leng, D.S. (2024). scRGCL: a cell type annotation method for single-cell RNA-
L., Xu, L., Dong, G., et al. (2019). Deconvolution of single-cell multi-omics seq data using residual graph convolutional neural network with contras-
layers reveals regulatory heterogeneity. Nat. Commun. 10, 470. tive learning. Brief. Bioinform. 26, bbae662.
25.Ma, S., Zhang, B., LaFave, L.M., Earl, A.S., Chiang, Z., Hu, Y., Ding, J., 41.Badia-I-Mompel, P., Wessels, L., Mu¨ller-Dott, S., Trimbour, R., Ramirez
Brack, A., Kartha, V.K., Tay, T., et al. (2020). Chromatin Potential Identified Flores, R.O., Argelaguet, R., and Saez-Rodriguez, J. (2023). Gene
iScience 29, 115323, April 17, 2026 9

iScience
ll
OPEN ACCESS Article
regulatory network inference in the era of single-cell multi-omics. Nat. Rev. Seq and chromatin accessibility data in human and mouse. Nucleic Acids
Genet. 24, 739–754. Res. 45, D658–D662.
42. Stock, M., Losert, C., Zambon, M., Popp, N., Lubatti, G., Ho¨rmanseder, E.,
48.Raghav, P., Kumar, R., Lathwal, A., and Sharma, N. (2024). Computational
Heinig, M., and Scialdone, A. (2025). Leveraging prior knowledge to infer
Biology for Stem Cell Research (Elsevier).
gene regulatory networks from single-cell RNA-sequencing data. Mol.
Syst. Biol. 21, 214–230. 49.Feng, C., Song, C., Song, S., Zhang, G., Yin, M., Zhang, Y., Qian, F., Wang,
43. Dua, D., and Graff, C. (2017). UCI Machine Learning. http://archive.ics.uci. Q., Guo, M., and Li, C. (2024). KnockTF 2.0: a comprehensive gene
edu/ml. expression profile database with knockdown/knockout of transcription
(co-)factors in multiple species. Nucleic Acids Res. 52, D183–D193.
44. Frankish, A., Carbonell-Sala, S., Diekhans, M., Jungreis, I., Loveland, J.E.,
Mudge, J.M., Sisu, C., Wright, J.C., Arnan, C., Barnes, I., et al. (2023).
50.Shen, W.-K., Chen, S.Y., Gan, Z.Q., Zhang, Y.Z., Yue, T., Chen, M.M., Xue,
GENCODE: reference annotation for the human and mouse genomes in
Y., Hu, H., and Guo, A.Y. (2023). AnimalTFDB 4.0: a comprehensive animal
2023. Nucleic Acids Res. 51, D942–D949.
transcription factor database updated with variation and expression anno-
45. Zou, Z., Ohta, T., and Oki, S. (2024). ChIP-Atlas 3.0: a data-mining suite to tations. Nucleic Acids Res. 51, D39–D45.
explore chromosome architecture together with large-scale regulome
data. Nucleic Acids Res. 52, W45–W53. 51.Parkinson, H., Kapushesky, M., Shojatalab, M., Abeygunawardena, N.,
Coulson, R., Farne, A., Holloway, E., Kolesnykov, N., Lilja, P., Lukk, M.,
46. Zheng, R., Wan, C., Mei, S., Qin, Q., Wu, Q., Sun, H., Chen, C.H., Brown,
et al. (2007). ArrayExpress–a public database of microarray experiments
M., Zhang, X., Meyer, C.A., and Liu, X.S. (2019). Cistrome Data Browser:
and gene expression profiles. Nucleic Acids Res. 35, D747–D750.
expanded datasets and new tools for gene regulatory analysis. Nucleic
Acids Res. 47, D729–D735. 52.Peidli, S., Green, T.D., Shen, C., Gross, T., Min, J., Garda, S., Yuan, B.,
47. Mei, S., Qin, Q., Wu, Q., Sun, H., Zheng, R., Zang, C., Zhu, M., Wu, J., Shi, Schumacher, L.J., Taylor-King, J.P., Marks, D.S., et al. (2024). scPerturb:
X., Taing, L., et al. (2017). Cistrome Data Browser: a data portal for ChIP- harmonized single-cell perturbation data. Nat. Methods 21, 531–540.
10 iScience 29, 115323, April 17, 2026

iScience
ll
Article OPEN ACCESS
STAR★METHODS
KEY RESOURCES TABLE
REAGENT or RESOURCE SOURCE IDENTIFIER
Deposited data
Constructed reference networks This paper https://scmogrndb.psu.edu/
Filtered single cell datasets This paper https://scmogrndb.psu.edu/
Software and algorithms
Code for constructing reference networks This paper https://github.com/UzunLab/SC-MO-
from raw ChIP-seq data GRN-DB
R version 4.2.3 R Core Team https://www.R-project.org/
GENCODE v38 human (hg38) mouse Frankish, A. et al.44 https://www.gencodegenes.org/
(mm39)
Other
Website This paper https://scmogrndb.psu.edu/
EXPERIMENTAL MODEL AND STUDY PARTICIPANT DETAILS
This study did not generate new experimental data or involve human participants, animals, cell lines, or primary cell cultures. All an-
alyses were performed using publicly available datasets obtained from previously published studies and community repositories.
Therefore, no experimental model or study participant information is applicable.
Sex, age, genotype, and sample allocation were not controlled by the present study and were reported as provided by the original
sources. Limitations related to incomplete metadata in original datasets are acknowledged.
METHOD DETAILS
Construction of reference networks
To establish a foundation of high-confidence GRNs, we systematically collected multiple complementary evidence types. Each line
of evidence provides unique strengths for capturing TF-TG interactions, and their intersection reduces bias and increases network
reliability. TF localization data were derived from ChIP-seq and ChIP-chip experiments, which map genome-wide TF binding sites.
Datasets were collected from BEELINE,8 ChIP-Atlas,45 Cistrome,46,47 and ESCAPE,27,48 as well as from individual publications
(Table S3). The binding peaks were assigned to the gene with the nearest transcription start site (TSS) based on the genomic dis-
tance, using the ‘‘nearestTSS’’ function in edgeR package (v4.0.16) with gene annotations from GENCODE44 v38 for human
(hg38) and mouse (mm39) genomes.
TF perturbation networks were constructed from loss-of-function and gain-of-function experiments, including TF knockdown,
knockout, and overexpression assays curated from KnockTF2.0,49ESCAPE,27,48and additional primary publications (Table S3).
The ground-truth networks were standardized into two-column lists of TF-to-TG interactions for consistency across networks. In
addition, protein-protein interaction data from STRING31(obtained through BEELINE8) were incorporated to capture regulatory com-
plexes, recognizing that such interactions represent indirect regulation rather than direct TF-DNA binding. Text-mined networks were
obtained from TRRUST,14which uses sentence-level text mining, along with individual studies that provided experimentally validated
TF-TG relationships. Literature-curated networks were obtained from BEELINE,8which selected published Boolean models from
four different cell types in both mouse and human. To further increase confidence, we generated dual-evidence networks by inter-
secting localization- and perturbation-based networks within the same cellular context.
Once networks were standardized, they were filtered using only known TFs50for both human and mouse to retain only TF-TG in-
teractions. All reference networks were annotated with metadata, including article identifier, organism, tissue or cell type, evidence
type, and summary statistics (number of TFs, TGs, and edges). The full repository of networks and their associated metadata is
accessible through the SC-MO-GRN-DB web interface.
Curation of single-cell datasets
To provide high-quality input data for benchmarking and inference of GRNs, we curated single-cell multiomic datasets spanning a
wide range of tissues, modalities, and organisms (Table S4). Datasets were obtained from BEELINE,8publicly available repositories
such as Gene Expression Omnibus (GEO) and ArrayExpress,51community resources including scPerturb,52and individual publica-
tions (Table S5). We retained only datasets that included wild-type or untreated control cells, as perturbed or genetically modified
iScience 29, 115323, April 17, 2026 e1

iScience
ll
OPEN ACCESS Article
cells may confound baseline regulatory relationships. When datasets included multiple annotated cell types, cells were stratified into
homogenous groups according to provided metadata, ensuring that each dataset corresponded to a single, unmodified cellular
state. For each dataset, we excluded cells with poor quality metrics to ensure consistency across studies. All curated datasets
were annotated with information about organism, tissue or cell type, molecular modality, number of cells, and primary literature refer-
ence (Table S5). The complete collection, including article identifiers and summary statistics, is made available through the SC-MO-
GRN-DB web interface.
Interactive browsing
The SC-MO-GRN-DB web interface was designed to facilitate intuitive exploration of reference networks and curated single-cell da-
tasets. Reference networks are indexed by evidence type, organism, and cell type, enabling users to rapidly identify networks of in-
terest based on experimental design or biological context. Each entry is annotated with standardized metadata, including the number
of TFs, TG, regulatory edges, and source publication, thereby allowing direct comparison across networks. Similarly, single-cell data-
sets are indexed by modality, organism, and cell type, with annotations for the number of cells and originating study. Individual entries
can be browsed to inspect specific datasets or perform batch selection and download for larger analyses. The browsing interface
was optimized for both usability and reproducibility, ensuring that researchers can efficiently identify and obtain the resources
most relevant to their work.
Interactive search
To support targeted exploration of regulatory relationships, we implemented a search module that allows queries by TF, TG, or TF-TG
pair. Searches by TF return all networks containing the specified TF, along with the number of unique targets it regulates in each
context. Searches by TG return all networks in which the gene appears as a target, while searches by TF-TG pair return only networks
in which the specified interaction is present. All results are displayed with network-level metadata, including evidence type, organism,
cell type, number of edges, and associated article identifier, ensuring that users can assess the context and reliability of each inter-
action. Search results can be exported individually or as batches, allowing seamless integration into downstream analyses. This func-
tionality supports both hypothesis-driven investigations, for example, querying a candidate TF for its regulatory footprint across tis-
sues, and broader network-level benchmarking tasks. By combining precise search capabilities with standardized metadata, the
search module enhances the flexibility and practical utility of SC-MO-GRN-DB for diverse research applications.
Analytical functions
To support quantitative benchmarking and downstream evaluation of GRNs, SC-MO-GRN-DB provides an analysis module that im-
plements standardized comparison, enrichment, and network-level control analyses. For GRN benchmarking, user-uploaded net-
works are compared against selected SC-MO-GRN-DB reference networks using edge-level overlap. True positives (TP) are defined
as inferred TF-TG interactions present in the reference network, false positives (FP) as inferred interactions absent from the reference,
and false negatives (FN) as reference interactions not recovered by the inferred network. From these quantities, performance metrics
including precision recall, F1 score, Jaccard index, and area under the precision-recall curve (AUC) are computed. Summary statis-
tics, including the total number of user-uploaded edges, reference edges, and TP, FP, and FN counts, are reported for each
comparison.
To enable functional interpretation of regulatory relationships, enrichment analysis is performed using a hypergeometric test. Given
a TF-specific target gene set, a pathway or gene set of interest, and a background gene set representing all expressed genes, the test
evaluates whether TF targets are significantly overrepresented in the pathway relative to the background. Enrichment statistics and
associated p-values are reported, allowing hypothesis-driven assessment of TF-pathway associations.
To characterize and compare reference network structure, control analyses based on graph-theoretic measures are implemented.
For a given reference network, centrality metrics, including degree centrality, betweenness centrality, closeness centrality, and
PageRank, are computed for each TF, providing complementary measures of regulatory influence within the network. In addition,
pairwise comparison of reference networks is supported by quantifying overlap between TF target sets. For each TF shared between
two networks, the number of common targets and targets unique to each network are calculated, enabling direct comparison of reg-
ulatory programs across networks.
QUANTIFICATION AND STATISTICAL ANALYSIS
All quantitative and statistical analyses in this study were performed using the comparison, enrichment, and control functions
described in the method detailssection. These analyses were used to evaluate overlap between networks, assess enrichment of
transcription factor target genes in gene sets of interest, and summarize properties of reference networks.
Exact values of n (defined as the number of transcription factors, target genes, regulatory edges, or networks, depending on the
analysis) are reported in the corresponding tables and figures. Because this study did not involve experimental group comparisons or
biological replicates, no sample size estimation or randomization procedures were applicable.
e2 iScience 29, 115323, April 17, 2026