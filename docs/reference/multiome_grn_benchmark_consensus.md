# Which benchmark should MeVD-GRN use? A literature consensus (2026-09-30)

**Question.** Across published work, which dataset, ground truth and protocol
is the most accepted benchmark for inferring TF→target gene regulatory
networks (GRNs) from *paired* scRNA-seq + scATAC-seq? High-impact journals
(Nature family, Genome Biology, NAR, and similar) count for more.

**Method.** Three independent web surveys, each covering different venues so
they didn't overlap:

| Survey | Papers | Per-paper log |
|---|---|---|
| High-impact journals (Nature family, Cell family, Science family, Genome Biol/Res, NAR, MSB) | 42, full text verified | [lit_survey/journals.md](lit_survey/journals.md) |
| ML venues, bioinformatics journals, preprints (NeurIPS/ICLR/ICML/KDD, ISMB, Bioinformatics, Brief Bioinform, bioRxiv, arXiv) | 59 | [lit_survey/ml_preprints.md](lit_survey/ml_preprints.md) |
| Benchmarks, reviews and ground-truth resources | 33 | [lit_survey/benchmarks.md](lit_survey/benchmarks.md) |

The 17 papers in `~/research/` are covered separately in
[literature_knowledge_base.md](literature_knowledge_base.md). The downloaded
full texts are in `~/.cache/lit_survey/` (not in git).

---

## 1. The answer

**Dataset: 10x Genomics "PBMC from a healthy donor, granulocytes removed
through cell sorting (10k)" Multiome (`pbmc_granulocyte_sorted_10k`).**
**Ground truth: blood-cell ChIP-seq from Cistrome, scored the way LINGER does.**

All three surveys reached this answer independently.

- **It is truly paired.** 10x Multiome (ATAC + Gene Expression) measures RNA
  and chromatin accessibility in the *same nucleus*. The two modalities share
  barcodes.
- **It is the most-used dataset by a wide margin.**
  - Journals survey: 17 papers use this exact file, and 22 if you count
    papers that use a 10x PBMC Multiome without naming the variant.
  - ML/preprint survey: 11 of the 15 paired papers that evaluate edges
    against a ground truth use it.
  - No other paired dataset appears in more than 8 papers (SHARE-seq mouse
    skin).
- **It is used in the top journals.** Users include SCENIC+ (Nat Methods
  2023), LINGER (Nat Biotechnol 2024), scGLUE (Nat Biotechnol 2022), DeepMAPS,
  scTFBridge, Epiregulon and regX (all Nat Commun), KEGNI (Genome Biol 2025),
  SCRIP (NAR 2022), DIRECT-NET (Sci Adv 2022), pgBoost (Nat Genet 2025), and
  the GRETA benchmark (bioRxiv 2024/25).
- **It is effectively a shared leaderboard.** LINGER's Cistrome evaluation was
  reused head-to-head by scTFBridge (Nat Commun 2025), KEGNI (Genome Biol 2025,
  using LINGER's own processed data) and regX (Nat Commun 2025). TRIPOD (Cell
  Syst) and Dictys (Nat Methods) independently chose Cistrome blood ChIP as
  well.

> Note: EpiAwareNet (KDD 2026) uses a *different* 10x PBMC file ("10k Human
> PBMCs Multiome v1.0, Chromium Controller") with DoRothEA labels. It is a
> useful second comparison, but it is not the standard file.

### Rejected alternatives

| Candidate | Why it's not the primary benchmark |
|---|---|
| scMultiomeGRN fetal lung (GSM4508936) | Used by 1 paper. The labels are motif scans of the same ATAC used as input, which is circular. The data is probably not truly paired (sci-ATAC3 + sci-RNA3). |
| BEELINE (7 scRNA sets) | The standard for *RNA-only* papers (about 25 of 29). It has no ATAC, so it cannot show what ATAC adds. Use it only if reviewers ask for GENELink-family comparisons. |
| SHARE-seq mouse skin (GSE140203) | Second most-used paired dataset (8 papers). Its ChIP ground truth is weak (Dictys: 33 ChIP sets / 16 TFs). A good secondary benchmark. |
| NeurIPS 2021 BMMC (GSE194122) | Large and paired, but scTFBridge reports no GRN ground truth for it. Its main use is geneRNIB's label-free evaluation. |
| DoRothEA / TRRUST as the main labels | Each paper uses a different curated set, so their numbers don't compare with each other. DoRothEA levels C–E include motif-derived edges, which is circular for an ATAC model. |

---

## 2. Numbers to beat (PBMC10k, LINGER-style Cistrome ChIP evaluation)

| Method | Venue | Metric | Value | Source |
|---|---|---|---|---|
| LINGER | Nat Biotechnol 2024 | average AUROC over the blood ChIP set | **0.714** | KEGNI's re-run (Genome Biol 2025) |
| LINGER | Nat Biotechnol 2024 | average AUPR ratio | about 2.25 | derived from "1.25 units above random"; figure-level, approximate |
| KEGNI | Genome Biol 2025 | average AUROC | 0.699 | KEGNI paper |
| scTFBridge | Nat Commun 2025 | STAT1 in CD14 monocytes: AUC / AUPR ratio | 0.693 / 2.221 | published version (the bioRxiv version says 0.776) |
| LINGER | — | STAT1 in CD14 monocytes: AUC | 0.719 | scTFBridge Fig. 5 |
| other methods (PCC, GENIE3, PIDC, SCENIC+, scMTNI, scMultiomeGRN) | — | AUPR ratio / AUC | about 1.17–1.29 / about 0.51–0.56 | LINGER, scTFBridge |

- **The real bar is LINGER.** It is pretrained on external ENCODE bulk data,
  so beating it without external pretraining is a strong result.
- **Scope of the ground truth:** it covers only **4 cell types**. KEGNI's cell
  counts are classical monocytes 1,848, naive CD4 T 1,373, naive B 282 and
  mDC 232. LINGER's set has 20 blood ChIP datasets.

**Secondary comparison (EpiAwareNet's protocol):** its own PBMC file with
DoRothEA labels and a random 80/20 edge split. EpiAwareNet reports AUROC
0.8218 (±5 kb) and 0.8708 (±100 kb). See the literature knowledge base,
section 2.9.

---

## 3. Protocol for a *supervised* model like MeVD-GRN

**Paired multiome has no accepted supervised split protocol.** The literature
is also now strongly critical of the random edge splits that most supervised
link-prediction papers use:

- **InfoSEM (ICML 2025):** a logistic regression that sees only one-hot gene
  IDs, with no expression, matches GENELink and scGREAT on the standard
  BEELINE random split (hESC AUPRC 0.600 vs 0.642). On unseen genes,
  supervised methods drop 42% (cell-type ChIP) or 79% (non-specific ChIP).
- **Stock et al. (bioRxiv, Dec 2025):** a simple out-degree sorter matches the
  GNNs. With degree-matched negatives, GENELink, GNNLink and GT-GRN all fall to
  near-random AUROC.
- **GLM-Prior:** on gene- or TF-disjoint splits, AUPRC is barely above chance
  (hESC 0.19 vs a 0.15 base rate).

For MeVD-GRN, this means:

1. **Labels.**
   - The Cistrome ChIP set is only about 20 TFs, too few to train on.
   - **Evaluate on the LINGER Cistrome set** for comparability, and keep those
     TFs out of training entirely.
   - **Train on a larger, non-circular label set.** Options: DoRothEA levels
     A–B only, CollecTRI, or ChIP-Atlas blood ChIP for the other TFs.
2. **Splits.** Report all three:
   - **TF-disjoint:** the headline, and the cleanest match to LINGER's setup.
   - **Target-gene-disjoint.**
   - **Random edge split:** only for comparability with EpiAwareNet-style
     numbers.
3. **Negatives.** Use degree-matched negatives, not uniform random ones.
   Stratify or weight results by TF degree.
4. **Trivial baselines are mandatory:**
   - a gene-ID-only model (InfoSEM-style)
   - a degree-only ranker
   - GRNBoost2, and Pearson correlation
5. **Other baselines to re-run on the same splits:** LINGER (bar), SCENIC+,
   scTFBridge, and ideally CellOracle, FigR and Pando.
6. **Rigor already in the repo:** 5 seeds, model selection on validation only,
   the negative-sampling leak fix (branch `leak-fix-and-benchmarks`), and
   label-free construction of the prior graphs.

---

## 4. Supporting and secondary benchmarks, in priority order

1. **Extra evidence on PBMC10k:**
   - KnockTF perturbations (LINGER, Dictys, TRIPOD and Epiregulon use them).
   - eQTLGen / GTEx eQTL and promoter-capture Hi-C (Javierre et al.) for the
     peak→gene layer.
2. **K562 multiome** with ENCODE ChIP labels and Replogle 2022 CRISPRi as a
   causal held-out test. This is the second pole in the benchmarks: BEAR-GRN,
   scE2G, and the SCENIC+ cell lines. It also overlaps with our existing
   SC-MO-GRN-DB K562 work.
3. **geneRNIB** (Open Problems `task_grn_inference`; bioRxiv, not yet
   peer-reviewed): a label-free, perturbation-based evaluation on OPSCA and
   NeurIPS 2021 BMMC.
   - There is no leakage risk.
   - The bar is **GRNBoost2 and Pearson**: GRNBoost2 ranked first on 10 of its
     11 datasets.
   - Adoption is still small (3 citations), so it is a bonus, not the main
     claim.
4. **SHARE-seq mouse skin** with Cistrome skin ChIP (Dictys protocol), for a
   second species or tissue.
5. **SC-MO-GRN-DB curriculum results.** These are the "cherry on the cake"
   (see below).

### BEAR-GRN matters for the SC-MO-GRN-DB results
**BEAR-GRN** (Karamveer…Uzun, Nat Commun, 18 Sep 2026,
10.1038/s41467-026-77838-w) comes from the same lab that built
SC-MO-GRN-DB.

- **What it benchmarks:** multiome GRN methods on K562, macrophage, iPSC and
  mouse embryo, using SC-MO-GRN-DB ground truths (ChIP-seq, KnockTF, STRING;
  "four complementary definitions").
- **Findings:** LINGER and DIRECT-NET came out on top. Ground-truth choice
  changes rankings a lot. Methods are "predominantly RNA-driven".
- **What it means for us:** it is now the published point of comparison for
  our SC-MO-GRN-DB numbers, so we should match its protocol when we present
  them.
- **Caveat:** only the abstract was verified. The exact ground truths and
  accessions still need checking, possibly from its Zenodo record.

---

## 5. Risks for MeVD-GRN's central claim

- **ATAC adds little under realistic evaluation.** geneRNIB, BEAR-GRN and
  PEREGGRN all find that motif/ATAC methods do no better than, or worse than,
  expression-only ones. ATAC clearly helps only in simulation (scMultiSim).
  MeVD-GRN's "dual modality helps" claim therefore needs a clean ablation on
  PBMC10k (with ATAC vs RNA only, same splits). The result may be null.
- **Motif circularity.** DoRothEA C–E and UniBind are partly motif-derived.
  Scoring an ATAC+motif model against them inflates results. MeVD-GRN's motif
  relation must be off whenever the labels are motif-derived.
- **Binding is not regulation.** ChIP labels measure binding. Perturbation
  data (KnockTF, CRISPRi) is the stronger causal check.
- **Ground-truth choice changes method rankings.** Choosing the reference
  network reverses 32% of pairwise rankings (Kendiukhov, arXiv 2026). Report
  more than one ground truth rather than picking the most favourable one.

---

## 6. Open items and uncertainties

- The exact 10x PBMC variant is unverified for TRIPOD, PRISM-GRN, SCARlink,
  scE2G and SCRIPT. regX's is inferred.
- LINGER's PBMC cell count is reported as 9,543 or 11,909 depending on the
  source. Use LINGER's released processed data (as KEGNI did) rather than
  re-deriving it.
- The LINGER average AUPR ratio (about 2.25) is derived, not printed.
- For BEAR-GRN and the Nat Rev Genet 2023 review, only the abstracts were
  verified. The Dictys details come from bioRxiv v1.
- The PubMed and bioRxiv claude.ai connectors were not authorised, so the
  surveys used the Europe PMC and NCBI APIs instead.
