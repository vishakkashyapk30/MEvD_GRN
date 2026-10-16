# MeVD-GRN against scMultiomeGRN, on scMultiomeGRN's own benchmark

Advisor meeting, Saturday 17 October 2026, 11:00. Covers the work of 11 to 17 October.

This document is written to be read. It explains what we compared, why the comparison is fair, what we found, and, just as important, what the findings do and do not mean. The headline result is positive but small, and a simple baseline that we added at the end changes how it should be read, so the later sections matter more than the summary alone.

## 1. The question and the short answer

Our model, MeVD-GRN, was built to predict which transcription factor (TF) regulates which gene. scMultiomeGRN, a recent paper in Nucleic Acids Research, solves a related but different problem: given single-cell RNA and ATAC data from a cell type, predict which pairs of TFs are linked to each other. This week we asked a simple question: if we take scMultiomeGRN's datasets, its data splits and its scoring code, and ask both models to solve its task, which one does better?

To answer it fairly we had to do three things. We rebuilt scMultiomeGRN's benchmark from raw data and checked that its own code, run unchanged on our rebuild, reproduces the numbers in its paper. We changed MeVD-GRN so that it solves the same task (undirected TF-to-TF links, where it used to predict directed TF-to-gene edges) and so that it can use the same inputs. And we ran both methods on all nine lung cell types and three hematopoiesis cell types, with ten official folds each.

The answer, on the held-out test links, is that MeVD-GRN scores higher than scMultiomeGRN, but only slightly. Averaged over the nine lung cell types, the area under the precision-recall curve (AUPR) is 0.935 against 0.918, a gain of 0.017. MeVD-GRN is ahead on eight of the nine lung cell types and on all three hematopoiesis cell types, and the gain is statistically significant on six of the nine lung cell types. Myeloid is the one cell type where scMultiomeGRN is slightly ahead.

Three further findings qualify that result, and we think they are the most useful part of the week. First, no single component of MeVD-GRN accounts for its score: removing any one of them (the pretrained gene embedding, the ATAC channel, the graph networks, the known links, the scMultiomeGRN inputs) changes AUPR by 0.001 or less, although removing all of them at once drops it to 0.63. Second, a baseline with no model at all, which scores a pair of TFs only by how many known links each already has, reaches 0.929 AUPR, beating scMultiomeGRN and recovering about two thirds of MeVD-GRN's advantage. Third, the benchmark's answer key is constructed from the same ATAC data that the models read. Taken together, these suggest the benchmark mostly measures how well a method can learn which TFs tend to have many links, and that MeVD-GRN's real advantage over that is about 0.006 AUPR.

![Where every method lands on held-out AUPR](figs_scmgrn/fig_ladder.png)

## 2. The benchmark: what task are both models solving?

scMultiomeGRN works on one cell type at a time. For each cell type it builds a graph whose nodes are transcription factors, and the task is to predict which pairs of TFs are connected. The graph is undirected: a link between TF A and TF B does not say which regulates which.

The answer key is not experimental. Following an earlier method called DeepTFni, a link between two TFs is defined by sequence and chromatin accessibility: TF A is linked to TF B if a binding-motif match for A (from the HOCOMOCO v11 motif library, with a FIMO significance of at most 1e-6) lies inside the promoter of B (within 2 kb of its transcription start site), in a stretch of DNA that is open in more than 10% of the cells of that cell type. We regenerated this key with the same tools and thresholds. The task is therefore, in effect, to recover a motif-in-promoter rule among TFs from single-cell data.

The paper uses two datasets. The first is human fetal lung (Domcke et al., Science 2020; GEO GSM4508936), with 72,662 cells in nine cell types. One detail matters for interpretation: the "RNA" part of this object is, according to its GEO description, gene-body accessibility rather than RNA expression, so both modalities are really derived from ATAC.

| Cell type | Cells | TFs (nodes) | Links (undirected) | Directed entries, ours / paper | Link density |
|---|---|---|---|---|---|
| Myeloid | 1,288 | 46 | 266 | 532 / 542 | 0.26 |
| Lymphoid | 1,860 | 41 | 301 | 602 / 613 | 0.37 |
| Megakaryocytes | 568 | 64 | 459 | 918 / 930 | 0.23 |
| Neuroendocrine | 544 | 69 | 674 | 1,348 / 1,361 | 0.29 |
| Lymphatic endothelial | 662 | 94 | 936 | 1,872 / 1,889 | 0.21 |
| Stromal | 43,289 | 94 | 1,044 | 2,088 / 2,104 | 0.24 |
| Bronchiolar and alveolar epithelial | 18,988 | 115 | 1,346 | 2,692 / 2,712 | 0.21 |
| Ciliated epithelial | 1,111 | 228 | 4,254 | 8,508 / 8,503 | 0.16 |
| Vascular endothelial | 4,352 | 271 | 5,026 | 10,052 / 9,609 | 0.14 |

The second dataset is human adult hematopoiesis (Buenrostro et al., Cell 2018), with three cell types. Here the ATAC data and the RNA data come from separate experiments, as in the scMultiomeGRN paper. We took the ATAC counts from the original GEO deposit (GSE96769) and the RNA from the paper's supplementary data, and for the CMP cell type we combined three RNA labels (cmp0, cmp1, cmp2) because the paper does not say which it used.

| Cell type | scATAC cells | scRNA cells | TFs (nodes) | Links (undirected) |
|---|---|---|---|---|
| HSC | 347 | 2,268 | 150 | 1,519 |
| CMP | 502 | 4,454 | 171 | 2,016 |
| GMP | 402 | 1,096 | 157 | 1,556 |

## 3. The two models

scMultiomeGRN (Xu et al., 2025) treats each TF as a node. For every TF it computes an accessibility score from the ATAC data and an importance profile towards all genes from the RNA data, and for every pair of TFs it computes a small joint-expression histogram. Two graph neural networks, one per modality, pass information between neighbouring TFs using attention, and a cross-modal attention layer then merges the two views into one representation per TF. The score of a pair is the inner product of the two representations, and the model is trained on the known links to reconstruct the graph. For the final network, the paper takes a majority vote over ten runs. We use the authors' code exactly as released.

MeVD-GRN is a different design. It first summarises each TF in a few numbers from the RNA and ATAC data, and adds a frozen, pretrained gene embedding from the Geneformer model, which gives every gene a fixed vector that is the same in every cell type. It then builds graphs between TFs without using any labels (one linking TFs that are co-expressed, and one linking each TF to the partners it is most likely to relate to, judged from chromatin accessibility together with expression), and runs two graph networks over them, one for the RNA view and one for the ATAC view. The two views are fused into a single embedding per TF, and a pair scorer turns the embeddings of two TFs into the probability that they are linked.

The diagram below shows the version used this week. Two things in it are new. The scorer is symmetric, meaning the score for the pair (u, v) is by construction the same as for (v, u), which is what an undirected task needs. And the model can optionally receive the same extra inputs that scMultiomeGRN uses, including the known training links, so that both methods have access to the same information.

![MeVD-GRN, version used for the scMultiomeGRN comparison](../figures/architecture_undirected_2026-10-17.png)

## 4. What we changed in MeVD-GRN to make the comparison fair

The original MeVD-GRN predicts a directed edge from a TF to a gene. It scores a pair using separate roles for the two ends: the target gene's RNA state is gated by how accessible that gene's DNA is, so the score depends on which end is the TF. For scMultiomeGRN's task this is the wrong shape, because the answer key is undirected. In an earlier attempt we simply averaged the two directions, but this week we replaced it with a scorer that is symmetric by design. It combines the two TF embeddings through terms that do not depend on order (their product, the absolute difference, their sum, and a diagonal bilinear term), adds an expression-correlation term and, when available, a pair feature, and is unit-tested to return identical scores for either order. We also fuse the RNA and ATAC graph networks into a single embedding per TF, where before the two views were kept apart until the very last step.

We then gave MeVD-GRN access to the information scMultiomeGRN uses, so that neither method has an informational head start. This is done by switching on optional inputs, none of which touches a held-out link:

| What scMultiomeGRN uses | How MeVD-GRN receives it |
|---|---|
| MAESTRO accessibility score of each TF | compressed to 64 numbers per TF (features only) and added through a small encoder |
| GRNBoost2 importance of each TF towards every gene | compressed to 64 numbers per TF, same encoder |
| 16 by 16 joint-expression histogram of each TF pair | symmetrised into 136 values, one extra term in the pair scorer |
| the known (training) links | fed to the graph networks as an extra graph |

We switched off MeVD-GRN's own motif relation, because the answer key is itself a motif scan and using it would be circular. Training follows the same budget as scMultiomeGRN (at most 2,000 epochs, stopping after 100 without improvement), uses one example per unordered pair, and draws its negative examples only from the fold's training negatives. The epoch to keep is chosen on a validation split, never on the test links. Everything else, including the two graph networks, the Geneformer embedding, the hidden size of 384 and the two layers, is MeVD-GRN's recommended configuration, fixed before any result existed.

To find out which of these changes matter, we ran four versions of MeVD-GRN. The headline variant H has everything. A1 removes the scMultiomeGRN input features. A2 removes those and also the known training links, so it never sees a known link. A3 is the old directed model with its two directions averaged.

| Variant | Symmetric scorer | scMultiomeGRN input features | Known links as a graph | Parameters |
|---|---|---|---|---|
| H (headline) | yes | yes | yes | 5.07 M |
| A1 | yes | no | yes | 4.88 M |
| A2 | yes | no | no | 3.70 M |
| A3 | no (old directed model) | no | yes | 4.44 M |

On the held-out links all four turn out to score within 0.001 AUPR of each other on average over the twelve cell types, and within about 0.01 on any single cell type. We return to what this means in section 9.

## 5. How the comparison was set up

Both methods are given the same data, the same splits and the same scoring. The splits are scMultiomeGRN's own: ten folds, generated with its random seed of 666. In each fold, one tenth of the known links, together with an equal number of randomly chosen non-links, is held out as the test set, and a further 18% of the links is set aside for validation. Our fold files are byte-identical to the output of the authors' own split code. The score is computed with the authors' own metric function: AUROC (the chance that a true link is ranked above a non-link) and AUPR (the area under the precision-recall curve), on the held-out test pairs. The official scMultiomeGRN code was run unmodified on our rebuilt data, with six random seeds on the first five lung cell types and two seeds on the rest; MeVD-GRN was run with five seeds on every cell type.

The scores we report as the headline are on the held-out test links only. This differs from how the paper itself reports results. The paper scores the full matrix of all TF pairs, which includes the links the model was trained on, so we also report that number, clearly labelled, in section 11 for comparison with the published figures.

We took several precautions so that nothing from the test set can reach training. The test links and non-links never appear in the training data, and this is checked in code for every fold. Negative training examples come only from the fold's training negatives. The graphs that MeVD-GRN reads are built without labels. Every MeVD-GRN result file records that it was run in this zero-leak mode, and the compiling step refuses to include any run that was not. No hyper-parameter was tuned on test results: the variants were decided before any result was available.

## 6. Did we rebuild the benchmark correctly?

Since the authors do not release the processed data, a comparison is only meaningful if our rebuild is faithful. We checked this in three ways. First, the number of links we regenerate for each lung cell type is within 5% of the number in the paper's Table S1 (the table in section 2 shows both; the largest gap is Vascular endothelial, 4.6% above). Second, our fold files are identical to those written by the official split code, in all 30 files for each of the twelve cell types. Third, and most convincingly, the official scMultiomeGRN code, run unchanged on our data, reproduces the paper's own scores on its own metric:

| Cell type | Paper, AUROC / AUPR | Our re-run of the official code |
|---|---|---|
| Myeloid | 0.95 / 0.85 | 0.943 / 0.849 |
| Lymphoid | 0.93 / 0.86 | 0.929 / 0.861 |
| Megakaryocytes | 0.92 / 0.78 | 0.917 / 0.788 |
| Neuroendocrine | 0.92 / 0.82 | 0.920 / 0.820 |
| Lymphatic endothelial | 0.93 / 0.81 | 0.915 / 0.781 |
| Bronchiolar and alveolar | 0.91 / 0.76 | 0.913 / 0.760 |
| Ciliated epithelial | 0.93 / 0.75 | 0.918 / 0.730 |
| Vascular endothelial | 0.90 / 0.67 | 0.899 / 0.668 |
| Stromal | 0.93 / 0.80 | 0.927 / 0.800 |
| HSC (hematopoiesis) | 0.92 / 0.67 | 0.927 / 0.696 |
| CMP (hematopoiesis) | 0.92 / 0.66 | 0.927 / 0.692 |
| GMP (hematopoiesis) | 0.93 / 0.71 | 0.928 / 0.681 |

Averaged over the nine lung cell types, the re-run scores 0.920 / 0.784, against the paper's reported average of 0.924 / 0.790. We are therefore comparing against a faithful reproduction of the published method, not a weakened one.

## 7. Main result: lung, held-out test links

The table gives AUROC / AUPR on the held-out test links, averaged over the ten folds and then over the seeds. The best score in each row is in bold.

| Cell type | **MeVD-GRN H** | A1 | A2 | A3 | scMultiomeGRN (official code) |
|---|---|---|---|---|---|
| Bronchiolar and alveolar | **0.914 / 0.930** | 0.913 / 0.930 | 0.913 / 0.929 | 0.915 / 0.931 | 0.901 / 0.911 |
| Ciliated epithelial | **0.928 / 0.940** | 0.928 / 0.940 | 0.927 / 0.939 | 0.928 / 0.940 | 0.912 / 0.915 |
| Vascular endothelial | **0.923 / 0.934** | 0.923 / 0.935 | 0.921 / 0.932 | 0.923 / 0.934 | 0.894 / 0.902 |
| Stromal | **0.929 / 0.941** | 0.929 / 0.941 | 0.929 / 0.941 | 0.928 / 0.940 | 0.920 / 0.923 |
| Lymphatic endothelial | **0.922 / 0.937** | 0.921 / 0.935 | 0.920 / 0.934 | 0.919 / 0.933 | 0.908 / 0.917 |
| Megakaryocytes | **0.905 / 0.923** | 0.907 / 0.925 | 0.902 / 0.923 | 0.901 / 0.920 | 0.898 / 0.915 |
| Myeloid | 0.927 / 0.933 | 0.928 / 0.936 | 0.928 / 0.937 | 0.930 / 0.939 | **0.930 / 0.943** |
| Neuroendocrine | **0.927 / 0.938** | 0.928 / 0.939 | 0.929 / 0.940 | 0.928 / 0.939 | 0.912 / 0.917 |
| Lymphoid | **0.926 / 0.935** | 0.922 / 0.932 | 0.927 / 0.937 | 0.921 / 0.926 | 0.917 / 0.918 |
| **Average of the nine** | **0.922 / 0.935** | 0.922 / 0.935 | 0.922 / 0.935 | 0.921 / 0.934 | 0.910 / 0.918 |

![Held-out test links, lung](figs_scmgrn/fig_test.png)

To judge whether the differences are real, we compared the two methods fold by fold within each cell type, since the same folds were used for both, using a Wilcoxon signed-rank test. The table gives the difference between MeVD-GRN H and the official code.

| Cell type | AUPR difference | p | AUROC difference | p |
|---|---|---|---|---|
| Vascular endothelial | +0.033 | 0.002 | +0.029 | 0.002 |
| Ciliated epithelial | +0.025 | 0.002 | +0.016 | 0.004 |
| Neuroendocrine | +0.021 | 0.002 | +0.015 | 0.004 |
| Lymphatic endothelial | +0.020 | 0.002 | +0.014 | 0.010 |
| Bronchiolar and alveolar | +0.019 | 0.004 | +0.013 | 0.014 |
| Stromal | +0.019 | 0.004 | +0.009 | 0.020 |
| Lymphoid | +0.017 | 0.105 | +0.009 | 0.375 |
| Megakaryocytes | +0.009 | 0.160 | +0.007 | 0.275 |
| Myeloid | -0.011 | 0.375 | -0.003 | 0.695 |
| **Pooled over all 90 fold pairs** | **+0.017** | **5e-11** | **+0.012** | **2e-9** |

MeVD-GRN is ahead on AUPR in 83% of the 90 cell-type-by-fold pairs, and on AUROC in 78%. These numbers should be read with some care. The gains are small, mostly one to three points of AUPR, while the fold-to-fold standard deviation is two to five points. With only ten folds per cell type, 0.002 is the smallest p-value this paired test can give. And the pooled p-value treats the folds of one cell type as independent although they share a graph, so it is optimistic. The fair summary is that MeVD-GRN is consistently a little better, not that it is dramatically better.

## 8. Main result: hematopoiesis

The hematopoiesis results show the same pattern on a second, independent dataset. MeVD-GRN H scores higher than the official code on all three cell types, and the four MeVD-GRN variants are again almost indistinguishable from each other.

| Cell type | **MeVD-GRN H** | A1 | A2 | A3 | scMultiomeGRN (official code) |
|---|---|---|---|---|---|
| HSC | **0.934 / 0.940** | 0.934 / 0.941 | 0.931 / 0.937 | 0.932 / 0.940 | 0.922 / 0.921 |
| CMP | **0.925 / 0.932** | 0.925 / 0.932 | 0.924 / 0.931 | 0.924 / 0.932 | 0.914 / 0.909 |
| GMP | **0.925 / 0.934** | 0.925 / 0.935 | 0.922 / 0.934 | 0.924 / 0.935 | 0.921 / 0.926 |
| Average | **0.928 / 0.936** | 0.928 / 0.936 | 0.926 / 0.934 | 0.927 / 0.936 | 0.919 / 0.919 |

![Held-out test links, hematopoiesis](figs_scmgrn_hema/fig_test.png)

Compared fold by fold with the official code, MeVD-GRN H gains +0.020 AUPR on HSC (p = 0.002), +0.023 on CMP (p = 0.002) and +0.008 on GMP (p = 0.049). Pooled over the 30 fold pairs the gain is +0.017 AUPR (p = 6e-8) and +0.009 AUROC (p = 8e-7), and MeVD-GRN is ahead in 90% of the pairs. Two features of the hematopoiesis data deserve a mention: the RNA and ATAC data come from separate experiments, so the two modalities are not measured in the same cells, and, as noted in section 2, the CMP RNA combines three labels because the paper is not explicit about which it used.

## 9. Why does MeVD-GRN score higher? Ablations

If MeVD-GRN's advantage came from one of its design choices, removing that choice should lower its score. We tested this with a set of ablations, each changing one component of the headline variant, and ran every one on all twelve cell types with five seeds. The table gives the mean held-out AUPR and AUROC over the twelve cell types, the difference from the headline variant H, and a fold-paired test against H over 120 cell-type-by-fold pairs.

| Variant | AUPR | diff vs H | p | AUROC | diff vs H | p |
|---|---|---|---|---|---|---|
| **H**, the headline model | 0.9349 | – | – | 0.9237 | – | – |
| only the scMultiomeGRN node features (no pair histograms) | 0.9347 | -0.0002 | 0.99 | 0.9237 | -0.0001 | 0.92 |
| only the pair histograms (no scMultiomeGRN node features) | 0.9347 | -0.0002 | 0.39 | 0.9237 | -0.0001 | 0.61 |
| A1: no scMultiomeGRN inputs | 0.9350 | +0.0001 | 0.33 | 0.9236 | -0.0001 | 0.94 |
| A2: A1 and no known training links | 0.9345 | -0.0004 | 0.19 | 0.9228 | -0.0009 | 0.002 |
| A3: old directed model, directions averaged | 0.9341 | -0.0008 | 0.43 | 0.9226 | -0.0011 | 0.07 |
| no Geneformer gene embedding | 0.9343 | -0.0006 | 0.30 | 0.9229 | -0.0008 | 0.03 |
| no ATAC channel | 0.9345 | -0.0003 | 0.40 | 0.9235 | -0.0002 | 0.27 |
| no graph network (features go straight to the scorer) | 0.9340 | -0.0009 | 0.005 | 0.9226 | -0.0011 | 0.02 |
| co-expression graph only (also drops the known links) | 0.9340 | -0.0009 | 0.07 | 0.9223 | -0.0014 | 0.0009 |
| TF-candidate graph only (also drops the known links) | 0.9339 | -0.0009 | 0.04 | 0.9225 | -0.0012 | 0.01 |
| base size: 128 wide, no Geneformer (0.53 M parameters instead of 5 M) | 0.9348 | -0.0000 | 0.71 | 0.9237 | -0.0001 | 0.92 |
| **floor**: base size, and also no ATAC channel, no graph network, no known links, no scMultiomeGRN inputs (0.10 M parameters) | 0.6326 | -0.3023 | 2e-21 | 0.6349 | -0.2888 | 2e-21 |
| *reference:* scMultiomeGRN, official code | 0.9180 | -0.0169 | 2e-15 | 0.9124 | -0.0114 | 7e-13 |

The first thing the table shows is how little any single change matters. Removing the pretrained embedding, the ATAC channel, the graph networks, either graph, the known links, or the scMultiomeGRN inputs, or shrinking the model by a factor of ten, changes the held-out score by 0.001 AUPR or less. A few of these tiny differences are statistically detectable (for instance, no graph network, p = 0.005), but they are about a twentieth of the gap to scMultiomeGRN and of no practical size. The gap to scMultiomeGRN, 0.017 AUPR, is an order of magnitude larger than any of them, so it is not something we can attribute to one design decision.

The second thing is the floor. When we remove everything optional at once, the model collapses to 0.633 AUPR, about 0.30 below the headline. So the information the model needs is not absent but redundant: several different routes carry it, any one of them is enough, and the model fails only when all of them are gone. To find out which routes are enough on their own, we added exactly one back to the floor model:

| Floor model plus only ... | AUPR | diff vs H | AUROC | diff vs H |
|---|---|---|---|---|
| nothing (the floor itself) | 0.633 | -0.302 | 0.635 | -0.289 |
| the ATAC channel | 0.836 | -0.099 | 0.823 | -0.101 |
| the graph network over label-free graphs | 0.875 | -0.060 | 0.862 | -0.062 |
| the scMultiomeGRN node features | 0.930 | -0.005 | 0.919 | -0.005 |
| **the Geneformer embedding** | **0.932** | **-0.003** | **0.921** | **-0.003** |

The result is striking. The Geneformer embedding is a fixed vector for each gene, identical in every cell type, so it carries no information about the cell data at all, only about which TFs are involved. A model that knows nothing but which two TFs form a pair scores 0.932 AUPR, within 0.003 of the full model and 0.014 above scMultiomeGRN. The scMultiomeGRN node features, which also act largely as an identity signature of each TF, do almost as well. This points to the identity of the TFs, rather than anything specific to the cell type, as the main source of signal, and the next section tests that idea directly.

![Where every method lands on held-out AUPR](figs_scmgrn/fig_ladder.png)

## 10. A trivial baseline: counting known links

The official folds split links, not TFs. Every TF therefore appears in the training links, and a TF that already has many known links is likely to have many among the held-out links too. This is a hub effect, and it would let a method score well without learning anything about regulation. We tested how much of the score it explains with a baseline that has no model, no features and no cell data. For a held-out pair of TFs (u, v), it simply counts how many training links u has and how many v has, using only the fold's training links, and scores the pair by the product of the two counts, or by a logistic regression on them. It is scored with exactly the same code as every other method.

| Method, held-out test links, 12 cell types | AUPR | AUROC |
|---|---|---|
| scMultiomeGRN, official code | 0.918 | 0.912 |
| trivial baseline: product of the two TFs' training-link counts | 0.924 | 0.912 |
| trivial baseline: logistic regression on the two counts | 0.929 | 0.917 |
| MeVD-GRN, headline variant H | 0.935 | 0.924 |

The counting baseline beats scMultiomeGRN on AUPR and matches or exceeds it on AUROC. On this split, then, scMultiomeGRN does not do better than a rule that only counts known links. For MeVD-GRN, of the 0.017 AUPR it gains over scMultiomeGRN, about 0.011 can be obtained by the counting baseline alone. What remains beyond hubness is about 0.006 AUPR and 0.007 AUROC. That remainder is small, but it is real: a fold-paired test over 120 pairs gives p = 1e-8 for AUPR and 7e-9 for AUROC.

This is in line with what we saw in the BEAR-GRN literature earlier, where simple in-degree rankers matched or beat published methods, and where reviewers expect splits that hold out whole TFs. The natural stricter version of this benchmark would hold out all the links of some TFs, which would remove the counting shortcut. We have not run it; it is the obvious next experiment.

## 11. The same benchmark under the paper's own evaluation

For comparison with the published figures, we also report the paper's own evaluation, which scores the full matrix of all TF pairs. This matrix includes the links each model was trained on, so it is not a leakage-free measure and we do not use it for any conclusion. It is useful for two things: the official code's numbers on it match the paper (shown in section 6), and it shows why the paper's evaluation flatters models that fit the training links well. AUROC / AUPR:

| Cell type | MeVD-GRN H | A3 | scMultiomeGRN (official code) | Paper |
|---|---|---|---|---|
| Bronchiolar and alveolar | 0.940 / 0.849 | 0.938 / 0.855 | 0.913 / 0.760 | 0.91 / 0.76 |
| Ciliated epithelial | 0.938 / 0.828 | 0.939 / 0.833 | 0.918 / 0.730 | 0.93 / 0.75 |
| Vascular endothelial | 0.935 / 0.800 | 0.935 / 0.802 | 0.899 / 0.668 | 0.90 / 0.67 |
| Stromal | 0.952 / 0.886 | 0.953 / 0.890 | 0.927 / 0.800 | 0.93 / 0.80 |
| Lymphatic endothelial | 0.939 / 0.850 | 0.942 / 0.862 | 0.915 / 0.781 | 0.93 / 0.81 |
| Megakaryocytes | 0.938 / 0.846 | 0.937 / 0.855 | 0.917 / 0.788 | 0.92 / 0.78 |
| Myeloid | 0.956 / 0.877 | 0.958 / 0.892 | 0.943 / 0.849 | 0.95 / 0.85 |
| Neuroendocrine | 0.951 / 0.893 | 0.952 / 0.900 | 0.920 / 0.820 | 0.92 / 0.82 |
| Lymphoid | 0.951 / 0.902 | 0.950 / 0.898 | 0.929 / 0.861 | 0.93 / 0.86 |
| **Average of the nine** | 0.945 / 0.859 | 0.945 / 0.865 | 0.920 / 0.784 | 0.924 / 0.79 |

![The paper's evaluation](figs_scmgrn/fig_paper.png)

On this evaluation MeVD-GRN's advantage looks much larger (about +0.075 AUPR pooled over the 90 fold pairs) than on the held-out links (+0.017). The difference is the training links: a higher-capacity model reconstructs the links it was trained on more completely, and the paper's evaluation counts that. It should not be read as better generalisation.

Hematopoiesis, same evaluation, AUROC / AUPR: HSC 0.945 / 0.790 for MeVD-GRN, 0.927 / 0.696 for the official code and 0.92 / 0.67 in the paper; CMP 0.945 / 0.800, 0.927 / 0.692 and 0.92 / 0.66; GMP 0.943 / 0.784, 0.928 / 0.681 and 0.93 / 0.71.

## 12. Limitations

The answer key is built from the same ATAC data that the models read. A link exists when a motif hit sits in an open promoter, and the models are given the accessibility of the TFs involved. A model can therefore score well by learning that construction rule without learning anything about regulation, and the same is true of scMultiomeGRN. This benchmark cannot distinguish a method that learns regulation from one that learns how the answer key is made, and everything in sections 9 and 10 is consistent with the second.

The folds split links, not TFs, so a TF's identity and its number of training links carry over to the test links. This is why a counting baseline matches the models, and a split that holds out whole TFs would be a stricter test.

The graphs are small. The lung cell types have between 41 and 271 TFs and between 266 and 5,026 links, so each held-out fold contains only 27 to 503 positive links. Differences of a point or two of AUPR are within the noise of a single fold, which is why we rely on paired comparisons over many folds, and why the pooled p-values, which treat the folds of one cell type as independent, are optimistic.

In the lung data the "RNA" modality is gene-body accessibility, not expression, so both modalities come from ATAC.

A few details differ from the official procedure. The Seurat object stores no variable-gene list, so the GRNBoost2 step uses genes detected in more than 10% of cells (plus all TFs), following the paper's own filter; this applies to both methods. The official code was run with six seeds on five cell types and two on the rest, because it is almost deterministic across seeds (on Myeloid, six seeds agree to within 0.001) and takes about five minutes per fold on the large cell types. For hematopoiesis, the ATAC counts come from the GEO matrix rather than per-cell sequence files, and the CMP RNA is our reading of three labels.

Finally, none of this shows how the methods compare on experimentally derived answer keys, which is what a benchmark such as BEAR-GRN offers.

## 13. What we conclude and what comes next

On scMultiomeGRN's own benchmark, with the same data, the same splits, the same scoring and a faithful reproduction of the published method as the baseline, MeVD-GRN scores consistently higher on held-out links: 0.935 against 0.918 AUPR, ahead on 8 of 9 lung and 3 of 3 hematopoiesis cell types. That is a fair claim and we can make it.

The more important finding is what the margin consists of. A rule that counts known links reproduces roughly two thirds of it, no single component of MeVD-GRN is responsible for it, and a fixed per-gene embedding with no cell information nearly matches the full model. What MeVD-GRN adds beyond the hub effect is about 0.6 points of AUPR. We think the honest description is that this benchmark mostly measures how well a method can learn which TFs tend to have links, that MeVD-GRN does so slightly better than scMultiomeGRN, and that the benchmark is a weak test of whether either method infers regulation.

The natural next step is a stricter version of the same comparison: split by TF, so that all links of some TFs are held out, which removes the counting shortcut and forces a method to use the cell data. Doing this on scMultiomeGRN's datasets would show whether MeVD-GRN's small advantage survives, and whether scMultiomeGRN does better or worse once hubness is taken away. The same idea applies to benchmarks with experimental answer keys, which we set aside this week.

## 14. Run record and where things are

Everything reported here is finished: all nine lung and three hematopoiesis cell types for both methods, the thirteen ablation variants on all twelve cell types, and the trivial baselines. Several problems came up while running on the real data and were fixed along the way. Ada's outbound internet turned out to be only about 10 to 140 KB/s per connection, so the large input files were streamed from the laptop over the VPN, which took about three and a half hours. The GEO lung file is compressed twice, the Seurat object has no variable-gene list, the official code ran out of memory on the largest cell types until the workers were given 30 GB, and a bug in our job queue made one worker stop early. All are recorded in `docs/experiments/scmultiomegrn_generalization.md` (section 11.8).

- Code: branch `scmgrn-tftf`. The model is in `src/models/mevd_grn.py` and `src/models/decoder.py`, the benchmark features in `src/benchmarks/scmgrn_features.py`, the scripts in `scripts/14b`, `16`, `17`, `18`, `19`, `20` and `21`, the cluster jobs in `slurm/scmgrn_*.sh`, and the tests in `tests/`.
- Compiled results: `docs/experiments/scmgrn_results/` (lung, hematopoiesis and the ablation table). Fold-level results are on Ada in `~/scmgrn/results/{lung,hema}/<cell type>/<method>/seed*/fold*.json`.
- Figures: `docs/presentation/make_scmgrn_figs.py`, `make_scmgrn_ladder_fig.py` and `docs/figures/make_architecture_undirected.py`.
