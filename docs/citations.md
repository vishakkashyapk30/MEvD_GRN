# MEvD-GRN design choices and citations

Gene Regulatory Network Inference from Integrative Multi-Omics Data  
By Vishak Kashyap K, UG4 CND  
Advisor: Dr Vinod PK  
Updated: 5 Sep 2026

This note lists every major design choice in MEvD-GRN, explains it in plain English, and points to the GRN or ML paper that inspired it.


## 1. Dataset and evidence tiers

**What we do.** We use SC-MO-GRN-DB. For each cell type we get paired single-cell RNA and ATAC data plus reference regulatory edges. Those edges are labeled by how they were discovered: localization (TF physically seen near the gene), perturbation (removing the TF changes the gene), or dual-evidence (both).

**Why.** We need ground truth that is experimental, not another algorithm's guesses, and we need RNA and ATAC from the same cell-type setting. The three tiers also give us a natural easy-to-hard training order.

**Citation.** Valensi et al., SC-MO-GRN-DB, iScience 2026.  
https://doi.org/10.1016/j.isci.2026.115323  
https://scmogrndb.psu.edu


## 2. Turning ATAC peaks into gene-level scores

**What we do.** ATAC data is peaks by cells (open DNA intervals), not genes. We map peaks that fall within about 100 kb of a gene's transcription start site, sum that signal into a gene activity score, then keep per-gene mean, variance, and an openness number.

**Why.** The model talks about genes and TF to gene edges. Peaks are not genes, so we need a simple way to say "how open is this gene's neighborhood?" That is the usual gene-activity idea in scATAC analysis.

**Citations.**  
Signac gene-activity workflow: https://stuartlab.org/signac/articles/pbmc_vignette  
Duren et al., modeling regulation from paired expression and accessibility, PNAS 2018: https://doi.org/10.1073/pnas.1802973115


## 3. RNA node features and co-expression signatures

**What we do.** From the RNA matrix we do not feed raw cells into the GNN. For each gene we keep mean expression, variance, and detection rate across cells. We also build a short co-expression signature vector so the dot product of two genes approximates how correlated they are.

**Why.** The model treats genes as nodes. A few stable per-gene numbers are enough for an embedding, and co-expression has long been used as a cue that two genes may share regulation.

**Citations.**  
GRNBoost2 / tree-ensemble GRN context: Moerman et al., Bioinformatics 2019, https://doi.org/10.1093/bioinformatics/bty916  
GENIE3 (co-expression / tree-based GRN inference): Huynh-Thu et al., PLoS ONE 2010, https://doi.org/10.1371/journal.pone.0012776


## 4. Task setup: directed TF to gene link prediction

**What we do.** We score candidate directed edges "does TF i regulate gene j?" rather than clustering cells or reconstructing undirected co-expression graphs alone.

**Why.** A GRN is a directed map from regulators to targets. Framing it as link prediction lets us use standard supervised GNN training and metrics.

**Citations.**  
Kipf and Welling, variational graph auto-encoders / link prediction, https://arxiv.org/abs/1611.07308  
Asymmetric bilinear scoring (knowledge-graph style): Yang et al., DistMult, ICLR 2015, https://arxiv.org/abs/1412.6575


## 5. Separate RNA and ATAC encoders

**What we do.** RNA features and ATAC features each go through their own small neural network. We do not glue the two modalities into one vector at the start.

**Why.** Expression and chromatin accessibility mean different things. Keeping them separate until later matches how recent multi-omic GRN models handle modalities, and it lets accessibility act as a gate rather than just another feature channel.

**Citations.**  
Xu et al., scMultiomeGRN, Nucleic Acids Research 2025: https://doi.org/10.1093/nar/gkaf138  
Sankar et al., GraFRank (multi-modal GNN pattern), WWW 2021: https://doi.org/10.1145/3442381.3450120


## 6. Role-aware decoder (ATAC gates the target)

**What we do.** At scoring time, the target's ATAC embedding becomes a gate on the target's RNA embedding. The TF side uses RNA. We also add explicit co-expression and openness terms for the pair.

**Why.** Biologically, a TF can only usefully regulate a gene if that locus is accessible. Mixing RNA and ATAC as equals throws that away. The gate says "accessibility first, then expression compatibility."

**Citations.**  
Duren et al., PNAS 2018: https://doi.org/10.1073/pnas.1802973115  
Co-expression as a regulatory cue: GENIE3 / GRNBoost2 papers above


## 7. GraphSAGE message passing

**What we do.** After encoding, each gene talks to its neighbors with a 2-layer GraphSAGE. The updated embedding mixes the gene's own features with its neighborhood.

**Why.** A gene's regulatory context is not only its own mean expression. Neighbors that co-vary or look like plausible TF targets carry useful structure. GraphSAGE is a standard inductive GNN for that.

**Citations.**  
Hamilton, Ying, Leskovec, GraphSAGE, NeurIPS 2017: https://proceedings.neurips.cc/paper/2017/hash/5dd9db5e033da9c6fb5ba83c7a7ebea9-Abstract.html  
GMFGRN (GNN / matrix-factorization GRN inference): Li et al., Briefings in Bioinformatics 2024: https://doi.org/10.1093/bib/bbad529


## 8. Two structural prior graphs

**What we do.** Message passing uses two graphs at once: (1) an undirected co-expression kNN graph among genes, and (2) a directed TF to candidate-target graph built from accessibility and co-expression. Known labeled positives are kept out of that second graph so labels cannot leak into embeddings through messages.

**Why.** Different biological relations should not be forced into one edge type. Relational GNNs sum updates over several relation types. Co-expression captures modules; the TF-candidate graph captures plausible regulation direction without copying the training labels.

**Citations.**  
Schlichtkrull et al., R-GCN, https://arxiv.org/abs/1703.06103  
Co-expression gene graphs in GRN GNNs: GMFGRN, https://doi.org/10.1093/bib/bbad529  
Accessibility-aware candidates: Duren / Signac (sections 2 above)  
Link-prediction hygiene (no label leakage into structure): see also OGB practices, https://ogb.stanford.edu


## 9. Curriculum over evidence tiers

**What we do.** We train first on localization, then fine-tune on perturbation. Dual-evidence is held out and used only as a zero-shot test.

**Why.** Curriculum learning says start with easier or broader signal, then specialize. Localization is large and noisy; perturbation is smaller and more functional; dual is the strictest exam set.

**Citations.**  
Bengio et al., Curriculum Learning, ICML 2009: https://dl.acm.org/doi/10.1145/1553374.1553380  
Tier definitions: SC-MO-GRN-DB, https://doi.org/10.1016/j.isci.2026.115323


## 10. Memory replay against catastrophic forgetting

**What we do.** While training on a later tier, we keep replaying a small fraction of earlier-tier positive edges at lower loss weight.

**Why.** Sequential fine-tuning often overwrites what the model learned first (catastrophic forgetting). Experience replay is a standard continual-learning fix. We measured that forgetting is real on localization after Stage 2, so replay is on by default.

**Citations.**  
McCloskey and Cohen, 1989, catastrophic interference (Psychology of Learning and Motivation, Vol. 24)  
Rolnick et al., Experience Replay for Continual Learning, NeurIPS 2019: https://proceedings.neurips.cc/paper/2019/hash/fa7cdfad1a5aaf8370ebeda47a1ff1c3-Abstract.html  
Related (not yet used): Kirkpatrick et al., EWC, PNAS 2017: https://doi.org/10.1073/pnas.1611835114


## 11. Global edge splits and random negatives

**What we do.** Every unique TF-gene pair is assigned to train, validation, or test once across all tiers. Negatives are random non-positive TF-gene pairs, not restricted to the accessibility prior.

**Why.** The tiers nest heavily. If you split each tier independently, a "test" edge in dual can still be a training edge in localization. A global split blocks that leak. Random negatives avoid a trivial "high expression equals positive" shortcut that prior-restricted negatives can create.

**Citations.**  
Link-prediction split protocol: Kipf and Welling VGAE, https://arxiv.org/abs/1611.07308  
Nesting motivation: SC-MO-GRN-DB, https://doi.org/10.1016/j.isci.2026.115323


## 12. Hard negatives across tiers

**What we do.** Half of the training negatives can be hard cases, for example pairs that appear in a weaker tier but not in the current stronger tier (bind but may not regulate).

**Why.** Random negatives are often too easy. Hard negatives force the model to separate "looks related" from "actually regulates," which matches the localization versus perturbation distinction.

**Citations.**  
Hard-negative mining: Schroff et al., FaceNet, CVPR 2015: https://doi.org/10.1109/CVPR.2015.7298682  
Biological reading of the tiers: SC-MO-GRN-DB, https://doi.org/10.1016/j.isci.2026.115323


## 13. Metrics

**What we do.** We report AUPR as the main metric, plus AUROC, early precision, and EPR.

**Why.** True regulatory edges are rare compared with non-edges. AUPR cares more about ranking true edges near the top under that imbalance. AUROC alone can look optimistic.

**Citation.**  
Saito and Rehmsmeier, precision-recall for imbalanced data, PLoS ONE 2015: https://doi.org/10.1371/journal.pone.0118432  
Also common in GRNBoost2 / GENIE3-style evaluations.


## 14. Baselines we compare against

**What we do.** We benchmark against a strong multi-omic GNN and several RNA-only methods on the same splits.

| Method | Paper | Link |
|---|---|---|
| scMultiomeGRN | Xu et al., NAR 2025 | https://doi.org/10.1093/nar/gkaf138 |
| scMultiomeGRN code | Zenodo | https://doi.org/10.5281/zenodo.14848389 |
| GRNBoost2 | Moerman et al., Bioinformatics 2019 | https://doi.org/10.1093/bioinformatics/bty916 |
| RegDiffusion | Zhu and Slonim, J. Comp. Biol. 2024 | https://doi.org/10.1089/cmb.2024.0607 |
| GMFGRN / GMF-GAE | Li et al., Brief. Bioinform. 2024 | https://doi.org/10.1093/bib/bbad529 |


## 15. Edge embeddings (new)

**What we do.** Instead of only adding linear weights on co-expression and openness, we can map those two numbers through a small MLP into an edge embedding, then add that to the TF-target score (`use_edge_mlp`).

**Why.** Message-passing literature treats edge features as first-class inputs, not just scalars. scMultiomeGRN already uses rich edge attributes. This is our lighter version of that idea inside the role-aware decoder.

**Citations.**  
Gilmer et al., Neural Message Passing (MPNN), ICML 2017: https://arxiv.org/abs/1704.01212  
Xu et al., scMultiomeGRN: https://doi.org/10.1093/nar/gkaf138  
Sankar et al., GraFRank: https://doi.org/10.1145/3442381.3450120  
Related: Simonovsky and Komodakis, edge-conditioned convolutions, https://arxiv.org/abs/1704.02901


## 16. Ablating the two prior graphs (new)

**What we do.** We train variants that message-pass on co-expression only, TF-candidate only, or both, plus the older no-GNN control.

**Why.** We already know message passing helps overall. We still need to know whether both priors are necessary or whether one dominates. That is a clean relational-GNN ablation.

**Citations.**  
Multi-relation motivation: R-GCN, https://arxiv.org/abs/1703.06103  
GraphSAGE necessity (no_gnn control): https://proceedings.neurips.cc/paper/2017/hash/5dd9db5e033da9c6fb5ba83c7a7ebea9-Abstract.html
