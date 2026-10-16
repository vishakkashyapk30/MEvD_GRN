#!/usr/bin/env python
"""Builds docs/presentation/presentation_11th_oct_17th_oct.md (bullet-point slides).

Tables that compare models mark the best value in bold and the second best in italics, separately for each metric.
Ties share the mark. Numbers are the final compiled results (docs/experiments/scmgrn_results/).
"""
from pathlib import Path

OUT = Path(__file__).with_name("presentation_11th_oct_17th_oct.md")


def mark(vals, nd):
    """vals: list of floats (or None). Returns formatted strings: best bold, second best italic (ties share)."""
    shown = [None if v is None else f"{v:.{nd}f}" for v in vals]
    uniq = sorted({float(s) for s in shown if s is not None}, reverse=True)
    best = uniq[0] if uniq else None
    second = uniq[1] if len(uniq) > 1 else None
    out = []
    for s in shown:
        if s is None:
            out.append("–")
        elif float(s) == best:
            out.append(f"**{s}**")
        elif float(s) == second:
            out.append(f"*{s}*")
        else:
            out.append(s)
    return out


def pair_table(header, rows, nd=3, extra=None):
    """Columns are methods, each cell AUROC / AUPR; marks are applied across a row, per metric.
    extra: optional unranked last column (list of strings, one per row)."""
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    for i, (name, cells) in enumerate(rows):
        ro = mark([c[0] for c in cells], nd)
        pr = mark([c[1] for c in cells], nd)
        row = [name] + [f"{a} / {b}" for a, b in zip(ro, pr)]
        if extra is not None:
            row.append(extra[i])
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def column_table(header, rows, rank_cols, nd):
    """Rows are methods; the columns listed in rank_cols are marked down the column."""
    cols = list(zip(*[r for r in rows]))
    marked = {}
    for c in rank_cols:
        marked[c] = mark([None if v in (None, "–") else float(v) for v in cols[c]], nd)
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    for i, r in enumerate(rows):
        cells = [marked[c][i] if c in marked else str(r[c]) for c in range(len(r))]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


# ----------------------------------------------------------------------------- data
M5 = ["Cell type", "MeVD-GRN, full model", "Without scMultiomeGRN inputs", "Without scMultiomeGRN inputs or known links", "Old directed model", "scMultiomeGRN (official code)"]
LUNG = [
    ("Bronchiolar and alveolar", [(0.914, 0.930), (0.913, 0.930), (0.913, 0.929), (0.915, 0.931), (0.901, 0.911)]),
    ("Ciliated epithelial", [(0.928, 0.940), (0.928, 0.940), (0.927, 0.939), (0.928, 0.940), (0.912, 0.915)]),
    ("Vascular endothelial", [(0.923, 0.934), (0.923, 0.935), (0.921, 0.932), (0.923, 0.934), (0.894, 0.902)]),
    ("Stromal", [(0.929, 0.941), (0.929, 0.941), (0.929, 0.941), (0.928, 0.940), (0.920, 0.923)]),
    ("Lymphatic endothelial", [(0.922, 0.937), (0.921, 0.935), (0.920, 0.934), (0.919, 0.933), (0.908, 0.917)]),
    ("Megakaryocytes", [(0.905, 0.923), (0.907, 0.925), (0.902, 0.923), (0.901, 0.920), (0.898, 0.915)]),
    ("Myeloid", [(0.927, 0.933), (0.928, 0.936), (0.928, 0.937), (0.930, 0.939), (0.930, 0.943)]),
    ("Neuroendocrine", [(0.927, 0.938), (0.928, 0.939), (0.929, 0.940), (0.928, 0.939), (0.912, 0.917)]),
    ("Lymphoid", [(0.926, 0.935), (0.922, 0.932), (0.927, 0.937), (0.921, 0.926), (0.917, 0.918)]),
    ("Average of the nine", [(0.922, 0.935), (0.922, 0.935), (0.922, 0.935), (0.921, 0.934), (0.910, 0.918)]),
]
HEMA = [
    ("HSC", [(0.934, 0.940), (0.934, 0.941), (0.931, 0.937), (0.932, 0.940), (0.922, 0.921)]),
    ("CMP", [(0.925, 0.932), (0.925, 0.932), (0.924, 0.931), (0.924, 0.932), (0.914, 0.909)]),
    ("GMP", [(0.925, 0.934), (0.925, 0.935), (0.922, 0.934), (0.924, 0.935), (0.921, 0.926)]),
    ("Average of the three", [(0.928, 0.936), (0.928, 0.936), (0.926, 0.934), (0.927, 0.936), (0.919, 0.919)]),
]
PAPERPROT = [  # full model, old directed model, official ; paper (unranked)
    ("Bronchiolar and alveolar", [(0.940, 0.849), (0.938, 0.855), (0.913, 0.760)], "0.91 / 0.76"),
    ("Ciliated epithelial", [(0.938, 0.828), (0.939, 0.833), (0.918, 0.730)], "0.93 / 0.75"),
    ("Vascular endothelial", [(0.935, 0.800), (0.935, 0.802), (0.899, 0.668)], "0.90 / 0.67"),
    ("Stromal", [(0.952, 0.886), (0.953, 0.890), (0.927, 0.800)], "0.93 / 0.80"),
    ("Lymphatic endothelial", [(0.939, 0.850), (0.942, 0.862), (0.915, 0.781)], "0.93 / 0.81"),
    ("Megakaryocytes", [(0.938, 0.846), (0.937, 0.855), (0.917, 0.788)], "0.92 / 0.78"),
    ("Myeloid", [(0.956, 0.877), (0.958, 0.892), (0.943, 0.849)], "0.95 / 0.85"),
    ("Neuroendocrine", [(0.951, 0.893), (0.952, 0.900), (0.920, 0.820)], "0.92 / 0.82"),
    ("Lymphoid", [(0.951, 0.902), (0.950, 0.898), (0.929, 0.861)], "0.93 / 0.86"),
    ("Average of the nine", [(0.945, 0.859), (0.945, 0.865), (0.920, 0.784)], "0.924 / 0.79"),
]
ABL = [  # label, AUPR, diff, p, AUROC, diff, p
    ("MeVD-GRN, full model", 0.9349, "–", "–", 0.9237, "–", "–"),
    ("Only the scMultiomeGRN node features (no pair histograms)", 0.9347, "-0.0002", "0.99", 0.9237, "-0.0001", "0.92"),
    ("Only the pair histograms (no scMultiomeGRN node features)", 0.9347, "-0.0002", "0.39", 0.9237, "-0.0001", "0.61"),
    ("Without scMultiomeGRN inputs", 0.9350, "+0.0001", "0.33", 0.9236, "-0.0001", "0.94"),
    ("Without scMultiomeGRN inputs or known links", 0.9345, "-0.0004", "0.19", 0.9228, "-0.0009", "0.002"),
    ("Old directed model, both directions averaged", 0.9341, "-0.0008", "0.43", 0.9226, "-0.0011", "0.07"),
    ("No Geneformer gene embedding", 0.9343, "-0.0006", "0.30", 0.9229, "-0.0008", "0.03"),
    ("No ATAC channel", 0.9345, "-0.0003", "0.40", 0.9235, "-0.0002", "0.27"),
    ("No graph network (features go straight to the scorer)", 0.9340, "-0.0009", "0.005", 0.9226, "-0.0011", "0.02"),
    ("Co-expression graph only (this also drops the known links)", 0.9340, "-0.0009", "0.07", 0.9223, "-0.0014", "0.0009"),
    ("TF-candidate graph only (this also drops the known links)", 0.9339, "-0.0009", "0.04", 0.9225, "-0.0012", "0.01"),
    ("Base size: 128 wide and no Geneformer (0.53 M parameters instead of 5 M)", 0.9348, "-0.0000", "0.71", 0.9237, "-0.0001", "0.92"),
    ("Floor: base size, plus no ATAC channel, no graph network, no known links and no scMultiomeGRN inputs (0.10 M parameters)", 0.6326, "-0.3023", "2e-21", 0.6349, "-0.2888", "2e-21"),
    ("scMultiomeGRN, official code (for reference)", 0.9180, "-0.0169", "2e-15", 0.9124, "-0.0114", "7e-13"),
]
ROUTE = [
    ("Nothing added (the floor itself)", 0.633, "-0.302", 0.635, "-0.289"),
    ("The ATAC channel", 0.836, "-0.099", 0.823, "-0.101"),
    ("The graph network over label-free graphs", 0.875, "-0.060", 0.862, "-0.062"),
    ("The scMultiomeGRN node features", 0.930, "-0.005", 0.919, "-0.005"),
    ("The Geneformer embedding", 0.932, "-0.003", 0.921, "-0.003"),
]
BASE = [
    ("scMultiomeGRN, official code", 0.918, 0.912),
    ("No-model baseline: product of the two TFs' training-link counts", 0.924, 0.912),
    ("No-model baseline: logistic regression on the two counts", 0.929, 0.917),
    ("MeVD-GRN, full model", 0.935, 0.924),
]

# ----------------------------------------------------------------------------- slides
S = []
add = S.append
add("""# MeVD-GRN against scMultiomeGRN, on scMultiomeGRN's own benchmark

Advisor meeting, Saturday 17 October 2026, 11:00. Work from 11 to 17 October.

In the result tables, **bold** marks the best score in a row and *italics* the second best, for AUROC and AUPR separately.

## 1. The question and the short answer

- MeVD-GRN was built to predict which transcription factor (TF) regulates which gene. scMultiomeGRN, published in Nucleic Acids Research, tackles a related task: it predicts which pairs of TFs are linked in a given cell type.
- This week we asked a simple question. If both models are trained and tested on scMultiomeGRN's own datasets, splits and scoring code, which one does better?
- To make that a fair contest we rebuilt the benchmark from raw data, checked that the authors' own code reproduces their published numbers on our rebuild, and changed MeVD-GRN so that it solves the same task with the same inputs.
- On the held-out test links MeVD-GRN comes out ahead, but only by a little. Its average AUPR over the nine lung cell types is 0.935, against 0.918 for scMultiomeGRN.
- MeVD-GRN wins on eight of the nine lung cell types and on all three blood cell types. The difference is statistically significant on six of the nine lung cell types, and Myeloid is the only one where scMultiomeGRN is slightly better.
- Three later findings change how this result should be read, and we think they matter more than the headline number.
- First, no single part of MeVD-GRN explains its score. Taking any one part out moves AUPR by 0.001 or less, although taking everything optional out at once drops it to 0.63.
- Second, a rule with no model at all, which just counts how many known links each TF already has, reaches 0.929. That beats scMultiomeGRN and covers about two thirds of MeVD-GRN's lead.
- Third, the benchmark's answer key is computed from the same ATAC data the models read. Put together, this suggests the benchmark mostly rewards knowing which TFs tend to have many links, and that MeVD-GRN's real edge beyond that is about 0.006 AUPR.

![Where every method lands on held-out AUPR](figs_scmgrn/fig_ladder.png)

## 2. The task both models solve

- scMultiomeGRN looks at one cell type at a time and builds a graph whose nodes are TFs. The job is to predict which pairs of TFs are connected.
- The links have no direction. A link between TF A and TF B does not say which one regulates the other.
- The answer key does not come from experiments. Following an earlier method, DeepTFni, A and B are linked when a binding motif of A (HOCOMOCO v11, FIMO p at most 1e-6) sits inside the promoter of B, within 2 kb of its start site, in DNA that is open in more than 10% of the cells.
- We rebuilt this answer key with the same tools and the same thresholds. In practice, the task is to recover a motif-in-promoter rule among TFs from single-cell data.

## 3. The datasets

- The main dataset is human fetal lung (Domcke et al., Science 2020; GEO GSM4508936), with 72,662 cells across nine cell types.
- One detail matters here. According to its GEO description, the "RNA" part of this object is gene-body accessibility rather than RNA expression, so both modalities really come from ATAC.
- The last two columns compare the size of our rebuilt graphs with the paper's own Table S1. They agree to within 5% on every cell type.

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

- The second dataset is human adult hematopoiesis (Buenrostro et al., Cell 2018), with three blood cell types.
- Here the ATAC and RNA data come from separate experiments, exactly as in the scMultiomeGRN paper. We took the ATAC counts from GEO (GSE96769) and the RNA from the paper's supplementary data.
- For CMP we merged three RNA labels (cmp0, cmp1 and cmp2), because the paper does not say which ones it used.

| Cell type | scATAC cells | scRNA cells | TFs (nodes) | Links (undirected) |
|---|---|---|---|---|
| HSC | 347 | 2,268 | 150 | 1,519 |
| CMP | 502 | 4,454 | 171 | 2,016 |
| GMP | 402 | 1,096 | 157 | 1,556 |

## 4. How scMultiomeGRN works

- Each TF is a node. For every TF the model computes an accessibility score from the ATAC data and an importance profile over all genes from the RNA data.
- For every pair of TFs it also computes a small histogram of how their expression varies together.
- Two graph neural networks, one for each kind of data, pass information between neighbouring TFs using attention. A cross-modal attention layer then merges the two views into one representation per TF.
- A pair is scored by the inner product of the two TF representations, and the model learns to rebuild the known links.
- The paper's final network is a majority vote over ten runs. We ran the authors' code exactly as they released it.

## 5. How MeVD-GRN works (the version used this week)

- MeVD-GRN starts by summarising each TF in a handful of numbers taken from the RNA and ATAC data.
- It adds a pretrained gene embedding from Geneformer. This gives every gene a fixed vector that is the same in every cell type.
- It then builds two graphs between TFs without looking at any labels. One connects TFs that are co-expressed, and the other connects each TF to the partners it is most likely to relate to, judged from accessibility and expression together.
- Two graph networks run over these graphs, one for the RNA view and one for the ATAC view, and their outputs are merged into a single embedding per TF.
- A pair scorer turns the embeddings of two TFs into the chance that they are linked. The scorer gives the same answer for (u, v) as for (v, u), which is what an undirected task needs.
- When we switch them on, the model can also use the same extra inputs as scMultiomeGRN, including the known training links.

![MeVD-GRN, version used for the scMultiomeGRN comparison](../figures/architecture_undirected_2026-10-17.png)

## 6. What we changed in MeVD-GRN for a fair comparison

- The original MeVD-GRN predicts a directed edge from a TF to a gene, and its scorer treats the two ends differently. That is the wrong shape for an undirected answer key.
- We replaced it with a scorer that is symmetric by design. It combines the two TF embeddings only through terms that ignore order, and a unit test checks that both orders give identical scores.
- We also merged the RNA and ATAC graph networks into one embedding per TF. Before, the two views were kept apart until the very last step.
- We switched off MeVD-GRN's own motif relation. The answer key is itself a motif scan, so using it would be circular.
- Training uses the same budget as scMultiomeGRN, with at most 2,000 epochs and a stop after 100 epochs without improvement. The epoch we keep is picked on a validation split, never on the test links.
- The rest of the model is our standard configuration, which was fixed before any result came in: hidden size 384, two graph layers and the Geneformer embedding.
- To give MeVD-GRN the same information scMultiomeGRN has, we added scMultiomeGRN's inputs as optional extras. None of them involves a held-out link.

| What scMultiomeGRN uses | How MeVD-GRN receives it |
|---|---|
| MAESTRO accessibility score of each TF | compressed to 64 numbers per TF and added through a small encoder |
| GRNBoost2 importance of each TF towards every gene | compressed to 64 numbers per TF, through the same encoder |
| 16 by 16 joint-expression histogram of each TF pair | made symmetric (136 values) and added as one term in the scorer |
| The known training links | passed to the graph networks as an extra graph |

## 7. The four versions of MeVD-GRN we compared

- To see which of these changes matter, we ran four versions of the model on every cell type.
- The full model has everything switched on.
- The version without scMultiomeGRN inputs drops those extra input features.
- The version without scMultiomeGRN inputs or known links also drops the known training links, so it never sees a single known link.
- The old directed model is the original MeVD-GRN, with its two directions averaged.

| Version | Symmetric scorer | scMultiomeGRN input features | Known links as a graph | Parameters |
|---|---|---|---|---|
| MeVD-GRN, full model | yes | yes | yes | 5.07 M |
| Without scMultiomeGRN inputs | yes | no | yes | 4.88 M |
| Without scMultiomeGRN inputs or known links | yes | no | no | 3.70 M |
| Old directed model | no | no | yes | 4.44 M |

- All four score within 0.001 AUPR of each other on average over the twelve cell types, and within about 0.01 on any single cell type. Slide 13 explains what that tells us.

## 8. How the comparison was set up

- Both methods get the same data, the same splits and the same scoring code.
- The splits are scMultiomeGRN's own: ten folds made with its random seed, 666. Our fold files are byte-for-byte identical to the ones the authors' split code writes.
- In each fold, a tenth of the known links plus the same number of random non-links is held out for testing. Another 18% of the links is kept aside for validation.
- Scores come from the authors' own metric function. AUROC is the chance that a true link is ranked above a non-link, and AUPR is the area under the precision-recall curve.
- Our headline numbers are on the held-out test links only. The paper scores the full matrix of TF pairs, which includes the training links, so we report that separately on slide 15 to compare with the published figures.
- We ran the official code with six seeds on the first five lung cell types and two seeds on the rest. MeVD-GRN ran with five seeds everywhere.

## 9. How we kept test data out of training

- The test links and non-links never appear in the training data, and a check in the code enforces this on every fold.
- Negative training examples are drawn only from each fold's training negatives.
- The graphs MeVD-GRN reads are built without any labels.
- Every MeVD-GRN result file records that it ran in this zero-leak mode, and the script that compiles the results refuses any file that did not.
- No setting was tuned on test results. All the versions were fixed before the first result came in.

## 10. Did we rebuild the benchmark correctly?

- The authors do not release their processed data, so the comparison only means something if our rebuild matches theirs. We checked this three ways.
- Our link counts are within 5% of the paper's Table S1 on all nine lung cell types (slide 3). The largest gap is Vascular endothelial, at 4.6% above.
- Our fold files are identical to the official split code's output, in all 30 files for each of the twelve cell types.
- Most convincingly, the official scMultiomeGRN code run unchanged on our data gives back the paper's own scores on the paper's own metric.

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

- Averaged over the nine lung cell types, the re-run scores 0.920 / 0.784, against 0.924 / 0.790 in the paper. We are comparing against a faithful copy of the published method, not a weakened one.
""")

add("""## 11. Main result: lung, held-out test links

- The table shows AUROC / AUPR on the held-out links, averaged over the ten folds and then over the seeds.
""")
add(pair_table(M5, LUNG))
add("""
![Held-out test links, lung](figs_scmgrn/fig_test.png)

- To check whether these gaps are real, we compared the full model with the official code fold by fold, since both used the same folds, with a Wilcoxon signed-rank test.

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
| Pooled over all 90 fold pairs | +0.017 | 5e-11 | +0.012 | 2e-9 |

- MeVD-GRN is ahead on AUPR in 83% of the 90 cell-type and fold pairs, and on AUROC in 78%.
- The gains are small, mostly one to three points of AUPR, while scores move by two to five points from one fold to the next.
- With only ten folds per cell type, 0.002 is the smallest p-value this test can produce.
- The pooled p-value treats the folds of one cell type as independent, even though they share a graph, so it is too optimistic.
- The fair summary is that MeVD-GRN is consistently a little better, not dramatically better.

## 12. Main result: hematopoiesis

- The blood data show the same pattern on a second, independent dataset.
""")
add(pair_table(M5, HEMA))
add("""
![Held-out test links, hematopoiesis](figs_scmgrn_hema/fig_test.png)

- Fold by fold against the official code, the full model gains 0.020 AUPR on HSC (p = 0.002), 0.023 on CMP (p = 0.002) and 0.008 on GMP (p = 0.049).
- Pooled over the 30 fold pairs, the gain is 0.017 AUPR (p = 6e-8) and 0.009 AUROC (p = 8e-7), and MeVD-GRN is ahead in 90% of the pairs.
- The four MeVD-GRN versions are again almost impossible to tell apart.
- Keep in mind that the RNA and ATAC data were measured in different cells here, and that the CMP RNA merges three labels.

## 13. Why does MeVD-GRN score higher? The ablations

- If MeVD-GRN's lead came from one design choice, removing that choice should lower its score. We tested this by changing one part of the full model at a time and running every version on all twelve cell types with five seeds.
- The table gives the average held-out score over the twelve cell types, the difference from the full model, and a fold-paired test against the full model over 120 pairs.
""")
add(column_table(["Version", "AUPR", "Difference from the full model", "p", "AUROC", "Difference from the full model", "p"], ABL, rank_cols=[1, 4], nd=4))
add("""
- No single change matters much. Removing Geneformer, the ATAC channel, the graph networks, either graph, the known links or the scMultiomeGRN inputs, or shrinking the model tenfold, moves AUPR by 0.001 or less.
- A few of these tiny shifts are statistically detectable (no graph network gives p = 0.005), but each is about a twentieth of the gap to scMultiomeGRN.
- That gap, 0.017 AUPR, is about twenty times bigger than any of them, so we cannot pin it on one design decision.
- The floor tells a different story. With everything optional removed at once, the score collapses to 0.633, about 0.30 below the full model.
- So the useful information is not missing, it is duplicated. Several routes carry it, any one of them is enough, and the model fails only when all of them are gone.
- To see which routes work on their own, we added exactly one back to the floor model.
""")
add(column_table(["Floor model plus only...", "AUPR", "Difference from the full model", "AUROC", "Difference from the full model"], ROUTE, rank_cols=[1, 3], nd=3))
add("""
- The Geneformer embedding is a fixed vector per gene and is the same in every cell type. It knows nothing about the cell data, only which TFs are in the pair.
- Even so, a model that only knows which two TFs it is looking at scores 0.932 AUPR. That is 0.003 below the full model and 0.014 above scMultiomeGRN.
- The scMultiomeGRN node features, which also act mostly as an identity tag for each TF, do almost as well.
- This points to the identity of the TFs, rather than anything specific to the cell type, as the main source of signal. The next slide tests that directly.

![Where every method lands on held-out AUPR](figs_scmgrn/fig_ladder.png)

## 14. A no-model baseline: counting known links

- The official folds split links, not TFs, so every TF shows up in the training links.
- A TF with many known links in training is likely to have many among the test links too. A method can exploit this hub effect without learning anything about regulation.
- To measure how much this explains, we built a baseline with no model, no features and no cell data.
- For a test pair of TFs, the baseline counts how many training links each of the two TFs has. It then scores the pair by the product of the two counts, or with a logistic regression on them.
- It only uses each fold's training links, and it is scored with the same code as every other method.
""")
add(column_table(["Method (held-out test links, 12 cell types)", "AUPR", "AUROC"], BASE, rank_cols=[1, 2], nd=3))
add("""
- The counting baseline beats scMultiomeGRN on AUPR and matches or beats it on AUROC. On this split, scMultiomeGRN does no better than a rule that counts known links.
- Of the 0.017 AUPR that MeVD-GRN gains over scMultiomeGRN, about 0.011 can be had from counting alone.
- What MeVD-GRN adds beyond the hub effect is about 0.006 AUPR and 0.007 AUROC. That is small but real: a fold-paired test over 120 pairs gives p = 1e-8 for AUPR and 7e-9 for AUROC.
- This matches what we saw with BEAR-GRN, where simple in-degree rankers matched or beat published methods and reviewers expected splits that hold out whole TFs.
- A stricter version of this benchmark would hold out all the links of some TFs, which removes the counting shortcut. We have not run it yet, and it is the obvious next experiment.

## 15. The paper's own evaluation (for reference only)

- The paper scores the full matrix of all TF pairs, including the links each model was trained on. We show it only so the numbers can be compared with the published ones, and we draw no conclusions from it.
- The last column holds the published figures for scMultiomeGRN and is not ranked.
""")
add(pair_table(["Cell type", "MeVD-GRN, full model", "Old directed model", "scMultiomeGRN (official code)", "Paper (published, not ranked)"],
               [(r[0], r[1]) for r in PAPERPROT], extra=[r[2] for r in PAPERPROT]))
add("""
![The paper's evaluation](figs_scmgrn/fig_paper.png)

- On this evaluation MeVD-GRN's lead looks much bigger, about 0.075 AUPR pooled over the 90 fold pairs, against 0.017 on the held-out links.
- The difference comes from the training links. A larger model rebuilds the links it was trained on more completely, and this evaluation rewards that, so it should not be read as better generalisation.
- For the blood cell types, on the same evaluation: HSC scores 0.945 / 0.790 for MeVD-GRN, 0.927 / 0.696 for the official code and 0.92 / 0.67 in the paper. CMP scores 0.945 / 0.800, 0.927 / 0.692 and 0.92 / 0.66, and GMP scores 0.943 / 0.784, 0.928 / 0.681 and 0.93 / 0.71.

## 16. Limitations

- The answer key is computed from the same ATAC data the models read. A model can score well by learning how the key was built, without learning anything about regulation, and the same holds for scMultiomeGRN.
- The folds split links, not TFs. A TF's identity and its number of training links carry over to the test links, which is why counting links works so well.
- The graphs are small. The lung cell types have 41 to 271 TFs and 266 to 5,026 links, so each test fold holds only 27 to 503 true links.
- In the lung data, the "RNA" is gene-body accessibility rather than expression, so both modalities come from ATAC.
- The Seurat object has no list of variable genes, so the GRNBoost2 step uses genes seen in more than 10% of cells plus all TFs, following the paper's own filter. This applies to both methods.
- The official code ran with six seeds on five cell types and two on the rest. It gives nearly the same result for every seed (six seeds on Myeloid agree to within 0.001), and it takes about five minutes per fold on the large cell types.
- For the blood data, the ATAC counts come from the GEO matrix rather than per-cell sequence files, and the CMP RNA is our reading of three labels.
- None of this tells us how the two methods compare on experimentally measured answer keys, which is what a benchmark like BEAR-GRN offers.

## 17. What we conclude and what comes next

- On scMultiomeGRN's own benchmark, with the same data, splits and scoring, and a faithful copy of the published method as the baseline, MeVD-GRN scores consistently higher on held-out links: 0.935 against 0.918 AUPR. This is a fair claim and we can make it.
- What matters more is what that lead is made of. Counting known links reproduces about two thirds of it, no single part of MeVD-GRN is responsible for it, and a fixed per-gene embedding with no cell data comes close to the full model.
- What MeVD-GRN adds beyond the hub effect is about 0.6 points of AUPR.
- Our reading is that this benchmark mostly measures how well a method learns which TFs tend to have links. MeVD-GRN does that slightly better than scMultiomeGRN, but the benchmark says little about whether either method infers regulation.
- The next step is to repeat the comparison with a split by TF, holding out every link of some TFs. That removes the counting shortcut and forces both methods to use the cell data.
- Running that on scMultiomeGRN's datasets would show whether MeVD-GRN's small lead survives, and how scMultiomeGRN fares once the hub effect is gone. The same idea applies to benchmarks with experimental answer keys.

## 18. Run record and where things are

- Everything in this document is finished: all nine lung and three blood cell types for both methods, the thirteen ablation versions on all twelve cell types, and the no-model baselines.
- Ada's outbound internet turned out to manage only 10 to 140 KB/s per connection, so we streamed the big input files from the laptop over the VPN, which took about three and a half hours.
- Along the way we found that the GEO lung file is compressed twice, that the Seurat object has no variable-gene list, and that the official code runs out of memory on the largest cell types unless each worker gets 30 GB. We also fixed a bug that made one of our workers stop early. All of this is written up in `docs/experiments/scmultiomegrn_generalization.md`, section 11.8.
- The code is on branch `scmgrn-tftf`. The model lives in `src/models/mevd_grn.py` and `src/models/decoder.py`, the benchmark features in `src/benchmarks/scmgrn_features.py`, the scripts in `scripts/14b`, `16` to `21`, the cluster jobs in `slurm/scmgrn_*.sh` and the tests in `tests/`.
- The compiled results are in `docs/experiments/scmgrn_results/`. The fold-level files are on Ada under `~/scmgrn/results/`.
- The figures are made by `docs/presentation/make_scmgrn_figs.py`, `docs/presentation/make_scmgrn_ladder_fig.py` and `docs/figures/make_architecture_undirected.py`, and this document by `docs/presentation/build_deck_11_17.py`.
""")

OUT.write_text("\n".join(S))
print("written", OUT)
