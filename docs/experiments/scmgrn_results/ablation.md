| Variant | Cell types | test AUPR | diff vs H | p | test AUROC | diff vs H | p |
|---|---|---|---|---|---|---|---|
| H: headline (FM + 384-wide GNN, undirected, training graph, scMultiomeGRN inputs) | 12 | 0.9349 | – | – | 0.9237 | – | – |
| only the scMultiomeGRN node features (no pair histograms) | 12 | 0.9347 | -0.0002 | 0.99 | 0.9237 | -0.0001 | 0.92 |
| only the pair histograms (no scMultiomeGRN node features) | 12 | 0.9347 | -0.0002 | 0.39 | 0.9237 | -0.0001 | 0.61 |
| A1: no scMultiomeGRN inputs | 12 | 0.9350 | +0.0001 | 0.33 | 0.9236 | -0.0001 | 0.94 |
| A2: no scMultiomeGRN inputs, no training graph | 12 | 0.9345 | -0.0004 | 0.19 | 0.9228 | -0.0009 | 0.0016 |
| A3: old directed model, scores averaged | 12 | 0.9341 | -0.0008 | 0.43 | 0.9226 | -0.0011 | 0.072 |
| no Geneformer gene embedding | 12 | 0.9343 | -0.0006 | 0.3 | 0.9229 | -0.0008 | 0.031 |
| no ATAC channel | 12 | 0.9345 | -0.0003 | 0.4 | 0.9235 | -0.0002 | 0.27 |
| no graph network (features straight to the decoder) | 12 | 0.9340 | -0.0009 | 0.0049 | 0.9226 | -0.0011 | 0.021 |
| co-expression graph only (also drops the training links) | 12 | 0.9340 | -0.0009 | 0.066 | 0.9223 | -0.0014 | 0.00089 |
| TF-candidate graph only (also drops the training links) | 12 | 0.9339 | -0.0009 | 0.037 | 0.9225 | -0.0012 | 0.01 |
| base size (128-wide, no Geneformer) | 12 | 0.9348 | -0.0000 | 0.71 | 0.9237 | -0.0001 | 0.92 |
| floor: 128-wide, no Geneformer, no ATAC channel, no graph network, no training graph, no scMultiomeGRN inputs | 12 | 0.6326 | -0.3023 | 2e-21 | 0.6349 | -0.2888 | 2e-21 |
| floor + only the ATAC channel | 12 | 0.8355 | -0.0993 | 3.2e-21 | 0.8232 | -0.1006 | 2e-21 |
| floor + only the Geneformer embedding | 12 | 0.9318 | -0.0030 | 1.3e-06 | 0.9208 | -0.0029 | 6.6e-05 |
| floor + only the graph network (label-free graphs) | 12 | 0.8754 | -0.0595 | 3.2e-20 | 0.8616 | -0.0621 | 4.2e-21 |
| floor + only the scMultiomeGRN node features | 12 | 0.9297 | -0.0052 | 8.3e-11 | 0.9189 | -0.0048 | 2.2e-09 |
| scMultiomeGRN, official code (reference) | 12 | 0.9180 | -0.0169 | 1.9e-15 | 0.9124 | -0.0114 | 7.2e-13 |
| trivial baseline: product of the two TFs' training-link counts (no model, no features) | 12 | 0.9239 | -0.0110 | 4.9e-15 | 0.9122 | -0.0115 | 3.1e-14 |
| trivial baseline: logistic regression on the two TFs' training-link counts | 12 | 0.9291 | -0.0058 | 9.7e-09 | 0.9172 | -0.0066 | 6.9e-09 |
