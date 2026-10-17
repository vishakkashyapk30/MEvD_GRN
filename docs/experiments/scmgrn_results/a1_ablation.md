| Variant | Cell types | test AUPR | diff vs H | p | test AUROC | diff vs H | p |
|---|---|---|---|---|---|---|---|
| H: headline (FM + 384-wide GNN, undirected, training graph, scMultiomeGRN inputs) | 12 | 0.9349 | -0.0001 | 0.33 | 0.9237 | +0.0001 | 0.94 |
| only the scMultiomeGRN node features (no pair histograms) | 12 | 0.9347 | -0.0003 | 0.3 | 0.9237 | +0.0001 | 0.41 |
| only the pair histograms (no scMultiomeGRN node features) | 12 | 0.9347 | -0.0003 | 0.81 | 0.9237 | +0.0001 | 0.38 |
| A1: no scMultiomeGRN inputs | 12 | 0.9350 | – | – | 0.9236 | – | – |
| A2: no scMultiomeGRN inputs, no training graph | 12 | 0.9345 | -0.0005 | 0.23 | 0.9228 | -0.0008 | 0.053 |
| A3: old directed model, scores averaged | 12 | 0.9341 | -0.0009 | 0.22 | 0.9226 | -0.0010 | 0.14 |
| no Geneformer gene embedding | 12 | 0.9343 | -0.0007 | 0.092 | 0.9229 | -0.0007 | 0.18 |
| no ATAC channel | 12 | 0.9345 | -0.0004 | 0.24 | 0.9235 | -0.0001 | 0.83 |
| no graph network (features straight to the decoder) | 12 | 0.9340 | -0.0010 | 0.021 | 0.9226 | -0.0010 | 0.033 |
| co-expression graph only (also drops the training links) | 12 | 0.9340 | -0.0010 | 0.1 | 0.9223 | -0.0013 | 0.015 |
| TF-candidate graph only (also drops the training links) | 12 | 0.9339 | -0.0010 | 0.031 | 0.9225 | -0.0011 | 0.025 |
| base size (128-wide, no Geneformer) | 12 | 0.9348 | -0.0002 | 0.16 | 0.9237 | +0.0001 | 0.9 |
| floor: 128-wide, no Geneformer, no ATAC channel, no graph network, no training graph, no scMultiomeGRN inputs | 12 | 0.6326 | -0.3024 | 2e-21 | 0.6349 | -0.2887 | 2e-21 |
| floor + only the ATAC channel | 12 | 0.8355 | -0.0994 | 2e-21 | 0.8232 | -0.1004 | 2e-21 |
| floor + only the Geneformer embedding | 12 | 0.9318 | -0.0032 | 1.1e-06 | 0.9208 | -0.0028 | 6.3e-05 |
| floor + only the graph network (label-free graphs) | 12 | 0.8754 | -0.0596 | 3.7e-20 | 0.8616 | -0.0620 | 3.4e-21 |
| floor + only the scMultiomeGRN node features | 12 | 0.9297 | -0.0053 | 4.3e-10 | 0.9189 | -0.0047 | 1.8e-08 |
| scMultiomeGRN, official code (reference) | 12 | 0.9180 | -0.0170 | 5.4e-16 | 0.9124 | -0.0112 | 2.5e-13 |
| trivial baseline: product of the two TFs' training-link counts (no model, no features) | 12 | 0.9239 | -0.0111 | 1.2e-14 | 0.9122 | -0.0114 | 3.8e-14 |
| trivial baseline: logistic regression on the two TFs' training-link counts | 12 | 0.9291 | -0.0059 | 6.1e-09 | 0.9172 | -0.0065 | 1.4e-08 |
