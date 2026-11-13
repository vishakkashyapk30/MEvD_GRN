# Head-to-head: mevd_fm_h384_hubmotif (s12.9 rule)

| dataset | GT | metric | MeVD-GRN (sd, n) | LINGER | best other | in-degree | coverage | random | WIN |
|---|---|---|---|---|---|---|---|---|---|
| K562 | CellTypeExclusive | AUROC | 0.647 (0.011, 3) | 0.525 | Pando-GLM 0.570 | 0.626 | 0.500 | 0.500 | **yes** |
| K562 | CellTypeExclusive | AUPRC | 0.358 (0.003, 3) | 0.348 | TRIPOD 0.327 | 0.345 | 0.305 | 0.253 | **yes** |
| K562 | ChIP | AUROC | 0.603 (0.009, 3) | 0.536 | Pando-GLM 0.576 | 0.652 | 0.500 | 0.500 | no |
| K562 | ChIP | AUPRC | 0.473 (0.003, 3) | 0.430 | TRIPOD 0.405 | 0.484 | 0.431 | 0.331 | no |
| K562 | Core | AUROC | 0.532 (0.051, 3) | 0.507 | Pando-GLM 0.519 | 0.577 | 0.500 | 0.500 | no |
| K562 | Core | AUPRC | 0.136 (0.007, 3) | 0.134 | TRIPOD 0.138 | 0.148 | 0.134 | 0.104 | no |
| K562 | Intersection | AUROC | 0.308 (0.051, 3) | 0.355 | SCENIC+ 0.578 | 0.550 | 0.500 | 0.500 | no |
| K562 | Intersection | AUPRC | 0.088 (0.001, 3) | 0.099 | GRaNIE 0.120 | 0.100 | 0.094 | 0.120 | no |
| K562 | KO | AUROC | 0.555 (0.006, 3) | 0.401 | SCENIC+ 0.580 | 0.569 | 0.500 | 0.500 | no |
| K562 | KO | AUPRC | 0.125 (0.001, 3) | 0.139 | SCENIC+ 0.159 | 0.134 | 0.122 | 0.160 | no |
| K562 | Union | AUROC | 0.647 (0.011, 3) | 0.525 | Pando-GLM 0.570 | 0.626 | 0.500 | 0.500 | **yes** |
| K562 | Union | AUPRC | 0.358 (0.003, 3) | 0.344 | TRIPOD 0.322 | 0.343 | 0.303 | 0.244 | **yes** |
| Macrophage_S1 | CellTypeExclusive | AUROC | 0.625 (0.005, 3) | 0.596 | SCENIC+ 0.593 | 0.609 | 0.500 | 0.500 | **yes** |
| Macrophage_S1 | CellTypeExclusive | AUPRC | 0.368 (0.008, 3) | 0.316 | TRIPOD 0.326 | 0.361 | 0.328 | 0.265 | **yes** |
| Macrophage_S1 | ChIP | AUROC | 0.625 (0.005, 3) | 0.596 | SCENIC+ 0.599 | 0.609 | 0.500 | 0.500 | **yes** |
| Macrophage_S1 | ChIP | AUPRC | 0.373 (0.008, 3) | 0.321 | TRIPOD 0.332 | 0.365 | 0.331 | 0.272 | **yes** |
| Macrophage_S1 | Core | AUROC | 0.566 (0.040, 3) | 0.529 | Pando-XGB 0.572 | 0.545 | 0.500 | 0.500 | **yes** |
| Macrophage_S1 | Core | AUPRC | 0.149 (0.012, 3) | 0.145 | TRIPOD 0.161 | 0.151 | 0.142 | 0.104 | no |
| Macrophage_S2 | CellTypeExclusive | AUROC | 0.588 (0.040, 3) | 0.548 | CellOracle 0.605 | 0.610 | 0.500 | 0.500 | no |
| Macrophage_S2 | CellTypeExclusive | AUPRC | 0.363 (0.019, 3) | 0.319 | TRIPOD 0.330 | 0.369 | 0.338 | 0.265 | no |
| Macrophage_S2 | ChIP | AUROC | 0.588 (0.040, 3) | 0.548 | CellOracle 0.605 | 0.610 | 0.500 | 0.500 | no |
| Macrophage_S2 | ChIP | AUPRC | 0.367 (0.019, 3) | 0.324 | TRIPOD 0.336 | 0.374 | 0.342 | 0.271 | no |
| Macrophage_S2 | Core | AUROC | 0.531 (0.026, 3) | 0.501 | Pando-GLM 0.564 | 0.548 | 0.500 | 0.500 | no |
| Macrophage_S2 | Core | AUPRC | 0.144 (0.003, 3) | 0.147 | TRIPOD 0.160 | 0.156 | 0.146 | 0.104 | no |
| Naive_mESC | ChIP | AUROC | 0.541 (0.002, 2) | 0.476 | SCENIC+ 0.556 | 0.598 | 0.500 | 0.500 | no |
| Naive_mESC | ChIP | AUPRC | 0.285 (0.002, 2) | 0.160 | TRIPOD 0.192 | 0.308 | 0.269 | 0.153 | no |
| Naive_mESC | Intersection | AUROC | 0.575 (0.013, 2) | 0.755 | CellOracle 0.607 | 0.477 | 0.500 | 0.500 | no |
| Naive_mESC | Intersection | AUPRC | 0.113 (0.008, 2) | 0.106 | TRIPOD 0.124 | 0.095 | 0.097 | 0.096 | **yes** |
| Naive_mESC | KO | AUROC | nan | nan | None nan | nan | nan | 0.500 | no |
| Naive_mESC | KO | AUPRC | nan | nan | None nan | nan | nan | nan | no |
| Naive_mESC | Union | AUROC | 0.545 (0.003, 2) | 0.477 | GRaNIE 0.552 | 0.594 | 0.500 | 0.500 | no |
| Naive_mESC | Union | AUPRC | 0.263 (0.002, 2) | 0.142 | TRIPOD 0.171 | 0.283 | 0.247 | 0.134 | no |
| iPS | CellTypeExclusive | AUROC | 0.617 (0.057, 3) | 0.496 | GRaNIE 0.569 | 0.604 | 0.500 | 0.500 | **yes** |
| iPS | CellTypeExclusive | AUPRC | 0.132 (0.007, 3) | 0.109 | TRIPOD 0.116 | 0.129 | 0.122 | 0.101 | **yes** |
| iPS | ChIP | AUROC | 0.617 (0.057, 3) | 0.496 | GRaNIE 0.568 | 0.604 | 0.500 | 0.500 | **yes** |
| iPS | ChIP | AUPRC | 0.156 (0.012, 3) | 0.108 | TRIPOD 0.119 | 0.150 | 0.137 | 0.091 | **yes** |
| iPS | Core | AUROC | 0.636 (0.045, 3) | 0.511 | Pando-GLM 0.606 | 0.543 | 0.500 | 0.500 | **yes** |
| iPS | Core | AUPRC | 0.136 (0.006, 3) | 0.124 | TRIPOD 0.143 | 0.127 | 0.123 | 0.104 | **yes** |
| mESC_E7.5_rep1 | ChIP | AUROC | 0.549 (0.011, 3) | 0.652 | GRaNIE 0.570 | 0.604 | 0.500 | 0.500 | no |
| mESC_E7.5_rep1 | ChIP | AUPRC | 0.296 (0.003, 3) | 0.241 | TRIPOD 0.193 | 0.323 | 0.274 | 0.155 | no |
| mESC_E7.5_rep1 | Intersection | AUROC | 0.504 (0.065, 3) | 0.510 | SCENIC+ 0.567 | 0.485 | 0.500 | 0.500 | no |
| mESC_E7.5_rep1 | Intersection | AUPRC | 0.103 (0.009, 3) | 0.113 | TRIPOD 0.140 | 0.103 | 0.105 | 0.091 | no |
| mESC_E7.5_rep1 | KO | AUROC | 0.513 (0.029, 3) | 0.461 | TRIPOD 0.571 | 0.415 | 0.500 | 0.500 | **yes** |
| mESC_E7.5_rep1 | KO | AUPRC | 0.093 (0.002, 3) | 0.095 | TRIPOD 0.114 | 0.088 | 0.093 | 0.091 | no |
| mESC_E7.5_rep1 | Union | AUROC | 0.548 (0.012, 3) | 0.633 | GRaNIE 0.555 | 0.599 | 0.500 | 0.500 | no |
| mESC_E7.5_rep1 | Union | AUPRC | 0.269 (0.004, 3) | 0.205 | TRIPOD 0.168 | 0.293 | 0.250 | 0.131 | no |
| mESC_E7.5_rep2 | ChIP | AUROC | 0.563 (0.024, 3) | 0.644 | GRaNIE 0.594 | 0.596 | 0.500 | 0.500 | no |
| mESC_E7.5_rep2 | ChIP | AUPRC | 0.294 (0.009, 3) | 0.225 | TRIPOD 0.186 | 0.306 | 0.268 | 0.153 | no |
| mESC_E7.5_rep2 | Intersection | AUROC | 0.587 (0.053, 3) | 0.501 | Pando-GLM 0.547 | 0.483 | 0.500 | 0.500 | **yes** |
| mESC_E7.5_rep2 | Intersection | AUPRC | 0.112 (0.008, 3) | 0.108 | TRIPOD 0.121 | 0.101 | 0.103 | 0.091 | **yes** |
| mESC_E7.5_rep2 | KO | AUROC | 0.549 (0.037, 3) | 0.453 | Pando-GLM 0.541 | 0.415 | 0.500 | 0.500 | **yes** |
| mESC_E7.5_rep2 | KO | AUPRC | 0.098 (0.003, 3) | 0.097 | TRIPOD 0.108 | 0.092 | 0.096 | 0.091 | **yes** |
| mESC_E7.5_rep2 | Union | AUROC | 0.562 (0.025, 3) | 0.628 | GRaNIE 0.585 | 0.592 | 0.500 | 0.500 | no |
| mESC_E7.5_rep2 | Union | AUPRC | 0.269 (0.009, 3) | 0.192 | TRIPOD 0.160 | 0.278 | 0.244 | 0.128 | no |
| mESC_E8.5_rep1 | ChIP | AUROC | 0.544 (0.017, 2) | 0.658 | GRaNIE 0.588 | 0.607 | 0.500 | 0.500 | no |
| mESC_E8.5_rep1 | ChIP | AUPRC | 0.283 (0.009, 2) | 0.227 | CellOracle 0.185 | 0.312 | 0.265 | 0.154 | no |
| mESC_E8.5_rep1 | Intersection | AUROC | 0.541 (0.027, 2) | 0.585 | SCENIC+ 0.618 | 0.483 | 0.500 | 0.500 | no |
| mESC_E8.5_rep1 | Intersection | AUPRC | 0.094 (0.002, 2) | 0.104 | TRIPOD 0.103 | 0.091 | 0.092 | 0.091 | no |
| mESC_E8.5_rep1 | KO | AUROC | 0.512 (0.015, 2) | 0.487 | SCENIC+ 0.602 | 0.413 | 0.500 | 0.500 | **yes** |
| mESC_E8.5_rep1 | KO | AUPRC | 0.094 (0.000, 2) | 0.095 | TRIPOD 0.106 | 0.089 | 0.093 | 0.091 | no |
| mESC_E8.5_rep1 | Union | AUROC | 0.543 (0.017, 2) | 0.643 | GRaNIE 0.593 | 0.602 | 0.500 | 0.500 | no |
| mESC_E8.5_rep1 | Union | AUPRC | 0.259 (0.009, 2) | 0.195 | CellOracle 0.159 | 0.284 | 0.243 | 0.132 | no |
| mESC_E8.5_rep2 | ChIP | AUROC | 0.539 (0.004, 2) | 0.656 | GRaNIE 0.621 | 0.609 | 0.500 | 0.500 | no |
| mESC_E8.5_rep2 | ChIP | AUPRC | 0.277 (0.004, 2) | 0.224 | CellOracle 0.188 | 0.313 | 0.263 | 0.155 | no |
| mESC_E8.5_rep2 | Intersection | AUROC | 0.537 (0.141, 2) | 0.560 | TRIPOD 0.594 | 0.497 | 0.500 | 0.500 | no |
| mESC_E8.5_rep2 | Intersection | AUPRC | 0.113 (0.030, 2) | 0.110 | TRIPOD 0.139 | 0.098 | 0.098 | 0.091 | **yes** |
| mESC_E8.5_rep2 | KO | AUROC | 0.562 (0.049, 2) | 0.475 | CellOracle 0.596 | 0.420 | 0.500 | 0.500 | **yes** |
| mESC_E8.5_rep2 | KO | AUPRC | 0.101 (0.004, 2) | 0.097 | TRIPOD 0.118 | 0.091 | 0.095 | 0.091 | **yes** |
| mESC_E8.5_rep2 | Union | AUROC | 0.544 (0.001, 2) | 0.636 | GRaNIE 0.611 | 0.605 | 0.500 | 0.500 | no |
| mESC_E8.5_rep2 | Union | AUPRC | 0.255 (0.001, 2) | 0.192 | CellOracle 0.161 | 0.286 | 0.242 | 0.132 | no |

WIN in 25/70 cells; beats LINGER in 48/70; beats every released method in 30/70.