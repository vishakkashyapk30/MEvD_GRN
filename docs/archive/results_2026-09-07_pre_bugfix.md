# Results catalog (K562)

All recorded runs of MEvD-GRN, ablations, and baselines on K562.
Numbers are copied from `results/K562_results.json`, `results/ablations/*.json`,
`results/baselines/*.json`, and `results/summary_table.csv`.

Cell type: K562  
Gene universe: 22,943 genes (225 TFs)  
Eval metrics: AUPR (primary), AUROC, early precision (EP), EPR  
Unless noted, dual-evidence is zero-shot (never trained on) using val+test edges.


## How to read this file

Two training budgets appear in this folder:

1. **Full schedule** (Aug 2026 paper runs): Stage1 localization 30 epochs, Stage2 perturbation 15 epochs (and Stage3 when used). Config: `configs/k562.yaml`.
2. **Fast schedule** (5 Sep 2026 graph / edge ablations): Stage1 8 epochs, Stage2 5 epochs. Config: `configs/k562_fast.yaml`.

The fast run **overwrote** `results/ablations/full_curriculum_K562.json`. The full-schedule full-curriculum numbers are still in `results/K562_results.json` and `results/summary_table.csv` (rows `MEvD-GRN` and `abl:full_curriculum`).


## 1. Main shipped model (full schedule)

Source: `results/K562_results.json`  
Script: `scripts/03_train.py --config configs/k562.yaml`  
Params: 306,691  
Train time: ~996 s (~16.6 min)  
Protocol: loc then pert, memory replay on, dual held out  
Curriculum val AUPR checkpoints: Stage1 0.940, Stage2 0.578

| Tier | Role | AUPR | AUROC | EP | EPR |
|---|---|---:|---:|---:|---:|
| Localization | trained (then forgotten) | 0.420 | 0.326 | 0.425 | 12.5 |
| Perturbation | trained | 0.575 | 0.876 | 0.536 | 80.2 |
| Dual-evidence | zero-shot | **0.558** | **0.838** | 0.419 | 272.5 |

Takeaway: strong zero-shot dual; localization AUROC below 0.5 after Stage2 (catastrophic forgetting).


## 2. Baselines (same splits)

Sources: `results/baselines/*_K562.json`  
Script: `scripts/05_run_baselines.py`

### 2.1 scMultiomeGRN (adapted wrapper)

Note: adapted TF→all-gene decoder; paper flags reduced epoch budget vs published 2000.

| Tier | AUPR | AUROC | EP | EPR |
|---|---:|---:|---:|---:|
| Localization | **0.831** | 0.865 | 0.804 | 23.7 |
| Perturbation | 0.423 | 0.828 | 0.443 | 66.4 |
| Dual-evidence | 0.384 | 0.825 | 0.402 | 261.0 |

### 2.2 GRNBoost2 (RNA-only)

| Tier | AUPR | AUROC | EP | EPR |
|---|---:|---:|---:|---:|
| Localization | 0.535 | 0.499 | 0.537 | 15.8 |
| Perturbation | 0.206 | 0.550 | 0.215 | 32.3 |
| Dual-evidence | 0.173 | 0.525 | 0.175 | 113.6 |

### 2.3 RegDiffusion (RNA-only)

| Tier | AUPR | AUROC | EP | EPR |
|---|---:|---:|---:|---:|
| Localization | 0.546 | 0.514 | 0.947 | 27.9 |
| Perturbation | 0.185 | 0.493 | 0.732 | 109.6 |
| Dual-evidence | 0.165 | 0.490 | 0.691 | 449.3 |

### 2.4 GMF-GAE (RNA-only graph AE)

| Tier | AUPR | AUROC | EP | EPR |
|---|---:|---:|---:|---:|
| Localization | 0.526 | 0.514 | 0.544 | 16.0 |
| Perturbation | 0.218 | 0.564 | 0.249 | 37.4 |
| Dual-evidence | 0.160 | 0.516 | 0.149 | 97.1 |

### 2.5 Baseline vs main model (AUPR only)

| Method | Loc | Pert | Dual |
|---|---:|---:|---:|
| MEvD-GRN (main, full) | 0.420 | **0.575** | **0.558** |
| scMultiomeGRN | **0.831** | 0.423 | 0.384 |
| GRNBoost2 | 0.535 | 0.206 | 0.173 |
| RegDiffusion | 0.546 | 0.185 | 0.165 |
| GMF-GAE | 0.526 | 0.218 | 0.160 |


## 3. Full-schedule ablations (Aug 2026)

Sources: `results/ablations/*.json` dated ~21–22 Aug 2026, and `summary_table.csv`.  
Script: `scripts/06_ablation.py --config configs/k562.yaml`  
Same 2-stage main protocol unless noted.

### 3.1 Training-protocol ablations

| Ablation | What changed | Trained on | Loc AUPR | Pert AUPR | Dual AUPR |
|---|---|---|---:|---:|---:|
| full_curriculum | reference 2-stage (matches main) | loc→pert | 0.410 | 0.569 | 0.546 |
| loc_only | only localization | loc | 0.940 | 0.187 | 0.663 |
| pert_only | only perturbation | pert | 0.407 | 0.496 | 0.485 |
| dual_only | only dual (exception: not zero-shot) | dual | 0.779 | 0.197 | 0.479 |
| all_at_once | loc+pert jointly, no sequence | loc+pert | 0.937 | 0.348 | **0.788** |
| with_replay | 3-stage + consolidation; trains dual | loc→pert→dual | 0.765 | 0.411 | **0.940** |

Notes:
- `with_replay` and `dual_only` train on dual, so dual is **not** a zero-shot claim.
- `with_replay` dual n_pos in JSON is 3979 (test-only style), others often use val+test (7943).

### 3.2 Architecture ablations (still loc→pert, dual zero-shot)

| Ablation | What changed | Loc AUPR | Pert AUPR | Dual AUPR |
|---|---|---:|---:|---:|
| rna_only | ATAC features and openness zeroed | 0.414 | 0.540 | 0.543 |
| gated_fusion | legacy gated fusion + bilinear | 0.419 | 0.588 | 0.546 |
| concat_fusion | concat fusion + bilinear | 0.421 | 0.595 | 0.555 |
| no_gnn | no message passing | 0.498 | 0.341 | **0.241** |

Takeaways from full-schedule ablations:
- GNN is critical (`no_gnn` dual 0.241).
- ATAC currently adds little (`rna_only` ≈ full).
- Fusion style is second-order.
- `all_at_once` and `with_replay` beat shipped curriculum on dual.


## 4. Fast-schedule graph and edge runs (5 Sep 2026)

Config: `configs/k562_fast.yaml` (8+5 epochs)  
Log: `results/logs/abl_batch_fast.log`  
Batch finished: 03:04 IST  
JSON files: `coexpr_only`, `tf_cand_only`, `edge_mlp`, and overwritten `full_curriculum`

These numbers are **not** directly comparable to the paper full-schedule table, but they are comparable to each other.

| Run | Graph / decoder | Loc AUPR | Pert AUPR | Dual AUPR | Dual AUROC |
|---|---|---:|---:|---:|---:|
| full_curriculum (fast) | both graphs, linear edge terms | 0.404 | 0.507 | **0.451** | 0.751 |
| coexpr_only | co-expression graph only | 0.418 | 0.441 | 0.383 | 0.751 |
| tf_cand_only | TF-candidate graph only | 0.418 | 0.374 | 0.258 | 0.660 |
| edge_mlp | both graphs, nonlinear edge MLP | 0.400 | 0.495 | 0.402 | 0.740 |

Wall times (approx from log):
- coexpr_only: ~4.3 min
- tf_cand_only: ~1.4 min
- edge_mlp: ~4.7 min
- full_curriculum (fast): ~4.7 min
- Total batch: ~15 min

### 4.1 Reading the fast graph ablation

On this short schedule:

1. **Both graphs beat either alone on dual** (0.451 vs 0.383 vs 0.258).
2. **Co-expression alone is much stronger than TF-candidate alone** on dual and perturbation.
3. TF-candidate-only is close to the old `no_gnn` dual level (0.258 vs 0.241 full-schedule), so that prior alone is weak under this setup.
4. **edge_mlp** (0.402 dual) did **not** beat fast both-graphs linear (0.451) on this budget. Needs a full-schedule rerun before calling it a win or loss.

### 4.2 Full detail for fast runs

**coexpr_only**

| Tier | AUPR | AUROC | EP | EPR |
|---|---:|---:|---:|---:|
| Localization | 0.418 | 0.307 | 0.393 | 11.6 |
| Perturbation | 0.441 | 0.800 | 0.447 | 66.9 |
| Dual-evidence | 0.383 | 0.751 | 0.371 | 241.0 |

**tf_cand_only**

| Tier | AUPR | AUROC | EP | EPR |
|---|---:|---:|---:|---:|
| Localization | 0.418 | 0.307 | 0.393 | 11.6 |
| Perturbation | 0.374 | 0.747 | 0.393 | 58.9 |
| Dual-evidence | 0.258 | 0.660 | 0.284 | 184.6 |

**edge_mlp**

| Tier | AUPR | AUROC | EP | EPR |
|---|---:|---:|---:|---:|
| Localization | 0.400 | 0.277 | 0.401 | 11.8 |
| Perturbation | 0.495 | 0.830 | 0.487 | 72.9 |
| Dual-evidence | 0.402 | 0.740 | 0.375 | 243.7 |

**full_curriculum (fast overwrite)**

| Tier | AUPR | AUROC | EP | EPR |
|---|---:|---:|---:|---:|
| Localization | 0.404 | 0.282 | 0.404 | 11.9 |
| Perturbation | 0.507 | 0.832 | 0.484 | 72.5 |
| Dual-evidence | 0.451 | 0.751 | 0.377 | 245.1 |


## 5. File index

| Path | Contents |
|---|---|
| `results/K562_results.json` | Main full-schedule MEvD-GRN |
| `results/summary_table.csv` | Compiled table (main + baselines + Aug ablations) |
| `results/baselines/*.json` | Per-baseline metrics |
| `results/ablations/*.json` | Per-ablation metrics (fast files dated 5 Sep) |
| `results/logs/abl_batch_fast.log` | Fast ablation batch log |
| `results/logs/abl_*_fast.log` | Per-run fast logs |
| `results/checkpoints/` | Saved weights per run |


## 6. Caveats

1. Prefer `results/K562_results.json` for the **paper headline** MEvD-GRN numbers, not the Sep 5 overwritten `full_curriculum` JSON.
2. scMultiomeGRN comparison is protocol-matched on splits but not yet a faithful full-budget original-code run.
3. `with_replay` dual AUPR 0.940 is impressive but that run **trains on dual**; do not mix it with zero-shot dual claims.
4. Fast graph/edge results need a full-schedule confirmation before changing the default model.
5. Only K562 is fully reported here. ESC / transfer runs are not in these JSON files yet.


| Method | Loc | Pert | Dual |
|---|---:|---:|---:|
| **MEvD-GRN** | 0.420 | 0.575 | 0.558 |
| scMultiomeGRN | 0.831 | 0.523 | 0.504 |
| GRNBoost2 | 0.535 | 0.206 | 0.173 |
| RegDiffusion | 0.546 | 0.185 | 0.165 |
| GMF-GAE | 0.526 | 0.218 | 0.160 |