# MEvD-GRN — Multi-Evidence Distillation for Gene Regulatory Network Inference

A compact dual-modality GNN that predicts TF→target regulatory edges from
scRNA-seq + scATAC-seq, trained via a **three-stage evidence-aware curriculum**
over SC-MO-GRN-DB's tiered reference networks
(localization → perturbation → dual-evidence). See `MEvD_GRN_plan_v2.md` for the
full theory/spec.

---

## 1. Which data to download (and why)

The full SC-MO-GRN-DB single-cell collection is **~24.7 GB** (17.5 GB human +
7.2 GB mouse). You do **not** need it. The curriculum requires a cell type with
**all three evidence tiers**, and only **two** cell types in the entire database
have that complete triple:

| Cell type | Localization | Perturbation | Dual-evidence | Notes |
|-----------|-------------|--------------|---------------|-------|
| **K562** (human) | RN117 ChIPseq | RN118 KO | RN119 ChIPseq&KO | dual is from a *different* study → nesting is verified, not assumed |
| **ESC** (mouse)  | RN114 ChIP-chip | RN115 LOGOF | RN116 ChIPchip&LOGOF | **self-consistent** triple (one study, dual = loc ∩ pert) |

Every other cell type (HepG2, BJ, GM12878, H1, Macrophage, MCF7, HSC, DC) has
**localization only** → cannot exercise the curriculum.

**Single-cell datasets** (paired RNA+ATAC, kept deliberately small):

| Cell type | Dataset | Modalities | Cells | Size | Role |
|-----------|---------|-----------|-------|------|------|
| K562 | **DS025** | scRNA + scATAC (joint) | 434 | ~522 MB | primary (both modalities) |
| K562 | DS019 | scRNA | 953 | ~17 MB | optional extra RNA cells |
| ESC  | **DS010** | scRNA + scATAC (+scHiC) | 9021 | ~37 MB | primary (small + many cells) |
| ESC  | DS012 | scRNA + scATAC | 930 | ~12 MB | lightweight fallback |

Total for the two-cell-type experiment: **≈575 MB** instead of 24.7 GB.

`scripts/01_download_data.sh` fetches exactly these (per-file URLs:
`RN_TSV/RN###.tsv`, `DS_ZIP/DS###.zip`).

> Note: gene-level features are per-gene mean/variance, so RNA and ATAC are
> aggregated **independently** — paired barcodes are not required, which is why
> RNA-only and ATAC-only datasets can be combined freely.

---

## 2. Environment (Ada / SLURM)

Custom envs must live in `$HOME` (Ada policy); `/home` quota is 25 GB.

```bash
conda create -y -n mevd-grn python=3.10
conda activate mevd-grn
# torch matched to the cluster CUDA (driver permitting):
pip install torch==2.1.0 --index-url https://download.pytorch.org/whl/cu118
pip install torch_geometric torch_scatter torch_sparse \
    -f https://data.pyg.org/whl/torch-2.1.0+cu118.html
pip install -r requirements.txt
pip install -e .
```

We use **1 GPU** (model is ~125K params, trains in minutes); request 10 CPUs to
respect Ada's 1:10 GPU:CPU rule, `research` account, `medium` QoS.

> Local dev note: this machine has multiple Python interpreters; the one with
> torch+PyG is `/home/vishak/miniforge3/bin/python` (base env). If `python3`
> resolves to linuxbrew and can't import `torch_geometric`, call that interpreter
> explicitly or `conda activate base` first. On Ada, use the `mevd-grn` env above.

### scATAC needs TSS coordinates (offline)
Ada compute nodes have **no internet**, so peak→gene mapping needs an offline
GTF. On a **login node**, download GENCODE and point the config at it:
```bash
# human (K562): gencode.v38.annotation.gtf.gz   mouse (ESC): gencode.vM27.annotation.gtf.gz
# then set in configs/{k562,esc}.yaml:  atac: { gtf_path: /path/to/....gtf.gz }
```
Without a GTF the pipeline still runs with **ATAC features = 0** (RNA-only); the
gated fusion learns to ignore the empty modality.

---

## 3. Run

```bash
# (login node) fetch data
bash scripts/01_download_data.sh

# validate the code with no data / no GPU
python scripts/00_sanity_check.py

# preprocess -> train -> evaluate (or just: sbatch slurm/train.sh)
python scripts/02_preprocess.py --config configs/k562.yaml
python scripts/03_train.py      --config configs/k562.yaml --device cuda:0
python scripts/04_evaluate.py   --config configs/k562.yaml \
    --checkpoint results/checkpoints/K562/final_model.pt --transfer_config configs/esc.yaml

# baselines (GRNBoost2 / RegDiffusion / self-contained GMF-GAE), same test splits
pip install arboreto regdiffusion        # optional external baselines
python scripts/05_run_baselines.py --config configs/k562.yaml \
    --baselines grnboost2,regdiffusion,gmf_gae --device cuda:0

# ablations (paper claims)
sbatch slurm/ablation.sh
```

**Baselines** (`src/baselines/`, plan Part 9), all evaluated on the *identical*
test splits via `evaluate_scorer_on_splits`:
- **GRNBoost2** (arboreto) — tree-ensemble co-expression, RNA-only.
- **RegDiffusion** — diffusion GRN, RNA-only, GPU; fastest DL baseline.
- **GMF-GAE** — self-contained graph auto-encoder over a co-expression kNN graph
  (the "GNN + matrix factorization" family GMFGRN belongs to). Always runnable;
  the official GMFGRN repo can be wired in via `run_gmfgrn()` after cloning.
  Uninstalled external baselines fail gracefully (error saved to JSON, run continues).

---

## 4. Repo layout

```
configs/      default + per-cell-type YAML (k562, esc); *_localtest are smoke tests
src/data/     preprocessing (real text-matrix loaders), graph_builder, dataset
src/models/   encoders, gated fusion, GraphSAGE backbone, bilinear decoder, MEvDGRN
src/training/ curriculum, sampler (random + hard negatives), losses, trainer
src/evaluation/ metrics (AUPR/AUROC/EP/EPR), benchmarker
scripts/      00_sanity_check, 01_download_data.sh, 02_preprocess, 03_train,
              04_evaluate, 06_ablation
slurm/        train.sh, ablation.sh
```

---

## 5. Design notes / deviations from the original plan

These were forced by the **real** data format (verified against the live DB):

1. **Single-cell files are tab-delimited text, oriented features×cells** (not
   `.h5ad`, not cells×genes). Loaders transpose and parse accordingly.
2. **Gene-symbol case mismatch**: networks are UPPERCASE, RNA is mixed-case → all
   symbols canonicalised to UPPERCASE for matching.
3. **Nesting is verified, not assumed.** Measured on ESC: `dual ⊆ loc` and
   `dual ⊆ pert` = 100% (RN116 = RN114∩RN115), but `pert ⊆ loc` ≈ 12%
   (knockouts hit indirectly-regulated genes). `dual_evidence_mode: reconstruct`
   can force dual = loc∩pert when a self-consistent triple isn't available.
4. **Negatives = random TF×gene non-positives** (field standard), **not** the
   accessibility-ranked prior — restricting negatives to the prior inflates AUPR
   via an expression-level artifact.
5. **Prior graph excludes known positives** — it is a pure message-passing
   structure, so val/test edges cannot leak into node embeddings. The bilinear
   decoder scores any (TF, gene) pair, so positives need not lie in the prior.
```
