#!/bin/bash
# Step 1: split the Domcke lung Seurat RDS into per-cell-type MatrixMarket dirs.
#SBATCH --job-name=smg-prep
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=6000
#SBATCH --time=12:00:00
#SBATCH --output=logs/smg_prep_%j.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn_common.sh"
activate
export WAIT_RAW=1
fetch_raw GSM4508936_lung_filtered.seurat.RDS.gz "$SCMGRN_ROOT_LOCAL/raw/lung/GSM4508936_lung_filtered.seurat.RDS.gz"
gzip -t "$SCMGRN_ROOT_LOCAL/raw/lung/GSM4508936_lung_filtered.seurat.RDS.gz" && echo "[prep] RDS gzip integrity OK"
"$RSCRIPT_BIN" scripts/15_scmgrn_split_seurat.R \
   "$SCMGRN_ROOT_LOCAL/raw/lung/GSM4508936_lung_filtered.seurat.RDS.gz" \
   "$SCMGRN_ROOT_LOCAL/sep_data/lung" auto 1
push sep_data/lung
echo "[prep] done"
