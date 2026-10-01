#!/bin/bash
# LINGER baseline re-run (LingerGRN==1.110, method='LINGER') on the benchmark's
# exact cells, then export LINGER's cell-type-specific TG sets (candidate space
# `linger_tg`) and eval-TF trans scores. Needs: `bash scripts/25_pbmc_download.sh
# linger linger-env` done on the gateway, and the build job's cells/ on /share1.
#SBATCH --job-name=pbmc-linger
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3000
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/pbmc_linger_%j.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/pbmc_common.sh"
source "$CONDA_SH"; conda activate linger; export PYTHONUNBUFFERED=1; cd "$REPO"
echo "[env] host=$(hostname) scratch=$SCR"
for r in raw/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5 cells linger/data_bulk.tar.gz; do pull "$r"; done
tar -xzf "$PBMC_ROOT_LOCAL/linger/data_bulk.tar.gz" -C "$PBMC_ROOT_LOCAL/linger/" && rm "$PBMC_ROOT_LOCAL/linger/data_bulk.tar.gz"
python -u src/benchmarks/pbmc_linger_run.py --root "$PBMC_ROOT_LOCAL" \
    --grndir "$PBMC_ROOT_LOCAL/linger/data_bulk/" --outdir "$PBMC_ROOT_LOCAL/linger/out/" \
    --steps "${STEPS:-prep,train,celltype,export}"
# keep LINGER's cell-type-specific trans matrices + exports (not the 40 GB data_bulk)
mkdir -p "$PBMC_ROOT_LOCAL/linger/keep"
cp "$PBMC_ROOT_LOCAL"/linger/out/cell_type_specific_trans_regulatory_*.txt "$PBMC_ROOT_LOCAL/linger/keep/" || true
cp "$PBMC_ROOT_LOCAL"/linger/out/cell_population_trans_regulatory.txt "$PBMC_ROOT_LOCAL/linger/keep/" || true
rm -rf "$PBMC_ROOT_LOCAL/linger/data_bulk" "$PBMC_ROOT_LOCAL/linger/out"
push linger
for CT in "${PBMC_CTS[@]}"; do push "results/$CT/unsupervised/none/linger"; done
echo "[linger] done"
