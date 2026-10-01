#!/bin/bash
# Step 1 (CPU work; 1 GPU only for Ada's 1 GPU : 10 CPU policy):
#   cells (LINGER QC) -> processed dirs (4 cell types, label-free graphs) -> FM
#   -> label splits (CollecTRI + DoRothEA A/B, 3 regimes) -> trivial baselines
#   (degree, gene-ID LR, |Pearson|, signed Pearson) for every cell type.
#SBATCH --job-name=pbmc-build
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3000
#SBATCH --time=08:00:00
#SBATCH --output=logs/pbmc_build_%j.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/pbmc_common.sh"
activate
for r in raw groundtruth labels; do pull "$r"; done
R="$PBMC_ROOT_LOCAL"
python -u scripts/26_pbmc_build.py --config configs/pbmc/benchmark.yaml --root "$R" \
    --steps cells,processed,fm,splits
push cells; push processed; push splits
for CT in "${PBMC_CTS[@]}"; do
  python -u scripts/28_pbmc_baselines.py --root "$R" --cell_type "$CT" \
      --methods degree,geneid --sources collectri,dorothea_ab --regimes tf,target,random
  python -u scripts/28_pbmc_baselines.py --root "$R" --cell_type "$CT" --methods pearson,pearson_signed
  push "results/$CT"
done
echo "[build] done"
