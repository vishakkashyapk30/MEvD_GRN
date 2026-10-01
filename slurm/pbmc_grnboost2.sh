#!/bin/bash
# GRNBoost2 baseline (CPU-heavy: SGBM per target over ~25k targets x ~940 regulators),
# all 4 cell types in ONE job (1 submit slot). The GPU is requested only for the 1 GPU : 10 CPU policy.
#SBATCH --job-name=pbmc-gb2
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3000
#SBATCH --time=3-00:00:00
#SBATCH --output=logs/pbmc_gb2_%j.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/pbmc_common.sh"
activate
for r in raw groundtruth cells; do pull "$r"; done
for CT in mdc naive_b naive_cd4_t classical_monocyte; do      # smallest first
  pull "processed/$CT"; pull "results/$CT"
  python -u scripts/28_pbmc_baselines.py --root "$PBMC_ROOT_LOCAL" --cell_type "$CT" \
      --methods grnboost2 --n_jobs "${SLURM_CPUS_PER_TASK:-10}" --skip_existing
  push "results/$CT"
done
