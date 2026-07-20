#!/bin/bash
# Baseline array job (plan Step 13). One array task per baseline.
#   sbatch slurm/baselines.sh
#SBATCH --job-name=MEvD-base
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=long
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G
#SBATCH --time=04:00:00
#SBATCH --output=logs/base_%A_%a.out
#SBATCH --array=0-2

set -euo pipefail
mkdir -p logs
cd "$SLURM_SUBMIT_DIR"
source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate mevd-grn

# arboreto (GRNBoost2) is optional: pip install arboreto
# regdiffusion is optional:       pip install regdiffusion
BASE=(grnboost2 regdiffusion gmf_gae)
NAME=${BASE[$SLURM_ARRAY_TASK_ID]}
for CT in configs/k562.yaml configs/esc.yaml; do
  echo "### Baseline $NAME on $CT ###"
  python scripts/05_run_baselines.py --config "$CT" --baselines "$NAME" --device cuda:0
done
