#!/bin/bash
# Ablation array job (plan Part 10). One array task per ablation variant.
#   sbatch slurm/ablation.sh
#SBATCH --job-name=MEvD-abl
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=long
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G
#SBATCH --time=06:00:00
#SBATCH --output=logs/abl_%A_%a.out
#SBATCH --array=0-8

set -euo pipefail
mkdir -p logs
cd "$SLURM_SUBMIT_DIR"
source "$HOME/miniconda3/etc/profile.d/conda.sh" 2>/dev/null || true
conda activate mevd-grn

ABL=(full_curriculum loc_only pert_only dual_only all_at_once rna_only concat_fusion no_gnn with_replay)
NAME=${ABL[$SLURM_ARRAY_TASK_ID]}
echo "### Ablation: $NAME (K562) ###"
python scripts/06_ablation.py --config configs/k562.yaml --ablation "$NAME" --device cuda:0
