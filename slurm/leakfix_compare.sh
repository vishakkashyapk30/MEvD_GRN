#!/bin/bash
# Before/after comparison for the 2026-09-30 negative-sampling leak fix on the
# K562 headline setting (FM + h384/l2): all_at_once (the paper's recommended
# model) and full_curriculum, 5 seeds each with the fixed code, plus seed 42 with
# the legacy (leaky) pool, which should reproduce the pre-fix numbers in
# results/ablations/*_seed42_K562.json.
# Submit from the repo root on ada-gw1:  sbatch --array=0-11%2 slurm/leakfix_compare.sh
# (%2 keeps at most 2 GPUs busy; PBMC jobs have priority on the 4-GPU QoS cap.)
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -p u22
#SBATCH --constraint=2080ti
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --gres=gpu:1
#SBATCH --mem-per-cpu=3000
#SBATCH --time=1-00:00:00
#SBATCH --job-name=leakfix
#SBATCH --output=logs/leakfix_%A_%a.out
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ~/miniforge3/etc/profile.d/conda.sh
conda activate mevd-grn
nvidia-smi --query-gpu=name --format=csv,noheader
run() {  # config ablation seed tag
  echo "=== $(date) $2 seed=$3 tag=$4"
  python -u scripts/06_ablation.py --config "$1" --ablation "$2" --seed "$3" --tag "$4" --device cuda:0
}
# Tasks 0-9: fixed code, (all_at_once, full_curriculum) x seeds 42-46.
# Tasks 10-11: legacy (leaky) pool, seed 42, all_at_once then full_curriculum.
i=${SLURM_ARRAY_TASK_ID:-0}
abl=(all_at_once full_curriculum)
if [ "$i" -lt 10 ]; then
  s=$((42 + i / 2)); a=${abl[$((i % 2))]}
  run configs/sweep/k562_fm_h384_l2.yaml "$a" "$s" "fm_h384l2_fixneg_seed$s"
else
  a=${abl[$((i - 10))]}
  run configs/sweep/k562_fm_h384_l2_legacyneg.yaml "$a" 42 fm_h384l2_legacyneg_seed42
fi
echo "=== done $(date)"
