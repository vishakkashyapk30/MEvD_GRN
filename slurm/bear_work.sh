#!/bin/bash
# BEAR step 2 (GPU): MeVD-GRN TF-disjoint cross-fitted runs, sharded over
# independent single-GPU array tasks. Units = dataset x variant x seed, headline
# variant first; unit u runs on task u % N. Every unit pushes its GRN to /share1
# when done and finished units are skipped, so a resubmission resumes.
# Each fold is a ~10k-gene graph that trains in minutes on one 2080 Ti, so
# independent runs (not DDP) use the GPUs fully.
#SBATCH --job-name=bear-work
#SBATCH -p u22
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --constraint=2080ti
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3000
#SBATCH --time=2-00:00:00
#SBATCH --array=0-3%2
#SBATCH --output=logs/bear_work_%A_%a.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/bear_common.sh"
activate
export HF_HUB_OFFLINE=1
NSH=${SLURM_ARRAY_TASK_COUNT:-4}
DATASETS=(${DATASETS:-${BEAR_PHASE1[@]}})
SEEDS=(${SEEDS:-42 43 44 45 46})
VARIANTS=(${VARIANTS:-fm_h384 base fm_h384_L2 fm_h384_rna_only fm_h384_uniformneg})
UNITS=()
for V in "${VARIANTS[@]}"; do for S in "${SEEDS[@]}"; do for DS in "${DATASETS[@]}"; do
  UNITS+=("$DS $V $S"); done; done; done
declare -A STAGED=()
for u in "${!UNITS[@]}"; do
  (( u % NSH == ${SLURM_ARRAY_TASK_ID:-0} )) || continue
  read -r DS V S <<<"${UNITS[$u]}"
  if [ -z "${STAGED[$DS]:-}" ]; then stage_ds "$DS"; STAGED[$DS]=1; fi
  PD="$DS"; [[ "$V" == *_L2 ]] && PD="${DS}__L2"
  if [ ! -e "$BEAR_ROOT/processed/$PD/summary.json" ]; then
    echo "[warn] $DS not prepared (run slurm/bear_prep.sh first)"; continue; fi
  echo "===== unit $u: $DS $V seed $S ($(date +%T)) ====="
  python -u scripts/22_bear_train_eval.py --root "$BEAR_ROOT" --dataset "$DS" \
      --model_config "configs/bear/mevd_$V.yaml" --seed "$S" --device cuda:0 \
      --skip_existing --tmp_dir "$SCR/tmp" || { echo "[warn] unit $u failed"; continue; }
  push "results/$DS/mevd_$V"
done
echo "[work] shard ${SLURM_ARRAY_TASK_ID:-0} done"
