#!/bin/bash
# s12.6 dev selection on K562 (seed 42): one GPU per candidate, array 0-3 = M0 M1 M2 M3.
#   sbatch slurm/bear_ada_dev.sh      (from ~/mevd_grn on ada-gw1)
#SBATCH --job-name=bear-dev
#SBATCH -p u22
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --constraint=2080ti
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3000
#SBATCH --time=0-08:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/bear_dev_%A_%a.out
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/bear_ada_common.sh"
VARS=(fm_h384 fm_h384_uniformneg fm_h384_hub fm_h384_hubmotif)
V=${VARS[${SLURM_ARRAY_TASK_ID:-0}]}; DS=K562; SEED=42
activate
stage_in processed/$DS || exit 1
echo "===== dev $DS $V seed $SEED $(date +%T) ====="
python -u scripts/22_bear_train_eval.py --root "$BEAR_ROOT" --dataset $DS \
    --model_config configs/bear/mevd_$V.yaml --seed $SEED --device cuda:0 --tmp_dir "$SCR/tmp" \
    > "$SCR/train.log" 2>&1 &
PID=$!
# stream the log to $B/logs so it can be watched from the gateway
( while kill -0 $PID 2>/dev/null; do cp -f "$SCR/train.log" "$B/logs/dev_${V}_s${SEED}.log" 2>/dev/null; sleep 60; done ) &
wait $PID; RC=$?
cp -f "$SCR/train.log" "$B/logs/dev_${V}_s${SEED}.log"
if [ $RC -eq 0 ]; then
  publish_dir "$BEAR_ROOT/results/$DS/mevd_$V" "$B/results/$DS/mevd_$V" && echo "[publish] OK" || { echo "[publish] FAILED"; exit 2; }
fi
echo "[dev] $V rc=$RC"; exit $RC
