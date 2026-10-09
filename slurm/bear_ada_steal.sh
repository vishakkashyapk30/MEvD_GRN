#!/bin/bash
# Work-stealing worker for the BEAR runs (arrays of identical workers). Units are lines of
# $B/units.txt:  "<dataset> <variant> <seed>"  (variant = configs/bear/mevd_<variant>.yaml).
# A worker claims a unit with an atomic mkdir in $B/claims (home is NFS; mkdir is atomic),
# trains it on scratch, publishes the result with verified copies, then writes $B/done/<unit>.
# A unit whose worker died keeps its claim without a done marker: `bash slurm/bear_ada_reset.sh`
# releases such stale claims (only when no worker is running).
#   sbatch --array=0-3 slurm/bear_ada_steal.sh
#SBATCH --job-name=bear-work
#SBATCH -p u22
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --constraint=2080ti
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3000
#SBATCH --time=1-18:00:00
#SBATCH --output=logs/bear_work_%A_%a.out
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/bear_ada_common.sh"
activate
START=$(date +%s); BUDGET=${BUDGET_SEC:-140000}     # stop claiming new units ~12 h before the walltime
while true; do
  [ $(( $(date +%s) - START )) -gt "$BUDGET" ] && { echo "[steal] time budget reached"; break; }
  UNIT=""
  while read -r DS V SEED; do
    [ -z "${DS:-}" ] && continue; case "$DS" in \#*) continue;; esac
    U="${DS}__${V}__seed${SEED}"
    [ -e "$B/done/$U" ] && continue
    if mkdir "$B/claims/$U" 2>/dev/null; then
      echo "$(hostname) job=${SLURM_JOB_ID:-} $(date +%F_%T)" > "$B/claims/$U/owner"; UNIT="$U"; break
    fi
  done < "$B/units.txt"
  [ -z "$UNIT" ] && { echo "[steal] no unclaimed units left"; break; }
  read -r DS V SEED <<<"$(echo "$UNIT" | sed 's/__/ /g; s/seed//')"
  PD="$DS"; [[ "$V" == *_L2 ]] && PD="${DS}__L2"
  echo "===== unit $UNIT $(date +%T) ====="
  rm -rf "$BEAR_ROOT/processed" "$BEAR_ROOT/results"; mkdir -p "$BEAR_ROOT/results"
  stage_in "processed/$PD" || { echo "[steal] stage-in failed"; mv "$B/claims/$UNIT" "$B/failed/$UNIT.$(date +%s)"; continue; }
  EXTRA=""; [ "${MODALITY:-1}" = 1 ] && [[ "$V" == fm_h384 || "$V" == *_hub* ]] && EXTRA="--modality_check"
  python -u scripts/22_bear_train_eval.py --root "$BEAR_ROOT" --dataset "$DS" \
      --model_config "configs/bear/mevd_$V.yaml" --seed "$SEED" --device cuda:0 --tmp_dir "$SCR/tmp" $EXTRA \
      > "$SCR/train.log" 2>&1
  RC=$?
  cp -f "$SCR/train.log" "$B/logs/${UNIT}.log" 2>/dev/null
  if [ $RC -eq 0 ] && publish_dir "$BEAR_ROOT/results/$DS/mevd_$V" "$B/results/$DS/mevd_$V" \
        && [ -s "$B/results/$DS/mevd_$V/seed$SEED/grn.tsv.gz" ]; then
    echo "$(date +%F_%T) rc=0" > "$B/done/$UNIT"; echo "[steal] DONE $UNIT"
  else
    echo "[steal] FAILED $UNIT rc=$RC"; mv "$B/claims/$UNIT" "$B/failed/$UNIT.$(date +%s)"
  fi
done
echo "[steal] worker finished"
