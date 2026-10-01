#!/bin/bash
# Step 2: MeVD-GRN training on up to 4 independent single-GPU array tasks.
# Unit = (variant, cell type, label source, regime, seed). Units are ordered by
# priority: the pre-registered headline (fm_h384, CollecTRI, TF-disjoint, 5 seeds)
# and its RNA-only ablation first, then the other regimes / variants / DoRothEA.
# Work stealing: each task walks the list and CLAIMS the next free unit with an
# atomic mkdir on /share1 (via the gateway), so however many GPUs the QoS grants
# (other users' jobs share the 4-GPU cap), the priority order is respected.
# Results are pushed to /share1 after every unit; a resubmission (new job id ->
# new claim dir) skips finished units via --skip_existing.
#SBATCH --job-name=pbmc-work
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3000
#SBATCH --time=3-00:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/pbmc_work_%A_%a.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/pbmc_common.sh"
activate
SEEDS=(${SEEDS:-42 43 44 45 46})
# tier 1 = headline + key ablation (CollecTRI, TF-disjoint)
GRID=(${GRID:-"fm_h384:atac:collectri:tf" "fm_h384:rna:collectri:tf"
  "fm_h384:atac:collectri:target" "fm_h384:rna:collectri:target"
  "fm_h384:atac:collectri:random" "fm_h384:rna:collectri:random"
  "base:atac:collectri:tf" "base:rna:collectri:tf" "fm_h384_uniformneg:atac:collectri:tf"
  "fm_h384:atac:dorothea_ab:tf" "fm_h384:rna:dorothea_ab:tf"})
UNITS=()
for G in "${GRID[@]}"; do for S in "${SEEDS[@]}"; do for CT in "${PBMC_CTS[@]}"; do
  UNITS+=("$G:$CT:$S"); done; done; done
for r in groundtruth splits; do pull "$r"; done
declare -A STAGED=()
CLAIMS="$SHARE/claims/${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local}}"
claim() {   # atomic: succeeds for exactly one task
  if [ -d "$SHARE" ]; then mkdir -p "$CLAIMS" && mkdir "$CLAIMS/$1" 2>/dev/null
  else ssh "$MASTER" "mkdir -p '$CLAIMS' && mkdir '$CLAIMS/$1'" 2>/dev/null; fi
}
for k in "${!UNITS[@]}"; do
  IFS=: read -r MC MOD SRC REG CT S <<<"${UNITS[$k]}"
  V=$MC; [ "$MOD" = rna ] && V="${MC}_rnaonly"
  [ -e "$PBMC_ROOT_LOCAL/results/$CT/$SRC/$REG/$V/seed$S/metrics.json" ] && continue
  claim "u$(printf %03d $k)" || continue
  [ -n "${STAGED[$CT]:-}" ] || { pull "processed/$CT"; pull "results/$CT"; STAGED[$CT]=1; }
  FLAG=""; [ "$MOD" = rna ] && FLAG="--rna_only"
  echo "===== unit $k/${#UNITS[@]}: ${UNITS[$k]} ($(date +%T)) ====="
  python -u scripts/27_pbmc_train_eval.py --config configs/pbmc/benchmark.yaml \
      --model_config "configs/pbmc/mevd_$MC.yaml" $FLAG --root "$PBMC_ROOT_LOCAL" \
      --cell_type "$CT" --source "$SRC" --regime "$REG" --seed "$S" --device cuda:0 \
      --skip_existing --tmp_dir "$SCR" || echo "[warn] unit failed: ${UNITS[$k]}"
  push "results/$CT"
done
echo "[work] task ${SLURM_ARRAY_TASK_ID:-0} done (no unclaimed units left)"
