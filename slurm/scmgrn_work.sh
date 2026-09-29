#!/bin/bash
# Step 3: all GPU work, sharded over 4 independent single-GPU tasks (%4).
# Units = lung cell type x {MeVD-GRN variant x seed 42-46 (10 folds each),
# official scMultiomeGRN reproduction (10 folds, seed 666)}; unit k runs on
# shard k % 4. Why not DDP: every fold is a ~50-500-node graph that trains in
# seconds-to-minutes on one GPU, so independent runs use the 4 GPUs fully with
# no communication; one array job also respects medium QoS's 8-submitted-jobs
# cap (each array task counts). Results are pushed to /share1 after every unit
# and existing fold JSONs are skipped, so a resubmission resumes.
#SBATCH --job-name=smg-work
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=long
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3000
#SBATCH --time=3-00:00:00
#SBATCH --array=0-3%4
#SBATCH --output=logs/smg_work_%A_%a.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn_common.sh"
activate
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1     # Geneformer cached by `14_... fm`
NSH=${SLURM_ARRAY_TASK_COUNT:-4}
SEEDS=(${SEEDS:-42 43 44 45 46})
# variant = <model config stem>:<obs|noobs>; first = pre-registered headline
VARIANTS=(${VARIANTS:-mevd_fm_h384:obs mevd_fm_h384:noobs mevd_base:obs})
UNITS=()
for CT in "${LUNG_CTS[@]}"; do
  UNITS+=("official $CT")
  for V in "${VARIANTS[@]}"; do for S in "${SEEDS[@]}"; do UNITS+=("mevd $CT $V $S"); done; done
done
pull official/ScmultiomeGRN-main
declare -A STAGED=()
stage_ct() {
  local CT=$1; [ -n "${STAGED[$CT]:-}" ] && return 0
  pull "graph/lung/$CT"; pull "splits/lung/$CT"; pull "processed/lung/$CT"; pull "results/lung/$CT"
  rm -rf "$SCMGRN_ROOT_LOCAL/graph/lung/$CT/work"
  STAGED[$CT]=1
}
for k in "${!UNITS[@]}"; do
  (( k % NSH == SLURM_ARRAY_TASK_ID )) || continue
  read -r KIND CT V S <<<"${UNITS[$k]}"
  stage_ct "$CT"
  echo "===== unit $k: ${UNITS[$k]} ($(date +%T)) ====="
  if [ "$KIND" = official ]; then
    python -u scripts/18_scmgrn_official_repro.py --config configs/scmgrn/lung.yaml \
      --root "$SCMGRN_ROOT_LOCAL" --cell_type "$CT" --folds 1-10 --skip_existing \
      --work_dir "$SCR/official_work/$CT" || echo "[warn] official repro failed for $CT"
  else
    MC=${V%%:*}; OBS=${V##*:}
    if grep -q "use_fm: true" "configs/scmgrn/$MC.yaml" && \
       [ ! -e "$SCMGRN_ROOT_LOCAL/processed/lung/$CT/fm_gene_embeddings.npy" ]; then
      python -u scripts/16_scmgrn_build_benchmark.py --config configs/scmgrn/lung.yaml \
        --root "$SCMGRN_ROOT_LOCAL" --cell_type "$CT" --steps fm
      push "processed/lung/$CT"
    fi
    FLAG=""; [ "$OBS" = obs ] && FLAG="--obs_graph"
    python -u scripts/17_scmgrn_train_eval.py --config configs/scmgrn/lung.yaml \
      --model_config "configs/scmgrn/$MC.yaml" $FLAG --root "$SCMGRN_ROOT_LOCAL" \
      --cell_type "$CT" --seed "$S" --folds 1-10 --device cuda:0 --skip_existing --tmp_dir "$SCR"
  fi
  push "results/lung/$CT"
done
echo "[work] shard $SLURM_ARRAY_TASK_ID done"
