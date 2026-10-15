#!/bin/bash
# Step 3: all GPU work, sharded over 4 independent single-GPU tasks (%4), 2080 Ti nodes.
# Units = lung cell type x { official scMultiomeGRN reproduction (seeds 666, 42-46; 10 folds)
#   | MeVD-GRN variant x seed 42-46 (10 folds) }; unit k runs on shard k % 4.
# MeVD-GRN variants (undirected TF-TF link prediction; frozen before any result):
#   H  = mevd_fm_h384 : obs : smg      headline (scMultiomeGRN's inputs + train adjacency)
#   A1 = mevd_fm_h384 : obs : nosmg    without the smg input features
#   A2 = mevd_fm_h384 : noobs : nosmg  label-free graphs only
#   A3 = mevd_fm_h384_directed : obs : nosmg   old directed TF->target model, scores symmetrised
# Every fold is a <=700-node graph that trains in seconds-to-minutes on one GPU, so independent
# runs use the 4 GPUs fully with no communication (no DDP). Results are written to $SHARE (home)
# after every unit and existing fold JSONs are skipped, so a resubmission resumes.
#SBATCH --job-name=smg-work
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3000
#SBATCH --time=3-12:00:00
#SBATCH --array=0-3%4
#SBATCH --output=logs/smg_work_%A_%a.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn_common.sh"
activate
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1     # Geneformer cached by `14_... fm`
NSH=${SLURM_ARRAY_TASK_COUNT:-4}
SEEDS=(${SEEDS:-42 43 44 45 46})
OFFICIAL_SEEDS=(${OFFICIAL_SEEDS:-666 42 43 44 45 46})
# variant = <model config stem>:<obs|noobs>:<smg|nosmg>
VARIANTS=(${VARIANTS:-mevd_fm_h384:obs:smg mevd_fm_h384:obs:nosmg mevd_fm_h384:noobs:nosmg mevd_fm_h384_directed:obs:nosmg})
UNITS=()
for CT in "${CTS[@]}"; do
  for S in "${OFFICIAL_SEEDS[@]}"; do UNITS+=("official $CT x $S"); done
  for V in "${VARIANTS[@]}"; do for S in "${SEEDS[@]}"; do UNITS+=("mevd $CT $V $S"); done; done
done
pull official/ScmultiomeGRN-main
declare -A STAGED=()
stage_ct() {
  local CT=$1; [ -n "${STAGED[$CT]:-}" ] && return 0
  pull "graph/$DS/$CT"; pull "splits/$DS/$CT"; pull "processed/$DS/$CT"; pull "results/$DS/$CT"
  rm -rf "$SCMGRN_ROOT_LOCAL/graph/$DS/$CT/work"
  STAGED[$CT]=1
}
for k in "${!UNITS[@]}"; do
  (( k % NSH == SLURM_ARRAY_TASK_ID )) || continue
  read -r KIND CT V S <<<"${UNITS[$k]}"
  stage_ct "$CT"
  echo "===== unit $k: ${UNITS[$k]} ($(date +%T)) ====="
  if [ "$KIND" = official ]; then
    python -u scripts/18_scmgrn_official_repro.py --config configs/scmgrn/$DS.yaml \
      --root "$SCMGRN_ROOT_LOCAL" --cell_type "$CT" --folds 1-10 --skip_existing --model_seed "$S" \
      --work_dir "$SCR/official_work/$CT" || echo "[warn] official repro failed for $CT seed $S"
  else
    IFS=: read -r MC OBS SMG <<<"$V"
    if [ ! -e "$SCMGRN_ROOT_LOCAL/processed/$DS/$CT/fm_gene_embeddings.npy" ]; then   # fallback; normally built by scmgrn_build.sh
      python -u scripts/16_scmgrn_build_benchmark.py --config configs/scmgrn/$DS.yaml \
        --root "$SCMGRN_ROOT_LOCAL" --cell_type "$CT" --steps fm
      push "processed/$DS/$CT"
    fi
    FLAG=""; [ "$OBS" = obs ] && FLAG="--obs_graph"; [ "$SMG" = smg ] && FLAG="$FLAG --smg_feats"
    python -u scripts/17_scmgrn_train_eval.py --config configs/scmgrn/$DS.yaml \
      --model_config "configs/scmgrn/$MC.yaml" $FLAG --root "$SCMGRN_ROOT_LOCAL" \
      --cell_type "$CT" --seed "$S" --folds 1-10 --device cuda:0 --skip_existing --tmp_dir "$SCR"
  fi
  push "results/$DS/$CT"
done
echo "[work] shard $SLURM_ARRAY_TASK_ID done"
