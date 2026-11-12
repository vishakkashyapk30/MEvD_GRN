#!/bin/bash
# Step 3: GPU work as a shared queue (slurm/scmgrn_worker_lib.sh). 4 identical workers, one per
# 2080 Ti: together they work on ONE cell type at a time (official scMultiomeGRN re-runs first, then
# MeVD-GRN H / A1 / A2 / A3 x 5 seeds), claiming units dynamically, then move to the next cell type.
# Workers start together with the build and poll (5 min) for finished cell types. Each unit = 10 folds. Resubmitting resumes (done markers + existing fold JSONs are skipped).
# DS=hema for the hematopoiesis benchmark.
#SBATCH --job-name=smg-work
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=30G                      # the official code loads the dense expression matrix (Bronchiolar 19k cells x 14k genes, Stromal 43k x 11.5k): 6 GB (2 CPUs x 3 GB) gets it OOM-killed with no traceback
#SBATCH --time=3-12:00:00
#SBATCH --array=0-3%4
#SBATCH --output=logs/smg_work_%A_%a.out
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn_common.sh"
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn_worker_lib.sh"
activate
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1     # Geneformer cached by scmgrn_setup.sh
mkdir -p "$SHARE/claims/$DS" "$SHARE/done/$DS" "$SHARE/failed/$DS"
q_exec() { "$@"; }
res_dir() { [ "$SPLIT_KIND" = official ] && echo "results/$DS" || echo "results/${DS}_$SPLIT_KIND"; }
stage_ct() { pull official/ScmultiomeGRN-main; pull official_ablation/ScmultiomeGRN-main; local CT=$1
  pull "graph/$DS/$CT"; pull "splits/$DS/$CT"; pull "processed/$DS/$CT"; pull "$(res_dir)/$CT"
  rm -rf "$SCMGRN_ROOT_LOCAL/graph/$DS/$CT/work"; }
publish_ct() { push "$(res_dir)/$1"; }
WAIT_FOR_BUILD=1 worker_loop    # start as soon as the first cell type is built; wait for the rest
echo "[work] worker $SLURM_ARRAY_TASK_ID finished"
