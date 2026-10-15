#!/bin/bash
# Step 3: GPU work as a shared queue (slurm/scmgrn_worker_lib.sh). 4 identical workers, one per
# 2080 Ti: together they work on ONE cell type at a time (official scMultiomeGRN re-runs first, then
# MeVD-GRN H / A1 / A2 / A3 x 5 seeds), claiming units dynamically, then move to the next cell type.
# Each unit = 10 folds. Resubmitting resumes (done markers + existing fold JSONs are skipped).
# DS=hema for the hematopoiesis benchmark.
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
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn_common.sh"
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn_worker_lib.sh"
activate
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1     # Geneformer cached by scmgrn_setup.sh
mkdir -p "$SHARE/claims/$DS" "$SHARE/done/$DS" "$SHARE/failed/$DS"
q_exec() { "$@"; }
stage_ct() { pull official/ScmultiomeGRN-main; local CT=$1
  pull "graph/$DS/$CT"; pull "splits/$DS/$CT"; pull "processed/$DS/$CT"; pull "results/$DS/$CT"
  rm -rf "$SCMGRN_ROOT_LOCAL/graph/$DS/$CT/work"; }
publish_ct() { push "results/$DS/$1"; }
WAIT_FOR_BUILD=0 worker_loop
echo "[work] worker $SLURM_ARRAY_TASK_ID finished"
