#!/bin/bash
# GPU workers that run a PLAN file of worker configurations in order (slurm/scmgrn/plans/*.plan). Every worker claims
# units from the shared queue (slurm/scmgrn/scmgrn_worker_lib.sh), so N workers split each configuration; a worker that
# finds nothing claimable moves to the next line. Resubmitting resumes (done markers + existing fold files).
#   sbatch --export=ALL,PLAN=slurm/scmgrn/plans/hard_tests_2026-11.plan slurm/scmgrn/scmgrn_work_plan.sh
#SBATCH --job-name=smg-plan
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=22G                      # ~9 CPUs: 4 workers fit under the 40-CPU cap with room for CPU jobs
#SBATCH --time=4-00:00:00
#SBATCH --array=0-3%4
#SBATCH --output=logs/smg_plan_%A_%a.out
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn/scmgrn_common.sh"
activate
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
: "${PLAN:?set PLAN=slurm/scmgrn/plans/<file>.plan}"
while IFS= read -r line; do
  [[ -z "$line" || "$line" == \#* ]] && continue
  echo "===== plan line: $line ($(date +%T))"
  (
    eval "export $line"
    source "$REPO/slurm/scmgrn/scmgrn_worker_lib.sh" || exit 1
    mkdir -p "$SHARE/claims/$DS" "$SHARE/done/$DS" "$SHARE/failed/$DS"
    q_exec() { "$@"; }
    res_dir() { [ "$SPLIT_KIND" = official ] && echo "results/$DS" || echo "results/${DS}_$SPLIT_KIND"; }
    stage_ct() { pull official/ScmultiomeGRN-main; pull official_ablation/ScmultiomeGRN-main; local CT=$1
      pull "graph/$DS/$CT"; pull "splits/$DS/$CT"; pull "processed/$DS/$CT"; pull "$(res_dir)/$CT"
      rm -rf "$SCMGRN_ROOT_LOCAL/graph/$DS/$CT/work"; }
    publish_ct() { push "$(res_dir)/$1"; }
    # WORKERS > 1: that many worker loops share this GPU (MeVD-GRN runs are CPU-bound, ~0-5% GPU), each with its own
    # scratch tree and an equal share of the job's CPUs; they claim from the same queue as every other worker.
    W=${WORKERS:-1}; NCPU=${SLURM_CPUS_ON_NODE:-$(nproc)}; BASE_SCR=$SCR
    for k in $(seq 1 "$W"); do
      (
        if [ "$W" -gt 1 ]; then
          export SCR="$BASE_SCR/w$k" SCMGRN_ROOT_LOCAL="$BASE_SCR/w$k/scmgrn"; mkdir -p "$SCMGRN_ROOT_LOCAL"
          T=$(( NCPU / W > 0 ? NCPU / W : 1 )); export OMP_NUM_THREADS=$T MKL_NUM_THREADS=$T OPENBLAS_NUM_THREADS=$T
        fi
        WAIT_FOR_BUILD=0 worker_loop 2>&1 | sed -u "s/^/[w$k] /"
      ) &
      sleep 20                                        # stagger the first claims/stage-ins
    done
    wait
  )
done < "$REPO/$PLAN"
echo "[plan] worker ${SLURM_ARRAY_TASK_ID:-0} finished all lines"
