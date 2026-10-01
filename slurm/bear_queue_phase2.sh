#!/bin/bash
# Gateway-side helper (run with nohup on ada-gw1): the user's QoS allows 20
# submitted jobs and the PBMC / leak-fix arrays fill it, so this waits for free
# slots and then submits the phase-2 chain (iPS + mouse embryo x4 + naive mESC):
#   prep (CPU) -> work array 0-1%1 (GPU) -> score (CPU)      = 4 job slots
# Phase 2 runs the headline (fm_h384) and L2 variants, seeds 42-44 (s9).
#   setsid nohup bash slurm/bear_queue_phase2.sh > ~/bear_tmp/queue_phase2.log 2>&1 < /dev/null &
set -uo pipefail
cd "$(dirname "$0")/.."
export DATASETS="${DATASETS:-iPS mESC_E7.5_rep1 mESC_E7.5_rep2 mESC_E8.5_rep1 mESC_E8.5_rep2 Naive_mESC}"
export VARIANTS="${VARIANTS:-fm_h384 fm_h384_L2}" SEEDS="${SEEDS:-42 43 44}"
need() { while [ $(( $(squeue -u $USER -r -h | wc -l) + $1 )) -gt 20 ]; do sleep 300; done; }
need 1; P=$(sbatch --parsable --export=ALL --time=1-00:00:00 slurm/bear_prep.sh); echo "phase2 prep:  $P ($(date))"
need 2; W=$(sbatch --parsable --export=ALL --dependency=afterok:$P --array=0-1%1 slurm/bear_work.sh); echo "phase2 work:  $W ($(date))"
need 1; S=$(sbatch --parsable --export=ALL --dependency=afterany:$W slurm/bear_score.sh); echo "phase2 score: $S ($(date))"
