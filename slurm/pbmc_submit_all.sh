#!/bin/bash
# One-command launch on the Ada gateway (ada-gw1), from the repo root, after
#   bash scripts/25_pbmc_download.sh data
# Submits: build (1 job) -> work (array 0-3, up to 4 GPUs, work-stealing) [+ GRNBoost2, 1 job].
# Submit-slot use: 1 + 4 (+ 1 with GB2=1) <= 10 (other agents keep the rest; QoS cap 20).
#   bash slurm/pbmc_submit_all.sh                 # build + work
#   SKIP_BUILD=1 bash slurm/pbmc_submit_all.sh    # processed dirs/splits already on /share1
#   GB2=1 bash slurm/pbmc_submit_all.sh           # also GRNBoost2 (waits for build)
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs
DEP=""
if [ -z "${SKIP_BUILD:-}" ]; then
  B=$(sbatch --parsable slurm/pbmc_build.sh); echo "build job: $B"; DEP="--dependency=afterok:$B"
fi
W=$(sbatch --parsable $DEP slurm/pbmc_work.sh); echo "work  job: $W (array 0-3)"
if [ -n "${GB2:-}" ]; then
  G=$(sbatch --parsable $DEP slurm/pbmc_grnboost2.sh); echo "gb2   job: $G"
fi
echo "monitor: squeue -u $USER ; tail -f logs/pbmc_*"
