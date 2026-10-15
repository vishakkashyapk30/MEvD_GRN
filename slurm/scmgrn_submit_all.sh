#!/bin/bash
# One-command launch of the scMultiomeGRN-benchmark pipeline on Ada (run from the repo root on
# ada-gw2; this only calls sbatch/squeue, no compute):
#   bash slurm/scmgrn_submit_all.sh                 # setup -> prep -> build (x3) -> work (x4) = 9 jobs
#   SKIP_SETUP=1 bash slurm/scmgrn_submit_all.sh    # tools env + raw data already in place
#   SKIP_PREP=1 / SKIP_BUILD=1 ...                  # sep_data / graphs already in $HOME/scmgrn
# All GPU work uses RTX 2080 Ti nodes (--constraint=2080ti in each script). Medium QoS: 4 GPUs,
# 40 CPUs, 20 submitted jobs per user. Nothing here ever runs on the login node.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs
DEP=""
if [ -n "${RUN_SETUP:-}" ]; then   # default OFF: Ada's own internet is ~10 KB/s, so setup was done from the laptop
  J=$(sbatch --parsable slurm/scmgrn_setup.sh); echo "setup job: $J"; DEP="--dependency=afterok:$J"
fi
if [ -z "${SKIP_PREP:-}" ]; then
  P=$(sbatch --parsable $DEP slurm/scmgrn_prep.sh); echo "prep  job: $P"; DEP="--dependency=afterok:$P"
fi
if [ -z "${SKIP_BUILD:-}" ]; then
  B=$(sbatch --parsable $DEP slurm/scmgrn_build.sh); echo "build job: $B (array 0-2)"; DEP="--dependency=afterok:$B"
fi
# no dependency: the GPU workers wait for built cell types themselves
W=$(sbatch --parsable slurm/scmgrn_work.sh); echo "work  job: $W (array 0-3)"
echo "monitor: squeue -u $USER ; tail -f logs/smg_*"
echo "when done (as a batch job): sbatch --wrap 'python scripts/19_scmgrn_compile.py --config configs/scmgrn/lung.yaml --root \$HOME/scmgrn' -p u22 -A research --qos=medium -t 00:30:00 -n 1"
