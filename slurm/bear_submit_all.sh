#!/bin/bash
# One-command launch of the BEAR-GRN benchmark on Ada (run on ada-gw1 from the repo root,
# after the data are on /share1: bash slurm/bear_fetch_gateway.sh, or
# python3 scripts/20_bear_download.py --root /share1/$USER/mevd_grn/bear --datasets ... --what gt input inferred).
#   bash slurm/bear_submit_all.sh                    # prep -> work[0-3]%2 -> score  (6 job slots)
#   DATASETS="iPS" SKIP_PREP= bash slurm/bear_submit_all.sh
#   SKIP_PREP=1 bash slurm/bear_submit_all.sh         # processed dirs already on /share1
# GPU tasks: --constraint=2080ti, 1 GPU : 10 CPUs. PBMC jobs have GPU priority, so
# the work array runs at most 2 tasks at once (override: WORK_CONC=4).
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs
export DATASETS="${DATASETS:-K562 Macrophage_S1 Macrophage_S2}"
DEP=""
if [ -z "${SKIP_PREP:-}" ]; then
  P=$(sbatch --parsable --export=ALL slurm/bear_prep.sh); echo "prep  job: $P"; DEP="--dependency=afterok:$P"
fi
W=$(sbatch --parsable --export=ALL $DEP --array=0-3%${WORK_CONC:-2} slurm/bear_work.sh); echo "work  job: $W (array 0-3)"
S=$(sbatch --parsable --export=ALL --dependency=afterany:$W slurm/bear_score.sh); echo "score job: $S"
echo "monitor: squeue -u $USER ; tail -f logs/bear_*"
