#!/bin/bash
# One-command launch of the scMultiomeGRN-benchmark pipeline on Ada (run on the
# login node from the repo root, after `bash scripts/14_download_scmultiomegrn_data.sh env|data|fm`):
#   bash slurm/scmgrn_submit_all.sh            # prep -> build (x3) -> work (x4)  = 8 jobs
#   SKIP_PREP=1 bash slurm/scmgrn_submit_all.sh   # sep_data already on /share1
#   SKIP_BUILD=1 ...                              # also skip the build
# Only RTX 2080 Ti nodes are used for the GPU work: the script looks for a
# 2080 gres type / node feature, and otherwise excludes the 1080 Ti nodes
# (gnode01-40). Override with WORK_GPU_OPTS="...".
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs
echo "[ada] $(sinfo --version 2>/dev/null)"; sacctmgr -n show assoc user=$USER format=Account,QOS,DefaultQOS 2>/dev/null | head -5 || true

if [ -z "${WORK_GPU_OPTS:-}" ]; then
  GT=$(sinfo -h -p long -o "%G" | tr ',' '\n' | grep -io 'gpu:[a-z0-9_]*2080[a-z0-9_]*' | head -1 | cut -d: -f2 || true)
  FT=$(sinfo -h -p long -o "%f" | tr ',' '\n' | grep -i '2080' | head -1 || true)
  if [ -n "$GT" ]; then WORK_GPU_OPTS="--gres=gpu:$GT:1"
  elif [ -n "$FT" ]; then WORK_GPU_OPTS="--gres=gpu:1 --constraint=$FT"
  else
    EXCL=$(sinfo -h -N -p long -o "%N" | sort -u | awk 'match($0,/[0-9]+$/){n=substr($0,RSTART)+0; if(n>=1&&n<=40) print}' | paste -sd, -)
    WORK_GPU_OPTS="--gres=gpu:1${EXCL:+ --exclude=$EXCL}"
  fi
fi
echo "[gpu] work jobs use: $WORK_GPU_OPTS"

DEP=""
if [ -z "${SKIP_PREP:-}" ]; then
  P=$(sbatch --parsable slurm/scmgrn_prep.sh); echo "prep  job: $P"; DEP="--dependency=afterok:$P"
fi
if [ -z "${SKIP_BUILD:-}" ]; then
  B=$(sbatch --parsable $DEP slurm/scmgrn_build.sh); echo "build job: $B (array 0-2)"; DEP="--dependency=afterok:$B"
fi
W=$(sbatch --parsable $DEP $WORK_GPU_OPTS slurm/scmgrn_work.sh); echo "work  job: $W (array 0-3)"
echo "monitor: squeue -u $USER ; tail -f logs/smg_*"
echo "when done: python scripts/19_scmgrn_compile.py --config configs/scmgrn/lung.yaml --root \${SCMGRN_ROOT:-/share1/$USER/mevd_grn/scmgrn}"
