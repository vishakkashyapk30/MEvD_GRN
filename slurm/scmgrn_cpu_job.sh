#!/bin/bash
# Run one Python command of the scMultiomeGRN track as a CPU-only batch job (no compute on the login node).
#   sbatch slurm/scmgrn_cpu_job.sh scripts/30_scmgrn_make_tfholdout_splits.py --root '$SHARE'
# '$SHARE' (single-quoted at submission) is expanded here to the benchmark root.
#SBATCH --job-name=smg-cpu
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=06:00:00
#SBATCH --output=logs/smg_cpu_%j.out
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn_common.sh"
activate
args=(); for x in "$@"; do args+=("${x//\$SHARE/$SHARE}"); done
echo "[cpu] python -u ${args[*]}"
python -u "${args[@]}" || { echo "[cpu] FAILED"; exit 1; }
echo "[cpu] done"
