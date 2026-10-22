#!/bin/bash
# Compile result tables in a short CPU job (no compute on the login node).
#   sbatch --export=ALL,REF=official_scmultiomegrn_full,TAG=vs_smg slurm/scmgrn_compile.sh
#SBATCH --job-name=smg-compile
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=logs/smg_compile_%j.out
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn_common.sh"
activate
mkdir -p "$SHARE/tables"
for DSS in lung,hema lung hema; do
  python scripts/25_scmgrn_baseline_table.py --root "$SHARE" --ref "${REF:-official_scmultiomegrn_full}" \
    --datasets "$DSS" --out "$SHARE/tables/${TAG:-vs_smg}_${DSS/,/_}" > /dev/null
done
cat "$SHARE/tables/${TAG:-vs_smg}_lung_hema.md"
