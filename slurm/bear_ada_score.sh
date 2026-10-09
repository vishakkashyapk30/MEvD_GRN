#!/bin/bash
# CPU scoring of MeVD-GRN results on Ada (laptop disk is too small to pull all GRNs back).
# Reads $B/results/<ds>/mevd_*/seed*/grn.tsv.gz + GTs under $B, writes per-(method,GT) JSONs to
# $B/results/<ds>/scores (cached; n_rep 20). Pull the JSONs back (tiny) and compile locally.
#   DATASETS="K562 Macrophage_S1" sbatch slurm/bear_ada_score.sh
#SBATCH --job-name=bear-score
#SBATCH -p u22
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6000
#SBATCH --time=0-10:00:00
#SBATCH --output=logs/bear_score_%j.out
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/bear_ada_common.sh"
activate
DATASETS=(${DATASETS:-K562 Macrophage_S1 Macrophage_S2 iPS mESC_E7.5_rep1 mESC_E7.5_rep2 mESC_E8.5_rep1 mESC_E8.5_rep2 Naive_mESC})
for DS in "${DATASETS[@]}"; do
  ls "$B"/results/$DS/mevd_*/seed*/grn.tsv.gz >/dev/null 2>&1 || { echo "[score] $DS: no MeVD results yet"; continue; }
  python -u scripts/24_bear_score.py --root "$B" --datasets "$DS" --ours "${OURS:-all}" --n_rep "${NREP:-20}" --sparse_match
done
echo BEAR_ADA_SCORE_DONE
