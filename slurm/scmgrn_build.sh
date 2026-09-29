#!/bin/bash
# Step 2: per cell type, regenerate scMultiomeGRN's ground truth (FIMO), both
# methods' node features, the official 10-fold splits (verified byte-identical
# against the official split code), and MeVD-GRN's processed inputs.
# 3 array tasks x 10 CPUs; each task builds every 3rd lung cell type (CPU work;
# the GPU is requested only to respect Ada's 1 GPU : 10 CPU policy).
#SBATCH --job-name=smg-build
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=long
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3000
#SBATCH --time=2-00:00:00
#SBATCH --array=0-2
#SBATCH --output=logs/smg_build_%A_%a.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn_common.sh"
activate
NSH=${SLURM_ARRAY_TASK_COUNT:-3}
pull official/ScmultiomeGRN-main
pull raw/maestro
pull raw/genome/hg19.fa
[ -e "$SHARE/raw/genome/hg19.fa.fai" ] && pull raw/genome/hg19.fa.fai || true
for i in "${!LUNG_CTS[@]}"; do
  (( i % NSH == SLURM_ARRAY_TASK_ID )) || continue
  CT=${LUNG_CTS[$i]}
  if [ -e "$SHARE/processed/lung/$CT/summary.json" ] && [ -e "$SHARE/graph/lung/$CT/x_scRNA_grnboost2_feature.txt" ]; then
    echo "[build] $CT already built, skipping"; continue; fi
  pull "sep_data/lung/$CT"
  pull "graph/lung/$CT"          # resumes a partial FIMO scan
  python -u scripts/16_scmgrn_build_benchmark.py --config configs/scmgrn/lung.yaml \
      --root "$SCMGRN_ROOT_LOCAL" --cell_type "$CT" --fimo_bin "$FIMO_BIN" \
      --n_jobs "${SLURM_CPUS_PER_TASK:-10}" --steps gt,official,splits,verify,processed
  rm -rf "$SCMGRN_ROOT_LOCAL/graph/lung/$CT/work/peak.fasta"
  push "graph/lung/$CT"; push "splits/lung/$CT"; push "processed/lung/$CT"
  rm -rf "$SCMGRN_ROOT_LOCAL/sep_data/lung/$CT"
done
push raw/genome/hg19.fa.fai
echo "[build] task $SLURM_ARRAY_TASK_ID done"
