#!/bin/bash
# Step 2: per cell type, regenerate scMultiomeGRN's ground truth (FIMO), both
# methods' node features, the official 10-fold splits (verified byte-identical
# against the official split code), and MeVD-GRN's processed inputs.
# 3 CPU-only array tasks x 9 CPUs (no GPU: the 4 GPUs stay free for the work queue, and medium QoS
# caps CPUs at 40 = 3x9 + 4x3); each task builds every 3rd cell type, small graphs first.
#SBATCH --job-name=smg-build
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH -n 1
#SBATCH --cpus-per-task=9
#SBATCH --mem-per-cpu=3000
#SBATCH --time=2-00:00:00
#SBATCH --array=0-2
#SBATCH --output=logs/smg_build_%A_%a.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn_common.sh"
activate
export WAIT_RAW=1
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1     # Geneformer cached by scmgrn_setup.sh
NSH=${SLURM_ARRAY_TASK_COUNT:-3}
pull official/ScmultiomeGRN-main
pull raw/maestro
pull raw/fm_tf_table.npz
fetch_raw hg19.fa.gz "$SCMGRN_ROOT_LOCAL/raw/genome/hg19.fa.gz"
gunzip -f "$SCMGRN_ROOT_LOCAL/raw/genome/hg19.fa.gz"
for i in "${!CTS[@]}"; do
  (( i % NSH == SLURM_ARRAY_TASK_ID )) || continue
  CT=${CTS[$i]}
  if [ -e "$SHARE/processed/$DS/$CT/smg_pair_feats.npy" ] && [ -e "$SHARE/graph/$DS/$CT/x_scRNA_grnboost2_feature.txt" ]; then
    echo "[build] $CT already built, skipping"; continue; fi
  pull "sep_data/$DS/$CT"
  pull "graph/$DS/$CT"          # resumes a partial FIMO scan
  python -u scripts/16_scmgrn_build_benchmark.py --config configs/scmgrn/$DS.yaml \
      --root "$SCMGRN_ROOT_LOCAL" --cell_type "$CT" --fimo_bin "$FIMO_BIN" \
      --n_jobs "${SLURM_CPUS_PER_TASK:-10}" --steps gt,official,splits,verify,processed,smg,fm
  rm -rf "$SCMGRN_ROOT_LOCAL/graph/$DS/$CT/work/peak.fasta" "$SCMGRN_ROOT_LOCAL/graph/$DS/$CT/work/fimo-res"
  push "graph/$DS/$CT"; push "splits/$DS/$CT"; push "processed/$DS/$CT"
  rm -rf "$SCMGRN_ROOT_LOCAL/sep_data/$DS/$CT"
done
echo "[build] task $SLURM_ARRAY_TASK_ID done"
