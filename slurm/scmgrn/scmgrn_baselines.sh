#!/bin/bash
# Statistical / simple baselines on the scMultiomeGRN benchmark (scripts/scmgrn/22_scmgrn_simple_baselines.py). CPU only.
# One array task per (dataset, cell type), all families A1..A4 = 12 tasks (Ada's 20-submitted-jobs cap counts
# every array element); small cell types first, Stromal last.
# Reads and writes $SHARE (~/scmgrn, visible on the nodes) directly; existing fold JSONs are skipped, so a
# resubmission resumes. SMOKE=1 runs fold 1 only (the full run later skips those identical files).
#SBATCH --job-name=smg-base
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH -n 1
#SBATCH --cpus-per-task=5
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --array=0-11%2
#SBATCH --output=logs/smg_base_%A_%a.out
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn/scmgrn_common.sh"
activate
export OMP_NUM_THREADS=5 OPENBLAS_NUM_THREADS=5 MKL_NUM_THREADS=5
UNITS=()
for ct in Myeloid_cells Lymphoid_cells Megakaryocytes Neuroendocrine_cells Lymphatic_endothelial_cells; do UNITS+=("lung $ct"); done
for ct in HSC CMP GMP; do UNITS+=("hema $ct"); done
for ct in Bronchiolar_and_alveolar_epithelial_cells Ciliated_epithelial_cells Vascular_endothelial_cells Stromal_cells; do UNITS+=("lung $ct"); done
i=${SLURM_ARRAY_TASK_ID:-0}
read -r ds ct <<< "${UNITS[$i]}"; fam=${FAM:-A1,A2,A3,A4}
extra=(); [ "${SMOKE:-0}" = 1 ] && extra=(--folds 1)
echo "[base] task $i: $ds/$ct family $fam ${extra[*]}"
for sk in ${SPLIT_KINDS:-${SPLIT_KIND:-official}}; do     # e.g. SPLIT_KINDS="official tfholdout"
  python -u scripts/scmgrn/22_scmgrn_simple_baselines.py --root "$SHARE" --ds "$ds" --ct "$ct" --family "$fam" \
    --n_jobs "${SLURM_CPUS_PER_TASK:-5}" --skip_existing --split_kind "$sk" "${extra[@]}" || { echo "[base] FAILED $ds/$ct $fam $sk"; exit 1; }
done
echo "[base] done $ds/$ct $fam"
