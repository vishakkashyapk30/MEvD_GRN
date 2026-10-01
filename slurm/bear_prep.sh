#!/bin/bash
# BEAR step 1+3 (CPU): label-free processed dirs (+ Geneformer embeddings),
# unsupervised/trivial baselines, and scoring of the 9 released methods + baselines.
#SBATCH --job-name=bear-prep
#SBATCH -p u22
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3000
#SBATCH --time=0-08:00:00
#SBATCH --output=logs/bear_prep_%j.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/bear_common.sh"
activate
DATASETS=(${DATASETS:-${BEAR_PHASE1[@]}})
for DS in "${DATASETS[@]}"; do
  echo "===== $DS ($(date +%T)) ====="
  stage_ds "$DS"
  python -u scripts/21_bear_prepare.py --root "$BEAR_ROOT" --datasets "$DS" --fm
  python -u scripts/21_bear_prepare.py --root "$BEAR_ROOT" --datasets "$DS" --fm --regime L2
  push "processed/$DS"; push "processed/${DS}__L2"
  python -u scripts/23_bear_baselines.py --root "$BEAR_ROOT" --dataset "$DS" \
      --baselines indegree pearson coverage grnboost2 --n_jobs 10
  push "results/$DS"
  pull "INFERRED.GRNS/$(python -c "import yaml;print(yaml.safe_load(open('configs/bear/datasets.yaml'))['inferred_dir']['$DS'])")"
  python -u scripts/24_bear_score.py --root "$BEAR_ROOT" --datasets "$DS" --released all --ours all \
      --n_rep 1 --paper_csv "$REPO/configs/bear/paper_numbers/bear_grn_auroc_auprc_tidy.csv"
  push "results/$DS"
  rm -rf "$BEAR_ROOT/INFERRED.GRNS"
done
push "results/summary_scores.csv"; push "results/port_validation.csv"
echo "BEAR_PREP_DONE"
