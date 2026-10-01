#!/bin/bash
# BEAR step 4 (CPU): score everything in results/<ds> against every GT and compile.
#SBATCH --job-name=bear-score
#SBATCH -p u22
#SBATCH -A research
#SBATCH --qos=medium
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6000
#SBATCH --time=0-06:00:00
#SBATCH --output=logs/bear_score_%j.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/bear_common.sh"
activate
DATASETS=(${DATASETS:-${BEAR_PHASE1[@]}})
for DS in "${DATASETS[@]}"; do
  pull "results/$DS"; for rel in $(python - "$DS" <<'PY'
import sys, yaml
d = yaml.safe_load(open("configs/bear/datasets.yaml"))["datasets"][sys.argv[1]]
print(" ".join(sorted(set(d["gts"].values()))))
PY
); do pull "$rel"; done
  python -u scripts/24_bear_score.py --root "$BEAR_ROOT" --datasets "$DS" --ours all --n_rep "${NREP:-20}" \
      --sparse_match
  push "results/$DS/scores"
done
pull "results/summary_scores.csv" || true
python -u scripts/24_bear_score.py --root "$BEAR_ROOT" --datasets "${DATASETS[@]}" --compile_only \
    --paper_csv "$REPO/configs/bear/paper_numbers/bear_grn_auroc_auprc_tidy.csv"
push "results/summary_scores.csv"; push "results/summary_tables.md"
echo "BEAR_SCORE_DONE"
