#!/bin/bash
# Smoke test of the held-out-TF split (Myeloid, fold 1) for MeVD-GRN (2 epochs), the official scMultiomeGRN
# code (30 epochs) and the statistical baselines. Runs on node scratch only; nothing is published.
#SBATCH --job-name=smg-tfsmoke
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/smg_tfsmoke_%j.out
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn/scmgrn_common.sh"
activate
export HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
DS=lung; CT=${CT:-Myeloid_cells}; R="$SCMGRN_ROOT_LOCAL"
pull official/ScmultiomeGRN-main; pull "graph/$DS/$CT"; pull "splits/$DS/$CT"; pull "processed/$DS/$CT"
echo "=== MeVD-GRN"
python -u scripts/scmgrn/17_scmgrn_train_eval.py --config configs/scmgrn/$DS.yaml --model_config configs/scmgrn/mevd_fm_h384.yaml \
  --obs_graph --root "$R" --cell_type "$CT" --seed 42 --folds 1 --max_epochs 2 --tmp_dir "$SCR" --split_kind tfholdout 2>&1 | tail -2
echo "=== official scMultiomeGRN"
python -u scripts/scmgrn/18_scmgrn_official_repro.py --config configs/scmgrn/$DS.yaml --root "$R" --cell_type "$CT" --folds 1 \
  --max_epochs 30 --model_seed 666 --split_kind tfholdout --save_scores --work_dir "$SCR/w" 2>&1 | grep -E "^\[official|Error|Traceback" | head -5
echo "=== baselines"
python -u scripts/scmgrn/22_scmgrn_simple_baselines.py --root "$R" --ds $DS --ct "$CT" --folds 1 --split_kind tfholdout --n_jobs 2 \
  2>&1 | grep -E "heur_degree_product|lr_degree_only|unsup_pearson|lr_geneformer:|Traceback|Error" | head -8
python - "$R" "$CT" <<'PY'
import json, glob, sys
R, ct = sys.argv[1:]
for f in sorted(glob.glob(f"{R}/results/lung_tfholdout/{ct}/*/seed*/fold01.json")):
    m = json.load(open(f))["metrics"]
    g = lambda k: f"{m[k]['auroc']:.3f}" if k in m else "  -  "
    print(f"{f.split('/')[-3]:40s} test {g('test')}  C2 {g('test_c2')}  C3 {g('test_c3')}  C2-degmatched {g('test_degmatched')}")
PY
echo "[smoke] done"
