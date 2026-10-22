#!/bin/bash
# Smoke test for the scMultiomeGRN ablations (Myeloid, fold 1, 30 epochs) and the MeVD-GRN degree-matched re-run.
# Builds the patched official copy into $SHARE (scripts/23), then runs every tag on node scratch only: nothing
# is published, so the short runs can never mix with real results. Prints one metric line per tag; the "full"
# tag must equal the unpatched code ("orig") on the same seed and epochs.
#SBATCH --job-name=smg-offsmoke
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/smg_offsmoke_%j.out
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn_common.sh"
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn_worker_lib.sh"
activate
DS=lung; CT=${CT:-Myeloid_cells}
python -u scripts/23_patch_official_ablation.py --src "$SHARE/official/ScmultiomeGRN-main" \
  --dst "$SHARE/official_ablation/ScmultiomeGRN-main" || exit 1
pull official/ScmultiomeGRN-main; pull official_ablation/ScmultiomeGRN-main
pull "graph/$DS/$CT"; pull "splits/$DS/$CT"; pull "processed/$DS/$CT"
R="$SCMGRN_ROOT_LOCAL"
run() {  # run <tag> <repo> <json>
  python -u scripts/18_scmgrn_official_repro.py --config configs/scmgrn/$DS.yaml --root "$R" --official_repo "$2" \
    --cell_type "$CT" --folds 1 --max_epochs 30 --model_seed 666 --ablation "smoke_$1" --extra_kw "$3" --save_scores \
    --work_dir "$SCR/w_$1" 2>&1 | grep -E "^\[official|edge_attr|Error|error|Traceback" | head -8
  grep -h "edge:" "$SCR/w_$1"/fold01/stdout.log 2>/dev/null | head -1
}
run orig official/ScmultiomeGRN-main '{}'
for T in "${OFFABL_TAGS[@]}"; do echo "=== $T"; run "$T" official_ablation/ScmultiomeGRN-main "$(offabl_kw "$T")"; done
python - "$R" "$CT" <<'PY'
import json, sys, glob
R, ct = sys.argv[1:]
for f in sorted(glob.glob(f"{R}/results/lung/{ct}/official_scmultiomegrn_smoke_*/seed666/fold01.json")):
    m = json.load(open(f))["metrics"]
    print(f"{f.split('/')[-3]:45s} test AUROC {m['test']['auroc']:.4f} AUPR {m['test']['aupr']:.4f} | degmatched AUROC {m['test_degmatched']['auroc']:.4f}")
PY
echo "=== mevddm (2 epochs)"
python -u scripts/17_scmgrn_train_eval.py --config configs/scmgrn/$DS.yaml --model_config configs/scmgrn/mevd_fm_h384.yaml \
  --obs_graph --variant smoke_dm --save_scores --root "$R" --cell_type "$CT" --seed 42 --folds 1 --max_epochs 2 --tmp_dir "$SCR" 2>&1 | tail -2
python -c "import json;m=json.load(open('$R/results/lung/$CT/mevd_smoke_dm/seed42/fold01.json'))['metrics'];print('mevd dm', m['test_degmatched'])"
ls "$R/results/lung/$CT/mevd_smoke_dm/seed42/"
echo "[smoke] done"
