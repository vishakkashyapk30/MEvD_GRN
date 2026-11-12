#!/bin/bash
# CPU smoke test of the DeepTFni (official Trainer, model_cls=DeepTFni) and GENELink (scripts/32) baselines on
# Myeloid fold 1, official and held-out-TF splits, few epochs, node scratch only (nothing published).
#SBATCH --job-name=smg-gnnsmoke
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/smg_gnnsmoke_%j.out
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn_common.sh"
activate
DS=lung; CT=Myeloid_cells; R="$SCMGRN_ROOT_LOCAL"
pull official/ScmultiomeGRN-main; pull official_ablation/ScmultiomeGRN-main
pull "graph/$DS/$CT"; pull "splits/$DS/$CT"
python -c "import scanpy; print('scanpy', scanpy.__version__)"
for SK in official tfholdout; do
  echo "=== DeepTFni $SK"
  python -u scripts/18_scmgrn_official_repro.py --config configs/scmgrn/$DS.yaml --root "$R" \
    --official_repo official_ablation/ScmultiomeGRN-main --cell_type "$CT" --folds 1 --max_epochs 20 --model_seed 666 \
    --method_name deeptfni --extra_kw '{"model_cls": "DeepTFni", "use_edge_attr": false}' --split_kind $SK \
    --work_dir "$SCR/d_$SK" 2>&1 | grep -E "^\[official|Error|Traceback|error" | head -6
  echo "=== GENELink $SK"
  python -u scripts/32_scmgrn_genelink.py --config configs/scmgrn/$DS.yaml --root "$R" --cell_type "$CT" --folds 1 \
    --seed 8 --split_kind $SK --work_dir "$SCR/g_$SK" 2>&1 | tail -4
done
echo "[smoke] done"
