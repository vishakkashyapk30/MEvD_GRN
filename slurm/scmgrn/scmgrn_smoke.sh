#!/bin/bash
# End-to-end smoke test of the scMultiomeGRN pipeline on Ada with the official demo cell type
# (PBMC B_naive, hg38, restricted to chr19+chr22 so no whole genome is needed). CPU-only, ~1 h.
# Exercises: FIMO ground truth, MAESTRO RP, GRNBoost2, official splits + byte-identity check, smg features,
# Geneformer table, MeVD-GRN (H and A3) training, the OFFICIAL scMultiomeGRN code (lightning), and compile.
#SBATCH --job-name=smg-smoke
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3000
#SBATCH --time=03:00:00
#SBATCH --output=logs/smg_smoke_%j.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn/scmgrn_common.sh"
activate
R="$SCMGRN_ROOT_LOCAL"
pull official/ScmultiomeGRN-main; pull raw/maestro; pull raw/fm_tf_table.npz
mkdir -p "$R/raw/genome"
zcat "$SHARE/raw/genome/chr19.fa.gz" "$SHARE/raw/genome/chr22.fa.gz" > "$R/raw/genome/hg38_subset.fa"
CFG=archive/configs/pbmc_demo.yaml; CT=B_naive
python -u scripts/scmgrn/16_scmgrn_build_benchmark.py --config $CFG --root "$R" --cell_type $CT --fimo_bin "$FIMO_BIN" \
   --n_jobs 4 --chroms chr19,chr22 --steps gt,official,splits,verify,processed,smg,fm
cat "$R/graph/pbmc_demo/$CT/build_summary.json"
for spec in "mevd_fm_h384:--obs_graph --smg_feats" "mevd_fm_h384:--obs_graph" "mevd_fm_h384:" "mevd_fm_h384_directed:--obs_graph"; do
  MC=${spec%%:*}; FL=${spec#*:}
  python -u scripts/scmgrn/17_scmgrn_train_eval.py --config $CFG --model_config configs/scmgrn/$MC.yaml $FL --root "$R" \
     --cell_type $CT --seed 42 --folds 1-3 --device cpu --tmp_dir "$SCR" --max_epochs 60
done
python -u scripts/scmgrn/18_scmgrn_official_repro.py --config $CFG --root "$R" --cell_type $CT --folds 1-2 \
   --max_epochs 40 --work_dir "$SCR/official_work" --model_seed 666
python -u scripts/scmgrn/19_scmgrn_compile.py --config $CFG --root "$R" --n_folds 2 --head mevd_fm_h384_obs_smg || true
mkdir -p "$SHARE/smoke"; rsync -a "$R/results" "$R/graph/pbmc_demo/$CT/build_summary.json" "$R/graph/pbmc_demo/$CT/split_verification.json" "$SHARE/smoke/"
echo "[smoke] done"
