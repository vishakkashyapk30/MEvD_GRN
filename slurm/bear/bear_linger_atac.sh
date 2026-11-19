#!/bin/bash
# LINGER re-run on one BEAR-GRN dataset under three ATAC conditions (full, peakshuffle, random), run in parallel on
# one node so LINGER's 20.8 GB bulk reference is staged once (streamed from /share1 on the gateway into /scratch).
# Results -> $HOME/bear_v2/results/<DS>/linger_rerun_<mode>/seed0/. Advisor item 3 (2026-11-19).
#   sbatch --export=ALL,DS=Macrophage_S1 slurm/bear/bear_linger_atac.sh
#SBATCH --job-name=bear-linger
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=1-12:00:00
#SBATCH --output=logs/bear_linger_%j.out
set -uo pipefail
DS="${DS:-Macrophage_S1}"
MASTER=ada-gw2
SH=/share1/$USER/mevd_grn/bear
SCR=/scratch/$USER/linger_${SLURM_JOB_ID:-local}
trap 'rm -rf "$SCR"' EXIT
mkdir -p "$SCR/in" "$SCR/bulk" "$HOME/bear_v2/logs"
REPO="${SLURM_SUBMIT_DIR:-$PWD}"
source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate linger; export PYTHONUNBUFFERED=1
echo "[env] host=$(hostname) job=${SLURM_JOB_ID:-} gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1) df=$(df -h /scratch | tail -1)"
case "$DS" in
  Macrophage_S1) RNA=subsampled_RNA_Macrophase_buffer1_filtered.csv; ATAC=subsampled_ATAC_Macrophase_buffer1_filtered.csv;;
  Macrophage_S2) RNA=subsampled_RNA_Macrophase_buffer2_filtered.csv; ATAC=subsampled_ATAC_Macrophase_buffer2_filtered.csv;;
  *) echo "add the file names for $DS"; exit 1;;
esac
rsync -a "$MASTER:$SH/INPUT.DATA/$DS/$RNA" "$MASTER:$SH/INPUT.DATA/$DS/$ATAC" "$SCR/in/" || { echo "input stage failed"; exit 1; }
ssh "$MASTER" "cat /share1/$USER/mevd_grn/pbmc/linger/data_bulk.tar.gz" | tar -xz -C "$SCR/bulk" || { echo "bulk stage failed"; exit 1; }
GRN="$SCR/bulk/data_bulk/"                         # the archive unpacks to data_bulk/ (as in slurm/pbmc/pbmc_linger.sh)
echo "[stage] grndir=$GRN $(du -sh "$SCR/bulk" | cut -f1)"
cd "$REPO"
for M in full peakshuffle random; do
  python -u src/benchmarks/bear_linger_run.py --rna "$SCR/in/$RNA" --atac "$SCR/in/$ATAC" --grndir "$GRN" \
    --outdir "$SCR/work_$M/" --export_dir "$SCR/out/linger_rerun_$M/seed0" --mode "$M" > "$HOME/bear_v2/logs/linger_${SLURM_JOB_ID:-local}_${DS}_$M.txt" 2>&1 &
done
wait
for M in full peakshuffle random; do                 # retry alone any condition that failed in parallel (e.g. memory)
  [ -f "$SCR/out/linger_rerun_$M/seed0/grn.tsv.gz" ] && continue
  echo "[retry] $M sequentially"; rm -rf "$SCR/work_$M"
  python -u src/benchmarks/bear_linger_run.py --rna "$SCR/in/$RNA" --atac "$SCR/in/$ATAC" --grndir "$GRN" \
    --outdir "$SCR/work_$M/" --export_dir "$SCR/out/linger_rerun_$M/seed0" --mode "$M" >> "$HOME/bear_v2/logs/linger_${SLURM_JOB_ID:-local}_${DS}_$M.txt" 2>&1
done
for M in full peakshuffle random; do echo "== $M"; tail -5 "$HOME/bear_v2/logs/linger_${SLURM_JOB_ID:-local}_${DS}_$M.txt"; done
mkdir -p "$HOME/bear_v2/results/$DS"
for M in full peakshuffle random; do
  [ -f "$SCR/out/linger_rerun_$M/seed0/grn.tsv.gz" ] && rsync -a "$SCR/out/linger_rerun_$M" "$HOME/bear_v2/results/$DS/" && echo "published $M"
done
