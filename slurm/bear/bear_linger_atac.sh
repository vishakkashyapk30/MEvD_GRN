#!/bin/bash
# LINGER re-run on one BEAR-GRN dataset under the ATAC conditions in $MODES (full, peakshuffle, random), one after
# another; LINGER's 20.8 GB bulk reference is streamed from /share1 on the gateway into /scratch once per job.
# LINGER trains on the CPU (torch.cuda is only used for seeding), so no GPU is requested.
# Results -> $HOME/bear_v2/results/<DS>/linger_rerun_<mode>/seed<SEED>/. Advisor item 3 (2026-11-19, 2026-11-20).
#   sbatch --export=ALL,DS=Macrophage_S1 slurm/bear/bear_linger_atac.sh                      # 3 conditions, seed 0
#   sbatch --mem=40G --export=ALL,DS=iPS,MODES=full,SEED=0 slurm/bear/bear_linger_atac.sh    # one condition
#SBATCH --job-name=bear-linger
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=27G
#SBATCH --time=1-12:00:00
#SBATCH --output=logs/bear_linger_%j.out
set -uo pipefail
DS="${DS:-Macrophage_S1}"
SEED="${SEED:-0}"
MASTER=ada-gw2
SH=/share1/$USER/mevd_grn/bear
SCR=/scratch/$USER/linger_${SLURM_JOB_ID:-local}
trap 'rm -rf "$SCR"' EXIT
mkdir -p "$SCR/in" "$SCR/bulk" "$HOME/bear_v2/logs"
REPO="${SLURM_SUBMIT_DIR:-$PWD}"
source "$HOME/miniforge3/etc/profile.d/conda.sh"; conda activate linger; export PYTHONUNBUFFERED=1
echo "[env] host=$(hostname) job=${SLURM_JOB_ID:-} DS=$DS SEED=$SEED MODES=${MODES:-all} df=$(df -h /scratch | tail -1)"
case "$DS" in
  Macrophage_S1) RNA=subsampled_RNA_Macrophase_buffer1_filtered.csv; ATAC=subsampled_ATAC_Macrophase_buffer1_filtered.csv;;
  Macrophage_S2) RNA=subsampled_RNA_Macrophase_buffer2_filtered.csv; ATAC=subsampled_ATAC_Macrophase_buffer2_filtered.csv;;
  K562) RNA=subsampled_RNA_K562_human_filtered.csv; ATAC=subsampled_ATAC_K562_human_filtered.csv;;
  iPS) RNA=filtered_multiomics_common_RNA_iPS_cell_L2.csv; ATAC=filtered_multiomics_common_ATAC_iPS_cell_L2.csv;;
  *) echo "add the file names for $DS"; exit 1;;
esac
rsync -a "$MASTER:$SH/INPUT.DATA/$DS/$RNA" "$MASTER:$SH/INPUT.DATA/$DS/$ATAC" "$SCR/in/" || { echo "input stage failed"; exit 1; }
ssh "$MASTER" "cat /share1/$USER/mevd_grn/pbmc/linger/data_bulk.tar.gz" | tar -xz -C "$SCR/bulk" || { echo "bulk stage failed"; exit 1; }
GRN="$SCR/bulk/data_bulk/"                         # the archive unpacks to data_bulk/ (as in slurm/pbmc/pbmc_linger.sh)
echo "[stage] grndir=$GRN $(du -sh "$SCR/bulk" | cut -f1)"
cd "$REPO"
# LingerGRN 1.110 bug: TF_RE_LINGER_chr uses the loop variable `ii` after `for ii in range(times)`, which is unset when a
# chromosome has fewer than 50 genes (times = 0; hit on iPS, 4,808 genes). A patched copy goes first on PYTHONPATH (the
# loky workers re-import the package, so an in-process monkeypatch would not reach them). `ii = times` equals the
# original `ii + 1` whenever the loop ran, so results are unchanged otherwise.
LP=$(python -c "import LingerGRN, os; print(os.path.dirname(LingerGRN.__path__[0]))")
mkdir -p "$SCR/pkgs" && cp -r "$LP/LingerGRN" "$SCR/pkgs/"
sed -i 's/^    ii=ii+1$/    ii=times/' "$SCR/pkgs/LingerGRN/LL_net.py"
grep -q "^    ii=times$" "$SCR/pkgs/LingerGRN/LL_net.py" || { echo "LINGER patch failed"; exit 1; }
export PYTHONPATH="$SCR/pkgs${PYTHONPATH:+:$PYTHONPATH}"
# LINGER sizes its worker pools from the NODE's free RAM (5 GB per worker) and SLURM_CPUS_PER_TASK, not from the job's
# memory limit: 8 workers OOM-killed a 27 GB job (93598). Cap the workers and run the conditions one after another,
# each with its own temp dir (pybedtools cleans up shared /tmp files of other processes). Each condition is published
# as soon as it finishes.
export SLURM_CPUS_PER_TASK="${LINGER_CPUS:-2}"
mkdir -p "$HOME/bear_v2/results/$DS"
for M in ${MODES:-full random peakshuffle}; do
  LOG="$HOME/bear_v2/logs/linger_${SLURM_JOB_ID:-local}_${DS}_${M}_seed$SEED.txt"
  mkdir -p "$SCR/tmp_$M"
  echo "[run] $M $(date +%H:%M) workers=$SLURM_CPUS_PER_TASK"
  TMPDIR="$SCR/tmp_$M" python -u src/benchmarks/bear_linger_run.py --rna "$SCR/in/$RNA" --atac "$SCR/in/$ATAC" \
    --grndir "$GRN" --outdir "$SCR/work_$M/" --export_dir "$SCR/out/linger_rerun_$M/seed$SEED" --mode "$M" \
    --seed "$SEED" --labels "${LABELS:-auto}" --max_cells "${MAX_CELLS:-0}" > "$LOG" 2>&1
  tail -3 "$LOG"
  [ -f "$SCR/out/linger_rerun_$M/seed$SEED/grn.tsv.gz" ] && rsync -a "$SCR/out/linger_rerun_$M" "$HOME/bear_v2/results/$DS/" \
    && echo "published $M $(date +%H:%M)"
  rm -rf "$SCR/work_$M" "$SCR/tmp_$M"
done
