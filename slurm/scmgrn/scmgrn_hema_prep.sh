#!/bin/bash
# Hematopoiesis (Supp. Fig. S1) inputs, as a batch job (compute nodes have internet; downloads and
# bedtools counting must not run on the login node). Writes $HOME/scmgrn/sep_data/hema/{HSC,CMP,GMP}.
# Then run the pipeline with DS=hema:  DS=hema SKIP_SETUP=1 SKIP_PREP=1 bash slurm/scmgrn/scmgrn_submit_all.sh
#SBATCH --job-name=smg-hema-prep
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3000
#SBATCH --time=12:00:00
#SBATCH --output=logs/smg_hema_prep_%j.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn/scmgrn_common.sh"
activate
source "$CONDA_SH"
conda run -n "$TOOLS_ENV" which bedtools >/dev/null 2>&1 || conda install -y -n "$TOOLS_ENV" -c conda-forge -c bioconda bedtools
export PATH="$(dirname "$(conda run -n "$TOOLS_ENV" which bedtools)"):$PATH"
W="$SCR/hema"; mkdir -p "$W"; cd "$W"
G=https://raw.githubusercontent.com/pinellolab/scATAC-benchmarking/master/Real_Data/Buenrostro_2018/input
dl() { curl -fL -C - --retry 8 --retry-delay 10 -o "$2" "$1"; }
dl "$G/metadata.tsv" metadata.tsv
dl "$G/combined.sorted.merged.bed" combined.sorted.merged.bed
dl "https://ars.els-cdn.com/content/image/1-s2.0-S009286741830446X-mmc4.zip" mmc4.tar.gz   # a tar.gz despite the name
dl "https://www.dropbox.com/sh/8o8f0xu6cvr46sm/AAB6FMIDvHqnG6h7athgcm5-a/Buenrostro_2018.tar.gz?dl=1" bams.tar.gz
mkdir -p bams && tar xzf bams.tar.gz -C bams
BD=$(dirname "$(find bams -name '*.bam' | head -1)")
echo "[hema] bam dir: $BD ($(ls "$BD"/*.bam | wc -l) BAMs)"
cd "$REPO"
python scripts/scmgrn/14b_scmgrn_hema_prepare.py rna --mmc4 "$W/mmc4.tar.gz" --out "$SCMGRN_ROOT_LOCAL/sep_data/hema"
python scripts/scmgrn/14b_scmgrn_hema_prepare.py atac --bam_dir "$W/$BD" --metadata "$W/metadata.tsv" \
   --peaks "$W/combined.sorted.merged.bed" --out "$SCMGRN_ROOT_LOCAL/sep_data/hema" --n_jobs "${SLURM_CPUS_PER_TASK:-10}"
push sep_data/hema
echo "[hema-prep] done"
