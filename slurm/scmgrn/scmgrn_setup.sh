#!/bin/bash
# Step 0 (batch job; compute nodes have internet, the login node must not run this):
#  * tools env `scmgrn-tools` (MEME 5.4.1 FIMO, R + Matrix) + pip lightning/arboreto in `mevd-grn`
#  * official scMultiomeGRN code (Zenodo) + MAESTRO refgenes -> $HOME/scmgrn (home)
#  * lung RDS.gz (3.01 GB) and hg19.fa.gz -> ada-gw2:/share1/$USER/mevd_grn/scmgrn/raw (2 big files)
#  * Geneformer cached in the HF cache for offline jobs
#SBATCH --job-name=smg-setup
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=3000
#SBATCH --time=06:00:00
#SBATCH --output=logs/smg_setup_%j.out
set -euo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn/scmgrn_common.sh"
source "$CONDA_SH"
dl() { echo ">> $2"; mkdir -p "$(dirname "$2")"; curl -fL -C - --retry 8 --retry-delay 10 -o "$2" "$1"; }
if ! conda env list | grep -q "^$TOOLS_ENV "; then
  conda create -y -n "$TOOLS_ENV" -c conda-forge -c bioconda meme=5.4.1 r-base r-matrix
fi
conda activate mevd-grn
pip install "lightning>=2.0,<2.4" arboreto
echo "tools: $(conda run -n "$TOOLS_ENV" which fimo) ; $(conda run -n "$TOOLS_ENV" which Rscript)"

W="$SCR/dl"; mkdir -p "$SHARE/official" "$SHARE/raw/maestro" "$W"
dl "https://zenodo.org/records/14848389/files/ScmultiomeGRN-main.zip?download=1" "$W/ScmultiomeGRN-main.zip"
[ -d "$SHARE/official/ScmultiomeGRN-main" ] || unzip -q -o "$W/ScmultiomeGRN-main.zip" -d "$SHARE/official"
dl "https://raw.githubusercontent.com/liulab-dfci/MAESTRO/v1.2.1/MAESTRO/annotations/GRCh38_refgenes.txt" \
   "$SHARE/raw/maestro/GRCh38_refgenes.txt"
dl "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4508nnn/GSM4508936/suppl/GSM4508936_lung_filtered.seurat.RDS.gz" "$W/lung.RDS.gz"
test "$(stat -c %s "$W/lung.RDS.gz")" = 3011500876 && echo "lung RDS size OK" || { echo "lung RDS size differs from GEO"; exit 1; }
dl "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz" "$W/hg19.fa.gz"
gzip -t "$W/hg19.fa.gz"
put_raw "$W/lung.RDS.gz" GSM4508936_lung_filtered.seurat.RDS.gz
put_raw "$W/hg19.fa.gz" hg19.fa.gz
python - <<'PY'
from huggingface_hub import hf_hub_download, snapshot_download
for f in ("geneformer/token_dictionary_gc104M.pkl", "geneformer/gene_name_id_dict_gc104M.pkl"):
    print(hf_hub_download("ctheodoris/Geneformer", f))
print(snapshot_download("ctheodoris/Geneformer", allow_patterns=["Geneformer-V2-104M/*"]))
PY
echo "[setup] done"
