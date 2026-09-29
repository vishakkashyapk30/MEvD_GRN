#!/bin/bash
# =============================================================================
# scMultiomeGRN benchmark: downloads + tool setup. Run on the Ada LOGIN node
# (compute nodes have no internet); data goes straight into /share1.
#
#   bash scripts/14_download_scmultiomegrn_data.sh env     # once: tools env + pip deps
#   bash scripts/14_download_scmultiomegrn_data.sh data    # ~6.5 GB into $SCMGRN_ROOT
#   bash scripts/14_download_scmultiomegrn_data.sh fm      # cache Geneformer for offline jobs
#
# SCMGRN_ROOT defaults to /share1/$USER/mevd_grn/scmgrn.
# Sources (see docs/experiments/scmultiomegrn_generalization.md):
#   official code  Zenodo 10.5281/zenodo.14848389 (ScmultiomeGRN-main.zip, 8.4 MB;
#                  HOCOMOCO v11 motifs + GENCODE v19/v38 promoters inside)
#   lung           GEO GSM4508936_lung_filtered.seurat.RDS.gz (3.01 GB)
#   genome         UCSC hg19.fa.gz (~0.9 GB -> 3.1 GB)
#   MAESTRO        v1.2.1 annotations/GRCh38_refgenes.txt (22.7 MB)
# =============================================================================
set -euo pipefail
ROOT="${SCMGRN_ROOT:-/share1/$USER/mevd_grn/scmgrn}"
REPO="$(cd "$(dirname "$0")/.." && pwd)"
CONDA_SH="${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}"
fetch() { echo ">> $2"; mkdir -p "$(dirname "$2")"; curl -fL -C - --retry 8 --retry-delay 10 -o "$2" "$1"; }

cmd="${1:-data}"
case "$cmd" in
env)
  source "$CONDA_SH"
  # FIMO (MEME 5.4.1, the version pinned by the official README) + R/Matrix for the RDS split
  conda create -y -n scmgrn-tools -c conda-forge -c bioconda meme=5.4.1 r-base r-matrix
  conda activate mevd-grn
  # official code needs lightning; GRNBoost2 per-target routine comes from arboreto
  pip install "lightning>=2.0,<2.4" arboreto
  echo "tools: $(conda run -n scmgrn-tools which fimo) ; $(conda run -n scmgrn-tools which Rscript)"
  ;;
data)
  mkdir -p "$ROOT"/{raw/lung,raw/genome,raw/maestro,official}
  z="$ROOT/raw/ScmultiomeGRN-main.zip"
  fetch "https://zenodo.org/records/14848389/files/ScmultiomeGRN-main.zip?download=1" "$z"
  [ -d "$ROOT/official/ScmultiomeGRN-main" ] || unzip -q -o "$z" -d "$ROOT/official"
  fetch "https://raw.githubusercontent.com/liulab-dfci/MAESTRO/v1.2.1/MAESTRO/annotations/GRCh38_refgenes.txt" \
        "$ROOT/raw/maestro/GRCh38_refgenes.txt"
  if [ ! -s "$ROOT/raw/genome/hg19.fa" ]; then
    fetch "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz" "$ROOT/raw/genome/hg19.fa.gz"
    gunzip -f "$ROOT/raw/genome/hg19.fa.gz"
  fi
  fetch "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4508nnn/GSM4508936/suppl/GSM4508936_lung_filtered.seurat.RDS.gz" \
        "$ROOT/raw/lung/GSM4508936_lung_filtered.seurat.RDS.gz"
  # sizes as a cheap integrity check (expected lung RDS: 3011500876 bytes)
  ls -l "$ROOT/raw/lung" "$ROOT/raw/genome" "$ROOT/raw/maestro"
  test "$(stat -c %s "$ROOT/raw/lung/GSM4508936_lung_filtered.seurat.RDS.gz")" = 3011500876 \
    && echo "lung RDS size OK" || echo "WARNING: lung RDS size differs from GEO (3011500876)"
  du -sh "$ROOT"
  ;;
fm)
  source "$CONDA_SH"; conda activate mevd-grn
  python - <<'EOF'
from huggingface_hub import hf_hub_download, snapshot_download
for f in ("geneformer/token_dictionary_gc104M.pkl", "geneformer/gene_name_id_dict_gc104M.pkl"):
    print(hf_hub_download("ctheodoris/Geneformer", f))
print(snapshot_download("ctheodoris/Geneformer", allow_patterns=["Geneformer-V2-104M/*"]))
EOF
  ;;
*) echo "usage: $0 {env|data|fm}"; exit 1 ;;
esac
