#!/bin/bash
# =============================================================================
# PBMC10k / LINGER-Cistrome benchmark: downloads. Run on the Ada GATEWAY
# (ada-gw1; /share1 is visible only there). Data goes into $PBMC_ROOT.
#
#   bash scripts/25_pbmc_download.sh data    # ~0.3 GB: 10x h5, labels, GT, GTF, label sets
#   bash scripts/25_pbmc_download.sh linger  # 20.8 GB LINGER data_bulk (LINGER re-run only)
#   bash scripts/25_pbmc_download.sh fm      # cache Geneformer (HF) for the FM step
#   bash scripts/25_pbmc_download.sh linger-env  # conda env `linger` (LingerGRN==1.110 + bedtools)
#
# PBMC_ROOT defaults to /share1/$USER/mevd_grn/pbmc.
# Sources: docs/experiments/pbmc10k_linger_benchmark.md section 1.
# =============================================================================
set -euo pipefail
ROOT="${PBMC_ROOT:-/share1/$USER/mevd_grn/pbmc}"
fetch() { echo ">> $2"; mkdir -p "$(dirname "$2")"; [ -s "$2" ] && { echo "   exists"; return 0; }
          curl -sS -fL --retry 8 --retry-delay 10 -o "$2.part" "$1" && mv "$2.part" "$2"; }
gdrive() { fetch "https://drive.usercontent.google.com/download?id=$1&export=download&confirm=t" "$2"; }

GT_IDS="5967_MYC 8481_RUNX1 40215_IRF4 41287_STAT1 41288_STAT1 41289_STAT1 41290_STAT1 41301_IRF1
41302_IRF1 41303_IRF1 44092_ETS1 44093_ETS1 44094_ETS1 44097_RUNX1 44098_FOXP3 45178_CTCF
45444_CTCF 47435_REST 81223_RUNX1 85986_SPI1"

case "${1:-data}" in
data)
  X=https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_10k
  fetch "$X/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5" \
        "$ROOT/raw/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5"
  gdrive 17PXkQJr8fk0h90dCkTi3RGPmFNtDqHO_ "$ROOT/raw/PBMC_label.txt"          # LINGER cell labels
  fetch https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz \
        "$ROOT/raw/gencode.v38.annotation.gtf.gz"
  # Cistrome DB gene-score files as redistributed by KEGNI (LINGER's 20 blood ChIP sets, Table S1)
  for f in $GT_IDS; do
    fetch "https://raw.githubusercontent.com/Lipxiao/KEGNI/main/data/CristromeDB/${f}_gene_score_5fold.txt" \
          "$ROOT/groundtruth/${f}_gene_score_5fold.txt"
  done
  fetch "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes&organisms=9606&fields=sources,references,curation_effort,extra_attrs" \
        "$ROOT/labels/collectri_omnipath.tsv"
  fetch "https://omnipathdb.org/interactions?datasets=dorothea&dorothea_levels=A,B&genesymbols=yes&organisms=9606&fields=dorothea_level,sources,dorothea_curated,dorothea_chipseq,dorothea_tfbs,dorothea_coexp" \
        "$ROOT/labels/dorothea_AB_omnipath.tsv"
  test "$(stat -c %s "$ROOT/raw/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")" = 192125528 \
    && echo "h5 size OK" || echo "WARNING: h5 size differs from 10x (192125528)"
  test "$(wc -l < "$ROOT/raw/PBMC_label.txt")" = 9544 && echo "label file OK (9,543 cells)" \
    || echo "WARNING: PBMC_label.txt is not 9,544 lines (Google Drive may have returned HTML)"
  ls "$ROOT/groundtruth" | wc -l | xargs echo "ground-truth files (expect 20):"
  du -sh "$ROOT"
  ;;
linger)
  gdrive 1jwRgRHPJrKABOk7wImKONTtUupV7yJ9b "$ROOT/linger/data_bulk.tar.gz"     # 20,812,483,490 bytes
  test "$(stat -c %s "$ROOT/linger/data_bulk.tar.gz")" = 20812483490 && echo "data_bulk size OK"
  ;;
fm)
  source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"; conda activate mevd-grn
  python - <<'EOF'
from huggingface_hub import hf_hub_download, snapshot_download
for f in ("geneformer/token_dictionary_gc104M.pkl", "geneformer/gene_name_id_dict_gc104M.pkl"):
    print(hf_hub_download("ctheodoris/Geneformer", f))
print(snapshot_download("ctheodoris/Geneformer", allow_patterns=["Geneformer-V2-104M/*"]))
EOF
  ;;
linger-env)
  source "${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
  conda env list | grep -q "^linger " || \
    conda create -y -n linger -c conda-forge -c bioconda python=3.10 bedtools pybedtools=0.12.0
  conda activate linger
  pip install LingerGRN==1.110 h5py
  python -c "import LingerGRN, pybedtools, scanpy; print('LingerGRN OK', scanpy.__version__)"; which bedtools
  ;;
*) echo "usage: $0 {data|linger|fm|linger-env}"; exit 1 ;;
esac
