#!/bin/bash
# =============================================================================
# Download SC-MO-GRN-DB data for the "partial-evidence" cell lines used for
# zero-shot cross-cell-type TRANSFER evaluation (localization-only ground
# truth; no perturbation/dual-evidence network exists for these cell types,
# so a curriculum-trained model is scored zero-shot -- never fine-tuned on
# these cell types at all).
#
# Catalog (from https://scmogrndb.psu.edu/browse.php), all Human:
#   HepG2      RN102 (84 TF,  342862 edges)  DS007 scRNA only (426 cells)
#   BJ         RN200 (24 TF,  205697 edges)  DS023 scRNA/scHiC (63c) + DS031 scATAC (96c)
#   GM12878    RN201 (178 TF, 1410483 edges) DS028 scRNA (8274c) + DS033 scATAC (384c)
#   Macrophage RN204 (24 TF,  169972 edges)  DS026 scRNA/scATAC paired (3114 cells)
#   MCF7       RN205 (262 TF, 1809170 edges) DS027 scRNA/scATAC paired (2995 cells)
#   HSC        RN103 (11 TF,  30 edges -- tiny, bonus only) DS030 scRNA (183 cells)
#
# NOT included (documented limitations, see README):
#   H1 -- only single-cell dataset (DS032) is scATAC-only; no scRNA exists for
#         H1 in the DB, so it can't go through our RNA-centric pipeline.
#   DC -- only network (RN110) and dataset (DS002) are MOUSE, not human; would
#         need the ESC (mouse) trained model, not K562. Deferred until ESC is
#         downloaded/trained.
#
# Usage:  bash scripts/09_download_transfer_data.sh
# =============================================================================
set -euo pipefail

BASE="https://scmogrndb.psu.edu"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
NET_DIR="$ROOT/data/raw/networks"
SC_DIR="$ROOT/data/raw/single_cell"
mkdir -p "$NET_DIR" "$SC_DIR"

fetch() { # fetch <url> <outpath>
  echo ">> $2"
  curl -fL -C - --retry 5 --retry-delay 5 -o "$2" "$1"
}

echo "=== Reference networks (localization only) ==="
declare -A NETS=(
  [RN102]="HepG2_localization"
  [RN200]="BJ_localization"
  [RN201]="GM12878_localization"
  [RN204]="Macrophage_localization"
  [RN205]="MCF7_localization"
  [RN103]="HSC_localization"
)
for rn in "${!NETS[@]}"; do
  fetch "$BASE/RN_TSV/${rn}.tsv" "$NET_DIR/${NETS[$rn]}.tsv"
done

echo ""
echo "=== Single-cell datasets ==="
declare -A DSETS=(
  [DS007]="HepG2_scRNA"          # scRNA         Human HepG2   426 cells
  [DS023]="BJ_scRNA_scHiC"       # scRNA/scHiC   Human BJ       63 cells
  [DS031]="BJ_scATAC"            # scATAC        Human BJ       96 cells
  [DS028]="GM12878_scRNA"        # scRNA         Human GM12878 8274 cells
  [DS033]="GM12878_scATAC"       # scATAC        Human GM12878  384 cells
  [DS026]="Macrophage_multiome"  # scRNA/scATAC  Human Macrophage 3114 cells (paired)
  [DS027]="MCF7_multiome"        # scRNA/scATAC  Human MCF7    2995 cells (paired)
  [DS030]="HSC_scRNA"            # scRNA         Human HSC      183 cells
)
for ds in "${!DSETS[@]}"; do
  out="$SC_DIR/${DSETS[$ds]}.zip"
  fetch "$BASE/DS_ZIP/${ds}.zip" "$out"
  echo "   unzip -> $SC_DIR/${DSETS[$ds]}/"
  unzip -o -q "$out" -d "$SC_DIR/${DSETS[$ds]}"
done

echo ""
echo "=== Done. Summary ==="
echo "Networks:"; ls -lh "$NET_DIR" | grep -E "HepG2|BJ|GM12878|Macrophage|MCF7|HSC" || true
echo "Single-cell:"; du -sh "$SC_DIR"/{HepG2,BJ,GM12878,Macrophage,MCF7,HSC}* 2>/dev/null || true
cat <<'EOF'

Next: python scripts/02_preprocess.py --config configs/transfer/<celltype>.yaml
      python scripts/08_transfer_eval.py --config configs/transfer/<celltype>.yaml
EOF
