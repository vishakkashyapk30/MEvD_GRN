#!/bin/bash
# =============================================================================
# SC-MO-GRN-DB selective download for MEvD-GRN
# =============================================================================
# The full single-cell collection is ~24.7 GB (17.5 GB human + 7.2 GB mouse),
# so we DO NOT pull the per-organism mega-ZIPs. Instead we fetch only the
# individual files needed for the two cell types that have a COMPLETE
# localization -> perturbation -> dual-evidence reference-network triple:
#
#   K562  (human) : RN117 (ChIPseq=loc)  RN118 (KO=pert)      RN119 (ChIPseq&KO=dual)
#   ESC   (mouse) : RN114 (ChIPchip=loc) RN115 (LOGOF=pert)   RN116 (ChIPchip&LOGOF=dual)
#
# These are the ONLY two cell types in the database with all three tiers.
# The ESC triple (RN114/115/116) is self-consistent (same study, PMID 23794736,
# RN116 = intersection of RN114 & RN115), guaranteeing the nesting the
# curriculum assumes. The K562 dual (RN119) comes from a different study than
# RN117/RN118 -- nesting is verified (not assumed) during preprocessing.
#
# Per-file URL patterns discovered from the live site:
#   networks : https://scmogrndb.psu.edu/RN_TSV/RN###.tsv         (3-col TSV: Source\tTarget\tRelationship, has header)
#   datasets : https://scmogrndb.psu.edu/DS_ZIP/DS###.zip
#
# Usage:   bash scripts/01_download_data.sh
# On Ada:  run this from a login node (internet access); data lands in data/raw/
# =============================================================================
set -euo pipefail

BASE="https://scmogrndb.psu.edu"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
NET_DIR="$ROOT/data/raw/networks"
SC_DIR="$ROOT/data/raw/single_cell"
mkdir -p "$NET_DIR" "$SC_DIR"

# curl with resume (-C -) so interrupted large transfers continue.
fetch() { # fetch <url> <outpath>
  echo ">> $2"
  curl -fL -C - --retry 5 --retry-delay 5 -o "$2" "$1"
}

echo "=== Reference networks (tiny; always download all six) ==="
declare -A NETS=(
  [RN117]="K562_localization"  [RN118]="K562_perturbation"  [RN119]="K562_dual_evidence"
  [RN114]="ESC_localization"   [RN115]="ESC_perturbation"   [RN116]="ESC_dual_evidence"
)
for rn in "${!NETS[@]}"; do
  fetch "$BASE/RN_TSV/${rn}.tsv" "$NET_DIR/${NETS[$rn]}.tsv"
done

echo ""
echo "=== Single-cell datasets (paired RNA+ATAC, kept small on purpose) ==="
# K562: DS025 is the only joint scRNA+scATAC dataset (434 cells, ~522 MB, barcodes matched).
#       DS019 adds 953 RNA-only cells (~17 MB) for richer gene-level RNA statistics.
# ESC : DS010 is joint scRNA+scATAC(+scHiC) (9021 cells, ~37 MB) -- small AND high cell count.
#       DS012 (paired, 930 cells, ~12 MB) kept as a lightweight fallback.
declare -A DSETS=(
  [DS025]="K562_multiome"     # scRNA/scATAC  Human K562  434 cells   ~522 MB
  [DS019]="K562_scRNA"        # scRNA         Human K562  953 cells   ~17 MB
  [DS010]="ESC_multiome"      # scRNA/scATAC/scHiC Mouse ESC 9021 cells ~37 MB
  [DS012]="ESC_multiome_alt"  # scRNA/scATAC  Mouse ESC   930 cells   ~12 MB
)
for ds in "${!DSETS[@]}"; do
  out="$SC_DIR/${DSETS[$ds]}.zip"
  fetch "$BASE/DS_ZIP/${ds}.zip" "$out"
  echo "   unzip -> $SC_DIR/${DSETS[$ds]}/"
  unzip -o -q "$out" -d "$SC_DIR/${DSETS[$ds]}"
done

echo ""
echo "=== Done. Summary ==="
echo "Networks:"; ls -lh "$NET_DIR"
echo "Single-cell:"; du -sh "$SC_DIR"/*/ 2>/dev/null || true
cat <<'EOF'

Next:
  python scripts/02_preprocess.py --config configs/k562.yaml --cell_type K562
  python scripts/02_preprocess.py --config configs/esc.yaml  --cell_type ESC

Optional larger datasets (only if you want more cells; edit this script to add):
  K562 scATAC (DS020, 1923 cells, ~2.7 GB), K562 multiome alt (none)
  ESC  scATAC-heavy (DS014, 68794 cells, ~2.0 GB; DS011, 36535 cells, ~626 MB)
EOF
