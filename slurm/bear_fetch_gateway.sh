#!/bin/bash
# BEAR-GRN data -> /share1 (run ON THE GATEWAY ada-gw1; /share1 is not mounted on
# compute nodes). Whole zips are fetched with scripts/_bear_pget.py (16 parallel
# range requests, resumable; Zenodo throttles one connection to ~50 KB/s) and unpacked with
# Info-ZIP unzip, which handles the Deflate64 members of INFERRED.GRNS.zip that
# Python's zipfile cannot read.
#   bash slurm/bear_fetch_gateway.sh            # inputs + all GT variants + released GRNs (~2.3 GB zipped)
#   WITH_STABILITY=1 bash slurm/bear_fetch_gateway.sh   # + INPUT.DATA.STABILITY (1.7 GB)
#   WITH_STABILITY_GRNS=1 ...                          # + STABILITY_GRNS (12.5 GB; only for re-scoring their JI)
set -euo pipefail
ROOT="${BEAR_ROOT:-/share1/$USER/mevd_grn/bear}"
REC=https://zenodo.org/api/records/20704929/files
PGET="${PGET:-$(cd "$(dirname "$0")/.." && pwd)/scripts/_bear_pget.py}"
mkdir -p "$ROOT/zips"
cd "$ROOT"
fetch() {  # fetch <zip name>
  local z="$1"
  echo ">> $z"
  python3 "$PGET" "$REC/$z/content" "zips/$z" --conns "${CONNS:-16}"
}
ZIPS=(GROUND.TRUTHS.zip GROUND.TRUTHS.KO.zip GROUND.TRUTHS.UNION.zip GROUND.TRUTHS.INETRSECTION.zip
      CORE_GROUND_TRUTH.zip CELL.TYPE.EXCLUSIVE.GROUND.TRUTH.zip INPUT.DATA.zip INFERRED.GRNS.zip)
[ -n "${WITH_STABILITY:-}" ] && ZIPS+=(INPUT.DATA.STABILITY.zip)
[ -n "${WITH_STABILITY_GRNS:-}" ] && ZIPS+=(STABILITY_GRNS.zip)
for z in "${ZIPS[@]}"; do
  fetch "$z"
  unzip -q -n "zips/$z" -d "$ROOT" -x '__MACOSX/*' '*/.DS_Store'
done
# integrity: sizes vs the Zenodo record
python3 - "$ROOT" <<'EOF'
import json, os, sys, urllib.request
root = sys.argv[1]
rec = json.load(urllib.request.urlopen("https://zenodo.org/api/records/20704929", timeout=60))
for f in rec["files"]:
    p = os.path.join(root, "zips", f["key"])
    if os.path.exists(p):
        ok = os.path.getsize(p) == f["size"]
        print(("OK  " if ok else "BAD ") + f["key"], os.path.getsize(p), f["size"])
EOF
du -sh "$ROOT"/* | sort -h
echo "BEAR_FETCH_DONE"
