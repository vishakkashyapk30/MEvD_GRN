#!/bin/bash
# Download STRING v12.0 (human) inside a batch job and build the label-free TF-TF prior per cell type
# (scripts/scmgrn/24_build_string_tf_graph.py). CPU only. Downloads resume across retries.
#SBATCH --job-name=smg-string
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH -n 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/smg_string_%j.out
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn/scmgrn_common.sh"
activate
D="$SHARE/raw/string"; mkdir -p "$D"
B=https://stringdb-downloads.org/download
get() {   # get <url> <file> <expected bytes>
  local i
  for i in $(seq 1 30); do
    [ "$(stat -c %s "$2" 2>/dev/null || echo 0)" = "$3" ] && return 0
    curl -sfL -C - --retry 5 --retry-delay 10 -o "$2" "$1"
    echo "[string] $(basename "$2"): $(stat -c %s "$2" 2>/dev/null || echo 0) of $3 bytes (try $i)"
    sleep 5
  done
  [ "$(stat -c %s "$2" 2>/dev/null || echo 0)" = "$3" ]
}
get "$B/protein.info.v12.0/9606.protein.info.v12.0.txt.gz" "$D/9606.protein.info.v12.0.txt.gz" 1970090 || { echo "[string] info download FAILED"; exit 1; }
get "$B/protein.links.v12.0/9606.protein.links.v12.0.txt.gz" "$D/9606.protein.links.v12.0.txt.gz" 83164437 || { echo "[string] links download FAILED"; exit 1; }
gzip -t "$D/9606.protein.links.v12.0.txt.gz" && gzip -t "$D/9606.protein.info.v12.0.txt.gz" || { echo "[string] corrupt download"; exit 1; }
python -u scripts/scmgrn/24_build_string_tf_graph.py --root "$SHARE" --links "$D/9606.protein.links.v12.0.txt.gz" \
   --info "$D/9606.protein.info.v12.0.txt.gz" --min_score 700
echo "[string] done"
