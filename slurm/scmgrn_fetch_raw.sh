#!/bin/bash
# Download the big raw inputs ON ADA with many parallel HTTP range connections. Ada's outbound speed is limited
# PER CONNECTION (measured 2026-10-16: 1 connection ~140 KB/s, 8 connections ~1.4 MB/s aggregate), so a
# parallel range download is far faster than streaming from the laptop over the VPN. CPU-only batch job (no
# downloads on the login node). Output: ada-gw2:/share1/$USER/mevd_grn/scmgrn/raw/{lung RDS.gz, hg19.fa.gz}.
#SBATCH --job-name=smg-fetch
#SBATCH --account=research
#SBATCH --qos=medium
#SBATCH --partition=u22
#SBATCH --constraint=2080ti
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000
#SBATCH --time=06:00:00
#SBATCH --output=logs/smg_fetch_%j.out
set -uo pipefail
source "${SLURM_SUBMIT_DIR:-$PWD}/slurm/scmgrn_common.sh"
NCONN=${NCONN:-16}
W="$SCR/dl"; mkdir -p "$W"

fetch_parallel() {   # fetch_parallel <url> <total_bytes> <out_file> <n_parts>
  local url="$1" size="$2" out="$3" n="$4" i part=$(( ($2 + $4 - 1) / $4 ))
  rm -rf "$W/parts_$(basename "$out")"; mkdir -p "$W/parts_$(basename "$out")"
  for i in $(seq 0 $((n-1))); do
    local s=$((i*part)) e=$((i*part+part-1)); [ $e -ge $size ] && e=$((size-1))
    echo "$i $s $e"
  done | xargs -P "$NCONN" -L 1 bash -c '
    i=$0; s=$1; e=$2; f="'"$W/parts_$(basename "$out")"'/p$(printf %04d $i)"; want=$((e-s+1))
    for try in 1 2 3 4 5 6 7 8; do
      [ "$(stat -c %s "$f" 2>/dev/null || echo 0)" = "$want" ] && exit 0
      curl -sfL --retry 3 --retry-delay 5 -r $((s+$(stat -c %s "$f" 2>/dev/null || echo 0)))-$e -o "$f.tmp" "'"$url"'" && cat "$f.tmp" >> "$f" ; rm -f "$f.tmp"
      [ "$(stat -c %s "$f" 2>/dev/null || echo 0)" = "$want" ] && exit 0
      sleep 5
    done; echo "part $i FAILED"; exit 1'
  local rc=$?
  [ $rc -eq 0 ] || return 1
  cat "$W/parts_$(basename "$out")"/p* > "$out" && rm -rf "$W/parts_$(basename "$out")"
  [ "$(stat -c %s "$out")" = "$size" ]
}

# publish <local_file> <remote_name> : private temp name on /share1, size-verified, then atomic rename
publish() {
  local src="$1" name="$2" i want got; want=$(stat -c %s "$src")
  ssh -o BatchMode=yes "$RAW_HOST" "mkdir -p '$RAW_SHARE'" || return 1
  for i in 1 2 3 4 5; do
    scp -q -o BatchMode=yes "$src" "$RAW_HOST:$RAW_SHARE/$name.ada.part" \
      && got=$(ssh -o BatchMode=yes "$RAW_HOST" "stat -c %s '$RAW_SHARE/$name.ada.part'") && [ "$got" = "$want" ] \
      && ssh -o BatchMode=yes "$RAW_HOST" "mv -f '$RAW_SHARE/$name.ada.part' '$RAW_SHARE/$name'" && return 0
    echo "[publish] attempt $i failed for $name"; sleep 20
  done
  return 1
}

RDS_URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4508nnn/GSM4508936/suppl/GSM4508936_lung_filtered.seurat.RDS.gz"
RDS_SIZE=3011500876
HG_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
HG_SIZE=948731419
t0=$(date +%s)
fetch_parallel "$RDS_URL" $RDS_SIZE "$W/lung.RDS.gz" 48 || { echo "[fetch] RDS download FAILED"; exit 1; }
echo "[fetch] RDS downloaded in $(( $(date +%s) - t0 )) s"
gzip -t "$W/lung.RDS.gz" && echo "[fetch] RDS gzip OK"
publish "$W/lung.RDS.gz" GSM4508936_lung_filtered.seurat.RDS.gz && echo "[fetch] RDS PUBLISHED $(date +%T)"
rm -f "$W/lung.RDS.gz"
fetch_parallel "$HG_URL" $HG_SIZE "$W/hg19.fa.gz" 24 || { echo "[fetch] hg19 download FAILED"; exit 1; }
gzip -t "$W/hg19.fa.gz" && echo "[fetch] hg19 gzip OK"
publish "$W/hg19.fa.gz" hg19.fa.gz && echo "[fetch] hg19 PUBLISHED $(date +%T)"
echo "[fetch] done in $(( $(date +%s) - t0 )) s"
