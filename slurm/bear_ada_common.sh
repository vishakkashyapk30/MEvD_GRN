#!/bin/bash
# Shared helpers for the BEAR-GRN jobs on Ada (sourced). Everything lives in $HOME/bear_v2
# (home is visible on compute nodes, so there is NO compute->gateway copy: that is what lost
# the results of the earlier BEAR jobs 2400_2/3 and 2401 with "Permission denied").
# Jobs train on node-local /scratch and publish results to $B with verified, retried copies.
export B="${BEAR_V2:-$HOME/bear_v2}"
export REPO="${REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
export CONDA_SH="${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
export SCR="/scratch/$USER/bear_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
export BEAR_ROOT="$SCR/bear"
mkdir -p "$BEAR_ROOT/results" "$B"/{claims,done,failed,logs,results,processed}
cleanup() { rm -rf "$SCR"; }
trap cleanup EXIT
activate() {
  source "$CONDA_SH"; conda activate mevd-grn
  export PYTHONUNBUFFERED=1 HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1
  cd "$REPO"
  echo "[env] host=$(hostname) job=${SLURM_JOB_ID:-} task=${SLURM_ARRAY_TASK_ID:-} gpu=$(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null | head -1) scratch=$SCR"
}
# copy_verified <src_file> <dst_file> : up to 5 tries with backoff, size verified before returning 0
copy_verified() {
  local src="$1" dst="$2" i s1 s2
  mkdir -p "$(dirname "$dst")"
  s1=$(stat -c %s "$src")
  for i in 1 2 3 4 5; do
    if cp -f "$src" "$dst.part" && s2=$(stat -c %s "$dst.part") && [ "$s1" = "$s2" ] && mv -f "$dst.part" "$dst"; then return 0; fi
    echo "[copy] attempt $i failed for $dst; retry in $((i*20))s"; sleep $((i*20))
  done
  return 1
}
# publish_dir <src_dir> <dst_dir> : every file copied with copy_verified; 0 only if ALL verified
publish_dir() {
  local src="$1" dst="$2" f rel rc=0
  while IFS= read -r -d '' f; do
    rel="${f#$src/}"; copy_verified "$f" "$dst/$rel" || rc=1
  done < <(find "$src" -type f -print0)
  return $rc
}
# stage_in <relpath under $B> : $B/<rel> -> $BEAR_ROOT/<rel> (retry; processed dirs are small)
stage_in() {
  local rel="$1" i
  for i in 1 2 3 4 5; do
    rsync -a "$B/$rel/" "$BEAR_ROOT/$rel/" 2>/dev/null && return 0
    sleep $((i*10))
  done
  echo "[stage_in] FAILED $rel"; return 1
}
