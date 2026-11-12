#!/bin/bash
# Shared helpers for the PBMC10k / LINGER-Cistrome jobs (sourced, not submitted).
# Storage plan (docs/reference/ada.md, post-upgrade section): code + conda env in
# /home2/$USER (visible on nodes); data + results on /share1, which is visible on
# the gateway ada-gw1 ONLY -> jobs rsync ada-gw1:/share1/... <-> node-local /scratch.
export SHARE="${PBMC_ROOT:-/share1/$USER/mevd_grn/pbmc}"
export MASTER="${MASTER:-ada-gw1}"             # the only host that mounts /share1
export REPO="${REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
export CONDA_SH="${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
SCR_BASE="${SCR_BASE:-/scratch}"
export SCR="$SCR_BASE/$USER/${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}_pbmc"
export PBMC_ROOT_LOCAL="$SCR/pbmc"
mkdir -p "$PBMC_ROOT_LOCAL"
cleanup() { rm -rf "$SCR"; }
trap cleanup EXIT

# pull <relpath> : $SHARE/<relpath> -> $PBMC_ROOT_LOCAL/<relpath>
pull() {
  local rel="$1"; mkdir -p "$(dirname "$PBMC_ROOT_LOCAL/$rel")"
  if [ -e "$SHARE/$rel" ]; then rsync -a "$SHARE/$rel" "$(dirname "$PBMC_ROOT_LOCAL/$rel")/"
  else rsync -a "$MASTER:$SHARE/$rel" "$(dirname "$PBMC_ROOT_LOCAL/$rel")/" 2>/dev/null || echo "[pull] (absent) $rel"; fi
}
# push <relpath> : local -> share (directory contents merged, never deleted)
push() {
  local rel="$1"; [ -e "$PBMC_ROOT_LOCAL/$rel" ] || return 0
  if [ -d "$SHARE" ]; then
    mkdir -p "$(dirname "$SHARE/$rel")"; rsync -a "$PBMC_ROOT_LOCAL/$rel" "$(dirname "$SHARE/$rel")/"
  else
    ssh "$MASTER" "mkdir -p '$(dirname "$SHARE/$rel")'"
    rsync -a "$PBMC_ROOT_LOCAL/$rel" "$MASTER:$(dirname "$SHARE/$rel")/"
  fi
}
activate() {
  source "$CONDA_SH"; conda activate mevd-grn
  export PYTHONUNBUFFERED=1
  cd "$REPO"
  echo "[env] host=$(hostname) gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1) scratch=$SCR"
}
# evaluated cell types (configs/pbmc/benchmark.yaml), largest first
PBMC_CTS=(classical_monocyte naive_cd4_t naive_b mdc)
