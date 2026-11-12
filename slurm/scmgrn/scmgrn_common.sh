#!/bin/bash
# Shared helpers for the scMultiomeGRN-benchmark jobs (sourced, not submitted).
# Post-upgrade Ada (docs/reference/ada.md); the gateway used is ada-gw2 (gw1 is down):
#   $SHARE     = $HOME/scmgrn : persistent derived data + results. Home IS visible on compute
#                nodes, so jobs read/write it directly (no compute->gateway copies of results).
#   $RAW_SHARE = /share1/$USER/mevd_grn/scmgrn/raw on ada-gw2 : a FEW big raw files only
#                (lung RDS.gz, hg19.fa.gz; /share1 is at its file-count quota, ~3190/3200).
#                /share1 is not mounted on compute nodes -> fetched with verified scp.
#   $SCR       = node-local /scratch/$USER/<job>: all compute happens here, wiped on exit.
export SHARE="${SCMGRN_ROOT:-$HOME/scmgrn}"
export RAW_HOST="${RAW_HOST:-ada-gw2}"
export RAW_SHARE="${RAW_SHARE:-/share1/$USER/mevd_grn/scmgrn/raw}"
export REPO="${REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
export CONDA_SH="${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
export TOOLS_ENV="${TOOLS_ENV:-scmgrn-tools}"
SCR_BASE="${SCR_BASE:-/scratch}"
export SCR="$SCR_BASE/$USER/smg_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
export SCMGRN_ROOT_LOCAL="$SCR/scmgrn"
mkdir -p "$SCMGRN_ROOT_LOCAL" "$SHARE"
cleanup() { rm -rf "$SCR"; }
trap cleanup EXIT

# pull <relpath> : $SHARE/<relpath> -> $SCMGRN_ROOT_LOCAL/<relpath> (retried)
pull() {
  local rel="$1" i; mkdir -p "$(dirname "$SCMGRN_ROOT_LOCAL/$rel")"
  [ -e "$SHARE/$rel" ] || return 0
  for i in 1 2 3; do
    rsync -a "$SHARE/$rel" "$(dirname "$SCMGRN_ROOT_LOCAL/$rel")/" && return 0; sleep $((i*5)); done
  echo "[pull] FAILED $rel"; return 1
}
# push <relpath> : local -> $SHARE (merged, never deleted; retried; size-verified by rsync)
push() {
  local rel="$1" i; [ -e "$SCMGRN_ROOT_LOCAL/$rel" ] || return 0
  mkdir -p "$(dirname "$SHARE/$rel")"
  for i in 1 2 3 4 5; do
    rsync -a "$SCMGRN_ROOT_LOCAL/$rel" "$(dirname "$SHARE/$rel")/" && return 0
    echo "[push] attempt $i failed for $rel"; sleep $((i*20)); done
  echo "[push] FAILED $rel"; return 1
}
# fetch_raw <name> <dest_file> : scp ada-gw2:$RAW_SHARE/<name> -> dest, size-verified, 5 tries
fetch_raw() {
  local name="$1" dst="$2" i want got; mkdir -p "$(dirname "$dst")"
  if [ "${WAIT_RAW:-0}" = 1 ]; then       # inputs are being uploaded from the laptop: wait (<= 8 h) for the final name
    for i in $(seq 1 240); do ssh -o BatchMode=yes "$RAW_HOST" "test -e '$RAW_SHARE/$name'" && break
      echo "[fetch_raw] waiting for $name ($(date +%T))"; sleep 120; done
  fi
  want=$(ssh -o BatchMode=yes "$RAW_HOST" "stat -c %s '$RAW_SHARE/$name'") || return 1
  for i in 1 2 3 4 5; do
    scp -q -o BatchMode=yes "$RAW_HOST:$RAW_SHARE/$name" "$dst.part" && got=$(stat -c %s "$dst.part") \
      && [ "$got" = "$want" ] && mv -f "$dst.part" "$dst" && return 0
    echo "[fetch_raw] attempt $i failed for $name"; sleep $((i*20)); done
  return 1
}
# put_raw <src_file> <name> : scp -> ada-gw2:$RAW_SHARE/<name>, size-verified, 5 tries
put_raw() {
  local src="$1" name="$2" i want got; want=$(stat -c %s "$src")
  ssh -o BatchMode=yes "$RAW_HOST" "mkdir -p '$RAW_SHARE'" || return 1
  for i in 1 2 3 4 5; do
    scp -q -o BatchMode=yes "$src" "$RAW_HOST:$RAW_SHARE/$name.part" \
      && got=$(ssh -o BatchMode=yes "$RAW_HOST" "stat -c %s '$RAW_SHARE/$name.part'") && [ "$got" = "$want" ] \
      && ssh -o BatchMode=yes "$RAW_HOST" "mv -f '$RAW_SHARE/$name.part' '$RAW_SHARE/$name'" && return 0
    echo "[put_raw] attempt $i failed for $name"; sleep $((i*20)); done
  return 1
}
activate() {
  source "$CONDA_SH"; conda activate mevd-grn
  export FIMO_BIN="${FIMO_BIN:-$HOME/miniforge3/envs/$TOOLS_ENV/bin/fimo}"
  export RSCRIPT_BIN="${RSCRIPT_BIN:-$HOME/miniforge3/envs/$TOOLS_ENV/bin/Rscript}"
  [ -x "$FIMO_BIN" ] && [ -x "$RSCRIPT_BIN" ] || { echo "[env] missing fimo/Rscript in env $TOOLS_ENV"; exit 1; }
  export PYTHONUNBUFFERED=1
  cd "$REPO"
  echo "[env] host=$(hostname) job=${SLURM_JOB_ID:-} task=${SLURM_ARRAY_TASK_ID:-} gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1) scratch=$SCR share=$SHARE"
}
# DS=lung (default, the paper's main benchmark) | hema (Supp. Fig. S1); cell types SMALL graphs first
export DS="${DS:-lung}"
LUNG_CTS=(Myeloid_cells Lymphoid_cells Megakaryocytes Neuroendocrine_cells Lymphatic_endothelial_cells
          Bronchiolar_and_alveolar_epithelial_cells Stromal_cells Ciliated_epithelial_cells Vascular_endothelial_cells)
HEMA_CTS=(CMP HSC GMP)
if [ "$DS" = hema ]; then CTS=("${HEMA_CTS[@]}"); else CTS=("${LUNG_CTS[@]}"); fi
