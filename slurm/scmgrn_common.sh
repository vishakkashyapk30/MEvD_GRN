#!/bin/bash
# Shared helpers for the scMultiomeGRN-benchmark jobs (sourced, not submitted).
# Storage plan (docs/reference/ada.md): code + conda envs in $HOME; data and
# results on /share1 (master node only) -> staged to node-local scratch.
export SHARE="${SCMGRN_ROOT:-/share1/$USER/mevd_grn/scmgrn}"
export MASTER="${MASTER:-ada}"                 # host that mounts /share1
export REPO="${REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
export CONDA_SH="${CONDA_SH:-$HOME/miniconda3/etc/profile.d/conda.sh}"
export TOOLS_ENV="${TOOLS_ENV:-scmgrn-tools}"
SCR_BASE="${SCR_BASE:-/scratch}"
export SCR="$SCR_BASE/$USER/${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
export SCMGRN_ROOT_LOCAL="$SCR/scmgrn"
mkdir -p "$SCMGRN_ROOT_LOCAL"
cleanup() { rm -rf "$SCR"; }
trap cleanup EXIT

# pull <relpath> : $SHARE/<relpath> -> $SCMGRN_ROOT_LOCAL/<relpath>
pull() {
  local rel="$1"; mkdir -p "$(dirname "$SCMGRN_ROOT_LOCAL/$rel")"
  if [ -e "$SHARE/$rel" ]; then rsync -a "$SHARE/$rel" "$(dirname "$SCMGRN_ROOT_LOCAL/$rel")/"
  else rsync -a "$MASTER:$SHARE/$rel" "$(dirname "$SCMGRN_ROOT_LOCAL/$rel")/" 2>/dev/null || true; fi
}
# push <relpath> : local -> share (directory contents merged, never deleted)
push() {
  local rel="$1"; [ -e "$SCMGRN_ROOT_LOCAL/$rel" ] || return 0
  if [ -d "$(dirname "$SHARE/$rel")" ] || [ -d "$SHARE" ]; then
    mkdir -p "$(dirname "$SHARE/$rel")"; rsync -a "$SCMGRN_ROOT_LOCAL/$rel" "$(dirname "$SHARE/$rel")/"
  else
    ssh "$MASTER" "mkdir -p '$(dirname "$SHARE/$rel")'"
    rsync -a "$SCMGRN_ROOT_LOCAL/$rel" "$MASTER:$(dirname "$SHARE/$rel")/"
  fi
}
activate() {
  source "$CONDA_SH"; conda activate mevd-grn
  export FIMO_BIN="${FIMO_BIN:-$(conda run -n "$TOOLS_ENV" which fimo)}"
  export RSCRIPT_BIN="${RSCRIPT_BIN:-$(conda run -n "$TOOLS_ENV" which Rscript)}"
  export PYTHONUNBUFFERED=1
  cd "$REPO"
  echo "[env] host=$(hostname) gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1) scratch=$SCR"
}
# Lung cell types (configs/scmgrn/lung.yaml), largest first for load balance
LUNG_CTS=(Stromal_cells Bronchiolar_and_alveolar_epithelial_cells Vascular_endothelial_cells
          Lymphoid_cells Myeloid_cells Ciliated_epithelial_cells Lymphatic_endothelial_cells
          Megakaryocytes Neuroendocrine_cells)
