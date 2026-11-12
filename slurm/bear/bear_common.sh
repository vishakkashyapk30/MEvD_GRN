#!/bin/bash
# Shared helpers for the BEAR-GRN jobs (sourced, not submitted).
# Post-upgrade Ada (docs/reference/ada.md top): /share1 is visible on the gateway
# ada-gw1 only, so jobs stage $SHARE <-> node-local /scratch with rsync over ssh.
export SHARE="${BEAR_SHARE:-/share1/$USER/mevd_grn/bear}"
export GW="${GW:-ada-gw1}"
export REPO="${REPO:-${SLURM_SUBMIT_DIR:-$PWD}}"
export CONDA_SH="${CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
export SCR="${SCR_BASE:-/scratch}/$USER/bear_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-0}"
export BEAR_ROOT="$SCR/bear"
mkdir -p "$BEAR_ROOT"
cleanup() { rm -rf "$SCR"; }
trap cleanup EXIT
_src() { if [ -d "$SHARE" ]; then echo "$SHARE"; else echo "$GW:$SHARE"; fi; }
# pull <relpath> : $SHARE/<relpath> -> $BEAR_ROOT/<relpath>  (missing = no-op)
pull() {
  local rel="$1"; mkdir -p "$(dirname "$BEAR_ROOT/$rel")"
  rsync -a "$(_src)/$rel" "$(dirname "$BEAR_ROOT/$rel")/" 2>/dev/null || echo "[pull] (none) $rel"
}
# push <relpath> : $BEAR_ROOT/<relpath> -> $SHARE/<relpath>  (merge, never delete)
push() {
  local rel="$1"; [ -e "$BEAR_ROOT/$rel" ] || return 0
  if [ -d "$SHARE" ]; then mkdir -p "$(dirname "$SHARE/$rel")"; rsync -a "$BEAR_ROOT/$rel" "$(dirname "$SHARE/$rel")/"
  else ssh "$GW" "mkdir -p '$(dirname "$SHARE/$rel")'"; rsync -a "$BEAR_ROOT/$rel" "$GW:$(dirname "$SHARE/$rel")/"; fi
}
# stage_ds <dataset> : inputs + GTs + processed dir + existing results of one dataset
stage_ds() {
  local ds="$1"
  python - "$ds" <<'PY' > "$SCR/stage_$ds.txt"
import sys, yaml
c = yaml.safe_load(open("configs/bear/datasets.yaml"))
d = c["datasets"][sys.argv[1]]
paths = [d["rna"], d["atac"]] + list(d.get("labels", {}).values()) + list(d.get("gts", {}).values())
print("\n".join(sorted(set(paths))))
PY
  while read -r rel; do pull "$rel"; done < "$SCR/stage_$ds.txt"
  pull "gtf"; pull "networks"; pull "processed/$ds"; pull "processed/${ds}__L2"; pull "results/$ds"
}
activate() {
  source "$CONDA_SH"; conda activate mevd-grn
  export PYTHONUNBUFFERED=1 HF_HUB_OFFLINE="${HF_HUB_OFFLINE:-0}"
  cd "$REPO"
  echo "[env] host=$(hostname) gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1) scratch=$SCR root=$BEAR_ROOT"
}
BEAR_PHASE1=(K562 Macrophage_S1 Macrophage_S2)
