#!/bin/bash
# Laptop GPU worker (RTX 4060, 8 GB): joins the SAME queue as the Ada workers, claiming MeVD-GRN units
# over ssh (KINDS=mevd; the official code stays on Ada), staging only a cell type's small processed
# data (graph/splits/processed, no raw data), and publishing result JSONs back to Ada's $HOME/scmgrn.
# Waits (polling every 5 min) until Ada's build has finished a cell type. Run in a terminal:
#   DS=lung bash slurm/scmgrn/scmgrn_local_worker.sh          # Ctrl-C to stop; claims of a killed unit
#                                                      # stay in claims/ -> release them with the reset note in docs
set -uo pipefail
export DS="${DS:-lung}" KINDS=mevd WAIT_FOR_BUILD=1 DEVICE=cuda:0
export RAW_HOST="${RAW_HOST:-${ADA_USER:-vishakkashyap.k}@ada-gw2}"
export SHARE="/home2/${ADA_USER:-vishakkashyap.k}/scmgrn"                 # path ON ADA
export REPO="${REPO:-$(cd "$(dirname "$0")/.." && pwd)}"
export SCR="${LOCAL_SCR:-$HOME/.cache/mevd_scmgrn_local}"; export SCMGRN_ROOT_LOCAL="$SCR/scmgrn"
mkdir -p "$SCMGRN_ROOT_LOCAL"; cd "$REPO"
source "$REPO/slurm/scmgrn/scmgrn_worker_lib.sh" || exit 1
SSH="ssh -n -o BatchMode=yes -o ConnectTimeout=20 $RAW_HOST"   # -n: ssh must not swallow the unit list the worker loop reads from stdin
q_exec() { $SSH "$(printf '%q ' "$@")"; }
stage_ct() { local CT=$1 r
  # script 17 reads only graph/x_graph.txt (adjacency), splits and processed/: skip the big feature files in graph/
  mkdir -p "$SCMGRN_ROOT_LOCAL/graph/$DS/$CT"
  rsync -a -e "ssh -o BatchMode=yes" "$RAW_HOST:$SHARE/graph/$DS/$CT/x_graph.txt" "$RAW_HOST:$SHARE/graph/$DS/$CT/x_node.txt" \
        "$SCMGRN_ROOT_LOCAL/graph/$DS/$CT/" || return 1
  for r in splits processed; do
    mkdir -p "$SCMGRN_ROOT_LOCAL/$r/$DS"
    rsync -a -e "ssh -o BatchMode=yes" "$RAW_HOST:$SHARE/$r/$DS/$CT" "$SCMGRN_ROOT_LOCAL/$r/$DS/" || return 1
  done; }
publish_ct() { rsync -a -e "ssh -o BatchMode=yes" "$SCMGRN_ROOT_LOCAL/results/$DS/$1" "$RAW_HOST:$SHARE/results/$DS/" ; }
worker_loop
