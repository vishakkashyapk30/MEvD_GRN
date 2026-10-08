#!/bin/bash
# GPU-aware sequential launcher for BEAR MeVD-GRN runs on the shared laptop GPU.
# Rules (coordinator, 2026-10-09): never more than one BEAR GPU job at a time; start a
# job only when free GPU memory >= its measured peak + 1 GB headroom; never touch the
# other queue. The memory check uses the shared, flock-serialised gate
# ~/.cache/local_runs/scripts/gpu_gate.sh (the other queue uses it too).
#   setsid nohup bash slurm/bear_gpu_queue.sh ~/.cache/bear/gpu_queue.txt > ~/.cache/bear/gpu_queue.out 2>&1 &
# Queue file lines: "<need_mb> <command ...>"; '#' lines are comments. Lines are read
# fresh before every job, so the file can be edited while the launcher runs. A line
# is done when its exact text is in <queue>.done (resume-safe).
set -u
Q="$1"; DONE="$Q.done"; touch "$DONE"
HEAD_MB="${HEAD_MB:-1024}"
GATE="${GATE:-$HOME/.cache/local_runs/scripts/gpu_gate.sh}"
echo $$ > "$Q.pid"
cd "$(dirname "$0")/.."
free_mb() { nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader,nounits | head -1 \
            | awk -F', ' '{print $2-$1}'; }
while true; do
  LINE=$(grep -v '^\s*#' "$Q" | grep -v '^\s*$' | while read -r l; do grep -qxF "$l" "$DONE" || { echo "$l"; break; }; done)
  [ -z "$LINE" ] && { echo "[$(date +%T)] queue empty"; break; }
  NEED=$(echo "$LINE" | awk '{print $1}'); CMD=$(echo "$LINE" | cut -d' ' -f2-)
  # optional hold file (touch <queue>.hold to pause before the next job)
  while [ -e "$Q.hold" ]; do sleep 60; done
  echo "[$(date '+%F %T')] QUEUED (need $NEED+$HEAD_MB MiB free): $CMD"
  "$GATE" $((NEED + HEAD_MB)) bash -c "$CMD"; RC=$?
  echo "[$(date '+%F %T')] END rc=$RC: $CMD"
  echo "$LINE" >> "$DONE"
  sleep 5
done
