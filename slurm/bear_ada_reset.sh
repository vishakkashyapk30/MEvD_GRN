#!/bin/bash
# Release stale claims (claimed, no done marker, no BEAR worker running). Run on the gateway
# (short command): bash slurm/bear_ada_reset.sh
B="${BEAR_V2:-$HOME/bear_v2}"
if squeue -u "$USER" -h -o %j | grep -q '^bear-work'; then echo "bear-work jobs are running; not resetting"; exit 1; fi
for c in "$B"/claims/*; do u=$(basename "$c"); [ -e "$B/done/$u" ] || { echo "release $u"; rm -rf "$c"; }; done
