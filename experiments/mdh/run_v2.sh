#!/usr/bin/env bash
# v2 grid (README "v2"), then the audit and pooled scores for everything. Logs in logs/.
set -e
cd "$(dirname "$0")"
for ds in K562 Macrophage_S1 Macrophage_S2 iPS mESC_E7.5_rep1 mESC_E7.5_rep2 mESC_E8.5_rep1 mESC_E8.5_rep2 Naive_mESC; do
  python 06_midgrn_v2.py --dataset $ds --seeds 42 43 44 45 46 >> logs/06_$ds.log 2>&1
  echo "[06 done] $ds $(date +%T)"
done
python 02_audit.py > logs/02.log 2>&1; echo "[02 done] $(date +%T)"
python 04_pooled.py > logs/04.log 2>&1; echo "[04 done] $(date +%T)"
