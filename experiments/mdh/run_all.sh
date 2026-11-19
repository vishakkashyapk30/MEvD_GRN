#!/usr/bin/env bash
# Full MDH grid (README): model -> A1 -> audit -> pooled BEAR scores. Logs in logs/.
set -e
cd "$(dirname "$0")"
for ds in K562 Macrophage_S1 Macrophage_S2 iPS mESC_E7.5_rep1 mESC_E7.5_rep2 mESC_E8.5_rep1 mESC_E8.5_rep2 Naive_mESC; do
  python 03_midgrn.py --dataset $ds --seeds 42 43 44 45 46 --variants mid_full mid_dataonly mid_bce >> logs/03_$ds.log 2>&1
  echo "[03 done] $ds $(date +%T)"
done
python 01_gt_structure.py > logs/01.log 2>&1; echo "[01 done] $(date +%T)"
python 02_audit.py > logs/02.log 2>&1; echo "[02 done] $(date +%T)"
python 04_pooled.py > logs/04.log 2>&1; echo "[04 done] $(date +%T)"
