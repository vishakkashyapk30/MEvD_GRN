#!/bin/bash
# Local (laptop) CPU pipeline for ONE BEAR dataset, disk-frugal (s12, Ada locked):
#   fetch inputs + GTs -> label-free processed dirs (L1, L2) + Geneformer -> delete raw ATAC
#   -> baselines (in-degree x 5 seeds, coverage, |Pearson|, GRNBoost2) -> delete raw RNA
#   -> released BEAR GRNs fetched ONE AT A TIME, scored (n_rep 20, all GTs), deleted.
# Motif pair features (M3) are a separate step (needs the genome store):
#   bash slurm/bear_local_dataset.sh <DS> motif
# Usage: bash slurm/bear_local_dataset.sh <DS> [all|prep|baselines|released|motif]
set -uo pipefail
DS="$1"; STEP="${2:-all}"
REPO="$(cd "$(dirname "$0")/.." && pwd)"
ROOT="${BEAR_ROOT:-$HOME/.cache/bear/data}"
PY="${PY:-/home/vishak/miniforge3/bin/python}"
NREP="${NREP:-20}"
cd "$REPO"
log() { echo "[$(date +%T)] $DS: $*"; }
freegb() { df -BG --output=avail / | tail -1 | tr -dc 0-9; }
guard_disk() { local f; f=$(freegb); if [ "$f" -lt 3 ]; then log "ABORT: only ${f} GB free (< 3 GB rule)"; exit 3; fi; }
SPECIES=$($PY -c "import yaml;print(yaml.safe_load(open('configs/bear/datasets.yaml'))['datasets']['$DS']['species'])")
INF=$($PY -c "import yaml;print(yaml.safe_load(open('configs/bear/datasets.yaml'))['inferred_dir']['$DS'])")

if [ "$STEP" = all ] || [ "$STEP" = prep ]; then
  guard_disk
  if [ ! -e "$ROOT/processed/$DS/summary.json" ] || [ ! -e "$ROOT/processed/${DS}__L2/summary.json" ]; then
    log "fetch inputs (kept compressed) + GTs"
    $PY scripts/20_bear_download.py --root "$ROOT" --datasets "$DS" --what input --no_extract --conns 24 || exit 1
    $PY scripts/20_bear_download.py --root "$ROOT" --datasets "$DS" --what gt --conns 24 || exit 1
    guard_disk
    log "prepare L1 + L2"
    $PY scripts/21_bear_prepare.py --root "$ROOT" --datasets "$DS" --fm || exit 1
    $PY scripts/21_bear_prepare.py --root "$ROOT" --datasets "$DS" --fm --regime L2 || exit 1
  fi
  for f in "$ROOT/INPUT.DATA/$DS"/*ATAC*.csv; do [ -e "$f" ] && { log "rm $(basename "$f")"; rm -f "$f"; }; done
fi

if [ "$STEP" = all ] || [ "$STEP" = baselines ]; then
  if ! ls "$ROOT/INPUT.DATA/$DS"/*RNA*.csv >/dev/null 2>&1 && [ ! -e "$ROOT/zips/INPUT.DATA.zip.sparse" ]; then
    $PY scripts/20_bear_download.py --root "$ROOT" --datasets "$DS" --what input --no_extract --conns 24
  fi
  log "baselines"
  $PY scripts/23_bear_baselines.py --root "$ROOT" --dataset "$DS" \
      --baselines indegree coverage pearson grnboost2 --n_jobs 12 || exit 1
  rm -f "$ROOT/INPUT.DATA/$DS"/*RNA*.csv "$ROOT/zips/INPUT.DATA.zip.sparse"
fi

if [ "$STEP" = all ] || [ "$STEP" = released ]; then
  log "score baselines (+ any MeVD runs present)"
  $PY scripts/24_bear_score.py --root "$ROOT" --datasets "$DS" --ours all --n_rep "$NREP"
  for M in LINGER CellOracle SCENIC+ Pando Pando_xgb FigR TRIPOD GRaNIE DIRECTNET; do
    guard_disk
    if ls "$ROOT/results/$DS/scores/${M}__"*.json >/dev/null 2>&1 && [ -z "${FORCE_RELEASED:-}" ]; then
      n=$(ls "$ROOT/results/$DS/scores/${M}__"*.json | wc -l); log "skip released $M ($n cached)"; continue; fi
    log "released $M"
    $PY scripts/20_bear_download.py --root "$ROOT" --datasets "$DS" --what inferred --methods "$M" --conns 24 || continue
    $PY scripts/24_bear_score.py --root "$ROOT" --datasets "$DS" --released "$M" --n_rep "$NREP" \
        --paper_csv configs/bear/paper_numbers/bear_grn_auroc_auprc_tidy.csv
    rm -rf "$ROOT/INFERRED.GRNS/$INF/$M"
  done
  $PY scripts/24_bear_score.py --root "$ROOT" --datasets "$DS" --compile_only \
      --paper_csv configs/bear/paper_numbers/bear_grn_auroc_auprc_tidy.csv
fi

if [ "$STEP" = motif ]; then
  STORE="${GENOME_STORE:-$HOME/.cache/bear/genome/$([ "$SPECIES" = human ] && echo hg38 || echo mm10)_2bit}"
  JASPAR="${JASPAR:-$HOME/.cache/bear/genome/JASPAR2024_CORE_vertebrates_redundant_pfms_jaspar.txt}"
  for R in L1 L2; do
    log "motif features $R"
    $PY scripts/21_bear_prepare.py --root "$ROOT" --datasets "$DS" --regime "$R" --motif \
        --genome_store "$STORE" --jaspar "$JASPAR" || exit 1
  done
fi
log "done ($STEP); $(freegb) GB free"
