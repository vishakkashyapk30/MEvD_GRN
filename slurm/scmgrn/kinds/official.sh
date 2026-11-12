# Unit kind of the scMultiomeGRN end-to-end reproduction (branch scmgrn-reproduction): the unmodified official code
# (scripts/scmgrn/18) on our fold files. The official code is near-deterministic across seeds (Myeloid, 6 seeds:
# all within 0.001) and costs ~5 min/fold on the big cell types.
OFFICIAL_SEEDS=(${OFFICIAL_SEEDS:-666 42})
units_official() { local S; for S in "${OFFICIAL_SEEDS[@]}"; do echo "official $1 x $S"; done; }
run_kind_official() {
  local CT=$1 S=$3
  python -u scripts/scmgrn/18_scmgrn_official_repro.py --config configs/scmgrn/$DS.yaml \
    --root "$SCMGRN_ROOT_LOCAL" --cell_type "$CT" --folds 1-10 --skip_existing --model_seed "$S" --save_scores \
    --split_kind "$SPLIT_KIND" --work_dir "$SCR/official_work/${CT}_$S"
}
