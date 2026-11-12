# Unit kind of the GENELink baseline (branch baselines): scripts/scmgrn/32_scmgrn_genelink.py on our folds.
GENELINK_SEEDS=(${GENELINK_SEEDS:-8})
units_genelink() { local S; for S in "${GENELINK_SEEDS[@]}"; do echo "genelink $1 x $S"; done; }
run_kind_genelink() {
  local CT=$1 S=$3
  python -u scripts/scmgrn/32_scmgrn_genelink.py --config configs/scmgrn/$DS.yaml --root "$SCMGRN_ROOT_LOCAL" \
    --cell_type "$CT" --folds 1-10 --skip_existing --seed "$S" --save_scores --split_kind "$SPLIT_KIND" \
    --work_dir "$SCR/genelink_work/${CT}_$S"
}
