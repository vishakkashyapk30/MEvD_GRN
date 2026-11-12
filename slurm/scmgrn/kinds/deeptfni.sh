# Unit kind of the DeepTFni baseline (branch baselines): the official scMultiomeGRN Trainer with model_cls=DeepTFni
# on our folds, through the patched official copy so the unused pair histograms are not computed.
DEEPTFNI_SEEDS=(${DEEPTFNI_SEEDS:-666})
units_deeptfni() { local S; for S in "${DEEPTFNI_SEEDS[@]}"; do echo "deeptfni $1 x $S"; done; }
run_kind_deeptfni() {
  local CT=$1 S=$3
  python -u scripts/scmgrn/18_scmgrn_official_repro.py --config configs/scmgrn/$DS.yaml \
    --root "$SCMGRN_ROOT_LOCAL" --official_repo official_ablation/ScmultiomeGRN-main --cell_type "$CT" \
    --folds 1-10 --skip_existing --model_seed "$S" --method_name deeptfni --save_scores \
    --extra_kw '{"model_cls": "DeepTFni", "use_edge_attr": false}' \
    --split_kind "$SPLIT_KIND" --work_dir "$SCR/official_work/${CT}_deeptfni_$S"
}
