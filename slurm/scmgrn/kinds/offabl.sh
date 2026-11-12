# Unit kind of the scMultiomeGRN ablations (branch ablations): the patched official copy
# (scripts/scmgrn/23_patch_official_ablation.py), one unit per (cell type, tag, seed). "full" = no switch (the
# patched copy's default equals the original code; it also stores score matrices for the degree-matched tests).
OFFABL_TAGS=(${OFFABL_TAGS:-full no_prior string_prior no_edge atac_only scrna_only plus_geneformer})
OFFABL_SEEDS=(${OFFABL_SEEDS:-666})
offabl_kw() {
  case "$1" in
    full) echo '{}';;
    no_prior) echo '{"mp_graph": "none"}';;
    string_prior) echo '{"mp_graph": "external", "mp_adj_file": "{graph_dir}/x_string_graph.txt"}';;
    no_edge) echo '{"use_edge_attr": false}';;
    atac_only) echo '{"modalities": ["atac"]}';;
    scrna_only) echo '{"modalities": ["scrna"]}';;
    plus_geneformer) echo '{"gf_from_processed": true}';;
  esac
}
units_offabl() { local V S; for V in "${OFFABL_TAGS[@]}"; do for S in "${OFFABL_SEEDS[@]}"; do echo "offabl $1 $V $S"; done; done; }
run_kind_offabl() {
  local CT=$1 V=$2 S=$3
  python -u scripts/scmgrn/18_scmgrn_official_repro.py --config configs/scmgrn/$DS.yaml \
    --root "$SCMGRN_ROOT_LOCAL" --official_repo official_ablation/ScmultiomeGRN-main --cell_type "$CT" \
    --folds 1-10 --skip_existing --model_seed "$S" --ablation "$V" --extra_kw "$(offabl_kw "$V")" --save_scores \
    --split_kind "$SPLIT_KIND" --work_dir "$SCR/official_work/${CT}_${V}_$S"
}
