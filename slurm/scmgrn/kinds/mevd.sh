# Unit kinds of MeVD-GRN (main). Sourced by slurm/scmgrn/scmgrn_worker_lib.sh.
#   mevd   : one unit per (cell type, variant, seed); variant = <model config>:<obs|noobs>:<smg|smgnode|smgpair|nosmg>
#            (obs = the fold's training links as a message-passing relation; smg = scMultiomeGRN's own inputs).
#            The presented model is mevd_fm_h384:obs:nosmg.
#   mevddm : re-run of the presented model into mevd_fm_h384_obs_dm with saved score matrices (the earliest fold
#            files predate the degree-matched metrics).
DM_SEEDS=(${DM_SEEDS:-42 43 44 45 46})
units_mevd() { local V S; for V in "${VARIANTS[@]}"; do for S in "${SEEDS[@]}"; do echo "mevd $1 $V $S"; done; done; }
units_mevddm() { local S; for S in "${DM_SEEDS[@]}"; do echo "mevddm $1 mevd_fm_h384:obs:nosmg $S"; done; }
run_kind_mevd() {
  local CT=$1 V=$2 S=$3 MC OBS SMG FLAG=""
  IFS=: read -r MC OBS SMG <<<"$V"
  [ "$OBS" = obs ] && FLAG="--obs_graph"
  case "$SMG" in smg) FLAG="$FLAG --smg_feats";; smgnode) FLAG="$FLAG --smg_feats --smg_parts node";; smgpair) FLAG="$FLAG --smg_feats --smg_parts pair";; esac
  python -u scripts/scmgrn/17_scmgrn_train_eval.py --config configs/scmgrn/$DS.yaml \
    --model_config "configs/scmgrn/$MC.yaml" $FLAG --root "$SCMGRN_ROOT_LOCAL" \
    --cell_type "$CT" --seed "$S" --folds 1-10 --device "${DEVICE:-cuda:0}" --skip_existing --tmp_dir "$SCR" \
    --split_kind "$SPLIT_KIND" ${EXTRA_ARGS:-}
}
run_kind_mevddm() {
  local CT=$1 S=$3
  python -u scripts/scmgrn/17_scmgrn_train_eval.py --config configs/scmgrn/$DS.yaml \
    --model_config configs/scmgrn/mevd_fm_h384.yaml --obs_graph --variant fm_h384_obs_dm --save_scores \
    --root "$SCMGRN_ROOT_LOCAL" --cell_type "$CT" --seed "$S" --folds 1-10 --device "${DEVICE:-cuda:0}" \
    --skip_existing --tmp_dir "$SCR" --split_kind "$SPLIT_KIND"
}
