#!/bin/bash
# Shared work-queue logic for the Ada GPU workers (slurm/scmgrn_work.sh) and the laptop worker
# (slurm/scmgrn_local_worker.sh). Sourced after scmgrn_common.sh.
#
# Units are ordered cell-type-major (small graphs first, so the fidelity gate is checked early):
# every worker claims the next unclaimed unit of the CURRENT cell type, so all 4 Ada GPUs (plus the
# laptop GPU) work on one cell type at a time with perfect load balance, then move on together.
# A claim = atomic mkdir in $SHARE/claims/$DS (home is shared across nodes; the laptop claims over
# ssh); the unit's result JSONs are published before $SHARE/done/$DS/<unit> is written.
#
# MeVD-GRN variants (undirected TF-TF link prediction; frozen before any result):
#   H  = mevd_fm_h384 : obs : smg      headline          A1 = mevd_fm_h384 : obs : nosmg
#   A2 = mevd_fm_h384 : noobs : nosmg  label-free graphs  A3 = mevd_fm_h384_directed : obs : nosmg
SEEDS=(${SEEDS:-42 43 44 45 46})
OFFICIAL_SEEDS=(${OFFICIAL_SEEDS:-666 42 43 44 45 46})
VARIANTS=(${VARIANTS:-mevd_fm_h384:obs:smg mevd_fm_h384:obs:nosmg mevd_fm_h384:noobs:nosmg mevd_fm_h384_directed:obs:nosmg})
KINDS=${KINDS:-official mevd}                         # the laptop worker sets KINDS=mevd
if [ "$DS" = hema ]; then WORK_CTS=(HSC GMP CMP)
else WORK_CTS=(Myeloid_cells Lymphoid_cells Megakaryocytes Neuroendocrine_cells Lymphatic_endothelial_cells
               Bronchiolar_and_alveolar_epithelial_cells Stromal_cells Ciliated_epithelial_cells Vascular_endothelial_cells)
fi
[ -n "${CTS_ONLY:-}" ] && WORK_CTS=($CTS_ONLY)         # tests only
gen_units() {   # "kind ct variant seed" per line; official runs first within a cell type
  local CT S V
  for CT in "${WORK_CTS[@]}"; do
    for S in "${OFFICIAL_SEEDS[@]}"; do echo "official $CT x $S"; done
    for V in "${VARIANTS[@]}"; do for S in "${SEEDS[@]}"; do echo "mevd $CT $V $S"; done; done
  done
}
unit_name() { echo "$1__$2__${3//:/-}__$4"; }

# --- backend: where the shared claim/done state and published data live. Ada: local fs ($SHARE);
# laptop: ssh to $RAW_HOST. Defined by the driver: q_exec <cmd...>, stage_ct, publish_ct.
q_claim() { q_exec mkdir "$SHARE/claims/$DS/$1" 2>/dev/null; }
q_done() { q_exec test -e "$SHARE/done/$DS/$1"; }
q_mark_done() { q_exec bash -c "echo '$(hostname) $(date +%F_%T)' > '$SHARE/done/$DS/$1'"; }
q_fail() { q_exec bash -c "mkdir -p '$SHARE/failed/$DS' && mv '$SHARE/claims/$DS/$1' '$SHARE/failed/$DS/$1.$(date +%s)'"; }
ct_ready() { q_exec test -e "$SHARE/processed/$DS/$1/smg_pair_feats.npy" \
                 -a -e "$SHARE/processed/$DS/$1/fm_gene_embeddings.npy" -a -e "$SHARE/splits/$DS/$1"; }

run_unit() {   # run_unit <kind> <ct> <variant> <seed>   (data for <ct> already staged under $SCMGRN_ROOT_LOCAL)
  local KIND=$1 CT=$2 V=$3 S=$4 MC OBS SMG FLAG=""
  if [ "$KIND" = official ]; then
    python -u scripts/18_scmgrn_official_repro.py --config configs/scmgrn/$DS.yaml \
      --root "$SCMGRN_ROOT_LOCAL" --cell_type "$CT" --folds 1-10 --skip_existing --model_seed "$S" \
      --work_dir "$SCR/official_work/$CT"
  else
    IFS=: read -r MC OBS SMG <<<"$V"
    [ "$OBS" = obs ] && FLAG="--obs_graph"; [ "$SMG" = smg ] && FLAG="$FLAG --smg_feats"
    python -u scripts/17_scmgrn_train_eval.py --config configs/scmgrn/$DS.yaml \
      --model_config "configs/scmgrn/$MC.yaml" $FLAG --root "$SCMGRN_ROOT_LOCAL" \
      --cell_type "$CT" --seed "$S" --folds 1-10 --device "${DEVICE:-cuda:0}" --skip_existing --tmp_dir "$SCR" ${EXTRA_ARGS:-}
  fi
}

worker_loop() {
  local STAGED_CT="" waited=0 progressed KIND CT V S U
  while true; do
    progressed=0; local pending=0
    while read -r KIND CT V S; do
      [[ " $KINDS " == *" $KIND "* ]] || continue
      U=$(unit_name "$KIND" "$CT" "$V" "$S")
      q_done "$U" && continue
      pending=1
      ct_ready "$CT" || continue                                  # cell type not built yet
      q_claim "$U" || continue                                    # another worker has it
      progressed=1
      if [ "$STAGED_CT" != "$CT" ]; then
        rm -rf "$SCMGRN_ROOT_LOCAL/graph" "$SCMGRN_ROOT_LOCAL/splits" "$SCMGRN_ROOT_LOCAL/processed" "$SCMGRN_ROOT_LOCAL/results"
        stage_ct "$CT" && STAGED_CT=$CT || { echo "[worker] stage failed for $CT"; q_fail "$U"; STAGED_CT=""; continue; }
      fi
      echo "===== unit $U ($(date +%T)) ====="
      if run_unit "$KIND" "$CT" "$V" "$S" && publish_ct "$CT"; then
        q_mark_done "$U"; echo "[worker] DONE $U"
      else
        echo "[worker] FAILED $U"; q_fail "$U"
      fi
      break                                                       # re-scan from the top: current cell type first
    done < <(gen_units)
    if [ "$progressed" = 0 ]; then
      [ "$pending" = 0 ] && { echo "[worker] queue empty"; return 0; }
      [ "${WAIT_FOR_BUILD:-0}" = 1 ] || { echo "[worker] nothing claimable (others hold the rest)"; return 0; }
      echo "[worker] waiting for builds/claims ($(date +%T))"; sleep 300
    fi
  done
}
