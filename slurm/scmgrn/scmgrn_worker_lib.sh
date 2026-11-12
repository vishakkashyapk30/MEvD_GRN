#!/bin/bash
# Shared work-queue logic for the Ada GPU workers (slurm/scmgrn/scmgrn_work.sh) and the laptop worker
# (slurm/scmgrn/scmgrn_local_worker.sh). Sourced after scmgrn_common.sh.
#
# Units are ordered cell-type-major (small graphs first, so the fidelity gate is checked early):
# every worker claims the next unclaimed unit of the CURRENT cell type, so all 4 Ada GPUs (plus the
# laptop GPU) work on one cell type at a time with perfect load balance, then move on together.
# A claim = atomic mkdir in $SHARE/claims/$DS (home is shared across nodes; the laptop claims over
# ssh); the unit's result JSONs are published before $SHARE/done/$DS/<unit> is written.
#
# Unit kinds are plug-ins: every slurm/scmgrn/kinds/<kind>.sh defines units_<kind> <ct> (one "kind ct variant seed"
# line per unit) and run_kind_<kind> <ct> <variant> <seed>. main ships kinds/mevd.sh; the scmgrn-reproduction,
# ablations and baselines branches add theirs. KINDS lists the kinds to run, in order, within each cell type.
SEEDS=(${SEEDS:-42 43 44 45 46})
VARIANTS=(${VARIANTS:-mevd_fm_h384:obs:nosmg})
KINDS=${KINDS:-mevd}
# SPLIT_KIND=official (scMultiomeGRN's folds) | tfholdout (ours: whole TFs held out; scripts/scmgrn/30). Non-official
# splits write to results/<ds>_<kind>/ and their units get a __<kind> tag, so done markers never collide.
# UNIT_TAG: extra tag for re-runs of finished units under new settings (e.g. abldm = ablations re-run for the
# degree-matched metric, with EXTRA_ARGS="--variant_suffix _dm").
export SPLIT_KIND=${SPLIT_KIND:-official}
UNIT_TAG=${UNIT_TAG:-}
[ -z "$UNIT_TAG" ] && [ "$SPLIT_KIND" != official ] && UNIT_TAG=$SPLIT_KIND
if [ "$DS" = hema ]; then WORK_CTS=(HSC GMP CMP)
else WORK_CTS=(Myeloid_cells Lymphoid_cells Megakaryocytes Neuroendocrine_cells Lymphatic_endothelial_cells
               Bronchiolar_and_alveolar_epithelial_cells Ciliated_epithelial_cells Vascular_endothelial_cells Stromal_cells)
fi
[ -n "${CTS_ONLY:-}" ] && WORK_CTS=($CTS_ONLY)         # tests only
for _k in "${REPO:-$PWD}"/slurm/scmgrn/kinds/*.sh; do [ -e "$_k" ] && source "$_k"; done
for _k in $KINDS; do
  declare -F "units_$_k" >/dev/null && declare -F "run_kind_$_k" >/dev/null || {
    echo "[worker] unknown unit kind '$_k': no slurm/scmgrn/kinds/$_k.sh on this branch"; return 1 2>/dev/null || exit 1; }
done
gen_units() {   # "kind ct variant seed" per line, cell-type-major, kinds in KINDS order
  local CT K
  for CT in "${WORK_CTS[@]}"; do for K in $KINDS; do "units_$K" "$CT"; done; done
}
unit_name() { echo "$1__$2__${3//:/-}__$4${UNIT_TAG:+__$UNIT_TAG}"; }

# --- backend: where the shared claim/done state and published data live. Ada: local fs ($SHARE);
# laptop: ssh to $RAW_HOST. Defined by the driver: q_exec <cmd...>, stage_ct, publish_ct.
q_claim() { q_exec mkdir "$SHARE/claims/$DS/$1" 2>/dev/null; }
q_done() { q_exec test -e "$SHARE/done/$DS/$1"; }
q_mark_done() { q_exec bash -c "echo '$(hostname) $(date +%F_%T)' > '$SHARE/done/$DS/$1'"; }
q_fail() { q_exec bash -c "mkdir -p '$SHARE/failed/$DS' && mv '$SHARE/claims/$DS/$1' '$SHARE/failed/$DS/$1.$(date +%s)'"; }
ct_ready() { q_exec test -e "$SHARE/processed/$DS/$1/READY"; }      # written last by scmgrn_build.sh

run_unit() {   # run_unit <kind> <ct> <variant> <seed>   (data for <ct> already staged under $SCMGRN_ROOT_LOCAL)
  local KIND=$1; shift
  "run_kind_$KIND" "$@"
}

worker_loop() {
  local STAGED_CT="" waited=0 progressed KIND CT V S U
  declare -A TRIES=()                                             # a unit that fails twice in this worker is skipped
  while true; do
    progressed=0; local pending=0
    while read -r KIND CT V S; do
      [[ " $KINDS " == *" $KIND "* ]] || continue
      U=$(unit_name "$KIND" "$CT" "$V" "$S")
      q_done "$U" && continue
      [ "${TRIES[$U]:-0}" -ge 2 ] && continue
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
        echo "[worker] FAILED $U"; q_fail "$U"; TRIES[$U]=$(( ${TRIES[$U]:-0} + 1 ))
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
