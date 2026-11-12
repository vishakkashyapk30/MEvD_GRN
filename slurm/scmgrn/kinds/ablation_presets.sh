# MeVD-GRN ablation presets (branch ablations): ABLATION=2 = single-component ablations of the presented model
# (mevd_fm_h384:obs:nosmg); ABLATION=1 = the earlier set on top of the model with scMultiomeGRN inputs.
if [ "${ABLATION:-}" = 2 ]; then
  VARIANTS=(mevd_abl_nofm:obs:nosmg mevd_abl_noatac:obs:nosmg mevd_abl_nognn:obs:nosmg mevd_abl_coexpr:obs:nosmg mevd_abl_tfcand:obs:nosmg mevd_abl_base:obs:nosmg)
  KINDS=mevd
elif [ -n "${ABLATION:-}" ]; then
  VARIANTS=(mevd_abl_nofm:obs:smg mevd_abl_noatac:obs:smg mevd_abl_nognn:obs:smg mevd_abl_coexpr:obs:smg mevd_abl_tfcand:obs:smg mevd_abl_base:obs:smg mevd_fm_h384:obs:smgnode mevd_fm_h384:obs:smgpair mevd_abl_min:noobs:nosmg mevd_abl_fl_atac:noobs:nosmg mevd_abl_fl_fm:noobs:nosmg mevd_abl_fl_gnn:noobs:nosmg mevd_abl_min:noobs:smgnode)
  KINDS=mevd
fi
