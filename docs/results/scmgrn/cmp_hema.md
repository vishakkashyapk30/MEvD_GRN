# scMultiomeGRN benchmark: hema

## Headline: held-out test fold (zero leakage, same for both methods)

1:1 positives/negatives of the official fold; training edges are never scored.

AUROC / AUPR per cell type, mean over 10 folds then over seeds.

| Cell type | MeVD-GRN H (undirected, +train adj, +smg inputs) | A1 (no smg inputs) | A2 (label-free graphs only) | A3 (old directed model, symmetrised) | scMultiomeGRN official code (re-run) | degree_logreg | degree_logsum | degree_product | mevd_abl_base_obs_smg | mevd_abl_coexpr_obs_smg | mevd_abl_fl_atac_noobs | mevd_abl_fl_fm_noobs | mevd_abl_fl_gnn_noobs | mevd_abl_min_noobs | mevd_abl_min_noobs_smgnode | mevd_abl_noatac_obs_smg | mevd_abl_nofm_obs_smg | mevd_abl_nognn_obs_smg | mevd_abl_tfcand_obs_smg | mevd_fm_h384_obs_smgnode | mevd_fm_h384_obs_smgpair |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HSC | 0.934 / 0.940 | 0.934 / 0.941 | 0.931 / 0.937 | 0.932 / 0.940 | 0.922 / 0.921 | 0.923 / 0.933 | 0.919 / 0.929 | 0.919 / 0.929 | 0.934 / 0.941 | 0.930 / 0.937 | 0.874 / 0.885 | 0.931 / 0.938 | 0.887 / 0.895 | 0.548 / 0.550 | 0.926 / 0.932 | 0.935 / 0.941 | 0.931 / 0.938 | 0.930 / 0.937 | 0.931 / 0.938 | 0.934 / 0.941 | 0.935 / 0.942 |
| CMP | 0.925 / 0.932 | 0.925 / 0.932 | 0.924 / 0.931 | 0.924 / 0.932 | 0.914 / 0.909 | 0.918 / 0.925 | 0.915 / 0.922 | 0.915 / 0.922 | 0.924 / 0.930 | 0.925 / 0.932 | 0.867 / 0.881 | 0.919 / 0.926 | 0.899 / 0.906 | 0.637 / 0.652 | 0.915 / 0.922 | 0.925 / 0.932 | 0.924 / 0.931 | 0.923 / 0.929 | 0.923 / 0.931 | 0.925 / 0.932 | 0.925 / 0.932 |
| GMP | 0.925 / 0.934 | 0.925 / 0.935 | 0.922 / 0.934 | 0.924 / 0.935 | 0.921 / 0.926 | 0.918 / 0.931 | 0.914 / 0.926 | 0.914 / 0.926 | 0.925 / 0.935 | 0.921 / 0.933 | 0.841 / 0.863 | 0.921 / 0.930 | 0.855 / 0.866 | 0.598 / 0.620 | 0.922 / 0.931 | 0.924 / 0.934 | 0.923 / 0.933 | 0.923 / 0.934 | 0.922 / 0.934 | 0.925 / 0.935 | 0.925 / 0.935 |
| average | 0.928 / 0.936 | 0.928 / 0.936 | 0.926 / 0.934 | 0.927 / 0.936 | 0.919 / 0.919 | 0.920 / 0.930 | 0.916 / 0.926 | 0.916 / 0.926 | 0.928 / 0.936 | 0.925 / 0.934 | 0.861 / 0.876 | 0.924 / 0.931 | 0.880 / 0.889 | 0.595 / 0.607 | 0.921 / 0.929 | 0.928 / 0.936 | 0.926 / 0.934 | 0.925 / 0.933 | 0.925 / 0.934 | 0.928 / 0.936 | 0.928 / 0.937 |

## Paper protocol: full N x N matrix (CONTAINS the training edges; not leakage-free)

Only for comparison with the published numbers (paper average AUROC None / AUPR None, scMultiomeGRN; baselines in lung.yaml).

AUROC / AUPR per cell type, mean over 10 folds then over seeds.

| Cell type | MeVD-GRN H (undirected, +train adj, +smg inputs) | A1 (no smg inputs) | A2 (label-free graphs only) | A3 (old directed model, symmetrised) | scMultiomeGRN official code (re-run) | degree_logreg | degree_logsum | degree_product | mevd_abl_base_obs_smg | mevd_abl_coexpr_obs_smg | mevd_abl_fl_atac_noobs | mevd_abl_fl_fm_noobs | mevd_abl_fl_gnn_noobs | mevd_abl_min_noobs | mevd_abl_min_noobs_smgnode | mevd_abl_noatac_obs_smg | mevd_abl_nofm_obs_smg | mevd_abl_nognn_obs_smg | mevd_abl_tfcand_obs_smg | mevd_fm_h384_obs_smgnode | mevd_fm_h384_obs_smgpair |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HSC | 0.945 / 0.790 | 0.946 / 0.793 | 0.943 / 0.781 | 0.945 / 0.793 | 0.927 / 0.696 | 0.932 / 0.764 | 0.929 / 0.753 | 0.930 / 0.753 | 0.946 / 0.795 | 0.943 / 0.780 | 0.911 / 0.704 | 0.948 / 0.793 | 0.916 / 0.712 | 0.546 / 0.164 | 0.949 / 0.792 | 0.945 / 0.793 | 0.944 / 0.787 | 0.946 / 0.793 | 0.944 / 0.783 | 0.945 / 0.791 | 0.946 / 0.795 |
| CMP | 0.945 / 0.800 | 0.946 / 0.805 | 0.946 / 0.802 | 0.947 / 0.811 | 0.927 / 0.692 | 0.930 / 0.767 | 0.928 / 0.757 | 0.928 / 0.757 | 0.949 / 0.810 | 0.946 / 0.804 | 0.925 / 0.743 | 0.952 / 0.821 | 0.940 / 0.769 | 0.658 / 0.289 | 0.957 / 0.829 | 0.945 / 0.802 | 0.947 / 0.809 | 0.946 / 0.804 | 0.946 / 0.802 | 0.946 / 0.803 | 0.946 / 0.804 |
| GMP | 0.943 / 0.784 | 0.943 / 0.784 | 0.942 / 0.779 | 0.945 / 0.797 | 0.928 / 0.681 | 0.928 / 0.758 | 0.926 / 0.751 | 0.926 / 0.751 | 0.945 / 0.791 | 0.942 / 0.778 | 0.898 / 0.680 | 0.949 / 0.800 | 0.881 / 0.677 | 0.598 / 0.219 | 0.948 / 0.785 | 0.943 / 0.787 | 0.944 / 0.788 | 0.945 / 0.795 | 0.943 / 0.780 | 0.943 / 0.785 | 0.943 / 0.786 |
| average | 0.944 / 0.791 | 0.945 / 0.794 | 0.944 / 0.787 | 0.946 / 0.800 | 0.927 / 0.690 | 0.930 / 0.763 | 0.928 / 0.753 | 0.928 / 0.753 | 0.947 / 0.799 | 0.944 / 0.787 | 0.911 / 0.709 | 0.950 / 0.805 | 0.912 / 0.719 | 0.601 / 0.224 | 0.951 / 0.802 | 0.945 / 0.794 | 0.945 / 0.795 | 0.946 / 0.797 | 0.944 / 0.788 | 0.945 / 0.793 | 0.945 / 0.795 |

## Fold-paired comparison: mevd_fm_h384_obs_smg minus official_scmultiomegrn

| set / metric | mean diff | frac. folds better | Wilcoxon p | n pairs |
|---|---|---|---|---|
| test_aupr | +0.0171 | 0.90 | 6.1e-08 | 30 |
| test_auroc | +0.0088 | 0.87 | 8.3e-07 | 30 |
| all_aupr | +0.1016 | 1.00 | 1.9e-09 | 30 |
| all_auroc | +0.0170 | 1.00 | 1.9e-09 | 30 |

