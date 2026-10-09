"""After the s12.6 selection: write the frozen-headline configs and append the run lines to the
local GPU queue (slurm/bear_gpu_queue.sh). Order (s12.9): (2) seeds 42-44 on all datasets, ChIP-only
datasets first; (3) seeds 45-46; (4) L2; (5) ablations (K562, Macrophage_S1, mESC_E7.5_rep1; 3 seeds).

  python -m src.benchmarks.bear_queue --variant fm_h384_hubmotif [--dry]
"""
import argparse
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CFG = REPO / "configs" / "bear"
LOGS = Path.home() / ".cache/bear/logs"
QUEUE = Path.home() / ".cache/bear/gpu_queue.txt"
ROOT = Path.home() / ".cache/bear/data"
ORDER_DS = ["Macrophage_S1", "Macrophage_S2", "iPS", "K562", "mESC_E7.5_rep1", "mESC_E7.5_rep2",
            "mESC_E8.5_rep1", "mESC_E8.5_rep2", "Naive_mESC"]
ABL_DS = ["K562", "Macrophage_S1", "mESC_E7.5_rep1"]
PY = "/home/vishak/miniforge3/bin/python -u scripts/22_bear_train_eval.py --root %s --device cuda:0 " \
     "--skip_existing --tmp_dir /home/vishak/.cache/bear/tmp" % ROOT


def write_cfgs(base: str):
    """H = the frozen headline (= configs/bear/mevd_<base>.yaml); derived variants inherit it."""
    (CFG / "mevd_H.yaml").write_text(f"# Frozen headline (s12.6 selection): {base}\ninherit: mevd_{base}.yaml\n")
    (CFG / "mevd_H_L2.yaml").write_text("# Headline in regime L2 (compendium labels only, s12.5)\ninherit: mevd_H.yaml\n"
                                        "data:\n  evidence_tiers: [compendium]\n  tier_hierarchy: [compendium]\n"
                                        "curriculum:\n  main_curriculum_tiers: [compendium]\nbear:\n  regime: L2\n")
    (CFG / "mevd_H_noatac.yaml").write_text("# Ablation: headline without ATAC (features + openness zeroed"
                                            "; motif pair features are peak-derived, so also zeroed)\n"
                                            "inherit: mevd_H.yaml\nbear:\n  zero_atac: true\n  zero_motif: true\n")
    (CFG / "mevd_H_nofm.yaml").write_text("# Ablation: headline without Geneformer (same h384/l2)\n"
                                          "inherit: mevd_H.yaml\nmodel:\n  use_fm: false\n")
    (CFG / "mevd_H_seq.yaml").write_text("# Ablation: sequential curriculum (ChIP stage -> KO stage) instead of all_at_once\n"
                                         "inherit: mevd_H.yaml\ncurriculum:\n  protocol: sequential\n  stages:\n"
                                         "    - {name: Stage1_ChIP, evidence_tier: localization, n_epochs: 40, lr: 1.0e-3, neg_ratio: 2, freeze_encoder: false}\n"
                                         "    - {name: Stage2_KO, evidence_tier: perturbation, n_epochs: 20, lr: 3.0e-4, neg_ratio: 2, freeze_encoder: false}\n")


def line(ds, variant, seed, need=3500, extra=""):
    log = LOGS / f"{ds}_{variant}_s{seed}.log"
    return f"{need} {PY} --dataset {ds} --model_config configs/bear/mevd_{variant}.yaml --seed {seed} {extra}> {log} 2>&1"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", required=True)
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    write_cfgs(a.variant)
    mod = "--modality_check "
    lines = ["# --- frozen headline runs (s12.9 order) ---"]
    for seeds in ([42, 43, 44], [45, 46]):
        for ds in ORDER_DS:
            for s in seeds:
                lines.append(line(ds, "H", s, extra=mod if (s == 42) else ""))
    lines.append("# --- L2 regime ---")
    for s in [42, 43, 44, 45, 46]:
        for ds in ORDER_DS:
            lines.append(line(ds, "H_L2", s))
    lines.append("# --- ablations (3 seeds) ---")
    for v in ["H_noatac", "H_nofm", "H_seq"]:
        for ds in ABL_DS:
            if v == "H_seq" and ds == "Macrophage_S1":
                continue                       # one label tier only: no curriculum to ablate
            for s in [42, 43, 44]:
                lines.append(line(ds, v, s))
    if a.dry:
        print("\n".join(lines))
        return
    with open(QUEUE, "a") as f:
        f.write("\n".join(lines) + "\n")
    print(f"appended {len(lines)} lines to {QUEUE}")


if __name__ == "__main__":
    main()
