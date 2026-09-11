#!/usr/bin/env python
"""Zero-shot cross-cell-type transfer evaluation for MEvD-GRN.

Extends plan Sec 8.5 (originally just K562<->ESC) to SC-MO-GRN-DB's
"partial-evidence" cell lines that only have a LOCALIZATION reference
network (HepG2, BJ, GM12878, Macrophage, MCF7, HSC -- see
scripts/09_download_transfer_data.sh for the full catalog and the
documented H1/DC exclusions).

A model trained on the source cell type (default: K562, curriculum
localization -> perturbation) is applied WITH NO RETRAINING to a target
cell type's own scRNA/scATAC-derived graph, and scored against that target
cell type's OWN localization edges (100% held out -- nothing from the
target cell type is EVER seen during training: different genes, different
accessibility landscape, different cells, different study). This is a
strictly harder generalization test than the within-cell-type dual_evidence
zero-shot inference result, since even the gene universe differs.

NOTE on baseline comparability: scMultiomeGRN's message-passing graph is,
by architectural design, derived from the reference network of whatever
cell type it is run on (see src/baselines/scmultiomegrn_wrapper.py
docstring) and the unsupervised baselines (GRNBoost2/RegDiffusion/GMF-GAE)
have no pretrained weights at all -- neither can be transferred zero-shot
in the sense used here. For a baseline comparison point at each transfer
cell type, run scripts/02_preprocess.py + scripts/05_run_baselines.py
directly against that cell type's own config; those baselines get to fit
fresh on the target cell type's own train split, which is an easier task
than MEvD-GRN's true zero-shot transfer, so report that asymmetry alongside
any comparison.

Usage:
  python scripts/08_transfer_eval.py --transfer_config configs/transfer/hepg2.yaml \
      --source_config configs/k562.yaml \
      --checkpoint results/checkpoints/K562/final_model.pt --device cuda:0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.data.dataset import load_celltype_data
from src.data.graph_builder import _sample_cols
from src.models.mevd_grn import MEvDGRN
from src.training.trainer import MEvDTrainer
from src.utils.io import ensure_dir, load_config, save_json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--transfer_config", required=True)
    ap.add_argument("--source_config", default="configs/k562.yaml",
                    help="Config the checkpoint was TRAINED with (for model architecture dims).")
    ap.add_argument("--checkpoint", default="results/checkpoints/K562/final_model.pt")
    ap.add_argument("--device", default=None)
    ap.add_argument("--neg_eval_ratio", type=int, default=5)
    args = ap.parse_args()

    tcfg = load_config(args.transfer_config)
    scfg = load_config(args.source_config)
    cell_type = tcfg["cell_type"]
    source_cell_type = scfg["cell_type"]
    device = args.device or scfg["hardware"]["device"]
    if device.startswith("cuda") and not torch.cuda.is_available():
        device = "cpu"

    data = load_celltype_data(tcfg["paths"]["processed_dir"], ["localization"])
    if "localization" not in data.evidence or data.evidence["localization"].shape[1] == 0:
        raise RuntimeError(f"{cell_type}: no localization evidence found "
                           "-- run scripts/02_preprocess.py on the transfer config first")

    pos = data.evidence["localization"]
    gen = torch.Generator().manual_seed(int(tcfg["data"]["seed"]))
    n_neg = args.neg_eval_ratio * pos.shape[1]
    neg = _sample_cols(data.negative_pool, n_neg, gen)

    m = scfg["model"]
    # NOTE: must mirror every architecture flag scripts/03_train.py passes,
    # or the source checkpoint will silently load into the WRONG architecture.
    model = MEvDGRN(rna_in_dim=m["rna_in_dim"], atac_in_dim=m["atac_in_dim"],
                    hidden_dim=m["hidden_dim"], n_gnn_layers=m["n_gnn_layers"],
                    dropout=m["dropout"], integration=m.get("integration", "role_aware"),
                    graph_mode=m.get("graph_mode", "both"),
                    use_edge_mlp=bool(m.get("use_edge_mlp", False)),
                    combine_mode=m.get("combine_mode", "sum"),
                    use_fm=bool(m.get("use_fm", False)),
                    fm_in_dim=int(m.get("fm_in_dim", 768)))

    trainer = MEvDTrainer(model, data, scfg, device)
    trainer.load_checkpoint(args.checkpoint)
    result = trainer.evaluate_split({"pos": pos, "neg": neg})

    print(f"[transfer] {source_cell_type} -> {cell_type} (localization, "
          f"{pos.shape[1]} positives + {neg.shape[1]} negatives, "
          "ZERO-SHOT, no retraining/fine-tuning on target)", flush=True)
    print(f"  AUPR={result['aupr']:.4f} AUROC={result['auroc']:.4f} "
          f"EP={result['early_precision']:.4f} EPR={result['epr']:.2f}", flush=True)

    ensure_dir("results/transfer")
    # Tag includes use_fm -- otherwise a with-FM and a without-FM source
    # checkpoint targeting the same cell type pair silently overwrite each
    # other's result file (see the same fix in scripts/12_joint_train.py).
    fm_tag = "_fm" if model.use_fm else "_nofm"
    save_json({"source_cell_type": source_cell_type, "target_cell_type": cell_type,
               "tier": "localization", "mode": "zero_shot_transfer", "use_fm": model.use_fm,
               "n_positives": int(pos.shape[1]), "n_negatives": int(neg.shape[1]),
               "result": result},
              f"results/transfer/{source_cell_type}_to_{cell_type}{fm_tag}.json")


if __name__ == "__main__":
    main()
