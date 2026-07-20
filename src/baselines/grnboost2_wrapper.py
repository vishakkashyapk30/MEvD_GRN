"""GRNBoost2 baseline (plan 9.1).

GRNBoost2 = for each target gene, fit a gradient-boosting regressor predicting its
expression from all TF expressions; feature importances are the TF->target edge
weights. The original `arboreto` package is unmaintained and incompatible with
modern dask (legacy DataFrame removed), so we implement the identical algorithm
directly with LightGBM (a gradient-boosting backend). RNA-only.
"""
from __future__ import annotations

from typing import Dict, Iterable, List, Optional

import numpy as np

from src.baselines.common import DictScorer, load_expression_df


def run_grnboost2(rna_path: str, gene_index: Dict[str, int], tf_names: List[str],
                  pcfg: dict, seed: int = 42, n_workers: int = 4,
                  target_symbols: Optional[Iterable[str]] = None,
                  n_estimators: int = 60) -> DictScorer:
    try:
        from lightgbm import LGBMRegressor
    except ImportError as e:  # pragma: no cover
        raise ImportError("lightgbm not installed — `pip install lightgbm`") from e

    expr = load_expression_df(rna_path, gene_index, pcfg)              # cells × genes
    tf_present = [t for t in tf_names if t in expr.columns]
    X_tf = expr[tf_present].to_numpy(np.float32)                       # cells × n_tf
    tf_pos = {t: k for k, t in enumerate(tf_present)}                  # tf symbol -> col

    targets = list(target_symbols) if target_symbols is not None else list(expr.columns)
    targets = [g for g in targets if g in expr.columns]
    print(f"[grnboost2] {expr.shape[0]} cells, {len(tf_present)} TFs, "
          f"{len(targets)} target genes (LightGBM backend)", flush=True)

    lookup: Dict[tuple, float] = {}
    for n, g in enumerate(targets):
        y = expr[g].to_numpy(np.float32)
        if y.std() == 0:
            continue
        # Regularization tree-GBM; drop the target's own column if it is a TF (no self-edge)
        cols = X_tf
        drop = tf_pos.get(g)
        feat_idx = [k for k in range(X_tf.shape[1]) if k != drop]
        if drop is not None:
            cols = X_tf[:, feat_idx]
        model = LGBMRegressor(n_estimators=n_estimators, num_leaves=15, learning_rate=0.05,
                              subsample=0.9, colsample_bytree=0.7, n_jobs=n_workers,
                              random_state=seed, verbosity=-1)
        model.fit(cols, y)
        imp = model.feature_importances_.astype(np.float64)
        s = imp.sum()
        if s > 0:
            imp = imp / s                                             # normalize per target
        gj = gene_index[g]
        for local_k, w in zip(feat_idx, imp):
            if w > 0:
                lookup[(gene_index[tf_present[local_k]], gj)] = float(w)
        if (n + 1) % 2000 == 0:
            print(f"[grnboost2] fit {n+1}/{len(targets)} targets", flush=True)
    return DictScorer(lookup)
