#!/usr/bin/env python
"""A1 (P1): how much of each ground truth is explained by node marginals alone?

For every dataset x GT: fit the two-way logistic (beta-model) logit P(Y_ij) = c + a_i + b_j on the GT matrix
(GT TFs x measured genes, no self pairs) and report full-universe AUROC of
  marginal_oracle_loo  logit(TF rate) + logit(target rate over the OTHER GT TFs) (README D2)
  row_only  TF positive rate   col_only_loo  target GT in-degree over the other TFs
  beta_model_insample  in-sample two-way logistic (self-inclusion; upper bound only)
next to full-universe AUROC of methods (unscored pairs = 0): in-degree (TF-disjoint, seed 42), MeVD-GRN M3
(seed 42), every released method. Output: results/a1_gt_structure.csv
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import torch

import mdh_lib as L


def beta_model(Y: np.ndarray, mask: np.ndarray, dev: str) -> np.ndarray:
    y = torch.tensor(Y, dtype=torch.float32, device=dev)
    m = torch.tensor(mask, dtype=torch.float32, device=dev)
    R, n = Y.shape
    a = torch.zeros(R, 1, device=dev, requires_grad=True)
    b = torch.zeros(1, n, device=dev, requires_grad=True)
    opt = torch.optim.LBFGS([a, b], lr=1, max_iter=200, line_search_fn="strong_wolfe")

    def closure():
        opt.zero_grad()
        z = a + b
        loss = (torch.nn.functional.binary_cross_entropy_with_logits(z, y, reduction="none") * m).sum() / m.sum()
        loss = loss + 1e-6 * (a.pow(2).sum() + b.pow(2).sum())
        loss.backward()
        return loss
    opt.step(closure)
    return (a + b).detach().cpu().numpy()


def main():
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    rows_out = []
    for ds in L.DATASETS:
        D = L.load_ds(ds, features=False)
        rel = {m: L.read_released(ds, m) for m in L.RELEASED}
        ours = {}
        for name, path in [("indegree", "baseline_indegree/seed42"), ("MeVD-GRN", "mevd_fm_h384_hubmotif/seed42"),
                           ("pearson", "baseline_pearson/seed0"), ("grnboost2", "baseline_grnboost2/seed0")]:
            f = L.ROOT / "results" / ds / path / "grn.tsv.gz"
            if f.exists():
                ours[name] = L.read_ours(f)
        for gt in L.ds_gts(D):
            rows, Y, Z, M = L.eval_block(D, gt)
            if Y.sum() == 0:
                continue
            rec = {"dataset": ds, "gt": gt, "n_tfs": len(rows), "pos_rate": float(Y[M].mean()),
                   # leave-one-TF-out marginal oracle (README D2): knows only the marginals, no self-inclusion
                   "marginal_oracle_loo": L.full_auroc(L.loo_marginal_oracle(Y, M), Y, M),
                   "row_only": L.full_auroc(np.repeat(Y.mean(1, keepdims=True), Y.shape[1], 1), Y, M),
                   "col_only_loo": L.full_auroc(Y.sum(0, keepdims=True) - Y.astype(float), Y, M),
                   # in-sample beta-model: self-inclusion upper bound only
                   "beta_model_insample": L.full_auroc(L.fit_beta_model(Y, M), Y, M)}
            for name, e in list(ours.items()) + [(m, e) for m, e in rel.items() if e is not None]:
                rec[name] = L.full_auroc(L.edges_to_matrix(D, rows, e), Y, M)
            relv = [rec[m] for m in L.RELEASED if m in rec]
            rec["best_released"] = max(relv) if relv else np.nan
            rec["best_released_name"] = [m for m in L.RELEASED if m in rec and rec[m] == rec["best_released"]][0] if relv else ""
            rows_out.append(rec)
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in rec.items()}, flush=True)
    df = pd.DataFrame(rows_out)
    L.OUT.mkdir(exist_ok=True, parents=True)
    df.to_csv(L.OUT / "a1_gt_structure.csv", index=False)


if __name__ == "__main__":
    main()
