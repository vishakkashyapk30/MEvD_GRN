#!/usr/bin/env python
"""Parameter ladder of simple, mathematically grounded link models on the scMultiomeGRN benchmark (2026-11-19,
advisor item 4). Same folds, same evaluation (official test, degree-matched, caliper-matched; held-out-TF classes)
and same result format as scripts/scmgrn/22. Everything is fitted on the fold's training pairs only.

Network models (no cell data), N = number of TFs:
  ladder_beta          beta-model / degree-corrected logistic model: logit p_uv = b + beta_u + beta_v
                       (N + 1 parameters; Chatterjee, Diaconis & Sly, Ann. Appl. Probab. 2011), L2-regularised
  ladder_lsm{k}        latent-space distance model (Hoff, Raftery & Handcock, JASA 2002):
                       logit p_uv = a - ||z_u - z_v||^2, z in R^k, k = 2, 4, 8  (k N + 1 parameters)
  ladder_dclsm{k}      degree-corrected latent-space model: logit p_uv = a + beta_u + beta_v - ||z_u - z_v||^2
                       ((k + 1) N + 1 parameters; Krivitsky et al., Social Networks 2009)
Gaussian models on cell-data pair features (no training-graph input): MAESTRO and GRNBoost2 per-TF vectors reduced
to 8 dimensions each by label-free SVD, pair feature [x_u * x_v, |x_u - x_v|] for both (32 numbers):
  ladder_lda           linear discriminant analysis (Gaussian classes, shared covariance)
  ladder_qda           quadratic discriminant analysis (Gaussian classes, own covariance; reg 0.1)
  ladder_gnb           Gaussian naive Bayes
  ladder_gmm           Gaussian-mixture classifier (4 full-covariance components per class; log-likelihood ratio)
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from src.benchmarks import scmgrn_protocol as proto            # noqa: E402
from src.benchmarks.scmgrn_features import svd_reduce          # noqa: E402
from src.utils.io import load_config, save_json                # noqa: E402

_spec = importlib.util.spec_from_file_location("s22", REPO / "scripts/scmgrn/22_scmgrn_simple_baselines.py")
s22 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s22)
VERSION = "ladder-2026-11-19"


def fit_beta(n, pairs, y, l2=1e-2):
    from scipy import sparse
    from sklearn.linear_model import LogisticRegression
    r = np.repeat(np.arange(len(pairs)), 2)
    X = sparse.csr_matrix((np.ones(2 * len(pairs)), (r, pairs.reshape(-1))), shape=(len(pairs), n))
    m = LogisticRegression(C=1.0 / l2, max_iter=3000).fit(X, y)
    b = m.coef_.ravel()
    S = b[:, None] + b[None, :] + m.intercept_[0]
    return 1 / (1 + np.exp(-S))


def fit_lsm(n, pairs, y, k, degree_corrected, seed, epochs=400, l2=1e-3):
    g = torch.Generator().manual_seed(seed)
    z = torch.nn.Parameter(0.1 * torch.randn(n, k, generator=g))
    a = torch.nn.Parameter(torch.zeros(1))
    beta = torch.nn.Parameter(torch.zeros(n)) if degree_corrected else None
    params = [z, a] + ([beta] if beta is not None else [])
    opt = torch.optim.Adam(params, lr=0.05)
    P = torch.as_tensor(pairs, dtype=torch.long)
    Y = torch.as_tensor(y, dtype=torch.float32)
    pw = torch.tensor((len(y) - y.sum()) / max(y.sum(), 1), dtype=torch.float32)   # class balance
    def logit(i, j):
        out = a - ((z[i] - z[j]) ** 2).sum(-1)
        return out + beta[i] + beta[j] if beta is not None else out
    for _ in range(epochs):
        opt.zero_grad()
        loss = torch.nn.functional.binary_cross_entropy_with_logits(logit(P[:, 0], P[:, 1]), Y, pos_weight=pw)
        loss = loss + l2 * sum((p ** 2).sum() for p in params)
        loss.backward()
        opt.step()
    with torch.no_grad():
        iu, ju = np.triu_indices(n, 1)
        S = np.zeros((n, n))
        S[iu, ju] = torch.sigmoid(logit(torch.as_tensor(iu), torch.as_tensor(ju))).numpy()
    return S + S.T


def fit_gauss(name, Xtr, y, Xall, seed):
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
    from sklearn.mixture import GaussianMixture
    from sklearn.naive_bayes import GaussianNB
    if name == "gmm":
        g1 = GaussianMixture(4, covariance_type="full", reg_covar=1e-3, random_state=seed).fit(Xtr[y == 1])
        g0 = GaussianMixture(4, covariance_type="full", reg_covar=1e-3, random_state=seed).fit(Xtr[y == 0])
        return g1.score_samples(Xall) - g0.score_samples(Xall)
    m = {"lda": LinearDiscriminantAnalysis(), "qda": QuadraticDiscriminantAnalysis(reg_param=0.1),
         "gnb": GaussianNB()}[name].fit(Xtr, y)
    return m.predict_proba(Xall)[:, 1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--ds", required=True, choices=["lung", "hema"])
    ap.add_argument("--ct", default=None)
    ap.add_argument("--split_kind", default="official", choices=proto.SPLIT_KINDS)
    ap.add_argument("--folds", default=None)
    ap.add_argument("--models", default="beta,lsm2,lsm4,lsm8,dclsm2,dclsm4,dclsm8,lda,qda,gnb,gmm")
    ap.add_argument("--skip_existing", action="store_true")
    a = ap.parse_args()
    root = Path(a.root).expanduser()
    cfg = load_config(str(REPO / f"configs/scmgrn/{a.ds}.yaml"))
    models = a.models.split(",")
    for c in ([a.ct] if a.ct else list(cfg["cell_types"])):
        ct = s22.CellType(root, cfg, c, a.split_kind)
        n = ct.n
        rdir = root / proto.results_subdir(cfg["results_dir"], a.split_kind) / c
        M8 = svd_reduce(ct.get("maestro"), 8)
        G8 = svd_reduce(ct.get("grnboost"), 8)
        iu, ju = np.triu_indices(n, 1)
        def phi(i, j):
            return np.concatenate([M8[i] * M8[j], np.abs(M8[i] - M8[j]), G8[i] * G8[j], np.abs(G8[i] - G8[j])], 1)
        Xall = phi(iu, ju)
        for fold in ([int(x) for x in a.folds.split(",")] if a.folds else range(1, ct.n_folds + 1)):
            f = proto.load_fold(ct.sdir, fold, n)
            tr_pos, tr_neg = s22.unordered(f["train"]["pos"]), s22.unordered(f["train"]["neg"])
            pairs = np.concatenate([tr_pos, tr_neg])
            y = np.r_[np.ones(len(tr_pos), int), np.zeros(len(tr_neg), int)]
            held = {tuple(p) for p in s22.unordered(f["test"]["pos"])} | {tuple(p) for p in s22.unordered(f["test"]["neg"])}
            assert not (held & {tuple(p) for p in pairs}), "held-out pair in training"
            mu, sd = Xall.mean(0), Xall.std(0) + 1e-8
            Xtr = (phi(pairs[:, 0], pairs[:, 1]) - mu) / sd
            Xa = (Xall - mu) / sd
            for m in models:
                name = f"ladder_{m}"
                out = rdir / name / "seed0" / f"fold{fold:02d}.json"
                if a.skip_existing and out.exists():
                    continue
                t0 = time.time()
                if m == "beta":
                    S = fit_beta(n, pairs, y)
                elif m.startswith("dclsm"):
                    S = fit_lsm(n, pairs, y, int(m[5:]), True, fold)
                elif m.startswith("lsm"):
                    S = fit_lsm(n, pairs, y, int(m[3:]), False, fold)
                else:
                    s = fit_gauss(m, Xtr, y, Xa, fold)
                    S = np.zeros((n, n))
                    S[iu, ju] = s
                    S = S + S.T
                met = proto.evaluate_all(S, ct.adj, f, ct.sdir, fold)
                out.parent.mkdir(parents=True, exist_ok=True)
                save_json({"method": name, "family": "ladder", "version": VERSION, "cell_type": c,
                           "dataset": cfg["dataset"], "fold": fold, "seed": 0,
                           "note": "fit on the fold's training part only", "metrics": met}, out)
                keys = [k for k in ("test", "test_degmatched", "test_degcaliper", "test_c2", "test_c3") if k in met]
                print(f"{c} fold {fold} {name}: " + " ".join(f"{k} {met[k]['auroc']:.3f}" for k in keys)
                      + f" ({time.time() - t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
