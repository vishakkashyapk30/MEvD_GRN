#!/usr/bin/env python
"""MID-GRN v2 (post-hoc revision of the frozen v1, README "v2"). Same TF-disjoint folds and leak rules.

Changes from v1, each from MDH or from v1's instability, none chosen on test folds:
  1. Target marginal from more label data: in-degree among the training TFs (as v1) PLUS in-degree in the
     SC-MO-GRN-DB ChIP compendium (RN005 human / RN006 mouse; the L2 label source of bear_grn.md s6.2) with EVERY
     source TF that appears in any GT or label file of the dataset removed, plus label-free target features
     (ATAC activity, openness, RNA mean / variance / detection, Geneformer token present).
  2. Out-of-fold combination: 3 TF-disjoint inner folds over the training TFs give out-of-fold marginal features,
     alpha and g for EVERY training TF; the combiner is fit on all of them (v1: ~3-25 validation TFs).
  3. Combiner: LightGBM (fixed hyperparameters below), monotone increasing in both in-degrees and in g.
  4. The interaction head g is v1's (same pair features, matched-set conditional logit, same MLP); held-out TFs
     get the mean of the 3 inner models.
Two outputs per run, both written like scripts/22:
  mid2_marginal  combiner on marginal features only (TF alpha + target features): the strongest marginal baseline
  mid2_full      the same + g
  (and mid2_full_g/score_matrix.npz = the interaction head alone)

  python 06_midgrn_v2.py --dataset K562 --seeds 42 43 44
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import time
from pathlib import Path

import lightgbm as lgb
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.linear_model import Ridge

import mdh_lib as L
from src.benchmarks import bear_protocol as bp

_spec = importlib.util.spec_from_file_location("v1", Path(__file__).parent / "03_midgrn.py")
V1 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V1)

LGB = dict(objective="binary", learning_rate=0.05, num_leaves=31, min_child_samples=200, subsample=0.8,
           subsample_freq=1, colsample_bytree=1.0, n_estimators=300, verbose=-1, n_jobs=12)
N_INNER = 3
# TF marginal alpha for an UNSEEN TF is not identifiable from its node features: v1's calibration weight on alpha
# flipped sign across K562 folds (-0.40 .. +0.60), and with it the combiner ranks whole TF rows by a noisy guess.
# v2 therefore scores held-out TFs with target marginals + the pair head only (decided on K562, the dev set).
USE_ALPHA = False


def compendium_indegree(D: dict) -> np.ndarray:
    """Per measured gene: number of distinct compendium source TFs targeting it, excluding every source that is a
    source in ANY GT or label file of this dataset (all evaluated / trainable TFs)."""
    cfg = L.DCFG["datasets"][D["ds"]]
    excl = set()
    for p in list(cfg.get("gts", {}).values()) + list(cfg.get("labels", {}).values()):
        f = L.ROOT / p
        if f.exists():
            excl |= set(L.bm.read_table(str(f))["Source"].astype(str).str.upper())
    comp = L.DCFG["compendium"][cfg["species"]]
    c = L.bm.read_table(str(L.ROOT / comp))
    # the Target field holds "SYMBOL  Relationship" (space-separated): keep the first token
    c = pd.DataFrame({"s": c["Source"].astype(str).str.strip().str.upper(),
                      "t": c["Target"].astype(str).str.split().str[0].str.upper()}).drop_duplicates()
    c = c[~c.s.isin(excl)]
    cnt = c.groupby("t").s.nunique()
    out = np.zeros(D["n"], dtype=np.float64)
    idx = cnt.index.map(D["upper"])
    ok = idx.notna()
    out[idx[ok].astype(int)] = cnt.to_numpy()[ok]
    return out, len(excl), int(c.s.nunique())


def train_g(X: np.ndarray, Ytr: np.ndarray, rows_fit: np.ndarray, selfmask, seed: int, dev: str):
    """v1 interaction head: matched-set conditional logit over (same TF, same in-degree stratum of rows_fit)."""
    m, n, f = X.shape
    indeg = Ytr[rows_fit].sum(0).astype(np.float64)
    in_fit = np.zeros(m, bool); in_fit[rows_fit] = True
    loo = indeg[None, :] - Ytr
    q = V1.bins_from(indeg)
    Zs = V1.stratum(np.where(in_fit[:, None], loo, indeg[None, :]), q)
    mu = X[rows_fit].reshape(-1, f).mean(0); sd = X[rows_fit].reshape(-1, f).std(0) + 1e-6
    Xt = torch.tensor(((X - mu) / sd).reshape(m * n, f), dtype=torch.float32, device=dev)
    fr, fc = np.nonzero(selfmask[rows_fit]); fr = rows_fit[fr]
    fy = Ytr[fr, fc]
    nb = int(Zs.max()) + 1
    gid = (np.searchsorted(rows_fit, fr) * nb + Zs[fr, fc]).astype(np.int64)
    neg = np.where(~fy)[0]
    o = neg[np.argsort(gid[neg], kind="stable")]
    cnt = np.bincount(gid[neg], minlength=len(rows_fit) * nb)
    start = np.r_[0, np.cumsum(cnt)[:-1]]
    P = np.where(fy)[0]; P = P[cnt[gid[P]] > 0]
    torch.manual_seed(seed)
    model = V1.MLP(f).to(dev)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    pf = torch.tensor(fr[P] * n + fc[P], device=dev)
    nf = torch.tensor(fr[o] * n + fc[o], device=dev)
    st_t = torch.tensor(start[gid[P]], device=dev); ct_t = torch.tensor(cnt[gid[P]], device=dev)
    for _ in range(4):                    # v1 early-stopped at 3-5 epochs; fixed 4 here (no validation TFs)
        model.train()
        perm = torch.randperm(len(P), device=dev)
        for b in range(0, len(P), 4096):
            ib = perm[b:b + 4096]
            u = torch.rand(len(ib), V1.N_NEG, device=dev)
            nidx = st_t[ib, None] + (u * ct_t[ib, None]).long().clamp(max=ct_t[ib, None] - 1)
            logits = torch.cat([model(Xt[pf[ib]])[:, None], model(Xt[nf[nidx]])], 1)
            loss = nn.functional.cross_entropy(logits, torch.zeros(len(ib), dtype=torch.long, device=dev))
            opt.zero_grad(); loss.backward(); opt.step()
    return model, mu, sd


def g_predict(model, mu, sd, X, rows, dev):
    m, n, f = X.shape
    Xt = torch.tensor(((X[rows] - mu) / sd).reshape(-1, f), dtype=torch.float32, device=dev)
    return V1.predict(model, Xt).reshape(len(rows), n)


def pair_X(D, st, Ytr_src: np.ndarray, src_rows: np.ndarray):
    """v1 pair features with look-alike neighbours restricted to src_rows (whose labels are Ytr_src)."""
    m = len(D["tf"])
    in_S = np.zeros(m, bool); in_S[src_rows] = True
    Ys = Ytr_src & in_S[:, None]
    indeg = Ys.sum(0).astype(np.float64)
    mats = [st["pair"][c] for c in V1.DATA_COLS]
    for c in V1.LA_COLS:
        sim, valid = st["sims"][c]
        mats.append(V1.look_alike(sim, valid, Ys, in_S, indeg))
    return np.stack(mats, -1)


def run(D, st, tfeat, seed, dev, log):
    m, n = len(D["tf"]), D["n"]
    pos = D["tfpos"]
    Yall = L.label_matrix(D, D["lab"])
    selfmask = np.ones((m, n), bool); selfmask[np.arange(m), D["tf"]] = False
    outdeg = {int(k): int(v) for k, v in zip(*np.unique(D["lab"][0], return_counts=True))}
    out = {"mid2_marginal": np.full((m, n), np.nan), "mid2_full": np.full((m, n), np.nan)}
    gout = np.full((m, n), np.nan, dtype=np.float32)
    foldv = np.full(m, -1)
    for fold in L.folds(D, seed):
        t0 = time.time()
        held = np.array([pos[int(t)] for t in fold["held"]])
        T = np.array([pos[int(t)] for t in fold["train"]])
        inT = np.zeros(m, bool); inT[T] = True
        Ytr = Yall & inT[:, None]
        assert not Ytr[held].any()
        inner = bp.tf_disjoint_folds(fold["train"], outdeg, N_INNER, seed + 2000 + fold["k"])

        def marg_feats(rows, src):           # in-degree over src TFs (leave-own-row-out), compendium, target feats
            k = len(src)
            ind = Ytr[src].sum(0).astype(np.float64)[None, :] - (Ytr[rows] * np.isin(rows, src)[:, None])
            kk = k - np.isin(rows, src).astype(np.float64)[:, None]
            lr = np.log((ind + 0.5) / (kk + 1.0))
            return [lr] + [np.broadcast_to(tfeat[:, c][None, :], lr.shape) for c in range(tfeat.shape[1])]

        def alpha_fit(src):
            nx = st["node"]; mu_, sd_ = nx[src].mean(0), nx[src].std(0) + 1e-6
            rate = (Ytr[src].sum(1) + 0.5) / (n + 1.0)
            r = Ridge(alpha=10.0).fit((nx[src] - mu_) / sd_, np.log(rate / (1 - rate)))
            return r.predict((nx - mu_) / sd_)

        # ---- out-of-fold features for every training TF
        oof_rows, oof_feats, oof_g = [], [], []
        models = []
        for v, Jg in enumerate(inner):
            J = np.array([pos[int(t)] for t in Jg])
            I = np.setdiff1d(T, J)
            YI = Ytr & np.isin(np.arange(m), I)[:, None]
            X = pair_X(D, st, YI, I)
            mdl, mu, sd = train_g(X, YI, I, selfmask, seed * 1000 + fold["k"] * 10 + v, dev)
            models.append((mdl, mu, sd))
            oof_g.append(g_predict(mdl, mu, sd, X, J, dev))
            a = alpha_fit(I)
            oof_feats.append(marg_feats(J, I) + ([np.repeat(a[J][:, None], n, 1)] if USE_ALPHA else []))
            oof_rows.append(J)
        J_all = np.concatenate(oof_rows)
        F_tr = [np.concatenate([f[c] for f in oof_feats]) for c in range(len(oof_feats[0]))]
        G_tr = np.concatenate(oof_g)
        Ytr_rows = Ytr[J_all]
        msk = selfmask[J_all]
        Xm = np.stack([f[msk] for f in F_tr], 1)
        y = Ytr_rows[msk]
        nm_ = Xm.shape[1]
        mono_m = [1, 1] + [0] * (nm_ - 2)
        comb_m = lgb.LGBMClassifier(**LGB, monotone_constraints=mono_m, random_state=seed).fit(Xm, y)
        Xf = np.c_[Xm, G_tr[msk]]
        comb_f = lgb.LGBMClassifier(**LGB, monotone_constraints=mono_m + [1], random_state=seed).fit(Xf, y)
        # ---- held-out TFs: features from all training TFs; g = mean of the inner models
        Xh = pair_X(D, st, Ytr, T)
        g_h = np.mean([g_predict(md, mu, sd, Xh, held, dev) for md, mu, sd in models], 0)
        a = alpha_fit(T)
        Fh = marg_feats(held, T) + ([np.repeat(a[held][:, None], n, 1)] if USE_ALPHA else [])
        hm = selfmask[held]
        Xhm = np.stack([f[hm] for f in Fh], 1)
        hr, hc = np.nonzero(hm)
        out["mid2_marginal"][held[hr], hc] = comb_m.predict_proba(Xhm)[:, 1]
        out["mid2_full"][held[hr], hc] = comb_f.predict_proba(np.c_[Xhm, g_h[hm]])[:, 1]
        gout[held] = g_h
        foldv[held] = fold["k"]
        log({"fold": fold["k"], "n_held": len(held), "n_train": len(T), "inner_sizes": [len(x) for x in inner],
             "n_combiner_rows": int(len(y)), "pos_rate": float(y.mean()),
             "gain_importance_full": dict(zip(["indeg", "compendium"] + [f"t{c}" for c in range(nm_ - 2 - USE_ALPHA)] + (["alpha"] if USE_ALPHA else []) + ["g"],
                                              comb_f.booster_.feature_importance("gain").round(0).tolist())),
             "seconds": round(time.time() - t0, 1)})
    assert (foldv >= 0).all()
    return out, gout, foldv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--seeds", nargs="+", type=int, default=[42, 43, 44])
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()
    D = L.load_ds(args.dataset)
    st = V1.static_features(D)
    comp, n_excl, n_src = compendium_indegree(D)
    has_tok = (np.abs(D["fm"]).sum(1) > 0).astype(np.float64)
    tfeat = np.c_[np.log1p(comp), D["atac"], D["rna"], has_tok]       # target features (label-free + compendium)
    print(f"[{args.dataset}] compendium: {n_src} source TFs kept, {n_excl} dataset TFs excluded; "
          f"genes with compendium in-degree > 0: {(comp > 0).mean():.3f}", flush=True)
    m, n = len(D["tf"]), D["n"]
    for seed in args.seeds:
        assert L.check_folds(D, seed)
        outs = [L.ROOT / "results" / args.dataset / v / f"seed{seed}" for v in ("mid2_marginal", "mid2_full")]
        if all((o / "grn.tsv.gz").exists() for o in outs):
            print("[skip]", seed, flush=True)
            continue
        t0 = time.time()
        logs = []
        res, g, foldv = run(D, st, tfeat, seed, args.device,
                            lambda r: (logs.append(r), print(f"[{args.dataset} s{seed}] " + json.dumps(r), flush=True)))
        src, dst = np.nonzero(np.ones((m, n), bool))
        keep = D["tf"][src] != dst
        src, dst = src[keep], dst[keep]
        nmz = np.asarray(D["names"], dtype=object)
        for v, o in zip(("mid2_marginal", "mid2_full"), outs):
            o.mkdir(parents=True, exist_ok=True)
            pd.DataFrame({"Source": nmz[D["tf"][src]], "Target": nmz[dst], "Score": res[v][src, dst],
                          "fold": foldv[src].astype(np.int16)}).to_csv(o / "grn.tsv.gz", sep="\t", index=False,
                                                                       float_format="%.7g")
            json.dump({"dataset": args.dataset, "variant": v, "seed": seed, "fold_logs": logs,
                       "compendium_sources": n_src, "excluded_dataset_tfs": n_excl, "seconds": time.time() - t0},
                      open(o / "meta.json", "w"), indent=1)
        go = L.ROOT / "results" / args.dataset / "mid2_full_g" / f"seed{seed}"
        go.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(go / "score_matrix.npz", S=g, tf=D["tf"])
        print(f"[done] {args.dataset} seed {seed} in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
