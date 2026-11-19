#!/usr/bin/env python
"""MID-GRN (Marginal-Interaction Decomposition), design frozen in README.md, on BEAR-GRN's TF-disjoint folds.

  python 03_midgrn.py --dataset K562 --seeds 42 --variants mid_full mid_dataonly mid_bce

Per (dataset, seed, fold): marginal head (training-TF in-degree z_j, ridge TF marginal alpha_i), interaction head
g_ij = MLP(pair features) trained with a matched-set conditional logit (same TF, same target in-degree stratum),
then s_ij = sigmoid(w . [z, z^2, alpha, g]) fitted on the inner-validation TFs. Held-out TFs never enter any
training input, feature or statistic (asserted). Outputs, under $BEAR_ROOT/results/<ds>/:
  <variant>/seed<s>/grn.tsv.gz        full score (Source, Target, Score, fold), like scripts/22
  <variant>_g/seed<s>/score_matrix.npz interaction-only score g (TF nodes x genes), for the matched AUROC
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score

import mdh_lib as L

DATA_COLS = ["motif_best", "motif_hits", "coexpr", "coexpr_resid", "fm_cos"]
LA_COLS = ["la_fm", "la_coexpr", "la_motif"]
N_BG = 5            # leading co-expression programs removed for coexpr_resid
K_NB, TAU = 10, 10.0
N_NEG = 8
NBINS = 20


def unit(x):
    return x / (np.linalg.norm(x, axis=1, keepdims=True) + 1e-8)


def static_features(D: dict) -> dict:
    """Fold-independent pair features (m TF nodes x n genes) and TF-TF similarities. Label-free."""
    tf = D["tf"]
    S, fm = D["sig"], D["fm"]
    mot = D["motif"]
    has = mot[:, :, 2].max(1) > 0
    fmn = unit(fm)
    sig_n = unit(S)
    m0 = mot[:, :, 0]
    mc = m0 - m0.mean(1, keepdims=True)
    mc = mc / (np.linalg.norm(mc, axis=1, keepdims=True) + 1e-8)
    feats = {"motif_best": mot[:, :, 0], "motif_hits": mot[:, :, 1],
             "coexpr": S[tf] @ S.T, "coexpr_resid": S[tf, N_BG:] @ S[:, N_BG:].T,
             "fm_cos": fmn[tf] @ fmn.T}
    sims = {"la_fm": (fmn[tf] @ fmn[tf].T, np.ones(len(tf), bool)),
            "la_coexpr": (sig_n[tf] @ sig_n[tf].T, np.ones(len(tf), bool)),
            "la_motif": (mc @ mc.T, has)}
    pca = PCA(n_components=8, random_state=0).fit(fm)
    node = np.c_[D["rna"][tf], D["atac"][tf], pca.transform(fm[tf]), has.astype(np.float32),
                 mot[:, :, 1].mean(1), np.linalg.norm(S[tf], axis=1)]
    return {"pair": {k: v.astype(np.float32) for k, v in feats.items()}, "sims": sims, "node": node}


def look_alike(sim: np.ndarray, valid: np.ndarray, Ytr: np.ndarray, in_T: np.ndarray, indeg: np.ndarray):
    """LA_ij = sum_k w_ik Y_kj - (indeg_j - Y_ij[i in T]) / (|T| - [i in T]), k in top-K of T \\ {i}."""
    m = sim.shape[0]
    s = sim.astype(np.float64).copy()
    s[:, ~(in_T & valid)] = -np.inf
    np.fill_diagonal(s, -np.inf)
    top = np.argsort(-s, axis=1)[:, :K_NB]
    v = np.take_along_axis(s, top, 1)
    ok = np.isfinite(v)
    e = np.where(ok, np.exp(TAU * (np.where(ok, v, 0) - np.where(ok, v, -1e9).max(1, keepdims=True))), 0.0)
    w = e / np.maximum(e.sum(1, keepdims=True), 1e-12)
    W = np.zeros((m, m))
    np.put_along_axis(W, top, w, 1)
    assert not W[:, ~in_T].any(), "look-alike weight on a held-out TF"
    LA = W @ Ytr.astype(np.float64)
    cnt = in_T.sum() - in_T.astype(np.int64)
    mean = (indeg[None, :] - Ytr * in_T[:, None]) / np.maximum(cnt, 1)[:, None]
    LA = LA - mean
    has_nb = ok.any(1) & valid
    LA[~has_nb] = 0.0
    return LA.astype(np.float32)


def bins_from(indeg: np.ndarray):
    nz = indeg[indeg > 0]
    return np.unique(np.quantile(nz, np.linspace(0, 1, NBINS + 1)[1:-1])) if nz.size else np.zeros(0)


def stratum(d: np.ndarray, q: np.ndarray):
    return np.where(d > 0, 1 + np.searchsorted(q, d, side="right"), 0)


class MLP(nn.Module):
    def __init__(self, f: int):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(f, 64), nn.ReLU(), nn.Linear(64, 64), nn.ReLU(), nn.Linear(64, 1))

    def forward(self, x):
        return self.net(x).squeeze(-1)


def matched_val(g_rows: np.ndarray, Y: np.ndarray, Z: np.ndarray, mask: np.ndarray) -> float:
    U, Dn = L.matched_stats(g_rows, Y, Z, mask)
    return float(U.sum() / max(Dn.sum(), 1))


@torch.no_grad()
def predict(model, X: torch.Tensor, chunk=1 << 20) -> np.ndarray:
    model.eval()
    return torch.cat([model(X[i:i + chunk]) for i in range(0, X.shape[0], chunk)]).float().cpu().numpy()


def run_fold(D, st, fold, variant, seed, dev, log):
    tf, n, m = D["tf"], D["n"], len(D["tf"])
    pos = D["tfpos"]
    held = np.array([pos[int(t)] for t in fold["held"]])
    T = np.array([pos[int(t)] for t in fold["train"]])
    V = np.array([pos[int(t)] for t in fold["val"]])
    Fi = np.array([pos[int(t)] for t in fold["fit"]])
    in_T = np.zeros(m, bool); in_T[T] = True
    Yall = L.label_matrix(D, D["lab"])
    Ytr = Yall & in_T[:, None]                       # labels of training TFs only
    assert not Ytr[held].any()
    indeg = Ytr.sum(0).astype(np.float64)            # = the in-degree baseline of this fold
    loo = indeg[None, :] - Ytr                       # leave-own-TF-out for training rows
    lmu, lsd = np.log1p(indeg).mean(), np.log1p(indeg).std() + 1e-8
    z = (np.log1p(np.where(in_T[:, None], loo, indeg[None, :])) - lmu) / lsd   # (m, n)
    q = bins_from(indeg)
    Zs = stratum(np.where(in_T[:, None], loo, indeg[None, :]), q)
    selfmask = np.ones((m, n), bool); selfmask[np.arange(m), tf] = False

    # TF marginal alpha: ridge on TF node features, fit TFs only
    nodeX = st["node"]
    mu, sd = nodeX[Fi].mean(0), nodeX[Fi].std(0) + 1e-6
    rate = (Ytr[Fi].sum(1) + 0.5) / (n + 1.0)
    ridge = Ridge(alpha=10.0).fit((nodeX[Fi] - mu) / sd, np.log(rate / (1 - rate)))
    alpha = ridge.predict((nodeX - mu) / sd)

    cols = DATA_COLS + (LA_COLS if variant in ("mid_full", "mid_bce") else [])
    mats = [st["pair"][c] for c in DATA_COLS]
    if variant in ("mid_full", "mid_bce"):
        for c in LA_COLS:
            sim, valid = st["sims"][c]
            mats.append(look_alike(sim, valid, Ytr, in_T, indeg))
    if variant == "mid_bce":
        mats += [z.astype(np.float32), (z ** 2).astype(np.float32),
                 np.repeat(alpha[:, None], n, 1).astype(np.float32)]
        cols = cols + ["z", "z2", "alpha"]
    X = np.stack(mats, -1)                            # (m, n, f)
    fmu = X[Fi].reshape(-1, X.shape[-1]).mean(0)
    fsd = X[Fi].reshape(-1, X.shape[-1]).std(0) + 1e-6
    X = (X - fmu) / fsd
    Xt = torch.tensor(X.reshape(m * n, -1), dtype=torch.float32, device=dev)
    flat = lambda r, c: r * n + c                    # noqa: E731

    torch.manual_seed(seed * 100 + fold["k"])
    rng = np.random.default_rng(seed * 100 + fold["k"])
    model = MLP(X.shape[-1]).to(dev)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    fr, fc = np.nonzero(selfmask[Fi])
    fr = Fi[fr]
    fy = Ytr[fr, fc]
    vr, vc = np.nonzero(selfmask[V]); vr = V[vr]
    Yv, Zv, Mv = Ytr[V], Zs[V], selfmask[V]
    best, best_state, bad, hist = -1.0, None, 0, []

    if variant in ("mid_full", "mid_dataonly"):
        nb = int(Zs.max()) + 1
        gid = (np.searchsorted(Fi, fr) * nb + Zs[fr, fc]).astype(np.int64)
        neg = np.where(~fy)[0]
        o = neg[np.argsort(gid[neg], kind="stable")]
        G = len(Fi) * nb
        cnt = np.bincount(gid[neg], minlength=G)
        start = np.r_[0, np.cumsum(cnt)[:-1]]
        P = np.where(fy)[0]
        P = P[cnt[gid[P]] > 0]
        pos_flat = torch.tensor(flat(fr[P], fc[P]), device=dev)
        neg_flat_sorted = torch.tensor(flat(fr[o], fc[o]), device=dev)
        st_t = torch.tensor(start[gid[P]], device=dev)
        ct_t = torch.tensor(cnt[gid[P]], device=dev)
        for ep in range(8):
            model.train()
            perm = torch.randperm(len(P), device=dev)
            tot = 0.0
            for b in range(0, len(P), 4096):
                ib = perm[b:b + 4096]
                u = torch.rand(len(ib), N_NEG, device=dev)
                nidx = st_t[ib, None] + (u * ct_t[ib, None]).long().clamp(max=ct_t[ib, None] - 1)
                xp = Xt[pos_flat[ib]]
                xn = Xt[neg_flat_sorted[nidx]]
                logits = torch.cat([model(xp)[:, None], model(xn)], 1)
                loss = nn.functional.cross_entropy(logits, torch.zeros(len(ib), dtype=torch.long, device=dev))
                opt.zero_grad(); loss.backward(); opt.step()
                tot += float(loss) * len(ib)
            gv = predict(model, Xt[torch.tensor(flat(V[:, None], np.arange(n)[None, :]).ravel(), device=dev)]).reshape(len(V), n)
            mv = matched_val(gv, Yv, Zv, Mv)
            hist.append({"epoch": ep, "loss": tot / len(P), "val_matched": mv})
            if mv > best + 1e-4:
                best, bad = mv, 0
                best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
            else:
                bad += 1
                if bad >= 2:
                    break
    else:   # mid_bce: the usual recipe
        P = np.where(fy)[0]
        Nn = np.where(~fy)[0]
        yv_flat = Ytr[vr, vc]
        for ep in range(8):
            model.train()
            ns = rng.choice(Nn, size=min(len(Nn), 4 * len(P)), replace=False)
            idx = np.r_[P, ns]
            lab = torch.tensor(np.r_[np.ones(len(P)), np.zeros(len(ns))], dtype=torch.float32, device=dev)
            fl = torch.tensor(flat(fr[idx], fc[idx]), device=dev)
            perm = torch.randperm(len(idx), device=dev)
            tot = 0.0
            for b in range(0, len(idx), 4096 * 5):
                ib = perm[b:b + 4096 * 5]
                loss = nn.functional.binary_cross_entropy_with_logits(model(Xt[fl[ib]]), lab[ib])
                opt.zero_grad(); loss.backward(); opt.step()
                tot += float(loss) * len(ib)
            gvf = predict(model, Xt[torch.tensor(flat(vr, vc), device=dev)])
            ap = float(average_precision_score(yv_flat, gvf))
            gv = np.zeros((len(V), n)); gv[np.searchsorted(V, vr), vc] = gvf
            hist.append({"epoch": ep, "loss": tot / len(idx), "val_ap": ap,
                         "val_matched": matched_val(gv, Yv, Zv, Mv)})
            if ap > best + 1e-4:
                best, bad = ap, 0
                best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
            else:
                bad += 1
                if bad >= 2:
                    break
    model.load_state_dict(best_state)
    g_all = predict(model, Xt).reshape(m, n)

    if variant == "mid_bce":
        final = g_all
        w = None
    else:   # calibration on the inner-validation TFs, natural distribution
        A = np.c_[z[vr, vc], z[vr, vc] ** 2, alpha[vr], g_all[vr, vc]]
        lr = LogisticRegression(C=1.0, max_iter=1000).fit(A, Ytr[vr, vc])
        w = lr.coef_.ravel().tolist() + [float(lr.intercept_[0])]
        hr, hc = np.nonzero(selfmask[held]); hr = held[hr]
        final = np.zeros((m, n), dtype=np.float64)
        final[hr, hc] = lr.decision_function(np.c_[z[hr, hc], z[hr, hc] ** 2, alpha[hr], g_all[hr, hc]])
    rec = {"fold": fold["k"], "n_held": len(held), "n_fit": len(Fi), "n_val": len(V), "features": cols,
           "best_val": best, "history": hist, "calib_w[z,z2,alpha,g,b]": w}
    log(rec)
    return held, final, g_all


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--seeds", nargs="+", type=int, default=[42])
    ap.add_argument("--variants", nargs="+", default=["mid_full", "mid_dataonly", "mid_bce"])
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()
    D = L.load_ds(args.dataset)
    st = static_features(D)
    m, n = len(D["tf"]), D["n"]
    for seed in args.seeds:
        assert L.check_folds(D, seed), "stored in-degree folds missing"
        for variant in args.variants:
            out = L.ROOT / "results" / args.dataset / variant / f"seed{seed}"
            if (out / "grn.tsv.gz").exists():
                print(f"[skip] {out}", flush=True)
                continue
            t0 = time.time()
            logs = []
            final = np.full((m, n), np.nan)
            gmat = np.full((m, n), np.nan, dtype=np.float32)
            foldv = np.full(m, -1)
            for fold in L.folds(D, seed):
                held, fs, g = run_fold(D, st, fold, variant, seed, args.device,
                                       lambda r: (logs.append(r), print(f"[{args.dataset} s{seed} {variant}] "
                                                  + json.dumps({k: v for k, v in r.items() if k != 'history'}), flush=True)))
                final[held] = fs[held]
                gmat[held] = g[held]
                foldv[held] = fold["k"]
            assert (foldv >= 0).all()
            src, dst = np.nonzero(np.ones((m, n), bool))
            keep = D["tf"][src] != dst
            src, dst = src[keep], dst[keep]
            sc = final[src, dst]
            sc = 1.0 / (1.0 + np.exp(-sc))
            nm = np.asarray(D["names"], dtype=object)
            out.mkdir(parents=True, exist_ok=True)
            pd.DataFrame({"Source": nm[D["tf"][src]], "Target": nm[dst], "Score": sc,
                          "fold": foldv[src].astype(np.int16)}).to_csv(out / "grn.tsv.gz", sep="\t", index=False,
                                                                       float_format="%.7g")
            json.dump({"dataset": args.dataset, "variant": variant, "seed": seed, "folds_logs": logs,
                       "seconds": time.time() - t0, "protocol": "TF-disjoint cross-fit, folds = scripts/22,23"},
                      open(out / "meta.json", "w"), indent=1)
            if variant != "mid_bce":
                go = L.ROOT / "results" / args.dataset / f"{variant}_g" / f"seed{seed}"
                go.mkdir(parents=True, exist_ok=True)
                np.savez_compressed(go / "score_matrix.npz", S=gmat, tf=D["tf"])
            print(f"[done] {out} in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
