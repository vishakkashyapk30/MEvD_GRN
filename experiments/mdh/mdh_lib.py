"""Shared code for the MDH experiment (README.md): BEAR data, folds, matched AUROC.

Everything is indexed like the processed BEAR builds (scripts/bear/21):
  genes 0..n-1 (measured RNA genes), TF nodes = tf_indices.json (gene ids),
  label tiers / GTs = (2, E) gene-id edges restricted to measured genes.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from src.benchmarks import bear_metrics as bm      # noqa: E402
from src.benchmarks import bear_protocol as bp     # noqa: E402
from src.utils.io import load_config               # noqa: E402

ROOT = Path(os.environ.get("BEAR_ROOT", Path.home() / ".cache/bear/data"))
OUT = Path(__file__).resolve().parent / "results"
DATASETS = ["K562", "Macrophage_S1", "Macrophage_S2", "iPS", "mESC_E7.5_rep1", "mESC_E7.5_rep2",
            "mESC_E8.5_rep1", "mESC_E8.5_rep2", "Naive_mESC"]
GTS = ["ChIP", "KO", "Union", "Intersection", "Core"]
DCFG = load_config(REPO / "configs/bear/datasets.yaml")
BCFG = load_config(REPO / "configs/bear/mevd_base.yaml")["bear"]
RELEASED = ["CellOracle", "FigR", "GRaNIE", "LINGER", "Pando", "Pando_xgb", "SCENIC+", "TRIPOD"]


def _edges(p: Path) -> np.ndarray:
    return torch.load(p).numpy().astype(np.int64).reshape(2, -1)


def load_ds(ds: str, features: bool = True) -> dict:
    p = ROOT / "processed" / ds
    gi = json.load(open(p / "gene_index.json"))
    names = [None] * len(gi)
    for g, i in gi.items():
        names[i] = g
    n = len(names)
    upper = {}
    for i, g in enumerate(names):
        upper.setdefault(str(g).upper(), i)
    tf = np.asarray(json.load(open(p / "tf_indices.json")), dtype=np.int64)
    labels = {q.stem: _edges(q) for q in sorted((p / "labels").glob("*.pt"))}
    gts = {q.stem: _edges(q) for q in sorted((p / "gts").glob("*.pt"))}
    lab = bp.keys_to_edges(np.unique(np.concatenate([bp.edge_keys(e, n) for e in labels.values()])), n)
    D = dict(ds=ds, names=names, upper=upper, n=n, tf=tf, labels=labels, gts=gts, lab=lab,
             tfpos={int(t): k for k, t in enumerate(tf)})
    if features:
        z = np.load(p / "motif_pair_feats.npz")
        assert np.array_equal(z["tf_nodes"], tf)
        D["motif"] = z["F"].astype(np.float32)                      # (m, n, 3)
        sig = np.load(p / "rna_signature.npy").astype(np.float32)   # (n, d), s_i.s_j ~ corr
        D["sig"] = sig[:, np.argsort(-np.linalg.norm(sig, axis=0))]  # leading programs first
        D["fm"] = np.load(p / "fm_gene_embeddings.npy").astype(np.float32)
        D["rna"] = np.load(p / "rna_features_aligned.npy").astype(np.float32)
        D["atac"] = np.load(p / "atac_features_aligned.npy").astype(np.float32)
    return D


def label_matrix(D: dict, edges: np.ndarray) -> np.ndarray:
    """(m TF nodes, n genes) bool from (2, E) edges (sources outside TF nodes dropped)."""
    Y = np.zeros((len(D["tf"]), D["n"]), dtype=bool)
    pos = D["tfpos"]
    keep = np.array([int(s) in pos for s in edges[0]], dtype=bool)
    r = np.array([pos[int(s)] for s in edges[0][keep]], dtype=np.int64)
    Y[r, edges[1][keep]] = True
    return Y


def folds(D: dict, seed: int):
    """Exactly scripts/22 and scripts/23: outer TF-disjoint folds and the inner validation TFs."""
    outdeg = {int(k): int(v) for k, v in zip(*np.unique(D["lab"][0], return_counts=True))}
    outer = bp.tf_disjoint_folds(D["tf"], outdeg, int(BCFG["n_folds"]), seed)
    res = []
    for k, held in enumerate(outer):
        train = np.setdiff1d(D["tf"], held)
        inner = bp.tf_disjoint_folds(train, outdeg, max(2, int(round(1 / float(BCFG["val_tf_frac"])))),
                                     seed + 1000 + k)
        val = inner[0]
        res.append(dict(k=k, held=held, train=train, val=val, fit=np.setdiff1d(train, val)))
    return res


def check_folds(D: dict, seed: int) -> bool:
    m = ROOT / "results" / D["ds"] / "baseline_indegree" / f"seed{seed}" / "meta.json"
    if not m.exists():
        return False
    stored = json.load(open(m))["folds"]
    mine = [f["held"].tolist() for f in folds(D, seed)]
    assert [sorted(x) for x in stored] == [sorted(x) for x in mine], f"fold mismatch {D['ds']} seed {seed}"
    return True


# ------------------------------------------------------------------ evaluation universe
def eval_block(D: dict, gt: str):
    """Rows = GT TFs that are TF nodes; Y (R, n) bool; strata (R, n) by leave-one-TF-out GT in-degree;
    mask (R, n) = valid pairs (no self pairs)."""
    Yall = label_matrix(D, D["gts"][gt])
    rows = np.where(Yall.any(1))[0]
    Y = Yall[rows]
    loo = Y.sum(0)[None, :].astype(np.int64) - Y.astype(np.int64)
    Z = strata(loo)
    mask = np.ones_like(Y, dtype=bool)
    tf_gene = D["tf"][rows]
    mask[np.arange(len(rows)), tf_gene] = False
    return rows, Y, Z, mask


def strata(deg: np.ndarray, nbins: int = 20) -> np.ndarray:
    """0 = degree 0; then quantile bins of the non-zero values (global edges)."""
    out = np.zeros(deg.shape, dtype=np.int64)
    nz = deg > 0
    if nz.any():
        q = np.unique(np.quantile(deg[nz], np.linspace(0, 1, nbins + 1)[1:-1]))
        out[nz] = 1 + np.searchsorted(q, deg[nz], side="right")
    return out


def grouped_avg_rank(gid: np.ndarray, s: np.ndarray) -> np.ndarray:
    """1-based average ranks of s within each group gid (ties averaged)."""
    o = np.lexsort((s, gid))
    g, v = gid[o], s[o]
    n = len(o)
    starts = np.r_[0, np.flatnonzero(g[1:] != g[:-1]) + 1]
    gstart = np.repeat(starts, np.diff(np.r_[starts, n]))
    pos = np.arange(n) - gstart
    change = np.r_[True, (g[1:] != g[:-1]) | (v[1:] != v[:-1])]
    rid = np.cumsum(change) - 1
    avg = np.bincount(rid, weights=pos) / np.bincount(rid)
    r = np.empty(n, dtype=np.float64)
    r[o] = avg[rid] + 1.0
    return r


def matched_stats(S: np.ndarray, Y: np.ndarray, Z: np.ndarray, mask: np.ndarray):
    """Per-row U and n+n- summed over that row's strata. Matched AUROC = U.sum() / D.sum()."""
    R = Y.shape[0]
    nb = int(Z.max()) + 1
    rr = np.repeat(np.arange(R), Y.shape[1]).reshape(Y.shape)
    gid = (rr * nb + Z)[mask]
    s = S[mask].astype(np.float64)
    y = Y[mask]
    rk = grouped_avg_rank(gid, s)
    G = R * nb
    npos = np.bincount(gid, weights=y, minlength=G)
    ntot = np.bincount(gid, minlength=G).astype(np.float64)
    sum_rank_pos = np.bincount(gid, weights=rk * y, minlength=G)
    U = sum_rank_pos - npos * (npos + 1) / 2.0
    Dn = npos * (ntot - npos)
    U = np.where(Dn > 0, U, 0.0)
    rowU = np.bincount(np.arange(G) // nb, weights=U, minlength=R)
    rowD = np.bincount(np.arange(G) // nb, weights=Dn, minlength=R)
    return rowU, rowD


def row_percentiles(S: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Percentile of S_ij against a reference set common to all rows (genes valid in every row, i.e. not an
    evaluated TF's own gene), ties averaged. Rows that order targets identically get identical percentiles."""
    ref = mask.all(0)
    P = np.full(S.shape, np.nan)
    for r in range(S.shape[0]):
        rs = np.sort(S[r, ref])
        lo = np.searchsorted(rs, S[r], side="left")
        hi = np.searchsorted(rs, S[r], side="right")
        P[r] = (lo + hi) / (2.0 * rs.size)
        P[r, ~mask[r]] = np.nan
    return P


def dr_pairs(S: np.ndarray, Y: np.ndarray, mask: np.ndarray, dev: str | None = None):
    """Doubly-ranked AUROC. P = within-row percentile of S (removes any TF-level effect exactly); then within each
    target column, positives' P vs negatives' P (removes any target-level effect exactly: a score that orders the
    targets the same way in every row ties). Returns C[i,k] = sum_j over (Y_ij=1, Y_kj=0) of 1[P_ij>P_kj]+.5 1[=],
    N[i,k] = number of such j. DR-AUROC = sum C / sum N; with TF weights w: w'Cw / w'Nw (cluster bootstrap)."""
    dev = dev or ("cuda" if torch.cuda.is_available() else "cpu")
    P = torch.tensor(np.nan_to_num(row_percentiles(S, mask), nan=-1.0), dtype=torch.float32, device=dev)
    Yt = torch.tensor(Y, device=dev)
    Mt = torch.tensor(mask, device=dev)
    neg = (~Yt) & Mt
    R = Y.shape[0]
    C = torch.zeros(R, R, dtype=torch.float64, device=dev)
    N = torch.zeros(R, R, dtype=torch.float64, device=dev)
    for i in range(R):
        sel = (Yt[i] & Mt[i])[None, :] & neg
        cmp = (P[i][None, :] > P).double() + 0.5 * (P[i][None, :] == P).double()
        C[i] = (cmp * sel).sum(1)
        N[i] = sel.sum(1).double()
    return C.cpu().numpy(), N.cpu().numpy()


def dr_auroc(C: np.ndarray, N: np.ndarray, w: np.ndarray | None = None, F: np.ndarray | None = None) -> float:
    """F (R, R) bool: compare only TF pairs held out in the same fold. Cross-fitted outputs need it: another fold's
    model was trained on TF i's own labels, which biases i-vs-k comparisons across folds."""
    if F is not None:
        C, N = C * F, N * F
    if w is None:
        return float(C.sum() / max(N.sum(), 1))
    return float(w @ C @ w / max(w @ N @ w, 1e-12))


def fold_pairs(D: dict, rows: np.ndarray, seed: int) -> np.ndarray:
    f = np.full(len(D["tf"]), -1)
    for fd in folds(D, seed):
        f[[D["tfpos"][int(t)] for t in fd["held"]]] = fd["k"]
    fr = f[rows]
    return fr[:, None] == fr[None, :]


def fit_beta_model(Y: np.ndarray, mask: np.ndarray, dev: str | None = None):
    """In-sample two-way logistic logit P(Y_ij) = a_i + b_j (the marginal-only generator)."""
    dev = dev or ("cuda" if torch.cuda.is_available() else "cpu")
    y = torch.tensor(Y, dtype=torch.float32, device=dev)
    m = torch.tensor(mask, dtype=torch.float32, device=dev)
    a = torch.zeros(Y.shape[0], 1, device=dev, requires_grad=True)
    b = torch.zeros(1, Y.shape[1], device=dev, requires_grad=True)
    opt = torch.optim.LBFGS([a, b], lr=1, max_iter=300, line_search_fn="strong_wolfe")

    def closure():
        opt.zero_grad()
        loss = (torch.nn.functional.binary_cross_entropy_with_logits(a + b, y, reduction="none") * m).sum() / m.sum()
        loss = loss + 1e-6 * (a.pow(2).sum() + b.pow(2).sum())
        loss.backward()
        return loss
    opt.step(closure)
    return (a + b).detach().cpu().numpy()


def marginal_null_draws(Y: np.ndarray, mask: np.ndarray, n_draws: int = 3, seed: int = 0):
    """Synthetic GTs with the same TF and target propensities and no pair-specific structure."""
    p = 1.0 / (1.0 + np.exp(-fit_beta_model(Y, mask)))
    rng = np.random.default_rng(seed)
    return [(rng.random(Y.shape) < p) & mask for _ in range(n_draws)]


def loo_marginal_oracle(Y: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """logit(TF positive rate) + logit(target rate over the OTHER GT TFs): knows only marginals, no self-inclusion."""
    Yf = (Y & mask).astype(np.float64)
    nr = mask.sum(1, keepdims=True)
    rr = (Yf.sum(1, keepdims=True) + 0.5) / (nr + 1.0)
    nc = mask.sum(0, keepdims=True) - mask.astype(np.float64)
    cc = (Yf.sum(0, keepdims=True) - Yf + 0.5) / (nc + 1.0)
    lg = lambda x: np.log(x / (1 - x))      # noqa: E731
    return lg(rr) + lg(np.clip(cc, 1e-6, 1 - 1e-6))


def full_auroc(S, Y, mask) -> float:
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(Y[mask], S[mask]))


# ------------------------------------------------------------------ score matrices
def edges_to_matrix(D: dict, rows: np.ndarray, e: pd.DataFrame) -> np.ndarray:
    """e: standardised (tf, target, score) upper-case, |score|, max over duplicates. Unscored = 0."""
    S = np.zeros((len(rows), D["n"]), dtype=np.float64)
    rpos = {str(D["names"][D["tf"][r]]).upper(): k for k, r in enumerate(rows)}
    ri = e["tf"].map(rpos)
    ci = e["target"].map(D["upper"])
    ok = ri.notna() & ci.notna()
    sc = np.nan_to_num(e["score"].to_numpy(dtype=np.float64)[ok.to_numpy()], nan=0.0)
    S[ri[ok].astype(int).to_numpy(), ci[ok].astype(int).to_numpy()] = sc
    return S


def read_ours(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[""])
    return bm.standardise_edges(df, "Source", "Target", "Score")


def released_file(ds: str, method: str):
    d = ROOT / "INFERRED.GRNS" / DCFG["inferred_dir"][ds] / method
    if not d.exists():
        return None
    files = sorted(f for f in d.iterdir() if f.is_file())
    return files[0] if files else None


def read_released(ds: str, method: str):
    f = released_file(ds, method)
    if f is None or method not in bm.ROC_PR_COLUMNS:
        return None
    df = bm.read_table(str(f))
    a, b, c = bm.ROC_PR_COLUMNS[method]
    if c not in df.columns or a not in df.columns:
        return None
    return bm.standardise_edges(df, a, b, c)


def gt_full(ds: str, gt: str):
    p = DCFG["datasets"][ds]["gts"].get(gt)
    if p is None or not (ROOT / p).exists():
        return None
    return bm.load_gt(str(ROOT / p))


def pooled_score(e: pd.DataFrame, ds: str, gt: str, n_rep: int = 20) -> dict:
    g = gt_full(ds, gt)
    return bm.score_roc_pr(e, g, np.random.default_rng(43), n_rep=n_rep)


def ds_gts(D: dict):
    return [g for g in GTS if g in D["gts"]]


def boot_ci(diffs_fn, R_by_ds: dict, n_boot: int = 2000, seed: int = 0):
    """Cluster bootstrap over TFs within each dataset; diffs_fn(ds, idx) -> per-dataset statistic;
    returns (point, lo, hi) of the across-dataset mean."""
    rng = np.random.default_rng(seed)
    pt = np.mean([diffs_fn(ds, np.arange(R)) for ds, R in R_by_ds.items()])
    bs = []
    for _ in range(n_boot):
        bs.append(np.mean([diffs_fn(ds, rng.integers(0, R, R)) for ds, R in R_by_ds.items()]))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return float(pt), float(lo), float(hi), np.asarray(bs)
