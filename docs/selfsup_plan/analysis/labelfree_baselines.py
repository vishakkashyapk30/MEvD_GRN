#!/usr/bin/env python
"""Read-only CPU diagnostics for the self-supervised plan (docs/selfsup_plan/plan.md).

Mode `baselines` (one dataset per call): label-free, TF-agnostic or motif-based
edge scores over (TF node x measured gene), scored with BEAR-GRN's rules
(re-implemented vectorised here, semantics identical to
src/benchmarks/bear_metrics.score_roc_pr: AUROC on the method's own edges with
1:1 down-sampling; AUPRC on the GT-TF x GT-target universe, unpredicted = 0,
negatives capped at 10x positives, PRROC auc.integral). Also reports two
TF-specificity diagnostics on the measured sub-universe (GT TFs that are TF
nodes x GT targets that are measured genes):
  * within-target AUROC: for each target gene, rank the TFs by score, AUROC
    against that gene's GT column, averaged over genes with >= 1 pos and
    >= 1 neg TF. Any TF-agnostic score gives exactly 0.5.
  * within-TF AUROC: per TF across targets, averaged over TFs.
Rows prefixed `oracle_` are built FROM THE GROUND TRUTH and exist only to
measure how much of the benchmark is explained by node marginals (headroom
diagnostics). They are never candidate methods.

Mode `released` (K562 only; the only dataset whose released GRNs are on disk):
for each released BEAR method, BEAR AUROC/AUPRC of (a) the original scores,
(b) TF-identity permutations (3 draws), (c) the target-marginal projection
(each edge gets the mean score of its target over the method's TFs), plus the
within-target AUROC of the original scores.

Inputs: ~/.cache/bear/data/processed/<ds>/* (built label-free by scripts/21)
and the BEAR GT files. Nothing is written outside docs/selfsup_plan/analysis/.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score

REPO = Path("/home/vishak/research/mevd_grn")
sys.path.insert(0, str(REPO))
from src.benchmarks import bear_metrics as bm  # noqa: E402

ROOT = Path.home() / ".cache/bear/data"
OUT = REPO / "docs/selfsup_plan/analysis"
CFG = yaml.safe_load(open(REPO / "configs/bear/datasets.yaml"))


# --------------------------------------------------------------------------- fast BEAR scorer
class Scorer:
    """Pre-indexes one GT once; scores many score vectors over the same pairs."""

    def __init__(self, gt: pd.DataFrame):
        self.tfs = pd.unique(gt["Source"])
        self.tgs = pd.unique(gt["Target"])
        self.tf_ix = {t: i for i, t in enumerate(self.tfs)}
        self.tg_ix = {t: i for i, t in enumerate(self.tgs)}
        self.n_tf, self.n_tg = len(self.tfs), len(self.tgs)
        k = gt["Source"].map(self.tf_ix).to_numpy(np.int64) * self.n_tg + \
            gt["Target"].map(self.tg_ix).to_numpy(np.int64)
        self.is_pos = np.zeros(self.n_tf * self.n_tg, dtype=bool)
        self.is_pos[k] = True
        self.pos_idx = np.flatnonzero(self.is_pos)
        self.neg_idx = np.flatnonzero(~self.is_pos)
        self.n_neg10 = min(len(self.pos_idx) * 10, len(self.neg_idx))

    def keys(self, tf_names: np.ndarray, tg_names: np.ndarray):
        """-> (mask of pairs inside the tested space, universe keys of those pairs)."""
        a = pd.Series(tf_names).map(self.tf_ix)
        b = pd.Series(tg_names).map(self.tg_ix)
        ok = (a.notna() & b.notna()).to_numpy()
        key = a[ok].to_numpy(np.int64) * self.n_tg + b[ok].to_numpy(np.int64)
        return ok, key

    def score(self, key: np.ndarray, sc: np.ndarray, rng, n_rep: int = 3, rand: bool = False):
        lab = self.is_pos[key]
        pos_s, neg_s = sc[lab], sc[~lab]
        out = {"n_eval": int(len(key)), "n_pos": int(lab.sum())}
        if lab.sum() == 0 or (~lab).sum() == 0:
            out.update(AUROC=np.nan, AUPRC=np.nan)
            return out
        full = np.zeros(self.n_tf * self.n_tg)
        full[key] = np.nan_to_num(sc, nan=0.0)
        fixed = self.n_neg10 == len(self.neg_idx)
        au, pr = [], []
        for r in range(n_rep):
            if len(neg_s) < len(pos_s):
                ps, ns = rng.choice(pos_s, size=len(neg_s), replace=False), neg_s
            else:
                ps, ns = pos_s, rng.choice(neg_s, size=len(pos_s), replace=False)
            y = np.r_[np.ones(len(ps)), np.zeros(len(ns))]
            au.append(roc_auc_score(y, np.r_[ps, ns]))
            if fixed and r > 0:
                pr.append(pr[0])
            else:
                nsel = self.neg_idx if fixed else rng.choice(self.neg_idx, self.n_neg10, replace=False)
                pr.append(bm.prroc_auc_integral(full[self.pos_idx], full[nsel]))
        out.update(AUROC=float(np.mean(au)), AUPRC=float(np.mean(pr)))
        if rand:
            nsel = self.neg_idx if fixed else rng.choice(self.neg_idx, self.n_neg10, replace=False)
            u = rng.random(self.n_tf * self.n_tg)
            out["AUPRC_random"] = bm.prroc_auc_integral(u[self.pos_idx], u[nsel])
        return out


def mw_auc_axis(S: np.ndarray, Y: np.ndarray, axis: int) -> float:
    """Mean Mann-Whitney AUROC along `axis` (ties = 1/2), over slices with >=1 pos and >=1 neg."""
    R = rankdata(S, axis=axis, method="average")
    npos = Y.sum(axis=axis).astype(np.float64)
    n = Y.shape[axis]
    nneg = n - npos
    rs = (R * Y).sum(axis=axis)
    ok = (npos > 0) & (nneg > 0)
    auc = (rs[ok] - npos[ok] * (npos[ok] + 1) / 2) / (npos[ok] * nneg[ok])
    return float(auc.mean()) if ok.any() else float("nan")


def swap_gap(S: np.ndarray, Y: np.ndarray):
    """TF-identity swap test. A[u, t] = AUROC of TF u's score row against TF t's
    GT row (over targets). own = A[t, t]; other = mean_{u != t} A[u, t].
    Returns (mean own, mean other, mean gap, fraction of TFs whose own row beats
    every other row's 95th percentile). A TF-agnostic score gives gap = 0 exactly;
    target hubness and TF marginals cannot create a gap (shared by all rows /
    constant within a row)."""
    T, G = S.shape
    if T < 3:
        return (np.nan,) * 4
    R = rankdata(S, axis=1, method="average")
    Yf = Y.astype(np.float64)
    npos = Yf.sum(1)
    nneg = G - npos
    sums = R @ Yf.T                                    # [u, t]
    ok = (npos > 0) & (nneg > 0)
    A = (sums - (npos * (npos + 1) / 2)[None, :]) / (npos * nneg)[None, :]
    A = A[:, ok][ok, :] if ok.sum() >= 3 else None
    if A is None:
        return (np.nan,) * 4
    own = np.diag(A)
    oth = (A.sum(0) - own) / (A.shape[0] - 1)
    q95 = np.array([np.quantile(np.delete(A[:, t], t), 0.95) for t in range(A.shape[0])])
    return float(own.mean()), float(oth.mean()), float((own - oth).mean()), float((own > q95).mean())


def gts_of(ds):
    out = {}
    for g, p in CFG["datasets"][ds]["gts"].items():
        fp = ROOT / p
        if fp.exists():
            out[g] = bm.load_gt(str(fp))
    return out


def load_proc(ds):
    p = ROOT / "processed" / ds
    gi = json.load(open(p / "gene_index.json"))
    names = np.empty(len(gi), dtype=object)
    for g, i in gi.items():
        names[i] = str(g).upper()
    tf = np.asarray(json.load(open(p / "tf_indices.json")), dtype=np.int64)
    rna = np.load(p / "rna_features_aligned.npy")
    atac = np.load(p / "atac_features_aligned.npy")
    op = np.load(p / "openness.npy")
    m = np.load(p / "motif_pair_feats.npz")
    F, mt = m["F"].astype(np.float32), m["tf_nodes"]
    assert np.array_equal(mt, tf), "motif tf_nodes != tf_indices"
    return names, tf, rna, atac, op, F


def pct(x):
    return rankdata(x, method="average") / len(x)


# --------------------------------------------------------------------------- mode: baselines
def run_baselines(ds: str, n_rep: int, seed: int):
    t0 = time.time()
    names, tf, rna, atac, op, F = load_proc(ds)
    n, T = len(names), len(tf)
    # candidate block: every TF node x every measured gene, no self pairs
    src = np.repeat(np.arange(T), n)
    dst = np.tile(np.arange(n), T)
    keep = tf[src] != dst
    src, dst = src[keep], dst[keep]
    tfn, tgn = names[tf[src]], names[dst]
    has_motif = F[:, :, 2].max(axis=1) > 0                       # per TF node
    g_expr, g_det, g_atac = pct(rna[:, 0]), pct(rna[:, 2]), pct(atac[:, 0])
    g_rp, g_npk = pct(op[:, 0]), pct(op[:, 3])
    t_expr = pct(rna[tf, 0])
    motif = F[src, dst, 0]
    variants = {
        "coverage": np.ones(len(src)),
        "tgt_expr_mean": 1 + g_expr[dst],
        "tgt_detection": 1 + g_det[dst],
        "tgt_atac_activity": 1 + g_atac[dst],
        "tgt_rp_proximal": 1 + g_rp[dst],
        "tgt_n_peaks": 1 + g_npk[dst],
        "tfexpr_x_tgtatac": 1 + t_expr[src] * g_atac[dst],
        # motif x accessibility: TFs with a JASPAR motif use best relative PWM score in
        # the gene's top-5 RP peaks times target accessibility; TFs without a motif fall
        # back to target accessibility alone (so coverage is unchanged).
        "motif_x_atac": 1 + np.where(has_motif[src], motif * g_atac[dst], 0.5 * g_atac[dst]),
    }
    rows = []
    for gname, gt in gts_of(ds).items():
        sc = Scorer(gt)
        ok, key = sc.keys(tfn, tgn)
        # measured sub-universe matrix for TF-specificity diagnostics
        tf_in = np.array([names[t] in sc.tf_ix for t in tf])
        g_in = np.array([g in sc.tg_ix for g in names])
        rows_tf = np.flatnonzero(tf_in)
        cols_g = np.flatnonzero(g_in)
        Y = np.zeros((len(rows_tf), len(cols_g)), dtype=bool)
        rmap = {int(r): i for i, r in enumerate(rows_tf)}
        cmap = {int(c): j for j, c in enumerate(cols_g)}
        gtp = gt[gt["Source"].isin(set(names[tf[rows_tf]])) & gt["Target"].isin(set(names[cols_g]))]
        nm2tf = {names[tf[r]]: rmap[int(r)] for r in rows_tf}
        nm2g = {names[c]: cmap[int(c)] for c in cols_g}
        yi, yj = gtp["Source"].map(nm2tf), gtp["Target"].map(nm2g)
        okm = yi.notna() & yj.notna()
        if okm.sum() == 0:
            print(f"[{ds} {gname}] no GT edge inside the measured TF x gene block; skipped", flush=True)
            continue
        Y[yi[okm].astype(np.int64).to_numpy(), yj[okm].astype(np.int64).to_numpy()] = True
        # remove self pairs from diagnostics by masking (they are absent from candidates)
        def as_matrix(v):
            M = np.zeros((T, n))
            M[src, dst] = v
            return M[np.ix_(rows_tf, cols_g)]
        rng = np.random.default_rng(seed)
        first = True
        for vname, v in variants.items():
            r = sc.score(key, v[ok], rng, n_rep=n_rep, rand=first)
            if first:
                rand_ap = r["AUPRC_random"]
                first = False
            M = as_matrix(v)
            so, sx, sg, sf = swap_gap(M, Y)
            rows.append({"dataset": ds, "gt": gname, "variant": vname, **r, "AUPRC_random": rand_ap,
                         "within_target_AUROC": mw_auc_axis(M, Y, 0),
                         "within_tf_AUROC": mw_auc_axis(M, Y, 1),
                         "swap_own": so, "swap_other": sx, "swap_gap": sg, "swap_frac_tf_sig": sf,
                         "measured_density": float(Y.mean()), "n_tf_meas": int(Y.shape[0]),
                         "n_tg_meas": int(Y.shape[1])})
        # ---- oracle diagnostics (GT-derived; NOT methods) ----
        Yf = Y.astype(np.float64)
        indeg = Yf.sum(0, keepdims=True)
        outdeg = Yf.sum(1, keepdims=True)
        orc = {
            # leave-own-TF-out target in-degree: the strongest form of the supervised hub reference
            "oracle_loo_indegree": indeg - Yf,
            # TF out-degree (TF marginal): the within-target trivial reference
            "oracle_tf_outdegree": np.broadcast_to(outdeg - Yf, Yf.shape),
            # rank-1 marginal model (both marginals, own edge removed)
            "oracle_degree_product": (indeg - Yf + 1) * (outdeg - Yf + 1),
        }
        ii, jj = np.nonzero(np.ones_like(Y))
        o_tf, o_tg = names[tf[rows_tf]][ii], names[cols_g][jj]
        nonself = o_tf != o_tg
        ii, jj, o_tf, o_tg = ii[nonself], jj[nonself], o_tf[nonself], o_tg[nonself]
        ok2, key2 = sc.keys(o_tf, o_tg)
        for oname, M in orc.items():
            v = np.asarray(M, dtype=np.float64)[ii, jj] + 1.0
            r = sc.score(key2, v[ok2], rng, n_rep=n_rep)
            rows.append({"dataset": ds, "gt": gname, "variant": oname, **r, "AUPRC_random": rand_ap,
                         "within_target_AUROC": mw_auc_axis(np.asarray(M), Y, 0),
                         "within_tf_AUROC": mw_auc_axis(np.asarray(M), Y, 1),
                         "measured_density": float(Y.mean()), "n_tf_meas": int(Y.shape[0]),
                         "n_tg_meas": int(Y.shape[1])})
        print(f"[{ds} {gname}] done ({time.time()-t0:.0f}s)", flush=True)
    df = pd.DataFrame(rows)
    df["motif_tfs"] = int(has_motif.sum())
    df["tf_nodes"] = T
    df.to_csv(OUT / f"labelfree__{ds}.csv", index=False)
    print(f"[write] labelfree__{ds}.csv ({time.time()-t0:.0f}s)")


# --------------------------------------------------------------------------- mode: released (K562)
def run_released(n_rep: int, seed: int, methods=None, ds: str = "K562", inf_root: Path | None = None):
    t0 = time.time()
    rel_root = (inf_root or ROOT) / "INFERRED.GRNS" / CFG["inferred_dir"][ds]
    gts = gts_of(ds)
    rows = []
    for mdir, (tc, gc, scol) in bm.ROC_PR_COLUMNS.items():
        if methods and mdir not in methods:
            continue
        d = rel_root / mdir
        files = sorted(f for f in d.iterdir() if f.is_file()) if d.exists() else []
        if not files:
            continue
        df = bm.read_table(str(files[0]))
        if scol not in df.columns:
            continue
        e = bm.standardise_edges(df, tc, gc, scol)
        del df
        e = e[np.isfinite(e["score"].to_numpy(np.float64))]
        for gname in ["ChIP", "Union", "KO"]:
            gt = gts[gname]
            sc = Scorer(gt)
            ok, key = sc.keys(e["tf"].to_numpy(), e["target"].to_numpy())
            s = e["score"].to_numpy(np.float64)[ok]
            tfs_in = e["tf"].to_numpy()[ok]
            tgs_in = e["target"].to_numpy()[ok]
            rng = np.random.default_rng(seed)
            base = sc.score(key, s, rng, n_rep=n_rep)
            # within-target AUROC on the method's tested-space matrix (unpredicted = 0)
            ti = pd.Series(tfs_in).map(sc.tf_ix).to_numpy()
            gi_ = pd.Series(tgs_in).map(sc.tg_ix).to_numpy()
            used_t = np.unique(ti)
            used_g = np.unique(gi_)
            tpos = {t: i for i, t in enumerate(used_t)}
            gpos = {g: j for j, g in enumerate(used_g)}
            M = np.zeros((len(used_t), len(used_g)))
            M[[tpos[t] for t in ti], [gpos[g] for g in gi_]] = s
            Yb = sc.is_pos.reshape(sc.n_tf, sc.n_tg)[np.ix_(used_t, used_g)]
            wt = mw_auc_axis(M, Yb, 0)
            so, sx, sg, sf = swap_gap(M, Yb)
            # (b) TF-identity permutation among the method's tested TFs
            perm_au, perm_pr = [], []
            for p in range(3):
                pm = rng.permutation(used_t)
                remap = dict(zip(used_t, pm))
                ti_p = np.array([remap[t] for t in ti], dtype=np.int64)
                key_p = ti_p * sc.n_tg + gi_.astype(np.int64)
                # collapse duplicates created by the permutation (none: bijection on TFs)
                rp = sc.score(key_p, s, rng, n_rep=n_rep)
                perm_au.append(rp["AUROC"]); perm_pr.append(rp["AUPRC"])
            # (c) target-marginal projection: score := mean score of the target over TFs
            gmean = pd.Series(s).groupby(gi_).transform("mean").to_numpy()
            proj = sc.score(key, gmean, rng, n_rep=n_rep)
            # (d) TF-normalised scores: percentile rank within each TF (removes the TF marginal,
            # i.e. how the method spreads score mass across TFs); (0, 1], so coverage is unchanged
            tfn = pd.Series(s).groupby(ti).rank(pct=True).to_numpy()
            tnorm = sc.score(key, tfn, rng, n_rep=n_rep)
            rows.append({"dataset": ds, "gt": gname, "method": mdir, "n_eval": base["n_eval"],
                         "AUROC": base["AUROC"], "AUPRC": base["AUPRC"],
                         "AUROC_tfperm": float(np.mean(perm_au)), "AUPRC_tfperm": float(np.mean(perm_pr)),
                         "AUROC_tgtproj": proj["AUROC"], "AUPRC_tgtproj": proj["AUPRC"],
                         "AUROC_tfnorm": tnorm["AUROC"], "AUPRC_tfnorm": tnorm["AUPRC"],
                         "within_target_AUROC": wt, "swap_own": so, "swap_other": sx, "swap_gap": sg,
                         "swap_frac_tf_sig": sf, "n_tf_used": len(used_t), "n_tg_used": len(used_g)})
            print(f"[{ds} {gname}] {mdir}: AUROC {base['AUROC']:.3f} perm {np.mean(perm_au):.3f} "
                  f"proj {proj['AUROC']:.3f} | AUPRC {base['AUPRC']:.3f} perm {np.mean(perm_pr):.3f} "
                  f"proj {proj['AUPRC']:.3f} | within-target {wt:.3f} ({time.time()-t0:.0f}s)", flush=True)
    tag = "_".join(methods) if methods else "all"
    pd.DataFrame(rows).to_csv(OUT / f"released_tfspec__{ds}__{tag}.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["baselines", "released"])
    ap.add_argument("--dataset", default="K562")
    ap.add_argument("--methods", nargs="*", default=None)
    ap.add_argument("--n_rep", type=int, default=3)
    ap.add_argument("--seed", type=int, default=43)
    ap.add_argument("--inf_root", default=None, help="root holding INFERRED.GRNS/ (default ~/.cache/bear/data)")
    a = ap.parse_args()
    if a.mode == "baselines":
        run_baselines(a.dataset, a.n_rep, a.seed)
    else:
        run_released(a.n_rep, a.seed, a.methods, a.dataset, Path(a.inf_root) if a.inf_root else None)
