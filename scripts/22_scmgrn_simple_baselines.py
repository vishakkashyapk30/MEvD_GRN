#!/usr/bin/env python
"""Statistical and simple baselines on the scMultiomeGRN benchmark (no neural networks).

Four families, all scored per official fold on the held-out test pairs with the official metric_fn
(src/benchmarks/scmgrn_protocol.py). Anything fitted or counted uses the fold's TRAINING part only.

  A1  graph-topology heuristics on the training links (link-prediction literature): degree sum, common
      neighbours, Jaccard, Adamic-Adar, resource allocation, Katz (beta 0.05, paths <= 4), random walk with
      restart (alpha 0.15), truncated-SVD reconstruction (k 16), normalised-adjacency spectral reconstruction (k 16).
      (Degree product = preferential attachment is scripts/21's degree_product.)
  A2  unsupervised GRN-inference scores from the cell data (BEELINE / geneRNIB family; no labels at all):
      |Pearson|, |Spearman|, |partial correlation| (Ledoit-Wolf precision), mutual information and CLR,
      GRNBoost2 between the two TFs, GENIE3 (random forest, TFs as regulators and targets),
      per-target ridge and lasso regression (|coefficient|), |Pearson| of MAESTRO accessibility.
  A3  simple supervised models on features: logistic regression on MAESTRO / GRNBoost2 / Geneformer /
      pair histograms alone, MAESTRO + GRNBoost2 given independently (concatenated), multiplied across
      modalities (m_u * g_u per TF) and within the pair (x_u * x_v only); linear (ridge) regression on the 0/1
      label and a random forest on the concatenated features.
  A4  topology + features: logistic regression on degree / heuristic pair features alone, with MAESTRO +
      GRNBoost2, and with Geneformer. Graph-derived features of TRAINING pairs are computed edge-held-out
      (edge_heldout_feats): a test link is never in the training graph, so its own link must not inform its
      training counterpart's features either.
  Every logistic regression picks its L2 strength by 3-fold CV on the training pairs (LogisticRegressionCV).

Pair features are symmetric: [x_u * x_v, |x_u - x_v|] on per-TF vectors compressed to k = 32 dimensions with a
label-free SVD (src/benchmarks/scmgrn_features.svd_reduce). Every score matrix is symmetrised.

Besides the official test fold, each result also reports "test_degmatched": the same test positives with
negatives drawn to match their (training) degree distribution, following "Hidden sampling biases inflate
performance in gene regulatory network inference" (bioRxiv 2025). A pure degree rule is near chance there.

Writes <root>/results/<ds>/<ct>/<method>/seed0/foldNN.json (same format as scripts/21).
Usage: python scripts/22_scmgrn_simple_baselines.py --root ~/scmgrn --ds lung --ct Myeloid_cells --family A1
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.covariance import LedoitWolf
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import Lasso, LogisticRegression, LogisticRegressionCV, Ridge
from sklearn.preprocessing import StandardScaler

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from src.benchmarks import scmgrn_protocol as proto            # noqa: E402
from src.benchmarks.scmgrn_features import svd_reduce          # noqa: E402
from src.utils.io import load_config, save_json                # noqa: E402

K_SVD = 32
VERSION = 3                      # bumped when a method changes; fold files of other versions are recomputed
CHUNKS = 5                       # edge-held-out chunks for features computed from the training graph
MAX_CELLS_REGRESSION = 3000      # GENIE3 / lasso / ridge: cells subsampled (seeded) to keep runtime bounded


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


# ----------------------------------------------------------------------------- loading
def load_rows(path: Path, names, sep=",", chunksize=400):
    """Rows of a (large) feature-by-column text table whose index is in `names`, in `names` order."""
    want = set(names)
    parts = [ch[ch.index.isin(want)] for ch in pd.read_csv(path, sep=sep, index_col=0, chunksize=chunksize)]
    df = pd.concat(parts)
    df = df[~df.index.duplicated(keep="first")]
    return df.reindex(names).fillna(0.0)


class CellType:
    def __init__(self, root: Path, cfg: dict, ct: str, split_kind: str = "official"):
        self.root, self.cfg, self.ct, self.split_kind = root, cfg, ct, split_kind
        self.gdir = root / cfg["graph_dir"] / ct
        self.pdir = root / cfg["processed_dir"] / ct
        s = cfg["split"]
        self.sdir = proto.resolve_split_dir(str(root / cfg["splits_dir"] / ct), s, split_kind)
        self.n_folds = int(s["n_splits"])
        self.adj, self.nodes = proto.load_adj(str(self.gdir / "x_graph.txt"))
        self.n = len(self.nodes)
        self._cache = {}

    def get(self, key):
        if key in self._cache:
            return self._cache[key]
        if key == "expr":            # TF x cells, log1p counts
            v = np.log1p(load_rows(self.gdir / "tmp" / "scrna.csv", self.nodes).to_numpy(np.float64))
        elif key == "maestro":       # TF x cells, MAESTRO regulatory potential
            v = load_rows(self.gdir / "x_atac_rp_score_feature.txt", self.nodes, sep="\t").to_numpy(np.float64)
        elif key == "grnboost":      # TF x genes, GRNBoost2 importance TF -> gene
            df = pd.read_csv(self.gdir / "x_scRNA_grnboost2_feature.txt", sep="\t", index_col=0)
            self._cache["grnboost_df"] = df.reindex(self.nodes).fillna(0.0)
            v = self._cache["grnboost_df"].to_numpy(np.float64)
        elif key == "geneformer":    # TF x 768, rows in processed/gene_index order (= x_node.txt order)
            gi = json.load(open(self.pdir / "gene_index.json"))
            E = np.load(self.pdir / "fm_gene_embeddings.npy")
            order = [gi[n.upper()] for n in self.nodes]
            v = E[order].astype(np.float64)
        elif key == "openness":      # TF x 4 promoter-accessibility summaries (scripts/16 step processed)
            gi = json.load(open(self.pdir / "gene_index.json"))
            v = np.load(self.pdir / "openness.npy")[[gi[n.upper()] for n in self.nodes]].astype(np.float64)
        elif key == "histogram":     # N x N x 136 symmetric pair histograms (scripts/16 --steps smg)
            v = np.load(self.pdir / "smg_pair_feats.npy").astype(np.float32)
        else:
            raise KeyError(key)
        self._cache[key] = v
        return v

    def reduced(self, key):
        k = f"{key}_svd"
        if k not in self._cache:
            self._cache[k] = svd_reduce(self.get(key), K_SVD)
        return self._cache[k]


# ----------------------------------------------------------------------------- fold helpers
def unordered(a):
    a = np.asarray(a).reshape(-1, 2)
    a = a[a[:, 0] < a[:, 1]]
    return a


def train_adjacency(n, train_pos):
    A = np.zeros((n, n))
    A[train_pos[:, 0], train_pos[:, 1]] = 1
    A = np.maximum(A, A.T)
    np.fill_diagonal(A, 0)
    return A


def evaluate(S, ct: CellType, f, fold):
    return proto.evaluate_all(S, ct.adj, f, ct.sdir, fold)


def pair_block(X, i, j, products_only=False):
    a, b = X[i], X[j]
    return a * b if products_only else np.concatenate([a * b, np.abs(a - b)], axis=1)


def LR():
    """L2 logistic regression; strength picked by 3-fold CV (average precision) on the training pairs only."""
    return LogisticRegressionCV(Cs=[1e-3, 1e-2, 1e-1, 1.0, 10.0], cv=3, scoring="average_precision", max_iter=3000)


def edge_heldout_feats(build, A, tr_pos, tr_pairs, seed):
    """Training-pair features from graph-derived quantities without the pair's own link.

    A training positive is a link of the training graph, a test positive never is, so features computed on the
    full training graph (degree, embeddings, paths) differ between the two and a classifier learns the wrong
    pattern. Standard fix in link prediction: split the training positives into CHUNKS parts and compute the
    features of each part's pairs (negatives assigned to parts at random) on the graph without that part's
    links. build(A_sub) -> feats_fn(i, j). Test/all pairs use build(A) on the full training graph."""
    rng = np.random.default_rng(seed)
    pos_set = {(u, v): k for k, (u, v) in enumerate(tr_pos)}
    part_pos = rng.integers(CHUNKS, size=len(tr_pos))
    part = np.array([part_pos[pos_set[(u, v)]] if (u, v) in pos_set else rng.integers(CHUNKS) for u, v in tr_pairs])
    X = None
    for c in range(CHUNKS):
        Asub = A.copy()
        e = tr_pos[part_pos == c]
        Asub[e[:, 0], e[:, 1]] = 0
        Asub[e[:, 1], e[:, 0]] = 0
        idx = np.where(part == c)[0]
        Xc = build(Asub)(tr_pairs[idx, 0], tr_pairs[idx, 1])
        if X is None:
            X = np.zeros((len(tr_pairs), Xc.shape[1]))
        X[idx] = Xc
    return X


def fit_predict(make_model, feats_fn, n, tr_pairs, tr_lab, train_feats=None):
    """Fit on training pairs; return an N x N score matrix for all u<v pairs (diagonal 0)."""
    Xtr = feats_fn(tr_pairs[:, 0], tr_pairs[:, 1]) if train_feats is None else train_feats
    sc = StandardScaler().fit(Xtr)
    model = make_model().fit(sc.transform(Xtr), tr_lab)
    iu, ju = np.triu_indices(n, 1)
    Xall = sc.transform(feats_fn(iu, ju))
    if hasattr(model, "predict_proba"):
        p = model.predict_proba(Xall)[:, 1]
    else:
        p = model.predict(Xall)
    S = np.zeros((n, n))
    S[iu, ju] = p
    return S + S.T


# ----------------------------------------------------------------------------- A1 heuristics
def heuristics(A):
    n = A.shape[0]
    deg = A.sum(1)
    CN = A @ A
    union = deg[:, None] + deg[None, :] - CN
    jac = np.divide(CN, union, out=np.zeros_like(CN), where=union > 0)
    w_aa = np.where(deg > 1, 1.0 / np.log(np.maximum(deg, 2)), 0.0)
    AA = (A * w_aa[None, :]) @ A
    w_ra = np.where(deg > 0, 1.0 / np.maximum(deg, 1), 0.0)
    RA = (A * w_ra[None, :]) @ A
    beta, Ak, katz = 0.05, A.copy(), np.zeros_like(A)
    katz2 = np.zeros_like(A)                     # paths of length 2..4 only (no direct-edge term)
    for l in range(1, 5):
        katz += beta ** l * Ak
        if l >= 2:
            katz2 += beta ** l * Ak
        Ak = Ak @ A
    P = np.divide(A, deg[:, None], out=np.zeros_like(A), where=deg[:, None] > 0)
    alpha = 0.15
    R = alpha * np.linalg.inv(np.eye(n) - (1 - alpha) * P)
    U, s, Vt = np.linalg.svd(A)
    k = min(16, n - 1)
    svd = (U[:, :k] * s[:k]) @ Vt[:k]
    return {"degree_product": (deg[:, None] + 1) * (deg[None, :] + 1), "degree_sum": deg[:, None] + deg[None, :], "common_neighbours": CN, "jaccard": jac,
            "adamic_adar": AA, "resource_allocation": RA, "katz": katz, "katz_no_direct": katz2,
            "rwr": R + R.T, "svd_lowrank": svd}


def spectral_lowrank(A, k=16):
    """Low-rank reconstruction of the normalised training adjacency D^-1/2 A D^-1/2 (top-k eigenpairs)."""
    deg = A.sum(1)
    d = np.where(deg > 0, 1.0 / np.sqrt(np.maximum(deg, 1e-12)), 0.0)
    w, V = np.linalg.eigh(d[:, None] * A * d[None, :])
    top = np.argsort(-w)[:min(k, A.shape[0])]
    return (V[:, top] * w[top]) @ V[:, top].T


def topo_pair_feats(H, deg):
    ld = np.log1p(deg)

    def fn(i, j):
        return np.stack([ld[i] + ld[j], np.abs(ld[i] - ld[j]), ld[i] * ld[j],
                         np.log1p(H["common_neighbours"][i, j]), H["jaccard"][i, j], H["adamic_adar"][i, j],
                         H["resource_allocation"][i, j], H["katz_no_direct"][i, j]], axis=1)
    return fn


# ----------------------------------------------------------------------------- A2 unsupervised scores
def mutual_information(X, bins=8):
    n, c = X.shape
    B = np.zeros((n, c), dtype=np.int64)
    for i in range(n):
        lo, hi = X[i].min(), X[i].max()
        if hi == lo:
            continue
        B[i] = np.clip(((X[i] - lo) / (hi - lo) * bins).astype(int), 0, bins - 1)
    px = np.stack([np.bincount(B[i], minlength=bins) / c for i in range(n)])
    MI = np.zeros((n, n))
    off = (np.arange(n, dtype=np.int64) * bins * bins)[:, None]
    for i in range(n):
        code = B * bins + B[i][None, :] + off
        J = np.bincount(code.ravel(), minlength=n * bins * bins).reshape(n, bins, bins) / c   # J[j, b_j, b_i]
        denom = px[:, :, None] * px[i][None, None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            t = np.where(J > 0, J * np.log(J / np.where(denom > 0, denom, 1)), 0.0)
        MI[i] = t.sum(axis=(1, 2))
    np.fill_diagonal(MI, 0)
    return (MI + MI.T) / 2


def clr(M):
    mu, sd = M.mean(1, keepdims=True), M.std(1, keepdims=True) + 1e-12
    z = np.maximum(0, (M - mu) / sd)
    return np.sqrt(z ** 2 + z.T ** 2)


def per_target_regression(X, kind, rng, n_jobs):
    """W[r, t] = importance of TF r for TF t (other TFs as regulators)."""
    n, c = X.shape
    cols = rng.choice(c, size=min(c, MAX_CELLS_REGRESSION), replace=False)
    Z = StandardScaler().fit_transform(X[:, cols].T)          # cells x TFs
    W = np.zeros((n, n))
    for t in range(n):
        reg = np.array([r for r in range(n) if r != t])
        y = Z[:, t]
        if np.allclose(y.std(), 0):
            continue
        if kind == "ridge":
            w = np.abs(Ridge(alpha=1.0).fit(Z[:, reg], y).coef_)
        elif kind == "lasso":
            w = np.abs(Lasso(alpha=0.01, max_iter=3000).fit(Z[:, reg], y).coef_)
        else:  # genie3
            w = RandomForestRegressor(n_estimators=300, max_features="sqrt", n_jobs=n_jobs,
                                      random_state=0).fit(Z[:, reg], y).feature_importances_
        W[reg, t] = w
    return W + W.T


def cosine(M):
    Z = M / (np.linalg.norm(M, axis=1, keepdims=True) + 1e-12)
    return Z @ Z.T


def unsupervised_scores(ct: CellType, n_jobs):
    out = {}
    X = ct.get("expr")
    out["pearson_abs"] = np.abs(np.nan_to_num(np.corrcoef(X)))
    out["spearman_abs"] = np.abs(np.nan_to_num(stats.spearmanr(X, axis=1).correlation))
    prec = LedoitWolf().fit(StandardScaler().fit_transform(X.T)).precision_
    d = np.sqrt(np.diag(prec))
    out["partial_corr_abs"] = np.abs(-prec / np.outer(d, d))
    MI = mutual_information(X)
    out["mutual_information"] = MI
    out["clr"] = clr(MI)
    G = ct._cache.get("grnboost_df")
    if G is None:
        ct.get("grnboost")
        G = ct._cache["grnboost_df"]
    gtf = G.reindex(columns=ct.nodes).fillna(0.0).to_numpy(np.float64)       # TF -> TF importances
    out["grnboost_tf_tf"] = np.maximum(gtf, gtf.T)
    rng = np.random.default_rng(0)
    out["ridge_regression"] = per_target_regression(X, "ridge", rng, n_jobs)
    out["lasso_regression"] = per_target_regression(X, "lasso", np.random.default_rng(0), n_jobs)
    out["genie3"] = per_target_regression(X, "genie3", np.random.default_rng(0), n_jobs)
    out["grnboost_cosine"] = cosine(ct.get("grnboost"))
    out["geneformer_cosine"] = cosine(ct.get("geneformer"))
    out["maestro_pearson_abs"] = np.abs(np.nan_to_num(np.corrcoef(ct.get("maestro"))))
    for v in out.values():
        np.fill_diagonal(v, 0)
    return out


# ----------------------------------------------------------------------------- main loop
def run_cell_type(ct: CellType, families, n_jobs, skip_existing, folds=None):
    rdir = ct.root / proto.results_subdir(ct.cfg["results_dir"], ct.split_kind) / ct.ct
    n = ct.n

    def done(name, fold):
        fp = rdir / name / "seed0" / f"fold{fold:02d}.json"
        return skip_existing and fp.exists() and json.load(open(fp)).get("version") == VERSION

    def save(name, fam, fold, m):
        out = rdir / name / "seed0"
        out.mkdir(parents=True, exist_ok=True)
        save_json({"method": name, "family": fam, "version": VERSION, "cell_type": ct.ct, "dataset": ct.cfg["dataset"],
                   "fold": fold, "seed": 0, "note": "fit/count on the fold's training part only",
                   "metrics": m}, out / f"fold{fold:02d}.json")

    unsup = None
    if "A2" in families:
        log(f"{ct.ct}: computing unsupervised scores (A2)")
        unsup = unsupervised_scores(ct, n_jobs)

    for fold in (folds or range(1, ct.n_folds + 1)):
        f = proto.load_fold(ct.sdir, fold, n)
        tr_pos, tr_neg = unordered(f["train"]["pos"]), unordered(f["train"]["neg"])
        held = {tuple(p) for p in unordered(f["test"]["pos"])} | {tuple(p) for p in unordered(f["test"]["neg"])} \
            | {tuple(p) for p in unordered(f["val"]["pos"])} | {tuple(p) for p in unordered(f["val"]["neg"])}
        assert not (held & ({tuple(p) for p in tr_pos} | {tuple(p) for p in tr_neg})), "held-out pair in training"
        A = train_adjacency(n, tr_pos)
        deg = A.sum(1)
        tr_pairs = np.concatenate([tr_pos, tr_neg])
        tr_lab = np.concatenate([np.ones(len(tr_pos)), np.zeros(len(tr_neg))])
        results = {}

        if "A1" in families:
            H = heuristics(A)
            for k in ("degree_product", "degree_sum", "common_neighbours", "jaccard", "adamic_adar", "resource_allocation", "katz",
                      "rwr", "svd_lowrank"):
                results[("A1", f"heur_{k}")] = lambda H=H, k=k: H[k]
            results[("A1", "heur_spectral_lowrank")] = lambda: spectral_lowrank(A)

        if "A2" in families:
            for k, S in unsup.items():
                results[("A2", f"unsup_{k}")] = lambda S=S: S

        if "A3" in families:
            M, Gf, Fm = ct.reduced("maestro"), ct.reduced("grnboost"), ct.reduced("geneformer")
            Hst = ct.get("histogram")
            Op = ct.get("openness")
            blocks = {
                "lr_maestro": lambda i, j: pair_block(M, i, j),
                "lr_grnboost": lambda i, j: pair_block(Gf, i, j),
                "lr_geneformer": lambda i, j: pair_block(Fm, i, j),
                "lr_histogram": lambda i, j: Hst[i, j].astype(np.float64),
                "lr_promoter_openness": lambda i, j: np.concatenate([pair_block(Op, i, j), Op[i] + Op[j]], 1),
                "lr_maestro_grnboost_concat": lambda i, j: np.concatenate([pair_block(M, i, j), pair_block(Gf, i, j)], 1),
                "lr_maestro_x_grnboost": lambda i, j: pair_block(M * Gf, i, j),
                "lr_within_pair_product": lambda i, j: np.concatenate([pair_block(M, i, j, True), pair_block(Gf, i, j, True)], 1),
            }
            for name, fn in blocks.items():
                results[("A3", name)] = lambda fn=fn: fit_predict(LR, fn, n, tr_pairs, tr_lab)
            cat = blocks["lr_maestro_grnboost_concat"]
            results[("A3", "linreg_maestro_grnboost_concat")] = lambda: fit_predict(lambda: Ridge(alpha=1.0), cat, n, tr_pairs, tr_lab)
            results[("A3", "rf_maestro_grnboost_concat")] = lambda: fit_predict(
                lambda: RandomForestClassifier(n_estimators=300, min_samples_leaf=2, n_jobs=n_jobs, random_state=0),
                cat, n, tr_pairs, tr_lab)

        if "A4" in families:
            M, Gf, Fm = ct.reduced("maestro"), ct.reduced("grnboost"), ct.reduced("geneformer")

            def degf(Asub):
                ld = np.log1p(Asub.sum(1))
                return lambda i, j: np.stack([ld[i] + ld[j], np.abs(ld[i] - ld[j]), ld[i] * ld[j]], axis=1)

            def topof(Asub):
                return topo_pair_feats(heuristics(Asub), Asub.sum(1))

            def with_feats(graph_build, blocks):
                def build(Asub):
                    g = graph_build(Asub)
                    return lambda i, j: np.concatenate([g(i, j)] + [pair_block(B, i, j) for B in blocks], 1)
                return build
            combos = {
                "lr_degree_only": degf,
                "lr_degree_maestro_grnboost": with_feats(degf, [M, Gf]),
                "lr_degree_geneformer": with_feats(degf, [Fm]),
                "lr_topology_only": topof,
                "lr_topology_maestro_grnboost": with_feats(topof, [M, Gf]),
                "lr_topology_geneformer": with_feats(topof, [Fm]),
            }
            for name, build in combos.items():
                results[("A4", name)] = lambda build=build: fit_predict(
                    LR, build(A), n, tr_pairs, tr_lab, edge_heldout_feats(build, A, tr_pos, tr_pairs, fold))

        for (fam, name), make in results.items():
            if done(name, fold):
                continue
            t0 = time.time()
            m = evaluate(make(), ct, f, fold)
            save(name, fam, fold, m)
            log(f"{ct.ct} fold {fold} {name}: test AUROC {m['test']['auroc']:.3f} AUPR {m['test']['aupr']:.3f} | "
                f"degree-matched AUROC {m['test_degmatched']['auroc']:.3f} ({time.time() - t0:.1f}s)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--ds", required=True, choices=["lung", "hema"])
    ap.add_argument("--ct", default=None, help="one cell type (default: all of the dataset)")
    ap.add_argument("--family", default="A1,A2,A3,A4")
    ap.add_argument("--n_jobs", type=int, default=4)
    ap.add_argument("--skip_existing", action="store_true")
    ap.add_argument("--folds", default=None, help="comma list (smoke tests); default all")
    ap.add_argument("--split_kind", default="official", choices=proto.SPLIT_KINDS)
    a = ap.parse_args()
    root = Path(a.root).expanduser()
    cfg = load_config(str(REPO / f"configs/scmgrn/{a.ds}.yaml"))
    families = [x.strip() for x in a.family.split(",")]
    cts = [a.ct] if a.ct else list(cfg["cell_types"])
    for c in cts:
        run_cell_type(CellType(root, cfg, c, a.split_kind), families, a.n_jobs, a.skip_existing,
                      [int(x) for x in a.folds.split(",")] if a.folds else None)


if __name__ == "__main__":
    main()
