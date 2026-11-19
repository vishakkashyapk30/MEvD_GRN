#!/usr/bin/env python
"""BEAR-GRN: simple, Gaussian and small non-graph models (2026-11-19, advisor items 4-5).

Supervised models follow the SAME protocol as MeVD-GRN and the in-degree baseline (scripts/22, 23): TF-disjoint
K-fold cross-fitting with bear_protocol.tf_disjoint_folds(seed); a model is fitted on the labels of the training TFs
only and scores every (held-out TF, gene) pair, so no label of a scored TF is ever used (BEAR-GRN qualification).

Features of a pair (TF t, gene g), 20 numbers, every one available for every pair:
  indeg        log(1 + in-degree of g among the training TFs); for a TRAINING pair the pair's own link is left out
               (leave-one-out), so the feature never contains the label it is fitted on
  pearson      |Pearson r| of t and g across cells (results/<ds>/baseline_pearson)
  grnboost2    GRNBoost2 importance t -> g (results/<ds>/baseline_grnboost2; 0 if not emitted)
  rna_g(3)     mean, variance, detection rate of g          rna_t(3)  the same for t
  atac_g(3)    mean, variance, detection rate of g's MAESTRO-style regulatory potential
  open_g(4)    promoter-openness descriptors of g
  motif(3)     best relative PWM score of t in g's top peaks, number of hit peaks, has-motif
  gf_cos       cosine of t's and g's Geneformer embeddings
Feature sets: "data" = everything except indeg; "all" = everything.

Models (sklearn; features standardised on the training pairs):
  lr     logistic regression, L2 strength by 3-fold CV                     (~21 parameters)
  lda    linear discriminant analysis: Gaussian classes, shared covariance (Gaussian; ~250 parameters)
  qda    quadratic discriminant analysis: Gaussian classes, own covariance (reg 0.1)
  gnb    Gaussian naive Bayes                                              (~80 parameters)
  gmm    Gaussian-mixture classifier: 4-component full-covariance mixture per class, score = log-likelihood ratio
  rbf    kernel logistic regression with a Gaussian (RBF) kernel, Nystroem approximation (300 components)
  mlp    one hidden layer of 32 units                                      (~700 parameters)
  hgb    histogram gradient boosting (non-neural, non-graph reference)
  deg    no fit: score = the in-degree feature itself (set "all" only); with --label_sources it is the in-degree
         counted in those label files, i.e. the hub rule matching the training labels
Training pairs per fold: every training-TF link (subsampled to <= 60k) + 5 non-links of the SAME TF per link
(pairs in any label or ground-truth file excluded from the negatives), so TF identity carries no label information.

Label-free models (no labels at all) are in scripts/bear/26_bear_label_free.py.
Outputs: results/<ds>/simple_<model>_<set>[_<label sources>]/seed<s>/grn.tsv.gz, same format as scripts/23.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from src.benchmarks import bear_protocol as bp             # noqa: E402
from src.utils.io import load_config, load_json            # noqa: E402

FEATS = ["indeg", "pearson", "grnboost2", "rna_g_mean", "rna_g_var", "rna_g_det", "rna_t_mean", "rna_t_var",
         "rna_t_det", "atac_g_mean", "atac_g_var", "atac_g_det", "open_g1", "open_g2", "open_g3", "open_g4",
         "motif_best", "motif_n", "motif_has", "gf_cos"]
MODELS = ["lr", "lda", "qda", "gnb", "gmm", "rbf", "mlp", "hgb"]


def write(out_dir: Path, names, src, dst, score, fold=None, meta=None):
    out_dir.mkdir(parents=True, exist_ok=True)
    nm = np.asarray(names, dtype=object)
    df = pd.DataFrame({"Source": nm[src], "Target": nm[dst], "Score": score})
    if fold is not None:
        df["fold"] = fold
    df.to_csv(out_dir / "grn.tsv.gz", sep="\t", index=False, float_format="%.4g")
    with open(out_dir / "meta.json", "w") as f:
        json.dump(meta or {}, f, indent=1)
    print(f"[done] {out_dir} ({len(df):,} edges)", flush=True)


def pair_matrix(path: Path, gene_index: dict, tf_pos: dict, n_tf: int, n: int) -> np.ndarray:
    M = np.zeros((n_tf, n), dtype=np.float32)
    if not path.exists():
        return M
    d = pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[""])
    s = d["Source"].map(gene_index)
    t = d["Target"].map(gene_index)
    ok = s.notna() & t.notna()
    r = s[ok].astype(int).map(tf_pos)
    ok2 = r.notna().to_numpy()
    M[r[ok2].astype(int).to_numpy(), t[ok].astype(int).to_numpy()[ok2]] = \
        np.abs(pd.to_numeric(d["Score"][ok], errors="coerce").fillna(0).to_numpy()[ok2])
    return M


def make_model(name: str, seed: int):
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.kernel_approximation import Nystroem
    from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
    from sklearn.naive_bayes import GaussianNB
    from sklearn.neural_network import MLPClassifier
    from sklearn.pipeline import make_pipeline
    if name == "lr":
        return LogisticRegressionCV(Cs=[1e-3, 1e-2, 1e-1, 1, 10], cv=3, scoring="roc_auc", max_iter=2000)
    if name == "lda":
        return LinearDiscriminantAnalysis()
    if name == "qda":
        return QuadraticDiscriminantAnalysis(reg_param=0.1)
    if name == "gnb":
        return GaussianNB()
    if name == "gmm":
        return GMMClassifier(4, seed)
    if name == "rbf":
        return make_pipeline(Nystroem(kernel="rbf", gamma=1.0 / 20, n_components=300, random_state=seed),
                             LogisticRegression(C=1.0, max_iter=2000))
    if name == "mlp":
        return MLPClassifier(hidden_layer_sizes=(32,), early_stopping=True, max_iter=200, random_state=seed)
    if name == "hgb":
        return HistGradientBoostingClassifier(max_iter=200, random_state=seed)
    raise KeyError(name)


class GMMClassifier:
    """Class-conditional Gaussian mixtures; decision score = log p(x | link) - log p(x | non-link)."""

    def __init__(self, k: int, seed: int):
        self.k, self.seed = k, seed

    def fit(self, X, y):
        from sklearn.mixture import GaussianMixture
        self.g1 = GaussianMixture(self.k, covariance_type="full", reg_covar=1e-3, random_state=self.seed).fit(X[y == 1])
        self.g0 = GaussianMixture(self.k, covariance_type="full", reg_covar=1e-3, random_state=self.seed).fit(X[y == 0])
        return self

    def decision_function(self, X):
        return self.g1.score_samples(X) - self.g0.score_samples(X)


def score(model, X) -> np.ndarray:
    if hasattr(model, "predict_proba"):
        try:
            return model.predict_proba(X)[:, 1]
        except Exception:
            pass
    return model.decision_function(X)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--config", default=str(REPO / "configs/bear/datasets.yaml"))
    ap.add_argument("--model_config", default=str(REPO / "configs/bear/mevd_base.yaml"))
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--seeds", nargs="+", type=int, default=[42])
    ap.add_argument("--models", nargs="+", default=MODELS)
    ap.add_argument("--sets", nargs="+", default=["data", "all"])
    ap.add_argument("--max_pos", type=int, default=60000)
    ap.add_argument("--neg_ratio", type=int, default=5)
    ap.add_argument("--label_sources", nargs="*", default=None,
                    help="train on these label files only (e.g. perturbation); folds and negative exclusion still use "
                         "every label file, so the held-out TFs are identical. Output suffix _<sources>")
    args = ap.parse_args()
    root = Path(args.root)
    K = int(load_config(args.model_config)["bear"]["n_folds"])
    pdir = root / "processed" / args.dataset
    res = root / "results" / args.dataset
    gene_index = load_json(pdir / "gene_index.json")
    names = [None] * len(gene_index)
    for g, i in gene_index.items():
        names[i] = g
    n = len(names)
    tf_nodes = np.asarray(load_json(pdir / "tf_indices.json"), dtype=np.int64)
    tf_pos = {int(t): i for i, t in enumerate(tf_nodes)}
    nT = len(tf_nodes)
    lab_files = sorted((pdir / "labels").glob("*.pt"))
    labels = [torch.load(p).numpy().reshape(2, -1) for p in lab_files]
    lab = bp.keys_to_edges(np.unique(np.concatenate([bp.edge_keys(e, n) for e in labels])), n)
    suffix = ""
    if args.label_sources:
        sel = [e for p, e in zip(lab_files, labels) if p.stem in args.label_sources]
        assert sel, f"no label file among {args.label_sources}"
        suffix = "_" + "+".join(sorted(args.label_sources))
    gts = [torch.load(p).numpy().reshape(2, -1) for p in sorted((pdir / "gts").glob("*.pt"))] \
        if (pdir / "gts").exists() else []
    known = np.unique(np.concatenate([bp.edge_keys(e, n) for e in labels + gts]))
    outdeg = {int(k): int(v) for k, v in zip(*np.unique(lab[0], return_counts=True))}   # folds: every label file
    if args.label_sources:
        lab = bp.keys_to_edges(np.unique(np.concatenate([bp.edge_keys(e, n) for e in sel])), n)

    # ---- per-pair matrices (TF row x gene) and per-gene / per-TF vectors
    P = {"pearson": pair_matrix(res / "baseline_pearson/seed0/grn.tsv.gz", gene_index, tf_pos, nT, n),
         "grnboost2": pair_matrix(res / "baseline_grnboost2/seed0/grn.tsv.gz", gene_index, tf_pos, nT, n)}
    rna = np.load(pdir / "rna_features_aligned.npy").astype(np.float32)
    atac = np.load(pdir / "atac_features_aligned.npy").astype(np.float32)
    opn = np.load(pdir / "openness.npy").astype(np.float32)
    mz = np.load(pdir / "motif_pair_feats.npz")
    assert np.array_equal(mz["tf_nodes"], tf_nodes)
    MF = mz["F"].astype(np.float32)                                   # (TF, gene, 3)
    fm = np.load(pdir / "fm_gene_embeddings.npy").astype(np.float32)
    fmn = fm / (np.linalg.norm(fm, axis=1, keepdims=True) + 1e-8)
    GF = fmn[tf_nodes] @ fmn.T                                        # (TF, gene)

    def feats(rows: np.ndarray, genes: np.ndarray, indeg: np.ndarray) -> np.ndarray:
        t = tf_nodes[rows]
        return np.column_stack([
            np.log1p(indeg), P["pearson"][rows, genes], P["grnboost2"][rows, genes],
            rna[genes, :3], rna[t, :3], atac[genes, :3], opn[genes, :4], MF[rows, genes, :3], GF[rows, genes],
        ]).astype(np.float32)

    from sklearn.preprocessing import StandardScaler
    for s in args.seeds:
        folds = bp.tf_disjoint_folds(tf_nodes, outdeg, K, s)
        out = {(m, fs): ([], [], [], []) for m in args.models for fs in args.sets}
        rng = np.random.default_rng(s)
        for k, held in enumerate(folds):
            t0 = time.time()
            train = bp.drop_tfs(lab, held)
            bp.assert_tf_disjoint({"labels": train}, held, what=f"simple fold {k}")
            indeg = np.bincount(train[1], minlength=n).astype(np.float32)
            fit_tfs = np.setdiff1d(tf_nodes, np.asarray(held))
            pos = train
            if pos.shape[1] > args.max_pos:
                pos = pos[:, rng.choice(pos.shape[1], args.max_pos, replace=False)]
            # negatives per TF: every positive (t, g) gets neg_ratio non-links (t, g') of the SAME TF, so TF identity
            # carries no label information and the models must rank genes within a TF (what BEAR-GRN scores; a
            # TF-uniform negative draw lets TF-level features learn "which TF is a hub", which does not transfer
            # to held-out TFs and inverted two folds of Macrophage_S1).
            src = np.repeat(pos[0], 3 * args.neg_ratio)
            cand = np.stack([src, rng.integers(0, n, len(src))])
            ok = (cand[0] != cand[1]) & ~np.isin(bp.edge_keys(cand, n), known)
            cand = cand[:, ok]
            per_tf = {int(t): args.neg_ratio * int(c) for t, c in zip(*np.unique(pos[0], return_counts=True))}
            order = rng.permutation(cand.shape[1])
            cand = cand[:, order]
            rank = pd.Series(cand[0]).groupby(cand[0]).cumcount().to_numpy()
            cand = cand[:, rank < np.array([per_tf[int(t)] for t in cand[0]])]
            tr = np.concatenate([pos, cand], axis=1)
            y = np.r_[np.ones(pos.shape[1], int), np.zeros(cand.shape[1], int)]
            rows_tr = np.array([tf_pos[int(t)] for t in tr[0]])
            Xtr = feats(rows_tr, tr[1], indeg[tr[1]] - y)               # leave-one-out in-degree for links
            pairs = bp.candidate_block(np.asarray(held), np.arange(n))
            rows_te = np.array([tf_pos[int(t)] for t in pairs[0]])
            Xte = feats(rows_te, pairs[1], indeg[pairs[1]])
            for fs in args.sets:
                cols = slice(1, None) if fs == "data" else slice(None)
                sc = StandardScaler().fit(Xtr[:, cols])
                A, B = sc.transform(Xtr[:, cols]), sc.transform(Xte[:, cols])
                for m in args.models:
                    if m == "deg":
                        assert fs == "all", "deg needs the in-degree feature (--sets all)"
                        o = out[(m, fs)]
                        o[0].append(pairs[0]); o[1].append(pairs[1]); o[2].append(Xte[:, 0].astype(np.float32))
                        o[3].append(np.full(pairs.shape[1], k, dtype=np.int16))
                        continue
                    mdl = make_model(m, s)
                    if m in ("qda", "gmm", "rbf", "mlp") and len(y) > 120000:      # cap the slow fits
                        sub = rng.choice(len(y), 120000, replace=False)
                        mdl.fit(A[sub], y[sub])
                    else:
                        mdl.fit(A, y)
                    o = out[(m, fs)]
                    o[0].append(pairs[0]); o[1].append(pairs[1]); o[2].append(score(mdl, B).astype(np.float32))
                    o[3].append(np.full(pairs.shape[1], k, dtype=np.int16))
            print(f"[{args.dataset} seed {s} fold {k}] {pos.shape[1]} pos / {cand.shape[1]} neg, "
                  f"{pairs.shape[1]:,} scored, {time.time() - t0:.0f}s", flush=True)
        for (m, fs), o in out.items():
            write(res / f"simple_{m}_{fs}{suffix}" / f"seed{s}", names, np.concatenate(o[0]), np.concatenate(o[1]),
                  np.concatenate(o[2]), np.concatenate(o[3]),
                  {"model": m, "features": fs, "seed": s, "feature_names": FEATS if fs == "all" else FEATS[1:],
                   "label_sources": args.label_sources or "all", "folds": [f.tolist() for f in folds], "protocol": "TF-disjoint cross-fitting (scripts/22 folds)"})



if __name__ == "__main__":
    main()
