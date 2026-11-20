"""Degree-controlled evaluation on BEAR-GRN for every method (2026-11-19).

Question: does any method (released BEAR methods, simple baselines, MeVD-GRN) rank true TF -> gene links above
non-links once the target gene's in-degree (number of regulating TFs) cannot help?

Evaluation frame (same pairs for every method): every (TF node, measured gene) pair of the dataset, restricted per
ground truth to its TFs and targets ("in-space" pairs). Pairs a method does not output score 0 (BEAR convention).
Covariate = target in-degree among the TRAINING TFs of the TF-disjoint fold that holds the scored TF out
(results/<ds>/baseline_indegree/seed<COV_SEED>), so it never uses the scored TF's own labels; the in-degree baseline
of that seed is exactly the covariate.

Tests (all AUROC; ties count 1/2):
  inspace     plain AUROC on the in-space pairs (reference; mixes TFs, so a TF-level score offset matters)
  perTF       mean over TFs of the AUROC within each TF (no matching; removes TF-level offsets only)
  degmatched  each link (t, g) paired with a non-link (t, g') of the SAME TF, drawn among the 3 nearest in
              log(1 + in-degree); the scored TF's own out-degree is matched by construction
  caliper     same TF, |log(1+indeg g') - log(1+indeg g)| <= 0.1, links without a match dropped
  adjusted    covariate-adjusted AUROC (Janes & Pepe, Biometrika 2009) with covariates (scored TF, in-degree stratum):
              links compared only with non-links of the same TF in the same in-degree stratum (0 alone, then 20
              quantile bins), AUROC per cell weighted by its number of links
  adjusted_indeg_only  the same with in-degree strata only (pooled over TFs; a method can still rank whole TFs)
  aiyappa     degree-corrected negatives (Aiyappa et al., ICML 2025): non-link pairs (t', g') with t' drawn as the
              TF of a random link and g' as the target of a random link, so both endpoint-degree distributions match
              the links'; as many as links
Matched sets are drawn once per (dataset, ground truth) with a fixed seed and reused for every method.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from src.benchmarks import bear_metrics as bm          # noqa: E402
from src.benchmarks.bear_hub import strata             # noqa: E402
from src.utils.io import load_config                   # noqa: E402

PAPER_NAME = {"Pando": "Pando-GLM", "Pando_xgb": "Pando-XGB", "DIRECTNET": "DIRECT-NET"}


def _read(p: Path) -> pd.DataFrame:
    return pd.read_csv(p, sep="\t", keep_default_na=False, na_values=[""])


def frame(root: Path, ds: str, cov_seed: int, cov_dir: str = "baseline_indegree") -> pd.DataFrame:
    d = _read(root / "results" / ds / cov_dir / f"seed{cov_seed}" / "grn.tsv.gz")
    f = pd.DataFrame({"tf": d["Source"].astype(str).str.upper(), "target": d["Target"].astype(str).str.upper(),
                      "indeg": pd.to_numeric(d["Score"], errors="coerce").fillna(0.0)})
    return f.groupby(["tf", "target"], as_index=False, sort=True)["indeg"].max()


def scores_on(f: pd.DataFrame, e: pd.DataFrame) -> np.ndarray:
    m = f[["tf", "target"]].merge(e[["tf", "target", "score"]], on=["tf", "target"], how="left")
    return m["score"].to_numpy(dtype=float, na_value=0.0)


def _auc(y, s):
    y = np.asarray(y, bool)
    return float(roc_auc_score(y, s)) if 0 < y.sum() < len(y) else np.nan


def _ap(y, s):
    y = np.asarray(y, bool)
    return float(average_precision_score(y, s)) if 0 < y.sum() < len(y) else np.nan


def build_tests(f: pd.DataFrame, gt: pd.DataFrame, seed: int = 0) -> dict | None:
    """Index sets into f for one ground truth (method-independent)."""
    ins = np.flatnonzero(f["tf"].isin(set(gt["Source"])).to_numpy() & f["target"].isin(set(gt["Target"])).to_numpy())
    if len(ins) == 0:
        return None
    keys = set(zip(gt["Source"], gt["Target"]))
    tf = f["tf"].to_numpy()[ins]
    tg = f["target"].to_numpy()[ins]
    y = np.fromiter(((a, b) in keys for a, b in zip(tf, tg)), bool, len(ins))
    if y.sum() == 0 or (~y).sum() == 0:
        return None
    indeg = f["indeg"].to_numpy()[ins]
    ld = np.log1p(indeg)
    rng = np.random.default_rng(seed)
    dm_p, dm_n, cal_p, cal_n = [], [], [], []
    for _t, ii in pd.Series(np.arange(len(ins))).groupby(tf):
        ii = ii.to_numpy()
        pos, neg = ii[y[ii]], ii[~y[ii]]
        if len(pos) == 0 or len(neg) == 0:
            continue
        o = np.argsort(ld[neg], kind="stable")
        nl = ld[neg][o]
        lp = ld[pos]
        k = np.searchsorted(nl, lp)
        win = np.clip(k[:, None] + np.arange(-3, 3)[None, :], 0, len(nl) - 1)      # 6 neighbours by rank
        dist = np.abs(nl[win] - lp[:, None]) + 1e-9 * np.abs(np.arange(-3, 3))[None, :]
        near = np.take_along_axis(win, np.argsort(dist, axis=1, kind="stable")[:, :3], axis=1)
        pick = near[np.arange(len(pos)), rng.integers(0, near.shape[1], len(pos))]
        dm_p.append(pos)
        dm_n.append(neg[o[pick]])
        a = np.searchsorted(nl, lp - 0.1, side="left")
        b = np.searchsorted(nl, lp + 0.1, side="right")
        ok = b > a
        if ok.any():
            cal_p.append(pos[ok])
            cal_n.append(neg[o[rng.integers(a[ok], b[ok])]])
    # Aiyappa et al.: both endpoints drawn from the links' endpoint distributions
    tf_codes, tf_u = pd.factorize(tf)
    tg_codes, tg_u = pd.factorize(tg)
    key = tf_codes.astype(np.int64) * len(tg_u) + tg_codes
    order = np.argsort(key)
    skey = key[order]
    lp_idx = np.flatnonzero(y)
    need = len(lp_idx)
    got = []
    for _ in range(50):
        if sum(len(g) for g in got) >= need:
            break
        i1 = lp_idx[rng.integers(0, len(lp_idx), 4 * need)]
        i2 = lp_idx[rng.integers(0, len(lp_idx), 4 * need)]
        kk = tf_codes[i1].astype(np.int64) * len(tg_u) + tg_codes[i2]
        loc = np.searchsorted(skey, kk)
        hit = (loc < len(skey)) & (skey[np.minimum(loc, len(skey) - 1)] == kk)
        idx = order[np.minimum(loc, len(skey) - 1)][hit]
        got.append(idx[~y[idx]])
    ay_n = np.concatenate(got)[:need] if got else np.array([], int)
    st = strata(indeg)
    cell = pd.factorize(pd.Series(tf).astype(str) + "|" + pd.Series(st).astype(str))[0]
    return {"ins": ins, "y": y, "strata": st, "cell": cell, "tfc": tf_codes,
            "dm": (np.concatenate(dm_p), np.concatenate(dm_n)) if dm_p else None,
            "cal": (np.concatenate(cal_p), np.concatenate(cal_n)) if cal_p else None,
            "ay": (lp_idx, ay_n) if len(ay_n) else None}


def _pair_auc(s, pn):
    """Pooled AUROC and AP over the matched set, plus the pair-level win rate: the share of links that score above
    their own matched non-link (ties 1/2), which estimates the degree-adjusted AUC directly (Pepe et al. 2013)."""
    if pn is None:
        return np.nan, np.nan, 0, np.nan
    p, n = pn
    yy = np.r_[np.ones(len(p), bool), np.zeros(len(n), bool)]
    ss = np.r_[s[p], s[n]]
    win = float(np.mean((s[p] > s[n]) + 0.5 * (s[p] == s[n]))) if len(p) == len(n) else np.nan
    return _auc(yy, ss), _ap(yy, ss), int(len(p)), win


def evaluate(s_full: np.ndarray, T: dict) -> dict:
    s = s_full[T["ins"]]
    y = T["y"]
    out = {"n_inspace": int(len(y)), "n_links": int(y.sum()),
           "AUROC_inspace": _auc(y, s), "AP_inspace": _ap(y, s)}
    per = []                                    # per-TF AUROC (TFs with >= 10 links and >= 10 non-links), mean
    for _t, ii in pd.Series(np.arange(len(y))).groupby(T["tfc"]):
        ii = ii.to_numpy()
        if y[ii].sum() >= 10 and (~y[ii]).sum() >= 10:
            per.append(roc_auc_score(y[ii], s[ii]))
    out["AUROC_perTF"] = float(np.mean(per)) if per else np.nan
    for name in ("dm", "cal", "ay"):
        a, ap, n, w = _pair_auc(s, T[name])
        tag = {"dm": "degmatched", "cal": "caliper", "ay": "aiyappa"}[name]
        out[f"AUROC_{tag}"], out[f"AP_{tag}"], out[f"n_{tag}"] = a, ap, n
        if name != "ay":
            out[f"pairwin_{tag}"] = w
    for tag, grp in (("adjusted", T["cell"]), ("adjusted_indeg_only", T["strata"])):
        num = den = 0.0
        for _g, ii in pd.Series(np.arange(len(y))).groupby(grp):
            ii = ii.to_numpy()
            npos = int(y[ii].sum())
            if npos and npos < len(ii):
                num += npos * roc_auc_score(y[ii], s[ii])
                den += npos
        out[f"AUROC_{tag}"] = num / den if den else np.nan
    return out


def released_files(root: Path, inf_dir: str) -> dict:
    d = root / "INFERRED.GRNS" / inf_dir
    out = {}
    if d.exists():
        for m in sorted(p.name for p in d.iterdir() if p.is_dir()):
            files = sorted(x for x in (d / m).iterdir() if x.is_file())
            if files:
                out[m] = files[0]
    return out


def method_edges(root: Path, ds: str, inf_dir: str, include: str) -> dict:
    """name -> (kind, seed, loader) for released methods and our results dirs matching `include`."""
    out = {}
    for m, p in released_files(root, inf_dir).items():
        cols = bm.ROC_PR_COLUMNS.get(m)
        if cols is None:
            continue
        out[(PAPER_NAME.get(m, m), "released", 0)] = (p, cols)
    pat = re.compile(include)
    for g in sorted((root / "results" / ds).glob("*/seed*/grn.tsv.gz")):
        name = g.parent.parent.name
        if pat.match(name):
            out[(name, "ours", int(g.parent.name.replace("seed", "")))] = (g, ("Source", "Target", "Score"))
    return out


def load_edges(path: Path, cols) -> pd.DataFrame | None:
    df = bm.read_table(str(path)) if not str(path).endswith(".tsv.gz") else _read(path)
    if cols[2] not in df.columns or cols[0] not in df.columns:
        return None
    return bm.standardise_edges(df, cols[0], cols[1], cols[2])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--datasets", nargs="+", required=True)
    ap.add_argument("--config", default=str(REPO / "configs/bear/datasets.yaml"))
    ap.add_argument("--gts", nargs="*", default=["ChIP", "KO", "Union", "Intersection"])
    ap.add_argument("--cov_seed", type=int, default=42)
    ap.add_argument("--cov_dir", default="baseline_indegree",
                    help="results/<ds>/<cov_dir>/seed<cov_seed> whose Score is the covariate (default: in-degree over "
                         "every label file); e.g. simple_deg_all_perturbation = in-degree over knockout labels only")
    ap.add_argument("--include", default=r"^(baseline_|mevd_|simple_|lf_|nn_|linger_)",
                    help="regex on results/<ds>/<dir> names to score besides the released methods")
    ap.add_argument("--out", default=None, help="CSV (default results/matched/<datasets>.csv)")
    args = ap.parse_args()
    root = Path(args.root)
    dcfg = load_config(args.config)
    rows = []
    for ds in args.datasets:
        c = dcfg["datasets"][ds]
        f = frame(root, ds, args.cov_seed, args.cov_dir)
        tests = {}
        for g in args.gts:
            p = c["gts"].get(g)
            if p and (root / p).exists():
                T = build_tests(f, bm.load_gt(str(root / p)))
                if T is not None:
                    tests[g] = T
        print(f"[{ds}] frame {len(f):,} pairs; GTs {list(tests)}", flush=True)
        for (name, kind, seed), (path, cols) in method_edges(root, ds, dcfg["inferred_dir"][ds], args.include).items():
            e = load_edges(path, cols)
            if e is None:
                continue
            s = scores_on(f, e)
            for g, T in tests.items():
                r = {"dataset": ds, "gt": g, "method": name, "kind": kind, "seed": seed, "covariate": args.cov_dir,
                     **evaluate(s, T)}
                rows.append(r)
            print(f"  {name} seed {seed}: " + ", ".join(
                f"{g} cal {rows[-len(tests) + i]['AUROC_caliper']:.3f} adj {rows[-len(tests) + i]['AUROC_adjusted']:.3f}"
                for i, g in enumerate(tests)), flush=True)
    out = Path(args.out) if args.out else root / "results" / "matched" / f"{'_'.join(args.datasets)}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"wrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
