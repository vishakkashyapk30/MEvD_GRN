"""Apply the pre-registered s12.6 selection rule to the K562 seed-42 dev runs.

Metric: mean over the 5 TF-disjoint folds of the best inner-validation AUPR on the NATURAL
distribution (all validation-TF x gene pairs), from results/K562/mevd_<variant>/seed42/meta.json.
Rule (fixed before any run): candidates are ordered by complexity M0 < M1 < M2 < M3; the headline
is the MOST complex Mk whose mean inner-val AUPR exceeds that of its immediate simpler model
M(k-1) by >= 0.005; if none does, M0. BEAR scores are never used.

  python -m src.benchmarks.bear_select --root ~/.cache/bear/data
"""
import argparse
import json
from pathlib import Path

ORDER = [("M0", "fm_h384"), ("M1", "fm_h384_uniformneg"), ("M2", "fm_h384_hub"), ("M3", "fm_h384_hubmotif")]
MARGIN = 0.005


def val(root: Path, variant: str):
    p = root / "results" / "K562" / f"mevd_{variant}" / "seed42" / "meta.json"
    if not p.exists():
        return None
    m = json.load(open(p))
    if not m.get("complete", False):
        return None
    fl = m["fold_log"]
    return {"mean": sum(f["best_val_aupr"] for f in fl) / len(fl), "folds": [f["best_val_aupr"] for f in fl],
            "peak_mb": max(f.get("peak_gpu_mem_mb", 0) for f in fl), "seconds": m.get("seconds"),
            "zero_leak": all(f.get("leak_status", {}).get("zero_leak") for f in fl)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    a = ap.parse_args()
    root = Path(a.root)
    res = {m: val(root, v) for m, v in ORDER}
    for m, v in ORDER:
        r = res[m]
        print(f"{m} {v:22s}", "not finished" if r is None else
              f"inner-val AUPR mean {r['mean']:.4f}  folds {[round(x, 3) for x in r['folds']]}  "
              f"peak GPU {r['peak_mb']} MiB  {r['seconds']:.0f}s  zero_leak={r['zero_leak']}")
    if any(r is None for r in res.values()):
        print("selection pending: not all four candidates finished")
        return
    chosen = ORDER[0]
    for k in range(1, 4):
        d = res[ORDER[k][0]]["mean"] - res[ORDER[k - 1][0]]["mean"]
        ok = d >= MARGIN
        print(f"{ORDER[k][0]} vs {ORDER[k-1][0]}: delta {d:+.4f} -> {'passes' if ok else 'fails'} (>= {MARGIN})")
        if ok:
            chosen = ORDER[k]
    print(f"HEADLINE = {chosen[0]} ({chosen[1]})")
    (root / "results" / "K562" / "headline.json").write_text(json.dumps(
        {"headline": chosen[0], "variant": chosen[1], "val": {m: res[m]["mean"] for m, _ in ORDER}}, indent=1))


if __name__ == "__main__":
    main()
