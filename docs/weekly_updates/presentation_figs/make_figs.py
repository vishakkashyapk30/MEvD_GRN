"""Charts for docs/weekly_updates/presentation.md. Reads result JSONs and BEAR score tables; no number is typed in
except the PBMC10k values, which come from docs/experiments/pbmc10k_linger_benchmark.md (sections 11-12)."""
import glob, json, os, re, statistics as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
OUT = os.path.dirname(__file__)
BLUE, ORANGE, GREEN, PURPLE, GREY, RED, INK = "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#9AA3AE", "#D55E00", "#1F2937"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#6B7280", "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "figure.dpi": 150})
T = ["localization", "perturbation", "dual_evidence"]
TL = ["Localization", "Perturbation", "Dual evidence\n(zero-shot)"]

def runs(pat):
    return [json.load(open(f))["results"] for f in sorted(glob.glob(f"{REPO}/results/ablations/{pat}"))]
def ms(rs, t, m="aupr"):
    v = [r[t][m] for r in rs]; return st.mean(v), (st.stdev(v) if len(v) > 1 else 0.0)

cur = runs("full_curriculum_fm_h384l2_labelfree_seed*_K562.json")
aao = runs("all_at_once_fm_h384l2_labelfree_seed*_K562.json")
sc = json.load(open(f"{REPO}/results/baselines/scmultiomegrn_K562.json"))["results"]
gb = json.load(open(f"{REPO}/results/baselines/grnboost2_K562.json"))["results"]
# leaky (pre-fix) references
aao_leaky = runs("all_at_once_fm_h384l2_seed4[2-6]_K562.json")
cur_leaky = runs("full_curriculum_fm_h384l2_legacyneg_seed42_K562.json")

# ---- fig 1: zero-leak K562, 3 tiers
fig, ax = plt.subplots(figsize=(9.5, 5))
models = [("MeVD-GRN, evidence curriculum", cur, GREEN), ("MeVD-GRN, single stage", aao, BLUE)]
w = 0.2; x = np.arange(3)
for i, (n, rs, c) in enumerate(models):
    mu = [ms(rs, t)[0] for t in T]; sd = [ms(rs, t)[1] for t in T]
    ax.bar(x + (i - 1.5) * w, mu, w, yerr=sd, color=c, label=f"{n} ({len(rs)} seeds)", capsize=3)
ax.bar(x + 0.5 * w, [sc[t]["aupr"] for t in T], w, color=ORANGE, label="scMultiomeGRN (our re-run, 1 seed)")
ax.bar(x + 1.5 * w, [gb[t]["aupr"] for t in T], w, color=GREY, label="GRNBoost2 (1 seed)")
ax.set_xticks(x); ax.set_xticklabels(TL); ax.set_ylabel("AUPR"); ax.set_ylim(0, 1.05)
ax.set_title("K562, zero leakage: curriculum wins the harder tiers", loc="left", fontweight="bold")
ax.legend(frameon=False, fontsize=10, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_k562_zero_leak.png", bbox_inches="tight"); plt.close(fig)

# ---- fig 2: the leak
fig, ax = plt.subplots(figsize=(9.5, 5.4))
tiers = ["perturbation", "dual_evidence"]; lab = ["Perturbation", "Dual evidence (zero-shot)"]
x = np.arange(2); w = 0.2
sets = [("single stage, WITH leak (old)", aao_leaky, BLUE, "//"), ("single stage, zero leak", aao, BLUE, None),
        ("curriculum, WITH leak (old, seed 42)", cur_leaky, GREEN, "//"), ("curriculum, zero leak", cur, GREEN, None)]
for i, (n, rs, c, h) in enumerate(sets):
    mu = [ms(rs, t)[0] for t in tiers]
    b = ax.bar(x + (i - 1.5) * w, mu, w, color=c if h is None else "white", edgecolor=c, hatch=h, linewidth=1.8, label=n)
    for xx, m in zip(x + (i - 1.5) * w, mu): ax.text(xx, m + 0.012, f"{m:.2f}", ha="center", fontsize=9)
ax.set_xticks(x); ax.set_xticklabels(lab); ax.set_ylabel("AUPR"); ax.set_ylim(0, 1.08)
ax.set_title("The leak inflated the single-stage model far more than the curriculum", loc="left", fontweight="bold")
ax.legend(frameon=False, fontsize=10, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
fig.tight_layout(); fig.savefig(f"{OUT}/fig2_leak_effect.png", bbox_inches="tight"); plt.close(fig)

# ---- fig 3: PBMC10k
rows = [("ATAC gene activity alone\n(same ranking for every TF)", 0.8141, PURPLE),
        ("MeVD-GRN (Geneformer + ATAC)", 0.7343, GREEN), ("MeVD-GRN, RNA only", 0.7197, GREEN),
        ("LINGER (published)", 0.7143, ORANGE), ("MeVD-GRN, no Geneformer", 0.6701, GREEN),
        ("MeVD-GRN, uniform negatives", 0.5973, GREEN), ("Pearson |r|", 0.5839, GREY), ("Target count (degree)", 0.5791, GREY),
        ("Gene ID only", 0.5782, GREY), ("SCENIC+ (published)", 0.5481, ORANGE), ("MeVD-GRN, no Geneformer, RNA only", 0.5438, GREEN),
        ("GENIE3 (published)", 0.5387, ORANGE)]
fig, ax = plt.subplots(figsize=(9.5, 5.6))
yy = np.arange(len(rows))[::-1]
ax.barh(yy, [r[1] - 0.5 for r in rows], left=0.5, color=[r[2] for r in rows])
for y, r in zip(yy, rows): ax.text(r[1] + 0.006, y, f"{r[1]:.3f}", va="center", fontsize=10)
ax.set_yticks(yy); ax.set_yticklabels([r[0] for r in rows], fontsize=10); ax.set_xlim(0.5, 0.88)
ax.set_xlabel("Cistrome AUROC (mean of 19 ChIP datasets) (from chance, 0.5)")
ax.set_title("PBMC10k: we edge LINGER, but an ATAC-only ranker\n(same ranking for every TF) beats both", loc="left", fontweight="bold")
fig.tight_layout(); fig.savefig(f"{OUT}/fig3_pbmc10k.png"); plt.close(fig)

# ---- fig 4: BEAR-GRN ChIP AUPRC, 9 datasets
md = open(os.path.expanduser("~/.cache/bear/data/results/summary_tables.md")).read()
m = re.search(r"## ChIP - AUPRC\n\n(\|.*?)(?:\n\n|\Z)", md, re.S)
rws = [r.split("|")[1:-1] for r in m.group(1).split("\n")]
hdr = [c.strip() for c in rws[0]][1:]
tab = {r[0].strip(): {h: float(re.sub(r" \(.*", "", c.strip())) for h, c in zip(hdr, r[1:])} for r in rws[2:]}
pub = [k for k in tab if k not in ("*random*", "baseline_coverage", "baseline_grnboost2", "baseline_indegree", "baseline_pearson")]
x = np.arange(len(hdr)); w = 0.16
best = [max(tab[p][h] for p in pub) for h in hdr]
series = [("LINGER", [tab["LINGER"][h] for h in hdr], ORANGE), ("best published method", best, "#C98A00"),
          ("rank by target count (in-degree)", [tab["baseline_indegree"][h] for h in hdr], RED),
          ("predict every measured pair", [tab["baseline_coverage"][h] for h in hdr], GREY),
          ("random", [tab["*random*"][h] for h in hdr], "#D1D5DB")]
fig, ax = plt.subplots(figsize=(11.5, 5.2))
for i, (n, v, c) in enumerate(series): ax.bar(x + (i - 2) * w, v, w, color=c, label=n)
ax.set_xticks(x); ax.set_xticklabels([h.replace("mESC_", "mouse ").replace("_", " ") for h in hdr], rotation=25, ha="right", fontsize=10)
ax.set_ylabel("ChIP AUPRC (BEAR-GRN scoring)")
ax.set_title("BEAR-GRN: ranking targets by how many TFs regulate them beats every published method", loc="left", fontweight="bold", fontsize=12)
ax.legend(frameon=False, fontsize=9.5, ncol=3, loc="upper right"); ax.set_ylim(0, 0.62)
fig.tight_layout(); fig.savefig(f"{OUT}/fig4_bear_baselines.png"); plt.close(fig)

# ---- fig 5: K562 ChIP, MeVD-GRN M0 vs baselines (AUPRC / AUROC)
sc_path = os.path.expanduser("~/.cache/bear/data/results/K562/scores/mevd_fm_h384__seed42__ChIP.json")
if os.path.exists(sc_path):
    m0 = json.load(open(sc_path))
    names = ["random", "LINGER", "MeVD-GRN M0\n(dev set, seed 42)", "rank by\ntarget count"]
    pr = [tab["*random*"]["K562"], tab["LINGER"]["K562"], m0["AUPRC"], tab["baseline_indegree"]["K562"]]
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.bar(names, pr, color=[GREY, ORANGE, GREEN, RED])
    for i, v in enumerate(pr): ax.text(i, v + 0.006, f"{v:.3f}", ha="center")
    ax.set_ylim(0, 0.56); ax.set_ylabel("K562 ChIP AUPRC")
    ax.set_title("First MeVD-GRN number under BEAR's protocol (M0)", loc="left", fontweight="bold", fontsize=12)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig5_bear_k562_m0.png"); plt.close(fig)
print("figures:", sorted(f for f in os.listdir(OUT) if f.endswith(".png")))
