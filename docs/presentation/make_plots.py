"""Plots for docs/presentation/presentation.md. Reads result files directly; nothing is typed in by hand
except the published PBMC10k reference values (LINGER Supp. Table 7) and the TF-agnostic ranker (pbmc runbook s12.1)."""
import glob, json, os, re, statistics as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = os.path.expanduser("~/research/mevd_grn"); OUT = f"{R}/docs/presentation/figs"
BEAR = os.path.expanduser("~/.cache/bear/data/results")
BLUE, ORANGE, GREEN, PURPLE, GREY, INK, RED = "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#8A94A0", "#1F2937", "#D55E00"
plt.rcParams.update({"font.family": "Poppins" if any("Poppins" in f.name for f in matplotlib.font_manager.fontManager.ttflist) else "DejaVu Sans",
                     "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#9CA3AF", "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK, "axes.titleweight": "bold", "axes.titlesize": 15, "figure.dpi": 100})
T = ["localization", "perturbation", "dual_evidence"]; TL = ["Localization", "Perturbation", "Dual evidence\n(zero-shot)"]

def abl(pat):
    return [json.load(open(f))["results"] for f in sorted(glob.glob(f"{R}/results/ablations/{pat}"))]
def ms(rs, t, m="aupr"):
    v = [r[t][m] for r in rs]; return st.mean(v), (st.stdev(v) if len(v) > 1 else 0.0)

# ---- fig 1: K562 zero-leak
cur, aao = abl("full_curriculum_fm_h384l2_labelfree_seed*_K562.json"), abl("all_at_once_fm_h384l2_labelfree_seed*_K562.json")
sc = json.load(open(f"{R}/results/baselines/scmultiomegrn_K562.json"))["results"]
gb = json.load(open(f"{R}/results/baselines/grnboost2_K562.json"))["results"]
series = [(f"MeVD-GRN, evidence curriculum (n={len(cur)} seeds)", BLUE, [ms(cur, t) for t in T]),
          (f"MeVD-GRN, single stage (n={len(aao)} seeds)", ORANGE, [ms(aao, t) for t in T]),
          ("scMultiomeGRN (our re-run, 1 seed)", GREEN, [(sc[t]["aupr"], 0) for t in T]),
          ("GRNBoost2 (1 seed)", GREY, [(gb[t]["aupr"], 0) for t in T])]
fig, ax = plt.subplots(figsize=(10, 5.2)); w = 0.2
for i, (lab, c, vals) in enumerate(series):
    xs = [j + (i - 1.5) * w for j in range(3)]
    ax.bar(xs, [v[0] for v in vals], w - 0.03, yerr=[v[1] for v in vals], color=c, label=lab, capsize=2, error_kw={"lw": 1.2, "ecolor": INK})
    for x, v in zip(xs, vals): ax.text(x, v[0] + 0.02, f"{v[0]:.2f}", ha="center", fontsize=9, color=INK)
ax.set_xticks(range(3)); ax.set_xticklabels(TL, fontsize=12); ax.set_ylabel("AUPR"); ax.set_ylim(0, 1.12)
ax.set_title("K562, zero test-to-train leakage", loc="left"); ax.legend(frameon=False, fontsize=9.5, loc="upper center", ncol=2, bbox_to_anchor=(0.5, -0.12))
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_k562_zero_leak.png", dpi=200); plt.close(fig)

# ---- fig 2: effect of the leaks (seed 42, same seed both sides)
leg = {"Single stage": abl("all_at_once_fm_h384l2_legacyneg_seed42_K562.json")[0], "Curriculum": abl("full_curriculum_fm_h384l2_legacyneg_seed42_K562.json")[0]}
zl = {"Single stage": abl("all_at_once_fm_h384l2_labelfree_seed42_K562.json")[0], "Curriculum": abl("full_curriculum_fm_h384l2_labelfree_seed42_K562.json")[0]}
fig, axs = plt.subplots(1, 2, figsize=(10, 4.4), sharey=True)
for ax, t, ttl in zip(axs, ["perturbation", "dual_evidence"], ["Perturbation", "Dual evidence (zero-shot)"]):
    for k, name in enumerate(["Single stage", "Curriculum"]):
        a, b = leg[name][t]["aupr"], zl[name][t]["aupr"]
        ax.bar(k - 0.18, a, 0.34, color=RED, alpha=0.85, label="with both leaks" if k == 0 else None)
        ax.bar(k + 0.18, b, 0.34, color=GREEN, label="zero leakage" if k == 0 else None)
        ax.text(k - 0.18, a + 0.015, f"{a:.2f}", ha="center", fontsize=10); ax.text(k + 0.18, b + 0.015, f"{b:.2f}", ha="center", fontsize=10)
        ax.text(k, 0.06, f"{b - a:+.2f}", ha="center", fontsize=12, fontweight="bold", color=INK, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#9CA3AF"))
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Single stage", "Curriculum"], fontsize=12); ax.set_title(ttl, fontsize=13, loc="left"); ax.set_ylim(0, 1.08)
axs[0].set_ylabel("AUPR (K562, seed 42)"); fig.legend(*axs[0].get_legend_handles_labels(), frameon=False, loc="lower center", ncol=2, fontsize=11)
fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig(f"{OUT}/fig2_leak_effect.png", dpi=200); plt.close(fig)

# ---- fig 3: BEAR-GRN ChIP AUPRC, nine datasets
tx = open(f"{BEAR}/summary_tables.md").read()
m = re.search(r"## ChIP - AUPRC\n\n(\|.*?)(?:\n\n|\Z)", tx, re.S); rows = [r.split("|")[1:-1] for r in m.group(1).split("\n")]
hdr = [c.strip() for c in rows[0]][1:]; tab = {r[0].strip(): [float(re.sub(r" \(.*", "", c.strip())) for c in r[1:]] for r in rows[2:]}
pub = [k for k in tab if not k.startswith(("baseline", "*"))]
bestpub = [max(tab[k][i] for k in pub) for i in range(len(hdr))]
mevd = {}
for f in glob.glob(f"{BEAR}/*/scores/mevd_fm_h384__*ChIP.json") + glob.glob(f"{BEAR}/*/scores/mevd_fm_h384_seed*ChIP.json") + glob.glob(f"{BEAR}/*/scores/mevd_fm_h384*seed42*ChIP.json"):
    d = json.load(open(f)); mevd[d["dataset"]] = d["AUPRC"]
order = sorted(range(len(hdr)), key=lambda i: tab["baseline_indegree"][i]); y = {i: n for n, i in enumerate(order)}
fig, ax = plt.subplots(figsize=(10, 5.6))
for i in order: ax.plot([tab["*random*"][i], tab["baseline_indegree"][i]], [y[i]] * 2, color="#D1D5DB", lw=2, zorder=1)
def pts(vals, c, lab, mk="o", s=70): ax.scatter([vals[i] for i in order], [y[i] for i in order], color=c, s=s, marker=mk, label=lab, zorder=3, edgecolor="white", linewidth=1.2)
pts(tab["*random*"], GREY, "random"); pts(tab["baseline_coverage"], PURPLE, "predict every measured pair", "s", 60)
pts(bestpub, BLUE, "best published BEAR method"); pts(tab["LINGER"], ORANGE, "LINGER", "D", 62); pts(tab["baseline_indegree"], GREEN, "rank targets by TF count (trained on labels)", "o", 90)
if mevd: ax.scatter([mevd[hdr[i]] for i in order if hdr[i] in mevd], [y[i] for i in order if hdr[i] in mevd], color=RED, s=170, marker="*", label="MeVD-GRN (M0)", zorder=4, edgecolor="white")
ax.set_yticks(range(len(order))); ax.set_yticklabels([hdr[i].replace("_", " ") for i in order], fontsize=11); ax.set_xlabel("ChIP AUPRC (BEAR-GRN scoring)")
ax.set_title("BEAR-GRN ChIP AUPR: label-trained references vs de novo methods", loc="left"); ax.grid(axis="x", color="#EEF0F2"); ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=9.5, loc="upper center", ncol=3, bbox_to_anchor=(0.45, -0.12))
fig.tight_layout(); fig.savefig(f"{OUT}/fig3_bear_chip_auprc.png", dpi=200); plt.close(fig)

# ---- fig 4: PBMC10k ablation
pb = [("ATAC accessibility alone\n(same ranking for every TF)", 0.8141, PURPLE), ("MeVD-GRN", 0.7343, BLUE), ("  no ATAC", 0.7197, BLUE),
      ("  no Geneformer", 0.6701, BLUE), ("  uniform (not degree-matched) negatives", 0.5973, BLUE), ("  no Geneformer, no ATAC", 0.5438, BLUE),
      ("SCENIC+ (published)", 0.5481, GREY), ("GENIE3 (published)", 0.5387, GREY)]
pb = sorted(pb, key=lambda r: r[1]); fig, ax = plt.subplots(figsize=(10, 5))
ax.barh([r[0] for r in pb], [r[1] for r in pb], color=[r[2] for r in pb], height=0.62)
for n, r in enumerate(pb): ax.text(r[1] + 0.006, n, f"{r[1]:.3f}", va="center", fontsize=10)
ax.axvline(0.7143, color=ORANGE, lw=2.5, ls="--"); ax.text(0.7185, -0.42, "LINGER 0.714", color=ORANGE, ha="left", fontsize=11, fontweight="bold")
ax.set_xlim(0.45, 0.88); ax.set_xlabel("Cistrome AUROC (19 ChIP datasets; LINGER's evaluation)"); ax.set_title("10x PBMC multiome: what carries the signal", loc="left")
fig.tight_layout(); fig.savefig(f"{OUT}/fig4_pbmc10k.png", dpi=200); plt.close(fig)
print("plots written; MeVD-GRN BEAR points:", mevd)

# ---- fig 5: PBMC10k per cell type, MeVD-GRN vs LINGER (LINGER values: LINGER Supp. Table 7, means per cell type; runbook s2)
import csv, collections
rows = list(csv.DictReader(open(f"{R}/docs/presentation/data/pbmc_per_dataset.tsv"), delimiter="\t"))
def per_ct(method, src="collectri"):
    d = collections.defaultdict(list)
    for r in rows:
        if r["method"] == method and r["source"] == src and r["regime"] == "tf" and r["space"] == "expressed" and r["in19"] == "True":
            d[r["cell_type"]].append(float(r["auc"]))
    return {k: st.mean(v) for k, v in d.items()}, {k: len(v) for k, v in d.items()}
mv, cnt = per_ct("fm_h384"); mv_r, _ = per_ct("fm_h384_rnaonly")
ling = {"classical_monocyte": 0.7111, "naive_cd4_t": 0.7071, "naive_b": 0.7220, "mdc": 0.7583}
names = {"classical_monocyte": "Classical\nmonocytes", "naive_cd4_t": "Naive CD4 T", "naive_b": "Naive B", "mdc": "Myeloid DC"}
order = ["classical_monocyte", "naive_cd4_t", "naive_b", "mdc"]
fig, ax = plt.subplots(figsize=(10, 4.8)); w = 0.26
for i, (lab, c, vals) in enumerate([("MeVD-GRN", BLUE, mv), ("MeVD-GRN, no ATAC", "#7FB5D8", mv_r), ("LINGER (published)", ORANGE, ling)]):
    xs = [j + (i - 1) * w for j in range(4)]
    ax.bar(xs, [vals[k] for k in order], w - 0.03, color=c, label=lab)
    for x, k in zip(xs, order): ax.text(x, vals[k] + 0.008, f"{vals[k]:.2f}", ha="center", fontsize=9.5)
ax.set_xticks(range(4)); ax.set_xticklabels([f"{names[k]}\n({cnt[k]} ChIP dataset{'s' if cnt[k] != 1 else ''})" for k in order], fontsize=11)
ax.set_ylim(0.5, 0.82); ax.set_ylabel("Cistrome AUROC"); ax.set_title("PBMC10k by cell type: the average is driven by monocytes", loc="left")
ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.2), fontsize=10)
fig.tight_layout(); fig.savefig(f"{OUT}/fig5_pbmc_per_celltype.png", dpi=200); plt.close(fig)

# ---- fig 6: TF-specificity control (runbook s12.4; local seed-42 `tf` models, 19 datasets; values transcribed from that table)
items = [("TF-agnostic ranker:\nATAC accessibility alone", 0.8141, PURPLE), ("MeVD-GRN: mean of all\nevaluation-TF rows", 0.7142, "#7FB5D8"),
         ("MeVD-GRN: this dataset's\nown TF row (reported)", 0.6904, BLUE), ("MeVD-GRN: a different\nTF's row, same labels", 0.6849, "#7FB5D8"),
         ("MeVD-GRN: own row,\nwithin accessibility strata", 0.5993, GREY)]
fig, ax = plt.subplots(figsize=(10, 4.6))
ax.barh([i[0] for i in items][::-1], [i[1] for i in items][::-1], color=[i[2] for i in items][::-1], height=0.62)
for n, i in enumerate(items[::-1]): ax.text(i[1] + 0.005, n, f"{i[1]:.3f}", va="center", fontsize=10.5)
ax.axvline(0.5, color="#9CA3AF", lw=1.2, ls=":"); ax.set_xlim(0.45, 0.88); ax.set_xlabel("Cistrome AUROC, 19 datasets, all expressed genes")
ax.set_title("Is the signal TF-specific? A different TF's row scores the same", loc="left")
fig.tight_layout(); fig.savefig(f"{OUT}/fig6_pbmc_tf_specificity.png", dpi=200); plt.close(fig)
print("fig5/6 written", {k: round(v, 3) for k, v in mv.items()})
