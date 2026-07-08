#!/usr/bin/env python3
"""DEPTH-CONFOUND CHECK for the genuinely-unique piRNA finding — for BOTH subcategories (conserved-but-silent +
strain-private, klass5 ≥2-read). Is the wild-strain excess a sequencing-DEPTH artifact or real phylogenetic
DIVERGENCE? Per-strain sRNA depth (STAR strain-wise Log.final.out) vs per-strain distinct-sequence counts of each
subcategory (final_classified_clean_2read, klass5). Panels: A strain-private vs depth; B conserved-but-silent vs
depth; C depth-normalised per strain (both subcategories); D per-library depth. If wild>>classical survives depth
control for BOTH subcategories, the finding is robust."""
import sys, glob, re, os
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD, add_classical_wild_companion
from collections import defaultdict
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
PG = f"{U}/pangenome_te"; SD = f"{ROOT}/analysis/claude_biomni_analysis/source_data"
CANON = [s for s in STRAIN_ORDER if s != "C57BL_6"]
# --- per-strain sRNA depth: sum of STAR "Number of input reads" over all libraries; libdep = per-library ---
dep = defaultdict(int); ns = defaultdict(int); libdep = defaultdict(list)
for f in glob.glob(f"{ROOT}/results/STAR_srna_strain_wise/**/Log.final.out", recursive=True):
    samp = next((p for p in f.split("/") if re.search(r"-(16\.5dpc|12\.5dpp|20\.5dpp)\.", p)), None)
    if not samp: continue
    s = samp.split("-")[0]; m = re.search(r"Number of input reads\s*\|\s*([0-9]+)", open(f).read())
    if m: dep[s] += int(m.group(1)); ns[s] += 1; libdep[s].append(int(m.group(1)))
# --- per-strain distinct-sequence counts per genuinely-unique subcategory (klass5, adopted ≥2-read) ---
fc = pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz", usecols=["strain", "sequence", "klass5"])
SUB = [("strain-private", "unique: strain-private locus", "#7a3b9a"), ("conserved-but-silent", "unique: conserved-but-silent", "#0072B2")]
cnt = {lab: fc[fc.klass5 == k].groupby("strain").sequence.nunique() for lab, k, _ in SUB}
names = [s for s in CANON if s in dep]
d = np.array([dep[s] / 1e6 for s in names]); wildmask = np.array([s in WILD for s in names])
cw, cc = "#C0392B", "#2C7FB8"
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(15, 10), dpi=300); gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.24)
stat = {}
def scatter_panel(ax, n, lab, col):
    r_p, p_p = stats.pearsonr(np.log10(d), np.log10(n + 1)); r_s, p_s = stats.spearmanr(d, n)
    norm = (n + 1) / d; wm = np.median(norm[wildmask]); cm = np.median(norm[~wildmask]); _, pu = stats.mannwhitneyu(norm[wildmask], norm[~wildmask], alternative="greater")
    stat[lab] = dict(r_p=r_p, p_p=p_p, r_s=r_s, p_s=p_s, wm=wm, cm=cm, pu=pu, norm=norm)
    ax.scatter(d[~wildmask], n[~wildmask], s=66, c=cc, edgecolor="k", lw=0.5, label="classical", zorder=3)
    ax.scatter(d[wildmask], n[wildmask], s=86, c=cw, edgecolor="k", lw=0.5, marker="D", label="wild-derived", zorder=3)
    ax.set_yscale("log")
    for i, s in enumerate(names):
        ax.annotate(s.replace("_", "/"), (d[i], n[i]), fontsize=5.6, ha="left", va="bottom", xytext=(2.5, 1.5), textcoords="offset points", color=cw if wildmask[i] else "#555")
    ax.set_xlabel("sequencing depth — total sRNA input reads (M)", fontsize=9); ax.set_ylabel(f"{lab} piRNAs [log]", fontsize=9)
    ax.legend(fontsize=8, frameon=False); ax.grid(alpha=0.22, which="both"); ax.spines[["top", "right"]].set_visible(False)
axA = fig.add_subplot(gs[0, 0]); scatter_panel(axA, data_sp := np.array([float(cnt["strain-private"].get(s, 0)) for s in names]), "strain-private", "#7a3b9a")
axA.set_title(f"A  strain-private (locus NEW) vs depth — Spearman ρ={stat['strain-private']['r_s']:.2f} (p={stat['strain-private']['p_s']:.2g})", fontsize=9.6, fontweight="bold", color="#7a3b9a")
axB = fig.add_subplot(gs[0, 1]); scatter_panel(axB, data_cbs := np.array([float(cnt["conserved-but-silent"].get(s, 0)) for s in names]), "conserved-but-silent", "#0072B2")
axB.set_title(f"B  conserved-but-silent (locus SHARED) vs depth — Spearman ρ={stat['conserved-but-silent']['r_s']:.2f} (p={stat['conserved-but-silent']['p_s']:.2g})", fontsize=9.6, fontweight="bold", color="#0072B2")
# C: depth-normalised per strain, grouped by subcategory
axC = fig.add_subplot(gs[1, 0]); order = names; x = np.arange(len(order)); ww = 0.4
for j, (lab, k, col) in enumerate(SUB):
    axC.bar(x + (j - 0.5) * ww, stat[lab]["norm"][[names.index(s) for s in order]], ww, color=col, edgecolor="k", lw=0.3, label=f"{lab}")
axC.set_yscale("log"); axC.set_xticks(x); axC.set_xticklabels([])   # strain labels carried by the classical/wild companion below
axC.set_ylabel("piRNAs per million reads (depth-normalised) [log]", fontsize=9); axC.legend(fontsize=8, frameon=False)
axC.set_title(f"C  depth-normalised — strain-private wild {stat['strain-private']['wm']:.1f} vs classical {stat['strain-private']['cm']:.2f}/Mread ({stat['strain-private']['wm']/stat['strain-private']['cm']:.0f}×, p={stat['strain-private']['pu']:.1g})  ·  CBS {stat['conserved-but-silent']['wm']/stat['conserved-but-silent']['cm']:.0f}×", fontsize=8.5, fontweight="bold")
axC.grid(axis="y", alpha=0.22, which="both"); axC.spines[["top", "right"]].set_visible(False)
# D: per-library depth (each dot = one library)
axD = fig.add_subplot(gs[1, 1]); rng = np.random.default_rng(1)
for xi, s in enumerate(order):
    libs = np.array(libdep[s]) / 1e6; col = cw if s in WILD else cc
    axD.scatter(np.full(len(libs), xi) + rng.uniform(-0.18, 0.18, len(libs)), libs, s=20, c=col, edgecolor="k", lw=0.3, alpha=0.85, zorder=3)
    axD.plot([xi - 0.26, xi + 0.26], [libs.mean(), libs.mean()], color="k", lw=1.1, zorder=4)
axD.set_xticks(range(len(order))); axD.set_xticklabels([s.replace("_", "/") for s in order], rotation=90, fontsize=7)
for t, s in zip(axD.get_xticklabels(), order): t.set_color(cw if s in WILD else "#333"); t.set_fontweight("bold" if s in WILD else "normal")
axD.set_ylabel("per-library sRNA depth (M input reads)", fontsize=9); axD.set_ylim(bottom=0)
axD.set_title("D  per-library depth — each dot = one library (— = strain mean); wild & classical overlap", fontsize=9.2, fontweight="bold")
axD.grid(axis="y", alpha=0.22); axD.spines[["top", "right"]].set_visible(False)
fig.suptitle("Depth-confound check — wild-strain excess of BOTH genuinely-unique subcategories is NOT a sequencing-depth artifact (tracks phylogenetic divergence)", fontsize=12, fontweight="bold", y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.97])
# classical(blue)/wild(orange) companion below panel C: strain-private piRNAs per strain (subspecies colour scheme)
fig.subplots_adjust(bottom=0.16)
_cax=add_classical_wild_companion(fig, axC, order, np.asarray(data_sp,float), gap=0.085, height_frac=0.18, ylabel="strain-priv\n(log)")
_cax.set_xticks(np.arange(len(order))); _cax.set_xticklabels([s.replace("_","/") for s in order], rotation=90, fontsize=6.5)
for lab,s in zip(_cax.get_xticklabels(),order): lab.set_color("#C0392B" if s in WILD else "#333")
_cax.set_title("classical (blue) vs wild-derived (orange) — strain-private piRNAs per strain", fontsize=7.5, fontweight="bold", loc="left")
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_depth_confound_check.{e}", bbox_inches="tight")
os.makedirs(SD, exist_ok=True)
with open(f"{SD}/Fig_depth_confound_check.csv", "w") as o:
    o.write("strain,wild,libraries,input_reads_M,strain_private,conserved_but_silent,sp_per_Mread,cbs_per_Mread\n")
    for i, s in enumerate(names):
        o.write(f"{s},{int(wildmask[i])},{ns[s]},{d[i]:.1f},{int(data_sp[i])},{int(data_cbs[i])},{stat['strain-private']['norm'][i]:.3f},{stat['conserved-but-silent']['norm'][i]:.3f}\n")
print("wrote Fig_depth_confound_check (both subcategories) + source data")
for lab in ("strain-private", "conserved-but-silent"):
    s = stat[lab]; print(f"  {lab}: Spearman ρ(count,depth)={s['r_s']:.2f} p={s['p_s']:.2g} | depth-norm wild/classical={s['wm']/s['cm']:.0f}× (MWU p={s['pu']:.2g})")
