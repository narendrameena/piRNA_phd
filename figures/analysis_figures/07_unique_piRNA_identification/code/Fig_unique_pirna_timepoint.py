#!/usr/bin/env python3
"""When are each strain's two kinds of UNIQUE piRNA born? BOTH genuinely-unique klass5 subcategories —
conserved-but-silent (locus shared but SILENT elsewhere = divergence-driven) and strain-private locus
(locus ABSENT in all other strains = insertion-driven) — counted per (strain, timepoint). Timepoints are
shown SIDE BY SIDE (grouped bars), log y (wild strains dominate). One row per subcategory (A = conserved-
but-silent, blue; B = strain-private locus, purple). Right column = pooled developmental fingerprint
(Σ over the 16 strains) for each subcategory. Input unique16/final_classified_clean_2read.csv.gz (klass5, ≥2-read)."""
import sys, os
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD, add_classical_wild_companion
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
SD = f"{ROOT}/analysis/claude_biomni_analysis/source_data"
CANON = [s for s in STRAIN_ORDER if s != "C57BL_6"]; WPOS = [i for i, s in enumerate(CANON) if s in WILD]
TPMAP = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}; TPO = ["E16.5", "P12.5", "P20.5"]
TPCOL = {"E16.5": "#4393C3", "P12.5": "#E8852B", "P20.5": "#B2182B"}   # developmental gradient blue->orange->red
SUB = [("unique: conserved-but-silent", "#0072B2", "Conserved-but-silent unique piRNAs  (locus shared but silent elsewhere — divergence-driven)"),
       ("unique: strain-private locus", "#7a3b9a", "Strain-private-locus unique piRNAs  (locus absent in all 15 other strains — insertion-driven)")]
d = pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz", usecols=["strain", "timepoint", "klass5"])
d["tp"] = d.timepoint.map(TPMAP)
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(15, 10.4), dpi=300)
gs = fig.add_gridspec(2, 2, width_ratios=[6, 1.05], height_ratios=[1, 1], hspace=0.46, wspace=0.11,
                      left=0.065, right=0.985, top=0.875, bottom=0.205)
x = np.arange(len(CANON)); bw = 0.27
src = {}
for r, (kl, accent, desc) in enumerate(SUB):
    sub = d[d.klass5 == kl]
    ct = pd.crosstab(sub.strain, sub.tp).reindex(index=CANON, columns=TPO).fillna(0).astype(int)
    src[kl] = ct
    axM = fig.add_subplot(gs[r, 0]); axS = fig.add_subplot(gs[r, 1])
    if WPOS: axM.axvspan(min(WPOS) - 0.5, max(WPOS) + 0.5, color="#C0392B", alpha=0.06, zorder=0)
    ymax = max(ct.values.max(), 10)
    for j, tp in enumerate(TPO):
        vals = ct[tp].values
        axM.bar(x + (j - 1) * bw, vals, bw, color=TPCOL[tp], edgecolor="white", linewidth=0.3, label=tp, zorder=3)
        for xi, v in zip(x + (j - 1) * bw, vals):
            if v > 0: axM.text(xi, v * 1.16, f"{v//1000}k" if v >= 1000 else f"{v}", ha="center", va="bottom", fontsize=4.3, rotation=90, color=TPCOL[tp], fontweight="bold")
            else: axM.text(xi, 1.16, "0", ha="center", va="bottom", fontsize=5, color=TPCOL[tp], fontweight="bold")
    axM.set_yscale("log"); axM.set_ylim(1, ymax * 3.2)
    axM.set_xticks(x); axM.set_xticklabels([] if r == 1 else [s.replace("_", "/") for s in CANON], rotation=45, ha="right", fontsize=7.5)   # bottom row labels carried by companion
    for t, X in zip(axM.get_xticklabels(), CANON):
        if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
    axM.set_ylabel("# unique piRNAs (log)", fontsize=8.5, color=accent, fontweight="bold")
    axM.spines[["top", "right"]].set_visible(False); axM.spines["left"].set_color(accent); axM.spines["left"].set_linewidth(1.8)
    axM.set_title(f"{'AB'[r]}   {desc}", fontsize=9.6, fontweight="bold", loc="left", color=accent)
    if r == 0:
        axM.legend(ncol=3, fontsize=8.5, frameon=False, loc="upper left", title="developmental window (bars side-by-side)", title_fontsize=8.5)
    # ---- pooled developmental fingerprint (Σ over 16 strains) ----
    pooled = ct.sum(axis=0)
    axS.bar(np.arange(3), [pooled[t] for t in TPO], 0.72, color=[TPCOL[t] for t in TPO], edgecolor="white", linewidth=0.4, zorder=3)
    for xi, t in enumerate(TPO):
        axS.text(xi, pooled[t] + pooled.max() * 0.02, f"{int(pooled[t]):,}", ha="center", va="bottom", fontsize=6.6, color=TPCOL[t], fontweight="bold")
    axS.set_xticks(np.arange(3)); axS.set_xticklabels(TPO, fontsize=7.2, rotation=45, ha="right")
    axS.set_ylim(0, pooled.max() * 1.20); axS.tick_params(labelsize=6.5)
    axS.set_title(f"pooled (Σ16)\nn = {int(pooled.sum()):,}", fontsize=8, fontweight="bold", color=accent)
    axS.spines[["top", "right"]].set_visible(False); axS.spines["left"].set_color(accent); axS.spines["left"].set_linewidth(1.8)
_tot=(src["unique: conserved-but-silent"].sum(1)+src["unique: strain-private locus"].sum(1)).reindex(CANON).values
_cax=add_classical_wild_companion(fig, axM, CANON, _tot, gap=0.085, height_frac=0.16, ylabel="total\nuniq (log)")
_cax.set_xticks(x); _cax.set_xticklabels([s.replace("_","/") for s in CANON], rotation=45, ha="right", fontsize=6.5)
for lab,s in zip(_cax.get_xticklabels(),CANON): lab.set_color("#C0392B" if s in WILD else "#333")
_cax.set_title("classical (blue) vs wild-derived (orange) — total genuinely-unique piRNAs per strain", fontsize=7.5, fontweight="bold", loc="left")
fig.suptitle("When are each strain's two kinds of UNIQUE piRNA born? — conserved-but-silent vs strain-private locus · 16 strains × 3 windows (timepoints side-by-side)",
             fontsize=11.2, fontweight="bold", y=0.96)
fig.text(0.5, 0.006, "Both genuinely-unique klass5 subcategories shown. A / blue = conserved-but-silent (locus present in other strains but silent there — divergence-driven). "
  "B / purple = strain-private locus (locus absent in all 15 other strains — insertion-driven; wild-dominated). Within each strain the three developmental windows are SIDE BY SIDE "
  "(E16.5 · P12.5 · P20.5, left→right). Log y (wild strains dominate). “0” = that strain × window had 0 such candidates. Right column = pooled developmental fingerprint (Σ over 16 strains). klass5, ≥2-read.",
  ha="center", fontsize=6.6, color="#555")
os.makedirs(SD, exist_ok=True)
for kl in src: src[kl].to_csv(f"{SD}/SourceData_unique_pirna_timepoint_{'CBS' if 'conserved' in kl else 'private'}.csv")
for e in ("pdf", "svg", "png"): fig.savefig(f"{U}/pangenome_te/Fig_unique_pirna_timepoint.{e}", bbox_inches="tight")
print("wrote Fig_unique_pirna_timepoint (BOTH unique subcategories; timepoints side-by-side) + source data")
for kl in src: print(f"\n== {kl} ==  (pooled by tp: {dict(src[kl].sum())})"); print(src[kl].to_string())
