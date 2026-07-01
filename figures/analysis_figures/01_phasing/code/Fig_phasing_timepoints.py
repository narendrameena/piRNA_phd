#!/usr/bin/env python3
"""
Nature-Genetics figure: C57BL/6NJ piRNA phasing across spermatogenesis.

Input : phasing_C57BL_6NJ_1random/ALL_summary.csv  (one row per sample)
Method: 1 random coordinate per read (STAR --outSAMmultNmax 1 --outMultimapperOrder
        Random), 25-32 nt, GenomicRanges::follow() 3'->5' adjacency; +1 nt = phased
        (Almeida et al. Genome Biol 2025, PMID 39844208).

Panels: A = +1 phasing fraction (%) by timepoint (bar=mean, dots=reps);
        B = phasing z-score(+1) by timepoint.
Style : Liberation Sans, Okabe-Ito colourblind-safe, vector + 300 dpi.
"""
import os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 7,
    "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "pdf.fonttype": 42, "svg.fonttype": "none",
})

BASE = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/phasing_C57BL_6NJ_1random"
df = pd.read_csv(f"{BASE}/ALL_summary.csv")

# map sample -> timepoint
def tp(s):
    if "16.5dpc" in s: return "E16.5"
    if "12.5dpp" in s: return "P12.5"
    if "20.5dpp" in s: return "P20.5"
    return "?"
df["timepoint"] = df["sample"].map(tp)
df["pct"] = df["frac_plus1"] * 100.0

ORDER = ["E16.5", "P12.5", "P20.5"]
LABEL = {"E16.5": "E16.5\n(fetal)", "P12.5": "P12.5\n(early postnatal)", "P20.5": "P20.5\n(pachytene)"}
# Okabe-Ito sequential warm gradient (light->dark) for developmental progression
COL = {"E16.5": "#F0C9A0", "P12.5": "#E69F00", "P20.5": "#B4500A"}

fig, (axA, axB) = plt.subplots(1, 2, figsize=(5.0, 2.6), dpi=300)

for ax, col, ylab, ymax in [(axA, "pct", "directly-adjacent piRNA pairs,\n+1 nt (% of all adjacent pairs)", None),
                            (axB, "zscore_plus1", "phasing z-score (+1 nt)", None)]:
    means, sds = [], []
    for i, t in enumerate(ORDER):
        v = df.loc[df.timepoint == t, col].values
        m, sd = v.mean(), v.std(ddof=1)
        means.append(m); sds.append(sd)
        ax.bar(i, m, width=0.62, color=COL[t], edgecolor="none", zorder=1)
        ax.errorbar(i, m, yerr=sd, fmt="none", ecolor="#444444", elinewidth=0.7, capsize=2.5, zorder=2)
        # individual replicate dots (jittered)
        jit = (np.arange(len(v)) - (len(v)-1)/2) * 0.11
        ax.scatter(np.full(len(v), i) + jit, v, s=11, color="#222222",
                   zorder=3, linewidths=0)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([LABEL[t] for t in ORDER], fontsize=6.3)
    ax.set_ylabel(ylab, fontsize=6.8)
    ax.margins(x=0.12)

axA.set_ylim(0, 70)
axA.set_title("A", loc="left", fontsize=9, fontweight="bold")
axB.set_ylim(0, 8)
axB.axhline(2, ls=(0,(3,2)), lw=0.6, color="#888888")  # z=2 significance guide
axB.text(2.4, 2.15, "z = 2", fontsize=5.6, color="#888888", va="bottom", ha="right")
axB.set_title("B", loc="left", fontsize=9, fontweight="bold")

# annotate the gradient on panel A
for i, t in enumerate(ORDER):
    m = df.loc[df.timepoint == t, "pct"].mean()
    axA.text(i, m + 2.5, f"{m:.0f}%", ha="center", va="bottom", fontsize=6.5, fontweight="bold", color=COL[t])

fig.suptitle("C57BL/6NJ piRNA phasing across spermatogenesis",
             fontsize=8.5, fontweight="bold", y=1.02)
fig.text(0.5, -0.06,
         "1 random coordinate/read (STAR --outSAMmultNmax 1 --outMultimapperOrder Random) · 25–32 nt · "
         "GenomicRanges::follow 3′→5′ adjacency · n=3 reps/timepoint",
         ha="center", fontsize=5.6, color="#555555")

fig.tight_layout(w_pad=1.6)
out = f"{BASE}/Fig_phasing_C57BL_6NJ_timepoints"
import os as _os; _SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/01_phasing/data/source_data"; _os.makedirs(_SD,exist_ok=True)
df[["sample","timepoint","frac_plus1","pct","zscore_plus1"]].to_csv(f"{_SD}/SourceData_Fig_phasing_timepoints.csv",index=False)   # per-replicate +1nt phasing (fraction / % / z-score) by timepoint
for ext in ("pdf", "svg", "png"):
    fig.savefig(f"{out}.{ext}", bbox_inches="tight")
print("wrote", out + ".{pdf,svg,png}")
print(df[["sample","timepoint","pct","zscore_plus1"]].to_string(index=False))
