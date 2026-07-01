#!/usr/bin/env python3
"""
All 16 strains together: piRNA phasing (+1 nt %) across E16.5 -> P12.5 -> P20.5.
One line per strain, strain order = project order (Strain_Snakefile).
Input: phasing_allstrains_1random/ALL_summary.csv
"""
import sys
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from strain_order import STRAIN_ORDER, WILD

plt.rcParams.update({"font.family":"Liberation Sans","font.size":8,"axes.linewidth":0.6,
    "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})

BASE="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/phasing_allstrains_1random"
df=pd.read_csv(f"{BASE}/ALL_summary.csv")
def tp(s): return "E16.5" if "16.5dpc" in s else "P12.5" if "12.5dpp" in s else "P20.5"
def strain(s):
    for tk in ("16.5dpc","12.5dpp","20.5dpp"):
        if f"-{tk}." in s: return s.split(f"-{tk}.")[0]
    return s
df["tp"]=df["sample"].map(tp); df["strain"]=df["sample"].map(strain); df["pct"]=df["frac_plus1"]*100

# CANONICAL FIGURE ORDER (thesis Fig 4.4, median P20.5 PC1) -- strain_order.py
PROJECT=[s for s in STRAIN_ORDER if s in set(df["strain"])]
TPO=["E16.5","P12.5","P20.5"]; X=np.arange(3)

mean=df.pivot_table(index="strain",columns="tp",values="pct",aggfunc="mean")[TPO]
cmap=plt.get_cmap("tab20")
colors={s:cmap(i % 20) for i,s in enumerate(PROJECT)}

fig,ax=plt.subplots(figsize=(5.6,4.2),dpi=300)
for s in PROJECT:
    y=mean.loc[s,TPO].values
    ls="--" if s in WILD else "-"; mk="D" if s in WILD else "o"
    lw=2.0 if s in WILD else 1.3
    ax.plot(X,y,ls=ls,marker=mk,ms=4.5,lw=lw,color=colors[s],label=s+(" ●" if s in WILD else ""),
            zorder=3 if s in WILD else 2, alpha=0.95)
    # individual rep dots (light)
    for j,t in enumerate(TPO):
        v=df[(df.strain==s)&(df.tp==t)]["pct"].values
        ax.scatter(np.full(len(v),j), v, s=5, color=colors[s], alpha=0.35, zorder=1, linewidths=0)

ax.set_xticks(X); ax.set_xticklabels(["E16.5\n(fetal)","P12.5\n(early postnatal)","P20.5\n(pachytene)"],fontsize=8)
ax.set_ylabel("+1-nt directly-adjacent piRNA pairs (% of all adjacent pairs)",fontsize=8)
ax.set_xlim(-0.15,2.15); ax.set_ylim(10,65)
ax.set_title("piRNA phasing across spermatogenesis — all 16 mouse strains",fontsize=9.5,fontweight="bold",pad=8)
# legend in PROJECT order, outside right, 1 column
leg=ax.legend(fontsize=6.3,frameon=False,loc="center left",bbox_to_anchor=(1.01,0.5),
              ncol=1,handlelength=2.0,title="strain (thesis order)\n● = wild-derived")
leg.get_title().set_fontsize(6.6)
fig.text(0.5,-0.02,"1 random coordinate/read (STAR --outSAMmultNmax 1 --outMultimapperOrder Random) · 25–32 nt · GenomicRanges::follow 3′→5′ adjacency",
         ha="center",fontsize=5.8,color="#666")
fig.tight_layout(rect=[0,0,0.80,1])
out=f"{BASE}/Fig_phasing_allstrains_lines"
import os as _os; _SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/01_phasing/data/source_data"; _os.makedirs(_SD,exist_ok=True)
df.to_csv(f"{_SD}/SourceData_Fig_phasing_allstrains_lines.csv",index=False)
for ext in ("pdf","svg","png"): fig.savefig(f"{out}.{ext}",bbox_inches="tight")
print("wrote",out+".{pdf,svg,png}")
print(mean.reindex(PROJECT).round(0).to_string())
