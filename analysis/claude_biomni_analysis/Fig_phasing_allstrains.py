#!/usr/bin/env python3
"""
Nature-Genetics figure: piRNA phasing across spermatogenesis in 16 mouse strains.
Input : phasing_allstrains_1random/ALL_summary.csv
Method: 1 random coordinate/read (STAR --outSAMmultNmax 1 --outMultimapperOrder Random),
        24-32 nt, GenomicRanges::follow 3'->5' adjacency; +1 nt = phased (Almeida GB2025).
Small-multiples: one panel/strain, x = E16.5 -> P12.5 -> P20.5 (FIXED developmental
order), y = +1 phasing %, bars=mean (orange gradient), dots=replicates. Wild strains
(CAST/PWK/SPRET/WSB) titled in red.
"""
import os, glob, sys
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from strain_order import STRAIN_ORDER, WILD

plt.rcParams.update({"font.family":"Liberation Sans","font.size":7,"axes.linewidth":0.5,
    "axes.spines.top":False,"axes.spines.right":False,"xtick.major.width":0.5,
    "ytick.major.width":0.5,"xtick.major.size":2,"ytick.major.size":2,
    "pdf.fonttype":42,"svg.fonttype":"none"})

BASE="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/phasing_allstrains_1random"
df=pd.read_csv(f"{BASE}/ALL_summary.csv")
def tp(s):
    return "E16.5" if "16.5dpc" in s else "P12.5" if "12.5dpp" in s else "P20.5" if "20.5dpp" in s else "?"
def strain(s):
    for tk in ("16.5dpc","12.5dpp","20.5dpp"):
        if f"-{tk}." in s: return s.split(f"-{tk}.")[0]
    return s
df["tp"]=df["sample"].map(tp); df["strain"]=df["sample"].map(strain); df["pct"]=df["frac_plus1"]*100

TPO=["E16.5","P12.5","P20.5"]
COL={"E16.5":"#F0C9A0","P12.5":"#E69F00","P20.5":"#B4500A"}
# CANONICAL FIGURE ORDER (thesis Fig 4.4, median P20.5 PC1) -- strain_order.py
order=[s for s in STRAIN_ORDER if s in set(df["strain"])]

ncol=4; nrow=int(np.ceil(len(order)/ncol))
fig,axes=plt.subplots(nrow,ncol,figsize=(7.2,7.6),dpi=300)
axes=axes.flatten()
for i,s in enumerate(order):
    ax=axes[i]; d=df[df["strain"]==s]
    for j,t in enumerate(TPO):
        v=d.loc[d.tp==t,"pct"].values
        if len(v)==0: continue
        ax.bar(j,v.mean(),width=0.66,color=COL[t],edgecolor="none",zorder=1)
        ax.scatter(np.full(len(v),j)+(np.arange(len(v))-(len(v)-1)/2)*0.13,v,s=7,color="#222",zorder=3,linewidths=0)
    ax.set_xticks(range(3)); ax.set_xticklabels(TPO,fontsize=5.6)
    ax.set_ylim(0,72); ax.set_yticks([0,25,50])
    tcol="#C0392B" if s in WILD else "#222222"
    nm=s.replace("_","/") if s=="C57BL_6NJ" else s.replace("_","/")
    ax.set_title(("● " if s in WILD else "")+s, fontsize=6.4, color=tcol, fontweight="bold", pad=2)
    if i%ncol==0: ax.set_ylabel("+1 phasing (%)",fontsize=6)
    ax.tick_params(labelsize=5.6)
for k in range(len(order),len(axes)): axes[k].axis("off")

fig.suptitle("piRNA phasing across spermatogenesis in 16 mouse strains",fontsize=10,fontweight="bold",y=1.060)
fig.text(0.5,1.022,"bars = mean +1-nt phasing fraction · dots = replicates · ● = wild-derived strain (CAST/PWK/SPRET/WSB)",
         ha="center",fontsize=6,color="#444")
fig.text(0.5,-0.018,"1 random coordinate/read (STAR --outSAMmultNmax 1 --outMultimapperOrder Random) · 24–32 nt · GenomicRanges::follow 3′→5′ adjacency · timepoints E16.5→P12.5→P20.5",
         ha="center",fontsize=5.6,color="#666")
fig.tight_layout(rect=[0,0,1,0.975],h_pad=1.25,w_pad=0.8)
out=f"{BASE}/Fig_phasing_allstrains"
for ext in ("pdf","svg","png"): fig.savefig(f"{out}.{ext}",bbox_inches="tight")
print("wrote",out+".{pdf,svg,png}")
# summary table
piv=df.pivot_table(index="strain",columns="tp",values="pct",aggfunc="mean").reindex(order)[TPO]
print(piv.round(0).to_string())
