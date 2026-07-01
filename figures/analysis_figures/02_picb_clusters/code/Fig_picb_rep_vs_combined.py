#!/usr/bin/env python3
"""
Effect of single-replicate vs combined-replicate PICB on cluster FINDING (count).
Panel A: each single replicate's cluster count (x) vs its group's combined count (y),
         log-log, y=x diagonal, coloured by timepoint -> combined ~= single replicate.
Panel B: combined / mean-replicate fold ratio by timepoint (centred on 1.0) -> no inflation.
Input : source_data/SourceData_PICB_cluster_counts.csv
Writes the plotted source data (per-group rep vs combined) in the same script.
"""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from strain_order import WILD

plt.rcParams.update({"font.family":"Liberation Sans","font.size":8,"axes.linewidth":0.6,
    "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})

BASE="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis"
d=pd.read_csv(f"{BASE}/source_data/SourceData_PICB_cluster_counts.csv")
TPMAP={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; TPO=["E16.5","P12.5","P20.5"]
COL={"E16.5":"#F0C9A0","P12.5":"#E69F00","P20.5":"#B4500A"}
d["tp"]=d.timepoint.map(TPMAP)
d=d[d.complete].copy()   # data-integrity: drop concurrency-corrupted PICB runs (complete=False; e.g. NZO_HlLtJ 12.5dpp rep3 — missing chromosomes)
rep=d[d.replicate!="combined"].copy()
comb=d[d.replicate=="combined"][["strain","timepoint","n_clusters"]].rename(columns={"n_clusters":"combined"})
rep=rep.merge(comb,on=["strain","timepoint"])          # each rep row gets its group's combined
grp=(rep.groupby(["strain","timepoint","tp"])
        .agg(rep_mean=("n_clusters","mean"),rep_min=("n_clusters","min"),
             rep_max=("n_clusters","max"),combined=("combined","first")).reset_index())
grp["fold"]=grp["combined"]/grp["rep_mean"]

fig,(axA,axB)=plt.subplots(1,2,figsize=(7.4,3.5),dpi=300,gridspec_kw={"width_ratios":[1.45,1]})

# Panel A — each replicate vs its combined (log-log)
for t in TPO:
    s=rep[rep.tp==t]
    axA.scatter(s["n_clusters"],s["combined"],s=16,color=COL[t],alpha=0.8,
                edgecolors="white",linewidths=0.3,label=t,zorder=3)
lo,hi=rep[["n_clusters","combined"]].min().min()*0.85, rep[["n_clusters","combined"]].max().max()*1.15
axA.plot([lo,hi],[lo,hi],ls="--",lw=0.8,color="#888",zorder=2)
axA.text(hi*0.9,hi*0.55,"y = x",fontsize=6.5,color="#888",rotation=0,ha="right")
axA.set_xscale("log"); axA.set_yscale("log"); axA.set_xlim(lo,hi); axA.set_ylim(lo,hi)
axA.set_xlabel("clusters in a single replicate",fontsize=8)
axA.set_ylabel("clusters in combined run",fontsize=8)
axA.legend(title="timepoint",fontsize=6.8,title_fontsize=7,frameon=False,loc="upper left")
axA.set_title("A  single replicate ≈ combined",fontsize=8.5,fontweight="bold",loc="left")

# Panel B — combined / mean-rep fold by timepoint
rng=np.random.default_rng(0)
for j,t in enumerate(TPO):
    v=grp[grp.tp==t]["fold"].values
    axB.scatter(np.full(len(v),j)+rng.uniform(-0.16,0.16,len(v)),v,s=16,color=COL[t],
                alpha=0.8,edgecolors="white",linewidths=0.3,zorder=3)
    axB.plot([j-0.25,j+0.25],[np.median(v)]*2,color="#222",lw=1.4,zorder=4)
axB.axhline(1.0,ls="--",lw=0.8,color="#888",zorder=1)
axB.set_xticks(range(3)); axB.set_xticklabels(TPO,fontsize=7.5)
axB.set_ylabel("combined / mean-replicate clusters",fontsize=8)
axB.set_ylim(0.9,1.12)
axB.set_title("B  no count gain from pooling",fontsize=8.5,fontweight="bold",loc="left")

fig.suptitle("PICB cluster finding: single replicate vs combined replicates (16 strains × 3 timepoints)",
             fontsize=9.5,fontweight="bold",y=1.04)
fig.text(0.5,-0.04,
   "each point = one strain×timepoint; PICB seed/core thresholds are depth-normalized (FPM/FPKM) → 3× reads does not add clusters. "
   "Combined never exceeds the max single replicate.",ha="center",fontsize=5.6,color="#666")
fig.tight_layout(rect=[0,0,1,0.98])
out=f"{BASE}/Fig_picb_rep_vs_combined"
for ext in ("pdf","svg","png"): fig.savefig(f"{out}.{ext}",bbox_inches="tight")
print("wrote",out)

# source data
sd=grp.copy(); sd["timepoint_label"]=sd["tp"]
sd=sd[["strain","timepoint","timepoint_label","rep_min","rep_mean","rep_max","combined","fold"]]
sd.to_csv(f"{BASE}/source_data/SourceData_Fig_picb_rep_vs_combined.csv",index=False)
print("fold combined/mean-rep by tp:")
print(grp.groupby("tp")["fold"].agg(["median","min","max"]).reindex(TPO).round(3).to_string())
print("combined > max single rep:", int((grp.combined>grp.rep_max).sum()), "/", len(grp))
