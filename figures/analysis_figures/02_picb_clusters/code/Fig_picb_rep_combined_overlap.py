#!/usr/bin/env python3
"""
Replicate-vs-combined PICB cluster OVERLAP (do single replicates find the SAME clusters?).
Panel A: composition of combined clusters by how many single replicates also contain them
         (3 / 2 / 1 / 0 reps) -> reproducible core + ~0% unique to combined.
Panel B: single-replicate clusters retained in combined, and replicate-replicate
         reproducibility (per strain×timepoint points + median), by timepoint.
Input: source_data/SourceData_PICB_rep_combined_overlap.csv
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams.update({"font.family":"Liberation Sans","font.size":8,"axes.linewidth":0.6,
    "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})
BASE="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis"
d=pd.read_csv(f"{BASE}/source_data/SourceData_PICB_rep_combined_overlap.csv")
TPMAP={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; TPO=["E16.5","P12.5","P20.5"]
d["tp"]=d.timepoint.map(TPMAP)
# data-integrity: drop groups containing a concurrency-corrupted PICB run (complete=False; e.g. NZO_HlLtJ 12.5dpp rep3 — missing chromosomes)
_cc=pd.read_csv(f"{BASE}/source_data/SourceData_PICB_cluster_counts.csv")
_gc=_cc.groupby(["strain","timepoint"]).complete.all().rename("group_complete").reset_index()
d=d.merge(_gc,on=["strain","timepoint"],how="left"); d=d[d.group_complete.fillna(True)].copy()

fig,(axA,axB)=plt.subplots(1,2,figsize=(7.6,3.6),dpi=300,gridspec_kw={"width_ratios":[1,1.1]})

# Panel A: stacked support composition (aggregate over strains, weighted by cluster number)
agg=d.groupby("tp")[["comb_support3","comb_support2","comb_support1","comb_support0"]].sum().reindex(TPO)
frac=agg.div(agg.sum(axis=1),axis=0)*100
cols={"comb_support3":"#1b7837","comb_support2":"#7fbf7b","comb_support1":"#d9f0d3","comb_support0":"#b2182b"}
labs={"comb_support3":"in all 3 reps","comb_support2":"in 2 reps","comb_support1":"in 1 rep","comb_support0":"unique to combined"}
x=np.arange(3); bottom=np.zeros(3)
for k in ["comb_support3","comb_support2","comb_support1","comb_support0"]:
    axA.bar(x,frac[k].values,bottom=bottom,width=0.62,color=cols[k],edgecolor="white",linewidth=0.4,label=labs[k],zorder=3)
    bottom+=frac[k].values
for xi in x:  # annotate the tiny 'unique' slice
    u=frac["comb_support0"].values[xi]
    axA.text(xi,101,f"{u:.1f}% new",ha="center",va="bottom",fontsize=5.6,color="#b2182b")
axA.set_xticks(x); axA.set_xticklabels(TPO,fontsize=8); axA.set_ylim(0,108)
axA.set_ylabel("% of combined clusters",fontsize=8)
axA.legend(fontsize=6.2,frameon=False,loc="lower center",bbox_to_anchor=(0.5,-0.30),ncol=2,handlelength=1.1)
axA.set_title("A  combined clusters = reproducible replicate clusters",fontsize=8,fontweight="bold",loc="left")

# Panel B: recovery + reproducibility per strain x tp
rng=np.random.default_rng(1)
series=[("frac_rep_in_combined","#0072B2","single-rep clusters kept in combined"),
        ("rep_rep_reproducibility","#E69F00","replicate–replicate reproducibility")]
off={0:-0.17,1:0.17}
for si,(col,c,lab) in enumerate(series):
    for j,t in enumerate(TPO):
        v=d[d.tp==t][col].values*100
        axB.scatter(np.full(len(v),j+off[si])+rng.uniform(-0.05,0.05,len(v)),v,s=12,color=c,
                    alpha=0.8,edgecolors="white",linewidths=0.3,zorder=3,label=lab if j==0 else "")
        axB.plot([j+off[si]-0.1,j+off[si]+0.1],[np.median(v)]*2,color="#222",lw=1.3,zorder=4)
axB.set_xticks(x); axB.set_xticklabels(TPO,fontsize=8); axB.set_ylim(60,100)
axB.set_ylabel("% clusters",fontsize=8)
axB.legend(fontsize=6.2,frameon=False,loc="lower center",bbox_to_anchor=(0.5,-0.30),ncol=1)
axB.set_title("B  a single replicate already recovers ~90%",fontsize=8,fontweight="bold",loc="left")

fig.suptitle("Combining replicates does not change which PICB clusters are found (16 strains × 3 timepoints)",
             fontsize=9,fontweight="bold",y=1.05)
fig.text(0.5,-0.13,">99% of combined clusters are present in ≥1 single replicate · same-strand interval overlap · each point = one strain×timepoint",
         ha="center",fontsize=5.6,color="#666")
fig.tight_layout(rect=[0,0.02,1,0.98])
out=f"{BASE}/Fig_picb_rep_combined_overlap"
for ext in ("pdf","svg","png"): fig.savefig(f"{out}.{ext}",bbox_inches="tight")
print("wrote",out)
# plotted-aggregate source data
frac.assign().to_csv(f"{BASE}/source_data/SourceData_Fig_picb_rep_combined_overlap_supportcomposition.csv")
print(frac.round(1).to_string())
