#!/usr/bin/env python3
"""Zamore stage (pre-pachytene/hybrid/pachytene) vs cross-strain conservation in our 16-strain PICB
cluster PAV. Each of the 214 Zamore loci -> coverage-weighted mean #strains carrying a cluster there.
Pachytene clusters are more conserved (nearly core) than pre-pachytene/hybrid (Ozata et al. 2020)."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"
z=pd.read_csv(f"{U}/zamore_loci_stage_conservation.csv").dropna(subset=["mean_nstrains"])
ORD=["prepachytene","hybrid","pachytene"]; LAB=["pre-pachytene","hybrid","pachytene"]
COL={"prepachytene":"#0072B2","hybrid":"#009E73","pachytene":"#D55E00"}
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(figsize=(5.6,4.2),dpi=300)
data=[z.loc[z.stage==s,"mean_nstrains"].values for s in ORD]
bp=ax.boxplot(data,patch_artist=True,widths=0.6,showfliers=False,medianprops=dict(color="black",lw=1.2))
for patch,s in zip(bp["boxes"],ORD): patch.set_facecolor(COL[s]); patch.set_alpha(0.75)
# jitter points
for i,s in enumerate(ORD):
    y=z.loc[z.stage==s,"mean_nstrains"].values; x=np.random.default_rng(0).normal(i+1,0.06,len(y))
    ax.scatter(x,y,s=8,color=COL[s],edgecolor="white",linewidth=0.2,zorder=3,alpha=0.6)
for i,s in enumerate(ORD):
    n=len(z[z.stage==s]); m=z.loc[z.stage==s,"mean_nstrains"].mean()
    ax.text(i+1,16.4,f"n={n}\nmean {m:.1f}",ha="center",va="bottom",fontsize=7,color=COL[s],fontweight="bold")
ax.set_xticks([1,2,3]); ax.set_xticklabels(LAB,fontsize=9)
ax.set_ylabel("# strains carrying the cluster (/16)",fontsize=9); ax.set_ylim(0,18)
ax.axhline(16,ls=":",lw=0.6,color="#888"); ax.text(0.6,16,"core (16)",fontsize=6,color="#888",va="bottom")
ax.set_title("Pachytene piRNA clusters are more cross-strain conserved\nthan pre-pachytene (Zamore stages × 16-strain PAV)",fontsize=9,fontweight="bold")
ax.spines[['top','right']].set_visible(False)
fig.tight_layout()
import os as _os; _SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/10_cluster_pangenome_PAV/data/source_data"; _os.makedirs(_SD,exist_ok=True)
z.to_csv(f"{_SD}/SourceData_Fig_zamore_stage_conservation.csv",index=False)   # boxplot+jitter data: per-locus Zamore stage + coverage-weighted mean #strains carrying
for e in ("pdf","svg","png"): fig.savefig(f"{U}/Fig_zamore_stage_conservation.{e}",bbox_inches="tight")
print("wrote Fig_zamore_stage_conservation.{png,pdf,svg}")
print(z.groupby("stage").mean_nstrains.agg(["size","mean","median"]).round(2).to_string())
