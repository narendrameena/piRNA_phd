#!/usr/bin/env python3
"""Strain-specific piRNA counts (edgeR QL-DA FDR<0.05 & logFC>0, intersected with presence/absence:
present >=2/3 reps in the strain, absent <2/3 in each other strain). Grouped bars, canonical strain
order/colours, value labels. Pilot = C57BL_6NJ, CAST_EiJ, SPRET_EiJ x {E16.5,P12.5,P20.5}.
NOTE: this set still contains SNP-variants of conserved piRNAs; Step 4 (STAR genome-anchored) splits
expressed-elsewhere vs genuinely-novel. These are the strain-specific DA candidates feeding Step 4.
"""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger"
TPMAP={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; TPO=["E16.5","P12.5","P20.5"]
SCOL={"C57BL_6NJ":"#0072B2","CAST_EiJ":"#009E73","SPRET_EiJ":"#D55E00"}

rows=[]
for f,tpk in [("16.5dpc","E16.5"),("12.5dpp","P12.5"),("20.5dpp","P20.5")]:
    d=pd.read_csv(f"{U}/{f}.strain_specific_DA.csv.gz")
    for s,n in d.strain.value_counts().items(): rows.append(dict(timepoint=tpk,strain=s,n=n))
tab=pd.DataFrame(rows).pivot(index="strain",columns="timepoint",values="n").reindex(columns=TPO)
PILOT=[s for s in STRAIN_ORDER if s in tab.index]; tab=tab.reindex(PILOT)
print(tab.to_string()); print("total:",int(tab.values.sum()))

plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(figsize=(6.6,4.0),dpi=300)
x=np.arange(len(TPO)); w=0.26
for i,s in enumerate(PILOT):
    vals=tab.loc[s,TPO].values.astype(float)
    bars=ax.bar(x+(i-1)*w,vals,w,color=SCOL[s],edgecolor="white",linewidth=0.5,label=s.replace("_","/"),zorder=3)
    for b,v in zip(bars,vals):
        ax.text(b.get_x()+b.get_width()/2,v+1500,f"{int(v/1000)}k",ha="center",va="bottom",fontsize=6.2,rotation=0,color=SCOL[s],fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(TPO,fontsize=9)
ax.set_ylabel("strain-specific piRNAs (edgeR DA ∩ presence/absence)",fontsize=8.5)
ax.set_ylim(0,tab.values.max()*1.18)
ax.set_title("Strain-specific unique piRNAs per timepoint (24–32 nt, FDR<0.05)",fontsize=9.2,fontweight="bold")
ax.legend(fontsize=7.5,frameon=False,loc="upper left",ncol=3)
ax.spines[['top','right']].set_visible(False)
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{U}/Fig_strain_specific_DA.{e}",bbox_inches="tight")
tab.to_csv(f"{U}/SourceData_strain_specific_DA.csv")
print("wrote Fig_strain_specific_DA.{png,pdf,svg} + SourceData_strain_specific_DA.csv")
