#!/usr/bin/env python3
"""PILOT figure: strain- & timepoint-specific unique-piRNA CANDIDATE counts (presence/absence
layer; pre-DA, pre-SNP/locus split; 3-strain pilot). Improved-method analogue of thesis Fig 5.26.
"""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from strain_order import STRAIN_ORDER, WILD, add_classical_wild_companion
plt.rcParams.update({"font.family":"Liberation Sans","font.size":9,"axes.linewidth":0.6,
    "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
d=pd.read_csv(f"{U}/strain_specific_presenceAbsence_candidates.csv.gz")
TPMAP={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; TPO=["E16.5","P12.5","P20.5"]
COL={"E16.5":"#F0C9A0","P12.5":"#E69F00","P20.5":"#B4500A"}
d["tp"]=d.timepoint.map(TPMAP)
order=[s for s in STRAIN_ORDER if s in set(d.strain)]
piv=d.groupby(["strain","tp"]).size().unstack().reindex(order)[TPO]
n=len(order); x=np.arange(n); bw=0.26
fig,ax=plt.subplots(figsize=(5.6,3.8),dpi=300)
for j,t in enumerate(TPO):
    xs=x+(j-1)*bw; v=piv[t].values
    ax.bar(xs,v,width=bw,color=COL[t],edgecolor="white",linewidth=0.4,label=t,zorder=3)
    for xi,val in zip(xs,v): ax.text(xi,val+400,f"{int(val/1000)}k",ha="center",va="bottom",fontsize=6)
ax.set_xticks(x)
labs=ax.set_xticklabels([])   # strain labels carried by the classical/wild companion below
for lab,s in zip(labs,order): lab.set_color("#C0392B" if s in WILD else "#222")
ax.set_ylabel("strain-specific unique piRNA\ncandidates (presence/absence)",fontsize=8.5)
ax.legend(title="timepoint",fontsize=7.5,title_fontsize=8,frameon=False,ncol=3,loc="upper center",bbox_to_anchor=(0.5,1.14))
ax.set_title("PILOT: strain-specific unique piRNAs (3-strain, presence/absence)",fontsize=9,fontweight="bold",pad=22)
fig.text(0.5,-0.06,"CANDIDATES only: ≥1 RPM in ≥2/3 reps of one strain, absent in the others (3-strain comparison) · "
  "pre-edgeR-DA, pre-SNP/locus split · counts will shrink after those",ha="center",fontsize=5.4,color="#666")
fig.tight_layout()
# classical(blue)/wild(orange) total-count companion per strain (subspecies colour scheme)
fig.subplots_adjust(bottom=0.34)
_cax=add_classical_wild_companion(fig,ax,order,piv.sum(1).reindex(order).values,gap=0.13,height_frac=0.20,ylabel="total\ncands")
_cax.set_xticks(np.arange(len(order))); _cax.set_xticklabels(order,rotation=45,ha="right",fontsize=7)
for lab,s in zip(_cax.get_xticklabels(),order): lab.set_color("#C0392B" if s in WILD else "#333")
_cax.set_title("classical (blue) vs wild-derived (orange) — total candidates per strain",fontsize=7.5,fontweight="bold",loc="left")
out=f"{U}/Fig_pilot_unique_pirna_counts"
for ext in ("pdf","svg","png"): fig.savefig(f"{out}.{ext}",bbox_inches="tight")
piv.to_csv(f"{U}/SourceData_pilot_unique_pirna_counts.csv")
print("wrote",out); print(piv.astype(int).to_string())
