#!/usr/bin/env python3
"""Strain-specific piRNA counts across ALL 16 strains x 3 timepoints (edgeR QL-DA FDR<0.05 & logFC>0,
intersected with presence/absence: present >=2/3 reps in the strain, absent <2/3 in each other strain).
The 16-strain generalisation of Fig_strain_specific_DA.py (pilot = C57/CAST/SPRET). Counts span 2 orders
of magnitude (classical ~10^2-10^3 vs wild ~10^4-10^5) so y is log; x = canonical strain order, wild=red,
bars coloured by developmental timepoint. Source = edger16/{tp}.strain_specific_DA.csv.gz.
NOTE: this set still contains SNP-variants of conserved piRNAs; Fig_step4_classification16 splits
expressed-elsewhere vs genuinely-novel. These are the strain-specific DA candidates feeding Step 4."""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/source_data"
TPO=["E16.5","P12.5","P20.5"]
TPC={"E16.5":"#0072B2","P12.5":"#009E73","P20.5":"#D55E00"}   # Okabe-Ito timepoints (theme-07 convention)
TPLAB={"E16.5":"E16.5 (prepachytene)","P12.5":"P12.5","P20.5":"P20.5 (pachytene)"}

rows=[]
for f,tpk in [("16.5dpc","E16.5"),("12.5dpp","P12.5"),("20.5dpp","P20.5")]:
    d=pd.read_csv(f"{U}/{f}.strain_specific_DA_2read.csv.gz")   # ADOPTED ≥2-read absence (edger16_2read.R)
    for s,n in d.strain.value_counts().items(): rows.append(dict(timepoint=tpk,strain=s,n=int(n)))
tab=(pd.DataFrame(rows).pivot(index="strain",columns="timepoint",values="n")
       .reindex(columns=TPO))
CANON=[s for s in STRAIN_ORDER if s in tab.index]; tab=tab.reindex(CANON).fillna(0)
WPOS=[i for i,s in enumerate(CANON) if s in WILD]
print(tab.astype(int).to_string()); print("total:",int(tab.values.sum()))
# per-replicate detection of the strain-specific set -> +-1 SD error bars (build_replicate_DA_detection.py)
repdet=pd.read_csv(f"{U}/replicate_detected_strain_specific_2read.csv")
REPSD=(repdet.pivot_table(index="strain",columns="timepoint",values="detected",aggfunc="std")
          .reindex(CANON).reindex(columns=TPO))

def klab(v):
    if v>=10000: return f"{v/1000:.0f}k"
    if v>=1000:  return f"{v/1000:.1f}k"
    return f"{int(v)}" if v>0 else ""

plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig,ax=plt.subplots(figsize=(13,4.8),dpi=300)
x=np.arange(len(CANON)); w=0.27
ymax=tab.values.max()
if WPOS: ax.axvspan(min(WPOS)-0.5,max(WPOS)+0.5,color="#C0392B",alpha=0.06,zorder=0)   # wild-derived block
for j,t in enumerate(TPO):
    xs=x+(j-1)*w; vals=tab[t].values.astype(float); errs=np.nan_to_num(REPSD[t].values)
    ax.bar(xs,vals,w,color=TPC[t],edgecolor="white",linewidth=0.3,label=TPLAB[t],zorder=3)
    m=vals>0
    ax.errorbar(xs[m],vals[m],yerr=errs[m],fmt="none",ecolor="#222",elinewidth=0.6,capsize=1.6,capthick=0.6,zorder=5)
    for xi,v,e in zip(xs,vals,errs):
        if v>0: ax.text(xi,(v+e)*1.14,klab(v),ha="center",va="bottom",fontsize=4.6,rotation=90,color=TPC[t])
ax.set_yscale("log"); ax.set_ylim(1,ymax*3.2)
ax.set_xticks(x); labs=ax.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=8)
for lab,s in zip(labs,CANON): lab.set_color("#C0392B" if s in WILD else "#333")
ax.set_ylabel("strain-specific piRNAs\n(edgeR DA ∩ presence/absence, log)",fontsize=8.5)
ax.text(np.mean(WPOS),ymax*2.4,"wild-derived",ha="center",va="top",fontsize=8,fontweight="bold",color="#C0392B")
ax.legend(fontsize=7.5,frameon=False,ncol=3,loc="lower left",bbox_to_anchor=(0,1.005),
          columnspacing=1.4,handlelength=1.3)
ax.set_title("Strain-specific unique piRNAs per timepoint, all 16 strains (24–32 nt, edgeR QL FDR<0.05 ∩ ≥2-read presence/absence — ADOPTED)",
             fontsize=9.6,fontweight="bold",loc="left",pad=28)
ax.spines[['top','right']].set_visible(False)
fig.text(0.5,-0.06,
    "x = 16 strains in canonical order (red = wild-derived) · bars coloured by developmental timepoint · y = log count of strain-specific DA piRNAs "
    "(edgeR QL FDR<0.05 & logFC>0 ∩ present ≥2/3 reps in strain, absent <2 reads in every other strain = ADOPTED ≥2-read absence) · "
    "error bar = ±1 SD across the 3 replicate libraries of the count detected (≥1 read) per replicate (small bars = strain-specific sets are highly replicate-reproducible) · "
    "wild strains carry 10–100× more · candidates still include SNP-variants — split into novel vs expressed-elsewhere in Fig_step4_classification16",
    ha="center",fontsize=5.8,color="#666")
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{U}/Fig_strain_specific_DA16.{e}",bbox_inches="tight")
tab.astype(int).to_csv(f"{SD}/Fig_strain_specific_DA16.csv")
print("wrote Fig_strain_specific_DA16.{png,pdf,svg} + source_data/Fig_strain_specific_DA16.csv")
