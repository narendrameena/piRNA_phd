#!/usr/bin/env python3
"""Circos #3 — STRAIN-PRIVATE / genuinely-unique piRNA density across all 16 strains x 3 timepoints, ONE GRCm39
circle. Shows WHERE each strain's genuinely-unique (novel-sequence: conserved-but-silent + strain-private-locus)
piRNAs map. Each strain = 3 nested timepoint sub-rings (E16.5->P12.5->P20.5, inward); each 2-Mb bin: bar HEIGHT
∝ log(# unique piRNAs in the bin), colour = same (magma). Data = unique16/loci/{X}.cand_GRCm39.bed (novel
candidates halLifted to GRCm39; name = strain|tp|seq, deduped per bin). Strain names (red = wild) at the spoke."""
import warnings; warnings.filterwarnings("ignore")
import sys,os,math; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
from circos_readkey import zoom_6nj
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.cm as cm, matplotlib.colors as mc
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG=f"{U}/pangenome_te"
FAI=f"{ROOT}/results/ref_genome/GRCm39.106.fasta.fai"
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]; TPS=[("E16.5","16.5dpc"),("P12.5","12.5dpp"),("P20.5","20.5dpp")]
CHROMS=[str(i) for i in range(1,20)]+["X"]; BIN=2_000_000
clen={}
for ln in open(FAI):
    f=ln.split("\t")
    if f[0] in CHROMS: clen[f[0]]=int(f[1])
gap=sum(clen.values())*0.006; LABELGAP=sum(clen.values())*0.085; off={}; cum=LABELGAP
for c in CHROMS: off[c]=cum; cum+=clen[c]+gap
TOT=cum
def theta(c,p): return 2*np.pi*(off[c]+p)/TOT
THLAB=2*np.pi*(LABELGAP*0.5)/TOT
bins=[]; binmap={}
for c in CHROMS:
    for b in range(int(np.ceil(clen[c]/BIN))):
        binmap[(c,b)]=len(bins); bins.append((c,b*BIN,min((b+1)*BIN,clen[c])))
NB=len(bins); bth=np.array([(theta(c,s)+theta(c,e))/2 for c,s,e in bins]); bwd=np.array([theta(c,e)-theta(c,s) for c,s,e in bins])
def density(X):
    seen=defaultdict(set); bed=f"{U}/unique16/loci/{X}.cand_GRCm39.bed"
    if os.path.exists(bed):
        for l in open(bed):
            f=l.split("\t")
            if f[0] in CHROMS and len(f)>3 and "|" in f[3]:
                tp=f[3].split("|")[1]; mid=(int(f[1])+int(f[2]))//2; b=mid//BIN
                if (f[0],b) in binmap: seen[(tp,binmap[(f[0],b)])].add(f[3])
    d=defaultdict(dict)
    for (tp,bi),nm in seen.items(): d[tp][bi]=len(nm)
    return d
DEN={X:density(X) for X in CANON}
allv=[v for X in CANON for tp in DEN[X] for v in DEN[X][tp].values()]
VMAX=np.percentile(allv,99) if allv else 1.0; LV=math.log10(VMAX+1)
PRIVCOL="#88419D"   # flat identity colour; bar HEIGHT carries density, so colour no longer duplicates height
print(f"strain-private density: VMAX(99pct)={VMAX:.0f} piRNAs/2Mb bin")
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig=plt.figure(figsize=(15.5,15.5),dpi=300); ax=fig.add_subplot(111,projection="polar")
ax.set_theta_direction(-1); ax.set_theta_offset(np.pi/2); ax.set_ylim(0,1.05); ax.axis("off")
for c in CHROMS:
    tt=np.linspace(theta(c,0),theta(c,clen[c]),60); ax.plot(tt,[1.0]*60,color="#222",lw=4,solid_capstyle="butt")
    ax.text((tt[0]+tt[-1])/2,1.04,c,ha="center",va="center",fontsize=8.5,fontweight="bold")
R_OUT,R_IN=0.95,0.18; grp_h=(R_OUT-R_IN)/len(CANON); gap_g=grp_h*0.18; sub_h=(grp_h-gap_g)/3
for k,X in enumerate(CANON):
    top=R_OUT-k*grp_h
    for j,(lab,tp) in enumerate(TPS):
        r=top-gap_g-(j+1)*sub_h; d=DEN[X].get(tp,{})
        if d:
            bb=sorted(d); hh=np.array([sub_h*min(1.0,max(0.12,math.log10(d[b]+1)/LV)) for b in bb])
            ax.bar(bth[bb],height=hh,width=bwd[bb],bottom=r,color=PRIVCOL,align="center",edgecolor="none",linewidth=0,rasterized=True)
    mid=top-gap_g-1.5*sub_h
    ax.text(THLAB,mid,X.replace("_","/"),fontsize=6.6,ha="center",va="center",fontweight="bold" if X in WILD else "normal",color="#C0392B" if X in WILD else "#222")
ax.annotate("",xy=(THLAB,R_IN-0.02),xytext=(THLAB,R_OUT+0.02),arrowprops=dict(arrowstyle="-|>",color="#888",lw=1.2))
ax.text(THLAB,R_OUT+0.05,"E16.5→P12.5→P20.5 (inward)",fontsize=6.3,ha="center",va="bottom",color="#888")
fig.text(0.5,0.072,"BAR HEIGHT ∝ log # genuinely-unique piRNAs per 2-Mb bin (colour = identity only, no longer duplicates height)",ha="center",fontsize=11,color="#333")
fig.suptitle("STRAIN-PRIVATE / genuinely-unique piRNA density circos — 16 strains × 3 timepoints in one GRCm39 circle\n"
             "each strain = 3 nested timepoint sub-rings (inward); BAR HEIGHT ∝ log # unique piRNAs per 2-Mb bin; strain names (red = wild) at the spoke",
             fontsize=12.5,fontweight="bold",y=0.99,linespacing=1.5)
_rt=R_OUT-gap_g
zoom_6nj(ax, rings=[("E16.5",_rt-0.5*sub_h),("P12.5",_rt-1.5*sub_h),("P20.5",_rt-2.5*sub_h)], theta_c=theta("1",0)-0.006)
import os as _os, pandas as _pd; _SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/14_circos_pangenome_TE/data/source_data"; _os.makedirs(_SD,exist_ok=True)
_pd.DataFrame([(X,tp,bins[bi][0],bins[bi][1],round(v,4)) for X,_td in DEN.items() for tp,_bd in _td.items() for bi,v in _bd.items()],columns=["strain","timepoint","chrom","bin_start","private_density"]).to_csv(f"{_SD}/SourceData_Fig_circos_private16.csv",index=False)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_circos_private16.{e}",bbox_inches="tight")
print("wrote Fig_circos_private16.{png,pdf,svg}")
