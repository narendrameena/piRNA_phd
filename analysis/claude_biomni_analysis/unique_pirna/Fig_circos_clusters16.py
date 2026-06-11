#!/usr/bin/env python3
"""Creative single-circos: piRNA-cluster developmental TIMELINE + SENSE/ANTISENSE across all 16 strains
(PICB-combined -> GRCm39). ONE GRCm39 circle. Each strain is a labelled track of THREE nested timepoint
sub-rings -- E16.5 -> P12.5 -> P20.5, read INWARD = developmental time. Within each sub-ring, every 2-Mb
cluster bin is coloured by STRAND in the GRCm39 frame: sense (+) blue, antisense (-) vermilion, bidirectional
(both strands at the bin) purple (CB-safe). So you read, per strain: cluster pattern across development (the 3
rings) AND sense/antisense direction (the colour); strain names are aligned to their track at the top spoke.
Data = cluster_pav/bytp/{X}.{tp}.stranded.in_GRCm39.bed (strand col6)."""
import warnings; warnings.filterwarnings("ignore")
import sys,os; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
PAV=f"{U}/cluster_pav/bytp"; PG=f"{U}/pangenome_te"; FAI=f"{ROOT}/results/ref_genome/GRCm39.106.fasta.fai"
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]; TPS=[("E16.5","16.5dpc"),("P12.5","12.5dpp"),("P20.5","20.5dpp")]
CHROMS=[str(i) for i in range(1,20)]+["X"]; BIN=2_000_000
SCOL={"sense":"#0072B2","antisense":"#D55E00","bidir":"#9467bd"}
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
def occ_str(bed):
    d=defaultdict(set)
    if os.path.exists(bed) and os.path.getsize(bed):
        for l in open(bed):
            f=l.split("\t")
            if f[0] in CHROMS:
                st=f[5].strip() if len(f)>5 and f[5].strip() in "+-." else "."
                mid=(int(f[1])+int(f[2]))//2; b=mid//BIN
                if (f[0],b) in binmap: d[binmap[(f[0],b)]].add(st)
    return d
def cat(sset): return "sense" if sset=={"+"} else ("antisense" if sset=={"-"} else "bidir")
CL={(X,tp):occ_str(f"{PAV}/{X}.{tp}.stranded.in_GRCm39.bed") for X in CANON for _,tp in TPS}
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig=plt.figure(figsize=(15.5,15.5),dpi=300); ax=fig.add_subplot(111,projection="polar")
ax.set_theta_direction(-1); ax.set_theta_offset(np.pi/2); ax.set_ylim(0,1.04); ax.axis("off")
for c in CHROMS:
    tt=np.linspace(theta(c,0),theta(c,clen[c]),60); ax.plot(tt,[1.0]*60,color="#222",lw=4,solid_capstyle="butt")
    ax.text((tt[0]+tt[-1])/2,1.035,c,ha="center",va="center",fontsize=8.5,fontweight="bold")
R_OUT,R_IN=0.95,0.18; grp_h=(R_OUT-R_IN)/len(CANON); gap_g=grp_h*0.18; sub_h=(grp_h-gap_g)/3
for k,X in enumerate(CANON):
    top=R_OUT-k*grp_h
    for j,(lab,tp) in enumerate(TPS):
        r=top-gap_g-(j+1)*sub_h; d=CL[(X,tp)]; grp=defaultdict(list)
        for b,sset in d.items(): grp[cat(sset)].append(b)
        for ck,bb in grp.items(): ax.bar(bth[bb],height=sub_h*0.92,width=bwd[bb],bottom=r,color=SCOL[ck],align="center",edgecolor="none",linewidth=0,rasterized=True)
    mid=top-gap_g-1.5*sub_h
    ax.text(THLAB,mid,X.replace("_","/"),fontsize=6.6,ha="center",va="center",fontweight="bold" if X in WILD else "normal",color="#C0392B" if X in WILD else "#222")
ax.annotate("",xy=(THLAB,R_IN-0.02),xytext=(THLAB,R_OUT+0.02),arrowprops=dict(arrowstyle="-|>",color="#888",lw=1.2))
ax.text(THLAB,R_OUT+0.05,"per strain: E16.5→P12.5→P20.5 (inward)",fontsize=6.5,ha="center",va="bottom",color="#888")
leg=[Line2D([0],[0],color=SCOL["sense"],lw=7,label="sense (+ strand)"),Line2D([0],[0],color=SCOL["antisense"],lw=7,label="antisense (− strand)"),Line2D([0],[0],color=SCOL["bidir"],lw=7,label="bidirectional (both strands)")]
fig.legend(handles=leg,loc="lower center",bbox_to_anchor=(0.5,0.05),ncol=3,fontsize=12,frameon=False,title="cluster strand (colour) — each strain's 3 nested sub-rings = E16.5/P12.5/P20.5 (outer→inner)",title_fontsize=11)
fig.suptitle("piRNA-cluster TIMELINE + SENSE/ANTISENSE circos — 16 strains × 3 timepoints (PICB-combined → GRCm39)\n"
             "each strain = 3 nested timepoint sub-rings (inward); colour = cluster strand; strain names (red = wild) at the top spoke; wedge = 2-Mb bin",
             fontsize=13.5,fontweight="bold",y=0.99,linespacing=1.5)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_circos_clusters16.{e}",bbox_inches="tight")
print("wrote Fig_circos_clusters16.{png,pdf,svg}")
