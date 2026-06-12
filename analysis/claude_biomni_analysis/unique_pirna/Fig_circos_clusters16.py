#!/usr/bin/env python3
"""Single piRNA-cluster circos: developmental TIMELINE + SENSE/ANTISENSE + EXPRESSION across all 16 strains in
ONE GRCm39 circle. Outer = CONSERVATION ring (viridis = # of 16 strains with a cluster in the 2-Mb bin,
private->core). Then each strain is a labelled track of THREE nested timepoint sub-rings (E16.5->P12.5->P20.5,
read inward); within each sub-ring every cluster bin is a bar whose HEIGHT is proportional to log(uniq-reads
FPM) — taller = higher expression — and whose COLOUR is the cluster STRAND in the GRCm39 frame (sense + blue,
antisense - vermilion, bidirectional purple, CB-safe). So one view shows location, conservation, development,
strand and expression level. Data = cluster_pav/bytp/{X}.{tp}.expr.in_GRCm39.bed (FPM col4, strand col6)."""
import warnings; warnings.filterwarnings("ignore")
import sys,os,math; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
from circos_readkey import zoom_6nj
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.cm as cm, matplotlib.colors as mc
from matplotlib.lines import Line2D
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
PAV=f"{U}/cluster_pav/bytp"; PG=f"{U}/pangenome_te"; FAI=f"{ROOT}/results/ref_genome/GRCm39.106.fasta.fai"
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]; TPS=[("E16.5","16.5dpc"),("P12.5","12.5dpp"),("P20.5","20.5dpp")]
CHROMS=[str(i) for i in range(1,20)]+["X"]; BIN=2_000_000; SCOL={"sense":"#0072B2","antisense":"#D55E00","bidir":"#9467bd"}
CMAP=plt.get_cmap("viridis"); CNORM=mc.Normalize(1,16)
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
def binexpr(bed):
    d=defaultdict(lambda:[0.0,0.0,0.0])
    if os.path.exists(bed) and os.path.getsize(bed):
        for l in open(bed):
            f=l.split("\t")
            if f[0] in CHROMS and len(f)>5:
                try: fpm=float(f[3])
                except: continue
                st=f[5].strip(); mid=(int(f[1])+int(f[2]))//2; b=mid//BIN
                if (f[0],b) in binmap: d[binmap[(f[0],b)]][0 if st=="+" else (1 if st=="-" else 2)]+=fpm
    return d
CL={(X,tp):binexpr(f"{PAV}/{X}.{tp}.expr.in_GRCm39.bed") for X in CANON for _,tp in TPS}
GMAX=max([sum(v) for d in CL.values() for v in d.values()]+[1.0]); LGM=math.log10(GMAX+1)
def cat(v):
    s,a,b=v
    if b>=max(s,a) and b>0: return "bidir"
    if s>0 and a>0 and min(s,a)>=0.3*max(s,a): return "bidir"
    return "sense" if s>=a else "antisense"
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig=plt.figure(figsize=(15.5,15.5),dpi=300); ax=fig.add_subplot(111,projection="polar")
ax.set_theta_direction(-1); ax.set_theta_offset(np.pi/2); ax.set_ylim(0,1.05); ax.axis("off")
for c in CHROMS:
    tt=np.linspace(theta(c,0),theta(c,clen[c]),60); ax.plot(tt,[1.0]*60,color="#222",lw=4,solid_capstyle="butt")
    ax.text((tt[0]+tt[-1])/2,1.04,c,ha="center",va="center",fontsize=8.5,fontweight="bold")
cons=np.zeros(NB)
for X in CANON:
    occ=set()
    for _,tp in TPS: occ|=set(CL[(X,tp)].keys())
    for b in occ: cons[b]+=1
R_OUT,R_IN=0.95,0.17; grp_h=(R_OUT-R_IN)/len(CANON); gap_g=grp_h*0.18; sub_h=(grp_h-gap_g)/3
for k,X in enumerate(CANON):
    top=R_OUT-k*grp_h
    for j,(lab,tp) in enumerate(TPS):
        r=top-gap_g-(j+1)*sub_h; d=CL[(X,tp)]
        if d:
            bb=sorted(d)
            hh=np.array([sub_h*min(1.0,max(0.14,math.log10(sum(d[b])+1)/LGM)) for b in bb])
            cols=[SCOL[cat(d[b])] for b in bb]
            ax.bar(bth[bb],height=hh,width=bwd[bb],bottom=r,color=cols,align="center",edgecolor="none",linewidth=0,rasterized=True)
    mid=top-gap_g-1.5*sub_h
    ax.text(THLAB,mid,X.replace("_","/"),fontsize=6.6,ha="center",va="center",fontweight="bold" if X in WILD else "normal",color="#C0392B" if X in WILD else "#222")
ax.annotate("",xy=(THLAB,R_IN-0.02),xytext=(THLAB,R_OUT+0.02),arrowprops=dict(arrowstyle="-|>",color="#888",lw=1.2))
ax.text(THLAB,R_OUT+0.05,"E16.5→P12.5→P20.5 (inward)",fontsize=6.3,ha="center",va="bottom",color="#888")
leg=[Line2D([0],[0],color=SCOL["sense"],lw=7,label="sense (+ strand)"),Line2D([0],[0],color=SCOL["antisense"],lw=7,label="antisense (− strand)"),Line2D([0],[0],color=SCOL["bidir"],lw=7,label="bidirectional (both strands)")]
fig.legend(handles=leg,loc="lower center",bbox_to_anchor=(0.5,0.05),ncol=3,fontsize=12,frameon=False,title="cluster strand (colour); BAR HEIGHT ∝ log expression (FPM); each strain = 3 nested sub-rings E16.5/P12.5/P20.5 (outer→inner)",title_fontsize=11)
# conservation outer ring removed (uninformative genome-wide aggregate); per-strain conservation lives in picb16
fig.suptitle("piRNA-cluster EXPRESSION circos — SENSE/ANTISENSE × EXPRESSION × TIMELINE, 16 strains in one circle (PICB-combined → GRCm39)\n"
             "each strain = 3 nested timepoint sub-rings (inward); colour = strand, BAR HEIGHT = log expression; strain names (red = wild) at the spoke",
             fontsize=12.5,fontweight="bold",y=0.99,linespacing=1.5)
_rt=R_OUT-gap_g
zoom_6nj(ax, rings=[("E16.5",_rt-0.5*sub_h),("P12.5",_rt-1.5*sub_h),("P20.5",_rt-2.5*sub_h)], theta_c=theta("1",0)-0.006)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_circos_clusters16.{e}",bbox_inches="tight")
print("wrote Fig_circos_clusters16.{png,pdf,svg}")
