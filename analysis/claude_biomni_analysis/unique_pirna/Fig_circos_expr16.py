#!/usr/bin/env python3
"""Information-rich circos: piRNA-cluster EXPRESSION + DIRECTIONALITY across 16 strains x 3 timepoints
(PICB-combined -> GRCm39). Three panels (E16.5/P12.5/P20.5). Each strain = a baseline ring; at every 2-Mb bin a
DIVERGENT bar shows expression by strand: SENSE (+) grows OUTWARD (blue), ANTISENSE (-) grows INWARD (vermilion),
bar length proportional to log(uniq-reads FPM). BIDIRECTIONAL clusters (PICB '*') split out+in -> they show bars
BOTH ways (the divergent pachytene-cluster signature); mono-directional clusters show one side only. So you read,
per strain and timepoint, where clusters are, how strongly they fire, and whether they are sense / antisense /
bidirectional. CB-safe (Okabe-Ito). Data = cluster_pav/bytp/{X}.{tp}.expr.in_GRCm39.bed (FPM col4, strand col6)."""
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
CHROMS=[str(i) for i in range(1,20)]+["X"]; BIN=2_000_000; SENSE="#0072B2"; ANTI="#D55E00"
clen={}
for ln in open(FAI):
    f=ln.split("\t")
    if f[0] in CHROMS: clen[f[0]]=int(f[1])
gap=sum(clen.values())*0.006; LABELGAP=sum(clen.values())*0.075; off={}; cum=LABELGAP
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
    d=defaultdict(lambda:[0.0,0.0,0.0])   # bin -> [sense,anti,bidir] FPM
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
GMAX=max([v[0]+v[2]/2 for d in CL.values() for v in d.values()]+[v[1]+v[2]/2 for d in CL.values() for v in d.values()]+[1.0])
import math; LGM=math.log10(GMAX+1); MAXH=0.018
def h(fpm): return MAXH*math.log10(fpm+1)/LGM
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig,axes=plt.subplots(1,3,figsize=(23.5,9.2),subplot_kw={"projection":"polar"},dpi=300)
base=np.linspace(0.92,0.30,len(CANON))
for ax,(lab,tp) in zip(axes,TPS):
    ax.set_theta_direction(-1); ax.set_theta_offset(np.pi/2); ax.set_ylim(0,1.06); ax.axis("off")
    for c in CHROMS:
        tt=np.linspace(theta(c,0),theta(c,clen[c]),60); ax.plot(tt,[1.0]*60,color="#222",lw=3.2,solid_capstyle="butt")
        ax.text((tt[0]+tt[-1])/2,1.04,c,ha="center",va="center",fontsize=6.5,fontweight="bold")
    for X,r in zip(CANON,base):
        d=CL[(X,tp)]
        if not d: continue
        bi=sorted(d); ang=bth[bi]; w=bwd[bi]
        outh=np.array([h(d[b][0]+d[b][2]/2) for b in bi]); inh=np.array([h(d[b][1]+d[b][2]/2) for b in bi])
        ax.bar(ang,outh,width=w,bottom=r,color=SENSE,align="center",edgecolor="none",linewidth=0,rasterized=True)
        ax.bar(ang,inh,width=w,bottom=r-inh,color=ANTI,align="center",edgecolor="none",linewidth=0,rasterized=True)
        ax.plot([0,2*np.pi],[r,r],color="#ddd",lw=0.2,zorder=0)
        ax.text(THLAB,r,X.replace("_","/"),fontsize=5.6,ha="center",va="center",fontweight="bold" if X in WILD else "normal",color="#C0392B" if X in WILD else "#222")
    ax.set_title(lab,fontsize=15,fontweight="bold",pad=12)
leg=[Line2D([0],[0],color=SENSE,lw=7,label="sense (+) — bar OUTWARD"),Line2D([0],[0],color=ANTI,lw=7,label="antisense (−) — bar INWARD"),Line2D([0],[0],color="#777",lw=7,label="bidirectional = bars BOTH ways")]
fig.legend(handles=leg,loc="lower center",bbox_to_anchor=(0.5,-0.02),ncol=3,fontsize=11,frameon=False,title="cluster expression by strand (bar length ∝ log FPM); strain name at the top spoke (red = wild-derived)",title_fontsize=10.5)
fig.suptitle("piRNA-cluster EXPRESSION + DIRECTIONALITY circos — 16 strains × 3 timepoints (PICB-combined → GRCm39): sense out / antisense in, bidirectional both ways",fontsize=13,fontweight="bold",y=1.02)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_circos_expr16.{e}",bbox_inches="tight")
print("wrote Fig_circos_expr16.{png,pdf,svg}")
