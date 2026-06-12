#!/usr/bin/env python3
"""Circos #2 — piRNA-read COVERAGE (sRNA 'bigwig' hit density) across all 16 strains x 3 timepoints, ONE circle.
Companion to the PICB-cluster circos: instead of called clusters, this shows the raw small-RNA read density
(where piRNAs actually map, how densely) from results/bigwig_tracks_parallel/{X}-{tp}.bw (own-genome PanSN,
placed on GRCm39 chromosomes by chrN+position — colinear at 2-Mb resolution). Each strain = 3 nested timepoint
sub-rings (E16.5->P12.5->P20.5, inward); each 2-Mb bin coloured by log mean coverage (magma heatmap). Strain
names (red = wild) at the top spoke; chromosome ideogram outer."""
import warnings; warnings.filterwarnings("ignore")
import sys,os,glob,math; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
from circos_readkey import zoom_6nj
import numpy as np, pyBigWig
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.cm as cm, matplotlib.colors as mc
from matplotlib.lines import Line2D
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG=f"{U}/pangenome_te"
BWDIR=f"{ROOT}/results/bamCoverageSrna/pacBio"; FAI=f"{ROOT}/results/ref_genome/GRCm39.106.fasta.fai"
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
def cov_bins(X,tp):                                  # per bin -> [plus-strand, minus-strand] mean coverage
    fp=sorted(glob.glob(f"{BWDIR}/{X}/{X}-{tp}.*_plusStrand.bw"))
    fm=sorted(glob.glob(f"{BWDIR}/{X}/{X}-{tp}.*_minusStrand.bw"))
    if not fp or not fm: return {}
    bwp=pyBigWig.open(fp[0]); bwm=pyBigWig.open(fm[0]); chl=bwp.chroms(); d={}
    for c in CHROMS:
        sc=f"{X}#1#chr{c}"
        if sc not in chl: continue
        L=chl[sc]
        for b in range(int(np.ceil(clen[c]/BIN))):
            s=b*BIN; e=min((b+1)*BIN,clen[c],L)
            if s>=L: break
            try: vp=bwp.stats(sc,s,e,type="mean")[0] or 0.0
            except: vp=0.0
            try: vm=bwm.stats(sc,s,e,type="mean")[0] or 0.0
            except: vm=0.0
            if vp>0 or vm>0: d[binmap[(c,b)]]=[vp,vm]
    bwp.close(); bwm.close(); return d
COV={(X,tp):cov_bins(X,tp) for X in CANON for _,tp in TPS}
allv=[v[0]+v[1] for d in COV.values() for v in d.values()]
VMAX=np.percentile(allv,99) if allv else 1.0; LV=math.log10(VMAX+1)
SCOL={"plus":"#0072B2","minus":"#D55E00","bidir":"#9467bd"}
def cat(p,m):   # dominant GENOMIC strand of coverage (plus/minus bigwig) = cluster architecture, NOT sense/antisense
    if p>0 and m>0 and min(p,m)>=0.3*max(p,m): return "bidir"
    return "plus" if p>=m else "minus"
print(f"coverage bins computed; VMAX(99pct)={VMAX:.1f}")
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
        r=top-gap_g-(j+1)*sub_h; d=COV[(X,tp)]
        if d:
            bb=sorted(d); hh=np.array([sub_h*min(1.0,max(0.12,math.log10(sum(d[b])+1)/LV)) for b in bb])
            ax.bar(bth[bb],height=hh,width=bwd[bb],bottom=r,color=[SCOL[cat(*d[b])] for b in bb],align="center",edgecolor="none",linewidth=0,rasterized=True)
    mid=top-gap_g-1.5*sub_h
    ax.text(THLAB,mid,X.replace("_","/"),fontsize=6.6,ha="center",va="center",fontweight="bold" if X in WILD else "normal",color="#C0392B" if X in WILD else "#222")
ax.annotate("",xy=(THLAB,R_IN-0.02),xytext=(THLAB,R_OUT+0.02),arrowprops=dict(arrowstyle="-|>",color="#888",lw=1.2))
ax.text(THLAB,R_OUT+0.05,"E16.5→P12.5→P20.5 (inward)",fontsize=6.3,ha="center",va="bottom",color="#888")
fig.legend(handles=[Line2D([0],[0],color=SCOL["plus"],lw=7,label="plus (+) strand dominant"),Line2D([0],[0],color=SCOL["minus"],lw=7,label="minus (−) strand dominant"),Line2D([0],[0],color=SCOL["bidir"],lw=7,label="dual-strand (both)")],loc="lower center",bbox_to_anchor=(0.5,0.04),ncol=3,fontsize=11,frameon=False,title="dominant GENOMIC strand of small-RNA coverage = cluster architecture (NOT sense/antisense); BAR HEIGHT ∝ log mean coverage per 2-Mb bin",title_fontsize=10.5)
fig.suptitle("piRNA-read COVERAGE circos (sRNA bigwig) — GENOMIC-STRAND-coloured (cluster architecture, NOT sense/antisense), 16 strains × 3 timepoints in one GRCm39 circle\n"
             "each strain = 3 nested timepoint sub-rings (inward); BAR HEIGHT ∝ log read density, COLOUR = dominant genomic strand (plus/minus bigwig); strain names (red = wild) at the spoke",
             fontsize=12.5,fontweight="bold",y=0.99,linespacing=1.5)
_rt=R_OUT-gap_g
zoom_6nj(ax, rings=[("E16.5",_rt-0.5*sub_h),("P12.5",_rt-1.5*sub_h),("P20.5",_rt-2.5*sub_h)], theta_c=theta("1",0)-0.006)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_circos_coverage16.{e}",bbox_inches="tight")
print("wrote Fig_circos_coverage16.{png,pdf,svg}")
