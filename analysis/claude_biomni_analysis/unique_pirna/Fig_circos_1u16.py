#!/usr/bin/env python3
"""Circos #9 — 1U SIGNATURE map: 1U (5'-Uridine) bias of piRNA-sized small RNAs across all 16 strains x 3
timepoints in ONE GRCm39 circle. 1U is the hallmark of PRIMARY piRNAs, so this shows WHERE (and when) the small
RNAs are genuinely piRNA-like. Each strain = 3 nested timepoint sub-rings (E16.5->P12.5->P20.5, inward); each 2-Mb
bin: bar HEIGHT = log read depth (how much small RNA), COLOUR = 1U fraction (viridis: low->high = purple->yellow).
A tall yellow bin = abundant, strongly-1U (bona-fide piRNA) production; tall purple = abundant but NOT 1U-biased
(degradation/other). Strand-aware (reverse reads use the complemented 5' base). Data = active_1u_bias_tp.tsv."""
import warnings; warnings.filterwarnings("ignore")
import sys,os,math; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
from circos_readkey import zoom_6nj
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.cm as cm, matplotlib.colors as mc
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
PG=f"{U}/pangenome_te"; FAI=f"{ROOT}/results/ref_genome/GRCm39.106.fasta.fai"
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]; TPS=["E16.5","P12.5","P20.5"]; CHROMS=[str(i) for i in range(1,20)]+["X"]; BIN=2_000_000
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
D=defaultdict(dict)   # (X,tp) -> bin -> (1Ufrac, depth)
for l in open(f"{U}/active_1u_bias_tp.tsv"):
    X,tp,c,b,frac,dep=l.rstrip("\n").split("\t")
    if X in CANON and (c,int(b)) in binmap: D[(X,tp)][binmap[(c,int(b))]]=(float(frac),int(dep))
deps=[v[1] for k in D for v in D[k].values()]; DMAX=np.percentile(deps,99) if deps else 1.0; LD=math.log10(DMAX+1)
CMAP=plt.get_cmap("viridis"); CN=mc.Normalize(0.40,0.95)   # 1U fraction range for contrast
print(f"1U circos: DMAX(99pct)={DMAX:.0f} reads; strains={len(set(k[0] for k in D))}; tps={sorted(set(k[1] for k in D))}")
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig=plt.figure(figsize=(15.5,15.5),dpi=300); ax=fig.add_subplot(111,projection="polar")
ax.set_theta_direction(-1); ax.set_theta_offset(np.pi/2); ax.set_ylim(0,1.05); ax.axis("off")
for c in CHROMS:
    tt=np.linspace(theta(c,0),theta(c,clen[c]),60); ax.plot(tt,[1.0]*60,color="#222",lw=4,solid_capstyle="butt")
    ax.text((tt[0]+tt[-1])/2,1.04,c,ha="center",va="center",fontsize=8.5,fontweight="bold")
R_OUT,R_IN=0.95,0.18; grp_h=(R_OUT-R_IN)/len(CANON); gap_g=grp_h*0.18; sub_h=(grp_h-gap_g)/3
for k,X in enumerate(CANON):
    top=R_OUT-k*grp_h
    for j,tp in enumerate(TPS):
        r=top-gap_g-(j+1)*sub_h; d=D.get((X,tp),{})
        if d:
            bb=sorted(d); hh=np.array([sub_h*min(1.0,max(0.12,math.log10(d[b][1]+1)/LD)) for b in bb])
            ax.bar(bth[bb],height=hh,width=bwd[bb],bottom=r,color=[CMAP(CN(d[b][0])) for b in bb],align="center",edgecolor="none",linewidth=0,rasterized=True)
    ax.text(THLAB,top-gap_g-1.5*sub_h,X.replace("_","/"),fontsize=6.6,ha="center",va="center",fontweight="bold" if X in WILD else "normal",color="#C0392B" if X in WILD else "#222")
ax.annotate("",xy=(THLAB,R_IN-0.02),xytext=(THLAB,R_OUT+0.02),arrowprops=dict(arrowstyle="-|>",color="#888",lw=1.2))
ax.text(THLAB,R_OUT+0.05,"E16.5→P12.5→P20.5 (inward)",fontsize=6.3,ha="center",va="bottom",color="#888")
_rt=R_OUT-gap_g
zoom_6nj(ax, rings=[("E16.5",_rt-0.5*sub_h),("P12.5",_rt-1.5*sub_h),("P20.5",_rt-2.5*sub_h)], theta_c=theta("1",0)-0.006)
sm=cm.ScalarMappable(norm=CN,cmap=CMAP); sm.set_array([])
cbar=fig.colorbar(sm,ax=ax,orientation="vertical",fraction=0.03,pad=0.02,shrink=0.55)
cbar.set_label("1U fraction (5′-Uridine) — piRNA signature: yellow = strongly 1U (bona-fide piRNA), purple = not 1U-biased",fontsize=8.5); cbar.ax.tick_params(labelsize=7)
fig.text(0.5,0.055,"BAR HEIGHT ∝ log small-RNA read depth per 2-Mb bin · COLOUR = 1U fraction (5′-U piRNA signature) — height and colour are independent",ha="center",fontsize=10.5,color="#333")
fig.suptitle("1U-SIGNATURE circos — where small RNAs are genuinely piRNA (5′-Uridine bias), 16 strains × 3 timepoints in one GRCm39 circle\n"
             "each strain = 3 nested timepoint sub-rings (inward); BAR HEIGHT = log read depth, COLOUR = 1U fraction (viridis); tall + yellow = abundant bona-fide piRNA production; strain names (red = wild) at the spoke",
             fontsize=12,fontweight="bold",y=0.99,linespacing=1.5)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_circos_1u16.{e}",bbox_inches="tight")
print("wrote Fig_circos_1u16.{png,pdf,svg}")
