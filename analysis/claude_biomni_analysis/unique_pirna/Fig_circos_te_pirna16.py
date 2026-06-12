#!/usr/bin/env python3
"""Circos #8 — TE-SILENCING ARMS-RACE, STRAIN- and TIMEPOINT-resolved: per strain, THREE timepoint groups
(E16.5/P12.5/P20.5); within each timepoint the active-TE THREAT vs the direct piRNA RESPONSE, both stacked by TE
family (L1/ERVK/IAP) so they compare family-by-family. For each strain × timepoint:
  * OUTER half-band = active-TE EXPRESSION (RNA-seq on TEs, grows OUTWARD) = threat
  * INNER half-band = piRNA-on-TE (sRNA/piRNA reads on the same TEs, grows INWARD) = direct response
Sense TE transcript × antisense piRNA-on-TE = ping-pong, so where a family's threat & response half-bands are both
tall, that family is being actively silenced at that timepoint. Direct piRNA-on-TE (not PICB clusters). Data:
active_te_expression_byfamily_tp.tsv + active_pirna_on_te_byfamily_tp.tsv (per strain, per timepoint, per family)."""
import warnings; warnings.filterwarnings("ignore")
import sys,os,math; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
from circos_readkey import zoom_6nj
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
PG=f"{U}/pangenome_te"; FAI=f"{ROOT}/results/ref_genome/GRCm39.106.fasta.fai"
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]; TPS=["E16.5","P12.5","P20.5"]; CHROMS=[str(i) for i in range(1,20)]+["X"]; BIN=2_000_000
TEFAM=[("L1","#1B7837"),("ERVK","#5AAE61"),("IAP","#A6DBA0")]   # TE threat = GREEN shades (dark→light = L1→ERVK→IAP)
PIFAM=[("L1","#762A83"),("ERVK","#9970AB"),("IAP","#C2A5CF")]   # piRNA response = PURPLE shades (same dark→light family order)
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
def load(fn):
    D=defaultdict(lambda:defaultdict(dict))
    if os.path.exists(fn):
        for l in open(fn):
            X,tp,c,b,g,v=l.rstrip("\n").split("\t")
            if X in CANON and (c,int(b)) in binmap: D[(X,tp)][binmap[(c,int(b))]][g]=float(v)
    return D
TE=load(f"{U}/active_te_expression_byfamily_tp.tsv"); PI=load(f"{U}/active_pirna_on_te_byfamily_tp.tsv")
def tot(D): return {(k[0],k[1],b):sum(D[k][b].values()) for k in D for b in D[k]}
TET=tot(TE); PIT=tot(PI)
TMAX=np.percentile(list(TET.values()),99) if TET else 1.0; LTE=math.log10(TMAX+1)
PMAX=np.percentile(list(PIT.values()),99) if PIT else 1.0; LPI=math.log10(PMAX+1)
print(f"te_pirna tp-resolved: TE99={TMAX:.0f} piRNA99={PMAX:.0f}; TE strains={len(set(k[0] for k in TE))} PI strains={len(set(k[0] for k in PI))}")
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig=plt.figure(figsize=(16,16),dpi=300); ax=fig.add_subplot(111,projection="polar")
ax.set_theta_direction(-1); ax.set_theta_offset(np.pi/2); ax.set_ylim(0,1.05); ax.axis("off")
for c in CHROMS:
    tt=np.linspace(theta(c,0),theta(c,clen[c]),60); ax.plot(tt,[1.0]*60,color="#222",lw=4,solid_capstyle="butt")
    ax.text((tt[0]+tt[-1])/2,1.04,c,ha="center",va="center",fontsize=8.5,fontweight="bold")
R_OUT,R_IN=0.95,0.17; grp_h=(R_OUT-R_IN)/len(CANON); gap_g=grp_h*0.12; tp_h=(grp_h-gap_g)/3; band_h=tp_h*0.46
def stack(key,B2T,D,base,maxh,LM,fams,outward=True):
    td=D.get(key,{})
    if not td: return
    bbA=np.array(sorted(td)); t=np.array([B2T[(key[0],key[1],bb)] for bb in bbA])
    H=np.array([maxh*min(1.0,max(0.10,math.log10(x+1)/LM)) for x in t])
    run=np.full(len(bbA),base,dtype=float) if outward else np.array([base-H[i] for i in range(len(bbA))])
    for g,gc in fams:
        seg=np.array([H[i]*(D[key][bbA[i]].get(g,0)/t[i] if t[i]>0 else 0.0) for i in range(len(bbA))])
        ax.bar(bth[bbA],height=seg,width=bwd[bbA],bottom=run,color=gc,align="center",edgecolor="none",linewidth=0,rasterized=True); run=run+seg
for k,X in enumerate(CANON):
    top=R_OUT-k*grp_h
    for j,tp in enumerate(TPS):
        tmid=top-gap_g-(j+0.5)*tp_h
        stack((X,tp),TET,TE,tmid+0.0008,band_h,LTE,TEFAM,outward=True)
        stack((X,tp),PIT,PI,tmid-0.0008,band_h,LPI,PIFAM,outward=False)
    ax.text(THLAB,top-gap_g-1.5*tp_h,X.replace("_","/"),fontsize=6.4,ha="center",va="center",fontweight="bold" if X in WILD else "normal",color="#C0392B" if X in WILD else "#222")
# 6NJ in-place timepoint labels (arrow style), at chr1 start
_t=R_OUT-gap_g
zoom_6nj(ax, rings=[("E16.5",_t-0.5*tp_h),("P12.5",_t-1.5*tp_h),("P20.5",_t-2.5*tp_h)], theta_c=theta("1",0)-0.006, fs=4.6, spread=True)
leg=[Line2D([0],[0],color="#1B7837",lw=7,label="TE threat: L1"),Line2D([0],[0],color="#5AAE61",lw=7,label="TE: ERVK"),Line2D([0],[0],color="#A6DBA0",lw=7,label="TE: IAP"),
     Line2D([0],[0],color="#762A83",lw=7,label="piRNA resp: L1"),Line2D([0],[0],color="#9970AB",lw=7,label="piRNA: ERVK"),Line2D([0],[0],color="#C2A5CF",lw=7,label="piRNA: IAP")]
fig.legend(handles=leg,loc="lower center",bbox_to_anchor=(0.5,0.045),ncol=3,fontsize=11,frameon=False,title="GREEN = active-TE EXPRESSION (threat, outward) · PURPLE = piRNA-on-TE (response, inward); shade = family (dark→light: L1→ERVK→IAP); per strain × 3 timepoints; HEIGHT ∝ log reads",title_fontsize=9.5)
fig.suptitle("TE-SILENCING ARMS-RACE circos — the transposon threat vs the piRNA counter-strike, family-by-family across development, 16 strains × 3 timepoints in one GRCm39 circle\n"
             "GREEN = active-TE EXPRESSION firing OUTWARD (threat) · PURPLE = piRNA-on-TE striking INWARD (response); shade = family (L1/ERVK/IAP); tall green meeting tall purple of the SAME family = the ping-pong cycle caught in the act",
             fontsize=11.5,fontweight="bold",y=0.99,linespacing=1.5)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_circos_te_pirna16.{e}",bbox_inches="tight")
print("wrote Fig_circos_te_pirna16.{png,pdf,svg}")
