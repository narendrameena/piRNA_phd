#!/usr/bin/env python3
"""Circos #10 — TE ARMS-RACE x PING-PONG BIOGENESIS combined, strain- and timepoint-resolved. Per strain x 3
timepoints, each timepoint shows the same threat-vs-response pair as te_pirna16, but the RESPONSE is recoloured by
the ping-pong signature so three signatures live in one view:
  * OUTER half = active-TE EXPRESSION (small-RNA sense-to-TE, threat) — STACKED by family (L1/ERVK/ERVL greens), grows outward
  * INNER half = piRNA-on-TE (sRNA response), HEIGHT = log piRNA reads, COLOUR = 10A fraction (PING-PONG secondary-
    piRNA signature; plasma: dark = phasing/low-ping-pong, bright = ping-pong-active), grows inward
1U (primary-piRNA hallmark) is near-uniformly high genome-wide (see 1u16), so 10A is the informative axis here:
bright inner band meeting a tall green threat = a TE family being silenced by an active ping-pong cycle. Data:
active_te_expression_byfamily_tp.tsv + active_pirna_on_te_byfamily_tp.tsv + active_1u_bias_tp.tsv (col7 = 10A)."""
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
PG=f"{U}/pangenome_te"; FAI=f"{ROOT}/results/ref_genome/GRCm39.106.fasta.fai"
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]; TPS=["E16.5","P12.5","P20.5"]; CHROMS=[str(i) for i in range(1,20)]+["X"]; BIN=2_000_000
TEFAM=[("L1","#1B7837"),("ERVK","#5AAE61"),("ERVL","#A6DBA0")]   # TE expression (SENSE-to-TE small RNA) = green family shades
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
def loadfam(fn):
    D=defaultdict(lambda:defaultdict(dict))
    if os.path.exists(fn):
        for l in open(fn):
            X,tp,c,b,g,v=l.rstrip("\n").split("\t")
            if X in CANON and (c,int(b)) in binmap: D[(X,tp)][binmap[(c,int(b))]][g]=float(v)
    return D
TE=loadfam(f"{U}/active_te_expression_sRNA_tp.tsv"); PI=loadfam(f"{U}/active_pirna_on_te_sRNA_tp.tsv")   # sRNA-only: sense-to-TE=expression, antisense-to-TE=piRNA (no small-RNA sense-to-TE, no liftover)
A10=defaultdict(dict)   # (X,tp) -> bin -> 10A fraction
if os.path.exists(f"{U}/active_1u_bias_tp.tsv"):
    for l in open(f"{U}/active_1u_bias_tp.tsv"):
        f=l.rstrip("\n").split("\t")
        if len(f)>=7 and f[0] in CANON and (f[2],int(f[3])) in binmap: A10[(f[0],f[1])][binmap[(f[2],int(f[3]))]]=float(f[6])
TET={(k[0],k[1],b):sum(TE[k][b].values()) for k in TE for b in TE[k]}
PIT={(k[0],k[1],b):sum(PI[k][b].values()) for k in PI for b in PI[k]}
TMAX=np.percentile(list(TET.values()),99) if TET else 1.0; LTE=math.log10(TMAX+1)
PMAX=np.percentile(list(PIT.values()),99) if PIT else 1.0; LPI=math.log10(PMAX+1)
PPMAP=plt.get_cmap("plasma"); a10v=[v for k in A10 for v in A10[k].values()]
PN=mc.Normalize(np.percentile(a10v,5) if a10v else 0.0, np.percentile(a10v,95) if a10v else 1.0)
print(f"combined: TE99={TMAX:.0f} piRNA99={PMAX:.0f} 10A range={PN.vmin:.2f}-{PN.vmax:.2f}")
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig=plt.figure(figsize=(16,16),dpi=300); ax=fig.add_subplot(111,projection="polar")
ax.set_theta_direction(-1); ax.set_theta_offset(np.pi/2); ax.set_ylim(0,1.05); ax.axis("off")
for c in CHROMS:
    tt=np.linspace(theta(c,0),theta(c,clen[c]),60); ax.plot(tt,[1.0]*60,color="#222",lw=4,solid_capstyle="butt")
    ax.text((tt[0]+tt[-1])/2,1.04,c,ha="center",va="center",fontsize=8.5,fontweight="bold")
R_OUT,R_IN=0.95,0.17; grp_h=(R_OUT-R_IN)/len(CANON); gap_g=grp_h*0.12; tp_h=(grp_h-gap_g)/3; band_h=tp_h*0.46
def te_stack(key,base):   # green family threat, grows OUT
    td=TE.get(key,{})
    if not td: return
    bbA=np.array(sorted(td)); t=np.array([TET[(key[0],key[1],bb)] for bb in bbA])
    H=np.array([band_h*min(1.0,max(0.10,math.log10(x+1)/LTE)) for x in t]); run=np.full(len(bbA),base,float)
    for g,gc in TEFAM:
        seg=np.array([H[i]*(TE[key][bbA[i]].get(g,0)/t[i] if t[i]>0 else 0) for i in range(len(bbA))])
        ax.bar(bth[bbA],height=seg,width=bwd[bbA],bottom=run,color=gc,align="center",edgecolor="none",linewidth=0,rasterized=True); run=run+seg
def pi_pp(key,base):      # piRNA-on-TE response, grows IN, colour = 10A (ping-pong)
    td=PI.get(key,{}); a=A10.get(key,{})
    if not td: return
    bb=sorted(td); H=np.array([band_h*min(1.0,max(0.10,math.log10(PIT[(key[0],key[1],b)]+1)/LPI)) for b in bb])
    cols=[PPMAP(PN(a.get(b,PN.vmin))) for b in bb]
    ax.bar(bth[bb],height=H,width=bwd[bb],bottom=np.array([base-H[i] for i in range(len(bb))]),color=cols,align="center",edgecolor="none",linewidth=0,rasterized=True)
for k,X in enumerate(CANON):
    top=R_OUT-k*grp_h
    for j,tp in enumerate(TPS):
        tmid=top-gap_g-(j+0.5)*tp_h
        te_stack((X,tp),tmid+0.0008); pi_pp((X,tp),tmid-0.0008)
    ax.text(THLAB,top-gap_g-1.5*tp_h,X.replace("_","/"),fontsize=6.4,ha="center",va="center",fontweight="bold" if X in WILD else "normal",color="#C0392B" if X in WILD else "#222")
_t=R_OUT-gap_g
zoom_6nj(ax, rings=[("E16.5",_t-0.5*tp_h),("P12.5",_t-1.5*tp_h),("P20.5",_t-2.5*tp_h)], theta_c=theta("1",0)-0.006, fs=4.6, spread=True)
sm=cm.ScalarMappable(norm=PN,cmap=PPMAP); sm.set_array([])
cbar=fig.colorbar(sm,ax=ax,orientation="vertical",fraction=0.03,pad=0.02,shrink=0.5)
cbar.set_label("piRNA-on-TE 10A fraction = PING-PONG signature (bright = ping-pong-active response)",fontsize=8.5); cbar.ax.tick_params(labelsize=7)
leg=[Line2D([0],[0],color="#1B7837",lw=7,label="TE threat: L1"),Line2D([0],[0],color="#5AAE61",lw=7,label="TE: ERVK"),Line2D([0],[0],color="#A6DBA0",lw=7,label="TE: ERVL"),Line2D([0],[0],color=PPMAP(0.85),lw=7,label="piRNA-on-TE response (colour = 10A ping-pong)")]
fig.legend(handles=leg,loc="lower center",bbox_to_anchor=(0.5,0.05),ncol=2,fontsize=10.5,frameon=False,title="per strain × 3 timepoints: GREEN TE-expression threat (out, family) + piRNA-on-TE response (in, HEIGHT = log reads, COLOUR = 10A ping-pong)",title_fontsize=9.5)
fig.suptitle("TE ARMS-RACE × PING-PONG circos — TE-expression threat vs piRNA-on-TE response, response coloured by 10A ping-pong signature, 16 strains × 3 timepoints\n"
             "outer = TE-family RNA expression (green, threat); inner = piRNA-on-TE (height = log reads, colour = 10A = ping-pong activity); bright inner meeting tall green = active ping-pong silencing of that TE family",
             fontsize=11,fontweight="bold",y=0.99,linespacing=1.5)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_circos_combined16.{e}",bbox_inches="tight")
print("wrote Fig_circos_combined16.{png,pdf,svg}")
