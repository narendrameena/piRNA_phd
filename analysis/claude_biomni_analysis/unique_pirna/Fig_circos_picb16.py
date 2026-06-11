#!/usr/bin/env python3
"""Circos of PICB piRNA clusters across all 16 strains x 3 timepoints, with SENSE/ANTISENSE (creative,
Nature-Genetics style, colourblind-safe). Clusters are small vs a 2.5 Gb genome, so the genome is BINNED
(2 Mb) and each occupied bin is a visible wedge. Three circular panels (E16.5/P12.5/P20.5): GRCm39 ideogram
(outer) + a CONSERVATION ring (viridis = # of 16 strains with a cluster in the bin, private->core) + 16
strain rings (canonical order; wild-derived bold) where each bin is coloured by the strand of that strain's
cluster(s) in the GRCm39 frame: sense (+) blue, antisense (-) vermilion, both/mixed grey (Okabe-Ito). Data =
per-timepoint stranded PICB-combined clusters halLifted to GRCm39."""
import warnings; warnings.filterwarnings("ignore")
import sys,os; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
from collections import defaultdict
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.cm as cm, matplotlib.colors as mc
from matplotlib.lines import Line2D
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
PAV=f"{U}/cluster_pav/bytp"; PG=f"{U}/pangenome_te"; FAI=f"{ROOT}/results/ref_genome/GRCm39.106.fasta.fai"
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]; TPS=[("E16.5","16.5dpc"),("P12.5","12.5dpp"),("P20.5","20.5dpp")]
CHROMS=[str(i) for i in range(1,20)]+["X"]; STRANDCOL={"+":"#0072B2","-":"#D55E00",".":"#9e9e9e"}; BIN=2_000_000
clen={}
for ln in open(FAI):
    f=ln.split("\t")
    if f[0] in CHROMS: clen[f[0]]=int(f[1])
gap=sum(clen.values())*0.006; off={}; cum=0
for c in CHROMS: off[c]=cum; cum+=clen[c]+gap
TOT=cum
def theta(c,p): return 2*np.pi*(off[c]+p)/TOT
bins=[]; binmap={}
for c in CHROMS:
    for b in range(int(np.ceil(clen[c]/BIN))):
        binmap[(c,b)]=len(bins); bins.append((c,b*BIN,min((b+1)*BIN,clen[c])))
NB=len(bins); bth=np.array([(theta(c,s)+theta(c,e))/2 for c,s,e in bins]); bwd=np.array([theta(c,e)-theta(c,s) for c,s,e in bins])
def stranded(bed):
    d=defaultdict(list)
    if os.path.exists(bed) and os.path.getsize(bed):
        for l in open(bed):
            f=l.rstrip("\n").split("\t")
            if f[0] in CHROMS:
                st=f[5].strip() if len(f)>5 and f[5].strip() in "+-." else "."
                d[f[0]].append((int(f[1]),int(f[2]),st))
    return d
def occ(sd):
    o=defaultdict(set)
    for c,L in sd.items():
        for s,e,st in L:
            for b in range(s//BIN,e//BIN+1):
                if (c,b) in binmap: o[binmap[(c,b)]].add(st)
    return o
CL={(X,tp):occ(stranded(f"{PAV}/{X}.{tp}.stranded.in_GRCm39.bed")) for X in CANON for _,tp in TPS}
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig,axes=plt.subplots(1,3,figsize=(23,9),subplot_kw={"projection":"polar"},dpi=300)
cmap=plt.get_cmap("viridis"); norm=mc.Normalize(1,16)
rings=np.linspace(0.90,0.30,len(CANON)); RH=0.030
for ax,(lab,tp) in zip(axes,TPS):
    ax.set_theta_direction(-1); ax.set_theta_offset(np.pi/2); ax.set_ylim(0,1.13); ax.axis("off")
    for c in CHROMS:
        tt=np.linspace(theta(c,0),theta(c,clen[c]),60); ax.plot(tt,[1.0]*60,color="#222",lw=3.2,solid_capstyle="butt")
        ax.text((tt[0]+tt[-1])/2,1.05,c,ha="center",va="center",fontsize=7,fontweight="bold")
    cons=np.zeros(NB)
    for X in CANON:
        for bi in CL[(X,tp)]: cons[bi]+=1
    ob=[bi for bi in range(NB) if cons[bi]>0]
    ax.bar(bth[ob],height=0.05,width=bwd[ob],bottom=0.945,color=[cmap(norm(cons[bi])) for bi in ob],align="center",edgecolor="none",linewidth=0)
    for X,r in zip(CANON,rings):
        grp=defaultdict(list)
        for bi,sts in CL[(X,tp)].items():
            col=STRANDCOL["+"] if sts=={"+"} else (STRANDCOL["-"] if sts=={"-"} else STRANDCOL["."])
            grp[col].append(bi)
        for col,bb in grp.items():
            ax.bar(bth[bb],height=RH if X in WILD else RH*0.8,width=bwd[bb],bottom=r,color=col,align="center",edgecolor="none",linewidth=0)
    ax.set_title(lab,fontsize=15,fontweight="bold",pad=12)
sm=cm.ScalarMappable(norm=norm,cmap=cmap); sm.set_array([])
cb=fig.colorbar(sm,ax=axes,orientation="horizontal",fraction=0.02,pad=0.02,aspect=55)
cb.set_label("conservation ring: # of 16 strains with a cluster in the 2-Mb bin (1 = strain-private → 16 = core)",fontsize=9); cb.ax.tick_params(labelsize=7.5)
sl=[Line2D([0],[0],color=STRANDCOL["+"],lw=4,label="sense (+ strand)"),Line2D([0],[0],color=STRANDCOL["-"],lw=4,label="antisense (− strand)"),Line2D([0],[0],color=STRANDCOL["."],lw=4,label="both strands / unstranded")]
fig.legend(handles=sl,loc="lower center",bbox_to_anchor=(0.5,-0.03),ncol=3,fontsize=10,frameon=False,title="strain-ring cluster strand (GRCm39 frame)",title_fontsize=10)
sleg=[Line2D([0],[0],color="#444",lw=2.6 if s in WILD else 1.3,label=s.replace("_","/")+(" ★" if s in WILD else "")) for s in CANON]
fig.legend(handles=sleg,loc="center left",bbox_to_anchor=(0.0,0.5),fontsize=6.4,frameon=False,title="strain rings\n(outer→inner)",title_fontsize=7.5)
fig.suptitle("piRNA-cluster CIRCOS — 16 strains × 3 timepoints (PICB-combined → GRCm39): conserved core + open strain-private accessory, sense/antisense resolved",fontsize=13.5,fontweight="bold",y=1.03)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_circos_picb16.{e}",bbox_inches="tight")
print("wrote Fig_circos_picb16.{png,pdf,svg}")
