#!/usr/bin/env python3
"""Circos #5 — INTEGRATED multi-omic circos, STRAIN- and TIMEPOINT-resolved, all 16 strains in ONE GRCm39 circle.
Each strain is a labelled group of SIX nested rings (outer→inner), every bar HEIGHT ∝ its metric (log-scaled):
  band 1  CONSERVATION  — this strain's cluster bins shaded + scaled by # of 16 strains sharing the bin (viridis)
  band 2  active-TE LOAD — young/active L1Md/ERVL/ERVK bp per 2-Mb bin (Oranges)   [genomic]
  band 3  STRAIN-PRIVATE — # genuinely-unique piRNAs per bin (teal)               [pooled over timepoints]
  rings 4-6  piRNA-CLUSTER EXPRESSION at E16.5 → P12.5 → P20.5 (height ∝ log FPM, colour = strand)
So for every strain you read its conservation context, its active-TE load, its strain-private innovation, and how
its piRNA cluster expression develops — all on the same genome. Strain names (red = wild) at the spoke."""
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
VMAP=plt.get_cmap("Greens"); VN=mc.Normalize(1,16); OMAP=plt.get_cmap("Greys"); TMAP=plt.get_cmap("PuRd")
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
# cluster expression per (strain,timepoint)
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
import pandas as _pd   # piRNA clusters from CURRENT pangenome projection (NOT past bytp liftover)
_PC=_pd.read_csv(f"{U}/cluster_pav/picb_pangenome_clusters.tsv",sep="\t",dtype={"g39_chrom":str})
CL=defaultdict(lambda:defaultdict(lambda:[0.0,0.0,0.0]))
for _r in _PC.itertuples(index=False):
    if _r.strain not in CANON or _r.g39_chrom not in CHROMS: continue
    _b=((int(_r.start)+int(_r.end))//2)//BIN
    if (_r.g39_chrom,_b) in binmap:
        CL[(_r.strain,_r.tp)][binmap[(_r.g39_chrom,_b)]][0 if _r.strand=="+" else (1 if _r.strand=="-" else 2)]+=float(_r.all_primary_FPM)
GMAX=max([sum(v) for d in CL.values() for v in d.values()]+[1.0]); LGM=math.log10(GMAX+1)
def cat(v):
    s,a,b=v
    if b>=max(s,a) and b>0: return "bidir"
    if s>0 and a>0 and min(s,a)>=0.3*max(s,a): return "bidir"
    return "sense" if s>=a else "antisense"
# conservation: # strains with a cluster in the bin (any timepoint)
cons=np.zeros(NB); OCC={}
for X in CANON:
    occ=set()
    for _,tp in TPS: occ|=set(CL[(X,tp)].keys())
    OCC[X]=occ
    for b in occ: cons[b]+=1
# active-TE load (genomic)
TE=defaultdict(lambda:defaultdict(dict))   # X -> bin -> {family: small-RNA sense-to-TE expression}
for l in open(f"{U}/active_te_expression_sRNA.tsv"):
    X,c,b,g,v=l.rstrip("\n").split("\t"); b=int(b)
    if X in CANON and (c,b) in binmap: TE[X][binmap[(c,b)]][g]=float(v)
TETOT={(X,b):sum(TE[X][b].values()) for X in CANON for b in TE[X]}
tev=list(TETOT.values()); TMX=np.percentile(tev,99) if tev else 1.0; LT=math.log10(TMX+1)
# strain-private piRNA density (pooled over timepoints)
PRIV=defaultdict(dict)
for X in CANON:
    seen=defaultdict(set); bed=f"{U}/unique16/loci/{X}.cand_GRCm39.bed"
    if os.path.exists(bed):
        for l in open(bed):
            f=l.split("\t")
            if f[0] in CHROMS and len(f)>3 and "|" in f[3]:
                mid=(int(f[1])+int(f[2]))//2; b=mid//BIN
                if (f[0],b) in binmap: seen[binmap[(f[0],b)]].add(f[3])
    for b,nm in seen.items(): PRIV[X][b]=len(nm)
pv=[v for X in CANON for v in PRIV[X].values()]; PMX=np.percentile(pv,99) if pv else 1.0; LP=math.log10(PMX+1)
print(f"integrated16 strain/tp-resolved: TMX={TMX:.0f} PMX={PMX:.0f} GMAX={GMAX:.0f}")
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig=plt.figure(figsize=(16.5,16.5),dpi=300); ax=fig.add_subplot(111,projection="polar")
ax.set_theta_direction(-1); ax.set_theta_offset(np.pi/2); ax.set_ylim(0,1.05); ax.axis("off")
for c in CHROMS:
    tt=np.linspace(theta(c,0),theta(c,clen[c]),60); ax.plot(tt,[1.0]*60,color="#222",lw=4,solid_capstyle="butt")
    ax.text((tt[0]+tt[-1])/2,1.04,c,ha="center",va="center",fontsize=8.5,fontweight="bold")
R_OUT,R_IN=0.95,0.15; grp_h=(R_OUT-R_IN)/len(CANON); gap_g=grp_h*0.10
bh=(grp_h-gap_g)*0.36/3            # 3 context bands take 36%
sh=(grp_h-gap_g)*0.64/3            # 3 timepoint rings take 64%
def bars(rb,h,binvals,colors,floor=0.12):
    bb=sorted(binvals);
    if not bb: return
    hh=np.array([h*min(1.0,max(floor,binvals[b])) for b in bb])
    ax.bar(bth[bb],height=hh,width=bwd[bb],bottom=rb,color=[colors(b) for b in bb],align="center",edgecolor="none",linewidth=0,rasterized=True)
for k,X in enumerate(CANON):
    top=R_OUT-k*grp_h; r=top-gap_g
    # band 1: conservation — flat GOLD (TE now uses green), bar HEIGHT ∝ # strains sharing the bin
    r-=bh; bars(r,bh,{b:cons[b]/16 for b in OCC[X]},lambda b:"#E6AB02",floor=0.18)
    # band 2: active-TE family EXPRESSION (small-RNA sense-to-TE) — STACKED L1/ERVK/ERVL, total HEIGHT ∝ log RNA reads
    r-=bh; td=TE.get(X,{})
    if td:
        bbA=np.array(sorted(td)); tot=np.array([TETOT[(X,b)] for b in bbA])
        Hh=np.array([bh*min(1.0,max(0.12,math.log10(t+1)/LT)) for t in tot]); bot=np.full(len(bbA),r,dtype=float)
        for g,gc in [("L1","#117733"),("ERVK","#44AA99"),("ERVL","#882255")]:
            seg=np.array([Hh[i]*(TE[X][bbA[i]].get(g,0)/tot[i] if tot[i]>0 else 0.0) for i in range(len(bbA))])
            ax.bar(bth[bbA],height=seg,width=bwd[bbA],bottom=bot,color=gc,align="center",edgecolor="none",linewidth=0,rasterized=True); bot=bot+seg
    # band 3: strain-private density — flat magenta, bar HEIGHT ∝ log count
    r-=bh; bars(r,bh,{b:math.log10(PRIV[X][b]+1)/LP for b in PRIV.get(X,{})},lambda b:"#C51B7D",floor=0.15)
    # rings 4-6: cluster expression at the 3 timepoints (height ∝ log FPM, strand colour)
    for j,(lab,tp) in enumerate(TPS):
        r-=sh; d=CL[(X,tp)]
        bars(r,sh,{b:math.log10(sum(d[b])+1)/LGM for b in d},lambda b,d=d:SCOL[cat(d[b])],floor=0.14)
    ax.text(THLAB,top-gap_g-(grp_h-gap_g)*0.5,X.replace("_","/"),fontsize=6.4,ha="center",va="center",fontweight="bold" if X in WILD else "normal",color="#C0392B" if X in WILD else "#222")
ax.annotate("",xy=(THLAB,R_IN-0.02),xytext=(THLAB,R_OUT+0.02),arrowprops=dict(arrowstyle="-|>",color="#888",lw=1.0))
ax.text(THLAB,R_OUT+0.05,"cons·TE·private bands, then E16.5→P12.5→P20.5",fontsize=5.4,ha="center",va="bottom",color="#888")
leg=[Line2D([0],[0],color="#E6AB02",lw=7,label="conservation (gold, height = # strains)"),Line2D([0],[0],color="#117733",lw=7,label="TE expr: L1"),Line2D([0],[0],color="#44AA99",lw=7,label="TE expr: ERVK"),Line2D([0],[0],color="#882255",lw=7,label="TE expr: ERVL"),Line2D([0],[0],color="#C51B7D",lw=7,label="strain-private piRNA (magenta)"),
     Line2D([0],[0],color=SCOL["sense"],lw=7,label="+ strand (uni-strand cluster)"),Line2D([0],[0],color=SCOL["antisense"],lw=7,label="− strand (uni-strand cluster)"),Line2D([0],[0],color=SCOL["bidir"],lw=7,label="dual-strand (bidirectional)")]
fig.legend(handles=leg,loc="lower center",bbox_to_anchor=(0.5,0.05),ncol=4,fontsize=10,frameon=False,title="per strain (outer→inner): conservation · active-TE family EXPRESSION (small-RNA sense-to-TE, L1/ERVK/ERVL stacked) · strain-private BANDS, then 3 timepoint cluster-expression rings — every bar HEIGHT ∝ log metric",title_fontsize=9.5)
fig.text(0.5,0.012,"every band/ring HEIGHT ∝ its metric (conservation = # strains; TE = log RNA expression, colour = family; private = log count; expression = log FPM, colour = strand) — colour never duplicates height",ha="center",fontsize=8.5,color="#333")
fig.suptitle("INTEGRATED multi-omic piRNA circos — STRAIN- and TIMEPOINT-resolved, 16 strains × 3 timepoints in one GRCm39 circle\n"
             "per strain: conservation · active-TE · strain-private BANDS + 3 timepoint cluster-expression sub-rings (strand colour); every bar HEIGHT ∝ log metric; strain names (red = wild) at the spoke",
             fontsize=12,fontweight="bold",y=0.995,linespacing=1.5)
_r0=R_OUT-gap_g
zoom_6nj(ax, rings=[("conserv",_r0-0.5*bh),("TE",_r0-1.5*bh),("private",_r0-2.5*bh),("E16.5",_r0-3*bh-0.5*sh),("P12.5",_r0-3*bh-1.5*sh),("P20.5",_r0-3*bh-2.5*sh)], theta_c=theta("1",0)-0.006, fs=4.0, spread=True)
import os as _os, pandas as _pd; _SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/14_circos_pangenome_TE/data/source_data"; _os.makedirs(_SD,exist_ok=True)
_pd.DataFrame([(X,tp,bins[bi][0],bins[bi][1],round(v[0],3),round(v[1],3),round(v[2],3)) for (X,tp),_dd in CL.items() for bi,v in _dd.items()],columns=["strain","timepoint","chrom","bin_start","plus_FPM","minus_FPM","dual_FPM"]).to_csv(f"{_SD}/SourceData_Fig_circos_integrated16_clusters.csv",index=False)
_pd.DataFrame([(X,bins[bi][0],bins[bi][1],fam,round(v,3)) for X,_bd in TE.items() for bi,_fd in _bd.items() for fam,v in _fd.items()],columns=["strain","chrom","bin_start","TE_family","value"]).to_csv(f"{_SD}/SourceData_Fig_circos_integrated16_TE.csv",index=False)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_circos_integrated16.{e}",bbox_inches="tight")
print("wrote Fig_circos_integrated16.{png,pdf,svg}")
