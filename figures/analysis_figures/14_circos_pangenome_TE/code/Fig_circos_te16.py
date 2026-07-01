#!/usr/bin/env python3
"""Circos #4 — TE × piRNA ARMS-RACE across all 16 strains × 3 timepoints in ONE GRCm39 circle. Each strain is a
labelled group: (top) an active-TE band = the genomic THREAT, drawn as a STACKED FAMILY COMPOSITION (per-strain
RepeatMasker bp of young/active families per 2-Mb bin: L1 = green, ERVK = teal, ERVL = wine; segment widths ≈
family share, total bar HEIGHT ∝ log total bp), then (inward) THREE timepoint sub-rings E16.5→P12.5→P20.5 = the
developing piRNA DEFENCE (PICB clusters, bar HEIGHT ∝ log FPM, colour = strand). A "HOW TO READ" key (top-left)
decodes one strain block. Colour carries family/strand (categorical), height carries amount — never duplicated.
TE data = active_te_expression_sRNA.tsv (small-RNA sense-to-TE, own-genome chr1-19,X). Cluster data =
cluster_pav/picb_pangenome_clusters.tsv (PICB-COMBINED clusters, pangenome-projected to GRCm39; FPM =
all_primary_FPM, genomic strand). Cluster colour = genomic-strand ARCHITECTURE (uni +/−, dual), NOT
sense/antisense (sense/antisense is defined only relative to a TE/gene feature)."""
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
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]; TPS=[("E16.5","16.5dpc"),("P12.5","12.5dpp"),("P20.5","20.5dpp")]
CHROMS=[str(i) for i in range(1,20)]+["X"]; BIN=2_000_000; SCOL={"sense":"#0072B2","antisense":"#D55E00","bidir":"#9467bd"}
TEFAM=[("L1","#117733"),("ERVK","#44AA99"),("ERVL","#882255")]; TEC=dict(TEFAM)   # stack order + CB-safe family colours
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
# --- active-TE density per FAMILY (the threat; genomic, one stacked band per strain) ---
TE=defaultdict(lambda:defaultdict(dict))   # X -> bin -> {family: small-RNA sense-to-TE expression}
for l in open(f"{U}/active_te_expression_sRNA.tsv"):
    X,c,b,g,v=l.rstrip("\n").split("\t"); b=int(b)
    if X in CANON and (c,b) in binmap: TE[X][binmap[(c,b)]][g]=float(v)
tetot={(X,b):sum(TE[X][b].values()) for X in CANON for b in TE[X]}
TMAX=np.percentile(list(tetot.values()),99) if tetot else 1.0; LT=math.log10(TMAX+1)
# --- cluster expression per timepoint (the defence) — PICB-COMBINED clusters, pangenome-projected (NOT past liftover/PAV) ---
import pandas as _pd
_PC=_pd.read_csv(f"{U}/cluster_pav/picb_pangenome_clusters.tsv",sep="\t",dtype={"g39_chrom":str})
CL=defaultdict(lambda:defaultdict(lambda:[0.0,0.0,0.0]))   # (X,tp) -> bin -> [+, -, dual] FPM (genomic-strand architecture)
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
print(f"TE×piRNA family-stacked EXPRESSION: TE-expr TMAX(99pct)={TMAX:.0f} RNA reads/2Mb, cluster GMAX={GMAX:.0f} FPM")
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig=plt.figure(figsize=(15.5,15.5),dpi=300); ax=fig.add_subplot(111,projection="polar")
ax.set_theta_direction(-1); ax.set_theta_offset(np.pi/2); ax.set_ylim(0,1.05); ax.axis("off")
for c in CHROMS:
    tt=np.linspace(theta(c,0),theta(c,clen[c]),60); ax.plot(tt,[1.0]*60,color="#222",lw=4,solid_capstyle="butt")
    ax.text((tt[0]+tt[-1])/2,1.04,c,ha="center",va="center",fontsize=8.5,fontweight="bold")
R_OUT,R_IN=0.95,0.17; grp_h=(R_OUT-R_IN)/len(CANON); gap_g=grp_h*0.12; te_h=grp_h*0.24; sub_h=(grp_h-gap_g-te_h)/3
for k,X in enumerate(CANON):
    top=R_OUT-k*grp_h
    # TE LOAD band (the threat) — STACKED family composition, total height ∝ log total bp
    rt=top-gap_g-te_h; td=TE.get(X,{})
    if td:
        bbA=np.array(sorted(td)); tot=np.array([tetot[(X,b)] for b in bbA])
        H=np.array([te_h*min(1.0,max(0.10,math.log10(t+1)/LT)) for t in tot])
        bottom=np.full(len(bbA),rt,dtype=float)
        for g,gc in TEFAM:
            seg=np.array([H[i]*(TE[X][bbA[i]].get(g,0)/tot[i] if tot[i]>0 else 0.0) for i in range(len(bbA))])
            ax.bar(bth[bbA],height=seg,width=bwd[bbA],bottom=bottom,color=gc,align="center",edgecolor="none",linewidth=0,rasterized=True)
            bottom=bottom+seg
    # 3 timepoint cluster sub-rings (the defence) — height ∝ log FPM, strand colour
    for j,(lab,tp) in enumerate(TPS):
        r=rt-(j+1)*sub_h; d=CL[(X,tp)]
        if d:
            bb=sorted(d); hh=np.array([sub_h*min(1.0,max(0.14,math.log10(sum(d[b])+1)/LGM)) for b in bb])
            ax.bar(bth[bb],height=hh,width=bwd[bb],bottom=r,color=[SCOL[cat(d[b])] for b in bb],align="center",edgecolor="none",linewidth=0,rasterized=True)
    ax.text(THLAB,rt-1.5*sub_h,X.replace("_","/"),fontsize=6.6,ha="center",va="center",fontweight="bold" if X in WILD else "normal",color="#C0392B" if X in WILD else "#222")
ax.annotate("",xy=(THLAB,R_IN-0.02),xytext=(THLAB,R_OUT+0.02),arrowprops=dict(arrowstyle="-|>",color="#888",lw=1.0))
ax.text(THLAB,R_OUT+0.05,"TE band, then E16.5→P12.5→P20.5 (inward)",fontsize=5.8,ha="center",va="bottom",color="#888")
# ---------- zoom callout: magnified C57BL/6NJ block (outermost, k=0) ----------
_rt=R_OUT-gap_g-te_h
zoom_6nj(ax, rings=[("TE band",R_OUT-gap_g-te_h/2),("E16.5",_rt-0.5*sub_h),("P12.5",_rt-1.5*sub_h),("P20.5",_rt-2.5*sub_h)], theta_c=theta("1",0)-0.006, fs=4.4, spread=True)
# ---------- legends ----------
leg=[Line2D([0],[0],color=TEC["L1"],lw=7,label="TE: L1 (LINE)"),Line2D([0],[0],color=TEC["ERVK"],lw=7,label="TE: ERVK/ERVB"),Line2D([0],[0],color=TEC["ERVL"],lw=7,label="TE: ERVL"),
     Line2D([0],[0],color=SCOL["sense"],lw=7,label="cluster + strand (uni)"),Line2D([0],[0],color=SCOL["antisense"],lw=7,label="cluster − strand (uni)"),Line2D([0],[0],color=SCOL["bidir"],lw=7,label="dual-strand (bidirectional)")]
fig.legend(handles=leg,loc="lower center",bbox_to_anchor=(0.5,0.045),ncol=3,fontsize=10.5,frameon=False,title="per strain: STACKED active-TE family EXPRESSION (small-RNA sense-to-TE) band (height ∝ log RNA reads, colour = family) + 3 timepoint PICB-cluster sub-rings E16.5/P12.5/P20.5 (height ∝ log FPM, colour = cluster genomic-strand architecture)",title_fontsize=10)
fig.suptitle("TE × piRNA ARMS-RACE circos — active-TE family EXPRESSION (small-RNA sense-to-TE) vs piRNA-cluster defence across development, 16 strains × 3 timepoints in one GRCm39 circle\n"
             "per strain: STACKED TE-family EXPRESSION band (L1/ERVK/ERVL) + 3 timepoint sub-rings (inward) = PICB-cluster expression (genomic-strand architecture colour); colour = TE family / cluster strand-architecture, height = log amount (never duplicated)",
             fontsize=12,fontweight="bold",y=0.99,linespacing=1.5)
import os as _os; _SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/14_circos_pangenome_TE/data/source_data"; _os.makedirs(_SD,exist_ok=True)
_pd.DataFrame([(X,bins[bi][0],bins[bi][1],fam,round(v,3)) for X in TE for bi,fd in TE[X].items() for fam,v in fd.items()],columns=["strain","chrom","bin_start","TE_family","sRNA_expression"]).to_csv(f"{_SD}/SourceData_Fig_circos_te16_TEband.csv",index=False)   # outer TE-expression band (per strain x 2Mb bin x family)
_pd.DataFrame([(X,tp,bins[bi][0],bins[bi][1],round(f[0],3),round(f[1],3),round(f[2],3)) for (X,tp),d in CL.items() for bi,f in d.items()],columns=["strain","timepoint","chrom","bin_start","plus_FPM","minus_FPM","dual_FPM"]).to_csv(f"{_SD}/SourceData_Fig_circos_te16_clusterrings.csv",index=False)   # 3 timepoint cluster sub-rings (per strain x tp x 2Mb bin, strand FPM)
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_circos_te16.{e}",bbox_inches="tight")
print("wrote Fig_circos_te16.{png,pdf,svg}")
