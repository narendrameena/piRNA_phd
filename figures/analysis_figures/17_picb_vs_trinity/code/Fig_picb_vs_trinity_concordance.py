#!/usr/bin/env python3
"""THEME 17 — PICB clusters vs Trinity precursors: genomic concordance (16 strains x 3 tp).
Both in the SAME strain REL-2205 assembly (verified chr1 len identical PanSN vs PICB refFasta). PICB = own-coord
clusters (merged); Trinity = 100/100 precursors (rpm>100 & rpkm>100, thesis Ch.6), mapped, union of 3 reps. Source:
analysis/claude_biomni_analysis/picb_vs_trinity/overlap_per_strain_tp.csv (compute_overlap.py).
Trinity intervals = EXON BLOCKS (bed12tobed6), not intron-inflated spans; counting unit = distinct contigs.
Panels: A locus counts per strain (PICB 3-7x more); B reciprocal overlap (most Trinity precursors sit ON a PICB
cluster, but Trinity recovers few clusters); C fragmentation + median length (Trinity over-segments, shorter)."""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from strain_order import STRAIN_ORDER, WILD
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
P=ROOT+"/analysis/claude_biomni_analysis/picb_vs_trinity"
d=pd.read_csv(f"{P}/overlap_per_strain_tp.csv"); TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
d["TP"]=d.tp.map(TPN); TPO=["E16.5","P12.5","P20.5"]
CP="#117733"  # PICB green   CT="#CC6677" Trinity rose
CT="#CC6677"; TCOL={"E16.5":"#F0C9A0","P12.5":"#E69F00","P20.5":"#B4500A"}
TPDOT={"E16.5":"#4393C3","P12.5":"#E8852B","P20.5":"#B2182B"}   # dot colour = timepoint (Panel A)
order=[s for s in STRAIN_ORDER if s in set(d.strain)]; WPOS=[i for i,s in enumerate(order) if s in WILD]
fig=plt.figure(figsize=(15,10.5),dpi=300); gs=fig.add_gridspec(2,2,hspace=0.42,wspace=0.22,left=0.07,right=0.985,top=0.91,bottom=0.07)
# ---- A: locus counts per strain (mean over tp, dots = 3 tp), log y ----
axA=fig.add_subplot(gs[0,:]); x=np.arange(len(order)); bw=0.38
if WPOS: axA.axvspan(min(WPOS)-0.5,max(WPOS)+0.5,color="#C0392B",alpha=0.06,zorder=0)
for k,(col,key,lab) in enumerate([(CP,"n_picb","PICB clusters"),(CT,"n_trin_contigs","Trinity precursors (200/200)")]):
    xs=x+(k-0.5)*bw
    for xi,s in zip(xs,order):
        sub=d[d.strain==s]; m=sub[key].mean()
        axA.bar(xi,m,bw,color=col,edgecolor="white",linewidth=0.3,zorder=2)
        for _,rr in sub.iterrows():
            axA.scatter(xi,rr[key],s=12,color=TPDOT[rr["TP"]],edgecolor="white",linewidth=0.25,zorder=3)
axA.set_yscale("log"); axA.set_ylim(top=axA.get_ylim()[1]*3.5); axA.set_ylabel("loci per strain (log; mean of 3 tp, dots = tp)",fontsize=10)  # headroom band for legends above the data
axA.set_xticks(x); axA.set_xticklabels([s.replace("_","/") for s in order],rotation=45,ha="right",fontsize=8)
for lab,s in zip(axA.get_xticklabels(),order): lab.set_color("#C0392B" if s in WILD else "#222")
leg1=axA.legend(handles=[Patch(facecolor=CP,label="PICB clusters"),Patch(facecolor=CT,label="Trinity precursors (100 RPM & 100 RPKM, thesis)")],fontsize=9,frameon=False,ncol=1,loc="upper right")
axA.add_artist(leg1)
axA.legend(handles=[Line2D([],[],marker='o',color='w',markerfacecolor=TPDOT[t],markeredgecolor='white',markersize=6,label=t) for t in TPO],fontsize=7.5,frameon=False,ncol=3,loc="upper left",title="dots = timepoint",title_fontsize=7.5)
axA.set_title("A  PICB calls ~2–4× more loci than Trinity (systematic across all 16 strains)",fontsize=11,fontweight="bold",loc="left")
axA.spines[["top","right"]].set_visible(False)
if WPOS: axA.text(np.mean(WPOS),axA.get_ylim()[1]*0.45,"wild-derived",ha="center",fontsize=8.5,fontweight="bold",color="#C0392B")
# ---- B: reciprocal overlap per tp ----
axB=fig.add_subplot(gs[1,0]); xb=np.arange(len(TPO)); bw2=0.38
g=d.groupby("TP")
m1=[g.get_group(t).pct_trin_in_picb.mean() for t in TPO]; s1=[g.get_group(t).pct_trin_in_picb.std() for t in TPO]
m2=[g.get_group(t).pct_picb_recovered.mean() for t in TPO]; s2=[g.get_group(t).pct_picb_recovered.std() for t in TPO]
axB.bar(xb-bw2/2,m1,bw2,yerr=s1,color="#CC6677",capsize=2,error_kw=dict(elinewidth=0.6),label="% Trinity precursors that overlap a PICB cluster")
axB.bar(xb+bw2/2,m2,bw2,yerr=s2,color="#117733",capsize=2,error_kw=dict(elinewidth=0.6),label="% PICB clusters recovered by Trinity")
for xi,v,s in zip(xb-bw2/2,m1,s1): axB.text(xi,v+s+2.2,f"{v:.0f}%",ha="center",fontsize=7,color="#CC6677",fontweight="bold")
for xi,v,s in zip(xb+bw2/2,m2,s2): axB.text(xi,v+s+2.2,f"{v:.0f}%",ha="center",fontsize=7,color="#117733",fontweight="bold")
axB.set_xticks(xb); axB.set_xticklabels(TPO,fontsize=9); axB.set_ylabel("% (mean±SD over 16 strains)",fontsize=9.5)
axB.set_ylim(0,85); axB.legend(fontsize=7.3,frameon=False,loc="upper right")
axB.set_title("B  Most Trinity precursors (~60–70%) sit ON a PICB cluster,\n     but Trinity recovers <18% of clusters (PICB far more comprehensive)",fontsize=10,fontweight="bold",loc="left")
axB.spines[["top","right"]].set_visible(False)
# ---- C: fragmentation + median length ----
axC=fig.add_subplot(gs[1,1])
mf=[g.get_group(t).frag_mean.mean() for t in TPO]; sf=[g.get_group(t).frag_mean.std() for t in TPO]
axC.bar(xb,mf,0.5,yerr=sf,color="#882255",capsize=2,error_kw=dict(elinewidth=0.6))
for xi,v in zip(xb,mf): axC.text(xi,v+0.1,f"{v:.1f}×",ha="center",fontsize=8,color="#882255",fontweight="bold")
axC.set_xticks(xb); axC.set_xticklabels(TPO,fontsize=9); axC.set_ylabel("Trinity contigs per recovered\nPICB cluster (mean±SD)",fontsize=9.5)
axC.set_ylim(0,max(mf)*1.5)
pl=[g.get_group(t).picb_med_len.mean() for t in TPO]; tl=[g.get_group(t).trin_med_len.mean() for t in TPO]
axC2=axC.twinx()
axC2.plot(xb,pl,"o-",color="#117733",lw=1.2,ms=5,label="PICB cluster median length")
axC2.plot(xb,tl,"s--",color="#CC6677",lw=1.2,ms=5,label="Trinity precursor median length")
axC2.set_ylabel("median locus length (bp)",fontsize=9.5); axC2.set_ylim(0,3400)
axC2.legend(fontsize=7.3,frameon=False,loc="lower center")
axC.set_title("C  Trinity FRAGMENTS clusters (over-segmentation) &\n     is shorter than PICB clusters",fontsize=10,fontweight="bold",loc="left")
axC.spines[["top"]].set_visible(False); axC2.spines[["top"]].set_visible(False)
fig.suptitle("PICB clusters vs Trinity precursors — genomic concordance (16 strains × 3 timepoints; same strain assembly, verified)",fontsize=13,fontweight="bold",y=0.965)
fig.text(0.5,0.015,"PICB = own-coord clusters (merged) · Trinity = de-novo contigs ≥500 bp, 100 RPM & 100 RPKM (thesis Ch.6), mapped, EXON BLOCKS (bed12tobed6), union of 3 reps · overlap by bedtools · "
  "most Trinity precursors fall on PICB clusters but recover only a small fraction of them — decisive test = piRNA read capture (companion figure).",ha="center",fontsize=7,color="#666")
out=f"{P}/Fig_picb_vs_trinity_concordance"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
d.to_csv(ROOT+"/analysis/claude_biomni_analysis/source_data/SourceData_Fig_picb_vs_trinity_concordance.csv",index=False)
print("wrote",out)
