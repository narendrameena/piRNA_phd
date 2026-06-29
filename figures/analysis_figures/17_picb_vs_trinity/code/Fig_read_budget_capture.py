#!/usr/bin/env python3
"""THEME 17 — piRNA READ-BUDGET + capture (16 strains x 3 tp). For each strain/tp two stacked views:
  (top, slim) READ BUDGET, absolute multimapper-weighted alignments (log): total small-RNA reads (all lengths)
     vs total 25-32 nt piRNA reads -> piRNA dominates the small-RNA pool (~86-92%).
  (bottom) CAPTURE, % of total piRNA: in PICB clusters vs Trinity EXON-blocks vs Trinity FULL-SPAN
     (intron-INCLUDED). The span bar shows the intron-inflation artifact (why exon blocks are used).
Data: picb_vs_trinity/read_budget_capture.csv (merge_read_budget.py). bar=mean over reps, dots=reps."""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from strain_order import STRAIN_ORDER, WILD, color_wild_labels
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
P=f"{ROOT}/analysis/claude_biomni_analysis/picb_vs_trinity"; T=f"{ROOT}/figures/analysis_figures/17_picb_vs_trinity"
d=pd.read_csv(f"{P}/read_budget_capture.csv")
TPS=["E16.5","P12.5","P20.5"]
CSR="#BBBBBB"; CPI="#E69F00"          # read budget: total sRNA grey, piRNA gold
CP="#117733"; CTE="#CC6677"; CTS="#882255"   # capture: PICB green, Trinity-exon rose, Trinity-span plum
order=[s for s in STRAIN_ORDER if s in set(d.strain)]; x=np.arange(len(order)); WPOS=[i for i,s in enumerate(order) if s in WILD]
fig=plt.figure(figsize=(15,13.5),dpi=300)
gs=fig.add_gridspec(6,1,height_ratios=[0.34,1,0.34,1,0.34,1],hspace=0.08,left=0.07,right=0.99,top=0.945,bottom=0.07)
for ti,tp in enumerate(TPS):
    axb=fig.add_subplot(gs[2*ti]); axc=fig.add_subplot(gs[2*ti+1],sharex=axb)
    sub=d[d.tp==tp]; g=sub.groupby("strain")
    # --- read budget (log absolute) ---
    tot=[g.get_group(s).total_all.mean() if s in g.groups else 0 for s in order]
    pir=[g.get_group(s).total_pirna.mean() if s in g.groups else 0 for s in order]
    if WPOS: axb.axvspan(min(WPOS)-0.5,max(WPOS)+0.5,color="#C0392B",alpha=0.06,zorder=0)
    axb.bar(x-0.2,tot,0.4,color=CSR,label="total small-RNA reads"); axb.bar(x+0.2,pir,0.4,color=CPI,label="total piRNA (25–32 nt)")
    axb.set_yscale("log"); axb.set_ylabel("reads (log)\nmm-weighted",fontsize=7.5); axb.tick_params(labelsize=6.5)
    plt.setp(axb.get_xticklabels(),visible=False); axb.spines[["top","right"]].set_visible(False)
    pf=100*np.nansum(pir)/max(np.nansum(tot),1)
    axb.set_title(f"{tp}  —  piRNA = {pf:.0f}% of small-RNA reads   ·   capture below = % of total piRNA",fontsize=10,fontweight="bold",loc="left")
    if ti==0: axb.legend(fontsize=7,frameon=False,ncol=2,loc="upper right")
    # --- capture (% of piRNA): PICB / Trinity-exon / Trinity-span(intron) ---
    if WPOS: axc.axvspan(min(WPOS)-0.5,max(WPOS)+0.5,color="#C0392B",alpha=0.06,zorder=0)
    bw=0.26
    for k,(col,key,lab) in enumerate([(CP,"pct_picb","in PICB clusters"),(CTE,"pct_trin_exon","in Trinity precursors (exon blocks)"),
                                      (CTS,"pct_trin_span","in Trinity precursors (FULL SPAN, intron-included)")]):
        xs=x+(k-1)*bw; m=[g.get_group(s)[key].mean() if s in g.groups else 0 for s in order]
        axc.bar(xs,m,bw,color=col,edgecolor="white",linewidth=0.3,label=lab,zorder=2,
                hatch=("///" if key=="pct_trin_span" else None))
        for xi,s in zip(xs,order):
            if s in g.groups:
                for v in g.get_group(s)[key].values: axc.scatter(xi,v,s=5,color="#222",zorder=3,linewidth=0)
    axc.set_ylabel("% of total piRNA",fontsize=8.5); axc.tick_params(labelsize=6.5); axc.spines[["top","right"]].set_visible(False)
    if ti==len(TPS)-1:
        axc.set_xticks(x); axc.set_xticklabels([s.replace("_","/") for s in order],rotation=45,ha="right",fontsize=8); color_wild_labels(axc,order)
    else: plt.setp(axc.get_xticklabels(),visible=False)
    if ti==0: axc.legend(fontsize=7.2,frameon=False,ncol=1,loc="upper right")
    me=sub.pct_trin_exon.mean(); ms=sub.pct_trin_span.mean()
    axc.text(0.005,0.93,f"intron inflation: Trinity {me:.1f}% → {ms:.1f}%  ({ms/me:.2f}×)",transform=axc.transAxes,fontsize=7.5,color=CTS,fontweight="bold",va="top")
fig.suptitle("piRNA read budget & capture — small-RNA pool → 25–32 nt piRNA → PICB vs Trinity (exon vs intron-included)",fontsize=13,fontweight="bold",y=0.975)
fig.text(0.5,0.012,"Per strain × tp · multimapper-weighted alignment counts from results/STAR_srna_strain_wise BAMs · capture = % of 25–32 nt piRNA inside each region set · "
  "Trinity FULL-SPAN (intron-included) over-counts vs EXON blocks (artifact) · PICB ≥ Trinity at E16.5/P12.5, ~tied at P20.5.",ha="center",fontsize=7,color="#666")
out=f"{T}/figures/Fig_read_budget_capture"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
d.to_csv(f"{T}/data/source_data/SourceData_Fig_read_budget_capture.csv",index=False)
print("wrote",out)
