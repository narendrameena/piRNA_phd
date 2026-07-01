#!/usr/bin/env python3
"""THEME 17 — piRNA read capture across ALL 16 strains x 3 timepoints, PER-REPLICATE (bar = mean, error = ±SD,
dots = the replicates). Of total piRNA (25-32 nt, multimapper-weighted), % mapping inside PICB clusters vs Trinity
precursor EXON blocks. Aggregates per-rep results from cap16_reps/{strain}-{tp}.{rep}.csv (cols strain,tpn,rep,
total,in_picb,in_trin) + rep1 (cap16_reps/*.1.csv copied from the rep1 array). 3 panels (E16.5/P12.5/P20.5),
per-strain (canonical order, wild shaded). DEVELOPMENTAL: PICB >> Trinity at E16.5 (~4x), > at P12.5 (~2x),
~tied at P20.5 (few dominant pachytene precursors well-assembled by Trinity)."""
import sys, glob; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from strain_order import STRAIN_ORDER, WILD
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
P="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity"
rows=[]
for f in glob.glob(f"{P}/cap16_reps/*.csv"):
    a=open(f).read().strip().split(",")
    if len(a)<6 or "NA" in a: continue
    rows.append(dict(strain=a[0],tp=a[1],rep=a[2],total=int(a[3]),in_picb=int(a[4]),in_trin=int(a[5])))
d=pd.DataFrame(rows); d["pct_picb"]=100*d.in_picb/d.total; d["pct_trin"]=100*d.in_trin/d.total
d.sort_values(["tp","strain","rep"]).to_csv(f"{P}/read_capture_pirna_all16_reps.csv",index=False)
TPO=["E16.5","P12.5","P20.5"]; CP="#117733"; CT="#CC6677"
order=[s for s in STRAIN_ORDER if s in set(d.strain)]; WPOS=[i for i,s in enumerate(order) if s in WILD]
rng=np.random.default_rng(0)
fig,axes=plt.subplots(3,1,figsize=(14,11),dpi=300,sharex=True)
for ax,tpn in zip(axes,TPO):
    sub=d[d.tp==tpn]; x=np.arange(len(order)); bw=0.4
    if WPOS: ax.axvspan(min(WPOS)-0.5,max(WPOS)+0.5,color="#C0392B",alpha=0.06,zorder=0)
    pm=[];tm=[]
    for xi,s in zip(x,order):
        for col,off,store,c in [("pct_picb",-bw/2,pm,CP),("pct_trin",bw/2,tm,CT)]:
            v=sub[sub.strain==s][col].values
            m=v.mean() if len(v) else np.nan; store.append(m)
            ax.bar(xi+off,m if len(v) else 0,bw,color=c,edgecolor="white",linewidth=0.3,zorder=2)
            if len(v)>1: ax.errorbar(xi+off,m,yerr=v.std(ddof=1),fmt="none",ecolor="#333",elinewidth=0.5,capsize=1.5,zorder=3)
            if len(v): ax.scatter(np.full(len(v),xi+off)+rng.uniform(-0.07,0.07,len(v)),v,s=6,color="#222",zorder=4,linewidths=0)
    mp=np.nanmean(pm); mt=np.nanmean(tm)
    ax.axhline(mp,color=CP,ls="--",lw=0.7,zorder=1); ax.axhline(mt,color=CT,ls="--",lw=0.7,zorder=1)
    miss=[s for s in order if not (sub.strain==s).any()]
    extra=f"  ·  missing: {', '.join(s.replace('_','/') for s in miss)}" if miss else ""
    ax.set_title(f"{tpn}  —  mean PICB {mp:.1f}%  vs  Trinity {mt:.1f}%  (ratio {mp/mt:.1f}×){extra}",fontsize=10.5,fontweight="bold",loc="left")
    ax.set_ylabel("% of total piRNA",fontsize=9.5); ax.spines[["top","right"]].set_visible(False)
    top=np.nanmax([np.nanmax(pm),np.nanmax(tm)]); ax.set_ylim(0,top*1.22)
axes[-1].set_xticks(np.arange(len(order))); axes[-1].set_xticklabels([s.replace("_","/") for s in order],rotation=45,ha="right",fontsize=9)
for lab,s in zip(axes[-1].get_xticklabels(),order): lab.set_color("#C0392B" if s in WILD else "#222")
axes[0].legend(handles=[Patch(facecolor=CP,label="in PICB clusters"),Patch(facecolor=CT,label="in Trinity precursors (exon blocks)"),
    plt.Line2D([],[],marker='o',color='w',markerfacecolor='#222',markersize=4,label='replicate')],fontsize=8.5,frameon=False,ncol=1,loc="upper left")
if WPOS: axes[0].text(np.mean(WPOS),axes[0].get_ylim()[1]*0.97,"wild-derived",ha="center",va="top",fontsize=8.5,fontweight="bold",color="#C0392B")
fig.suptitle("piRNA read capture across 16 strains — PICB vs Trinity (of total piRNA, 25–32 nt; multimapper-weighted; bar=mean, ±SD, dots=reps)",fontsize=12.5,fontweight="bold",y=0.995)
fig.text(0.5,0.005,"DEVELOPMENTAL: PICB captures ~4× more piRNA than Trinity at E16.5 and ~2× at P12.5, converging to ~tied at P20.5 — pachytene precursors are few + dominant, so de-novo Trinity catches up only there. "
  "dashed = panel mean · bar=mean of replicates, error=±SD, dots=reps · sets overlap (not additive) · Trinity = exon blocks.",ha="center",fontsize=7,color="#666")
fig.tight_layout(rect=[0,0.012,1,0.985])
out=f"{P}/Fig_picb_vs_trinity_readcapture_all16"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
d.sort_values(["tp","strain","rep"]).to_csv("/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/source_data/SourceData_Fig_picb_vs_trinity_readcapture_all16.csv",index=False)
print(f"n={d.groupby(['strain','tp']).ngroups} strain×tp cells, {len(d)} rep-rows"); print(d.groupby("tp")[["pct_picb","pct_trin"]].mean().round(2).reindex(TPO).to_string()); print("wrote",out)
