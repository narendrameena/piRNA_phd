#!/usr/bin/env python3
"""Pooled piRNA length distribution (count-weighted, all 27 pilot samples) and a DATA-DRIVEN
length window (no magic numbers): the contiguous length range around the mode whose per-length
read fraction is >= half the modal fraction (FWHM / 'major peak'); also report the contiguous
range capturing >=90% of in-peak reads.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
d=pd.read_csv(f"{U}/QC_length_distribution.csv")
prof=d.groupby("length")["reads"].sum()
prof=prof[(prof.index>=18)&(prof.index<=40)]
pct=(prof/prof.sum()*100)
mode=int(pct.idxmax()); half=pct.max()/2
# FWHM major peak: contiguous lengths around the mode with pct >= half-max
lo=mode
while lo-1 in pct.index and pct[lo-1]>=half: lo-=1
hi=mode
while hi+1 in pct.index and pct[hi+1]>=half: hi+=1
# broader: smallest contiguous window around mode capturing >=90% of total reads
def cover(a,b): return pct[(pct.index>=a)&(pct.index<=b)].sum()
a,b=mode,mode
while cover(a,b)<90:
    la=pct.get(a-1,0); rb=pct.get(b+1,0)
    if la>=rb and (a-1) in pct.index: a-=1
    elif (b+1) in pct.index: b+=1
    else: break
print("pooled count-weighted length distribution (% of reads):")
print(pct.round(2).to_string())
print(f"\nmode = {mode} nt ({pct.max():.1f}%)")
print(f"FWHM 'major peak' (>= half-max {half:.1f}%): {lo}-{hi} nt  (captures {cover(lo,hi):.1f}% of reads)")
print(f">=90%-coverage window: {a}-{b} nt")
# figure
fig,ax=plt.subplots(figsize=(6.0,3.6),dpi=300)
plt.rcParams.update({"font.family":"Liberation Sans"})
xs=pct.index.values
ax.bar(xs,pct.values,width=0.85,color=["#E69F00" if lo<=x<=hi else "#cccccc" for x in xs],
       edgecolor="white",linewidth=0.4,zorder=3)
ax.axhline(half,ls="--",lw=0.7,color="#888",zorder=2)
ax.text(xs.max(),half,f" half-max",va="bottom",ha="right",fontsize=6,color="#888")
ax.axvspan(lo-0.5,hi+0.5,color="#E69F00",alpha=0.10,zorder=0)
ax.set_xlabel("piRNA length (nt)",fontsize=9); ax.set_ylabel("% of reads (pooled, 27 samples)",fontsize=9)
ax.set_title(f"piRNA length distribution -> data-driven window {lo}-{hi} nt (FWHM, mode {mode})",fontsize=9,fontweight="bold")
ax.set_xticks(xs); ax.tick_params(labelsize=6.5)
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{U}/Fig_pirna_length_window.{e}",bbox_inches="tight")
pct.round(3).to_csv(f"{U}/SourceData_pirna_length_distribution.csv")
print("\nwrote Fig_pirna_length_window + SourceData_pirna_length_distribution.csv")
