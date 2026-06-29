#!/usr/bin/env python3
"""THEME 17 — DECISIVE test: of TOTAL piRNA (25-32 nt reads, multimapper-weighted), what fraction maps inside PICB
clusters vs Trinity precursor EXON blocks? A precursor method is 'better' if it explains more of the piRNA output.
Representative subset (pachytene P20.5 ×3 strains + P12.5 ×2): SPRET/CAST (wild) + C57BL_6NJ (classical).
BAM = results/STAR_srna_strain_wise (PanSN). Source: read_capture_pirna.csv (read_capture_pirna.py)."""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; P=ROOT+"/analysis/claude_biomni_analysis/picb_vs_trinity"
d=pd.read_csv(f"{P}/read_capture_pirna.csv")
d["lab"]=d.strain.str.replace("_","/")+"\n"+d.tp; CP="#117733"; CT="#CC6677"
x=np.arange(len(d)); bw=0.38
fig,ax=plt.subplots(figsize=(11,5.6),dpi=300)
ax.bar(x-bw/2,d.pct_picb,bw,color=CP,edgecolor="white",linewidth=0.4,label="in PICB clusters",zorder=2)
ax.bar(x+bw/2,d.pct_trin,bw,color=CT,edgecolor="white",linewidth=0.4,label="in Trinity precursors (RPM≥200 & RPKM≥200; exon blocks)",zorder=2)
for xi,v in zip(x-bw/2,d.pct_picb): ax.text(xi,v+0.8,f"{v:.0f}%",ha="center",fontsize=8,color=CP,fontweight="bold")
for xi,v in zip(x+bw/2,d.pct_trin): ax.text(xi,v+0.8,f"{v:.0f}%",ha="center",fontsize=8,color=CT,fontweight="bold")
# ratio annotation
for xi,(_,r) in zip(x,d.iterrows()):
    rr=r.pct_picb/max(r.pct_trin,1e-9)
    ax.text(xi,max(r.pct_picb,r.pct_trin)+6,f"{rr:.1f}×",ha="center",fontsize=7.5,color="#333",style="italic")
ax.set_xticks(x); ax.set_xticklabels(d.lab,fontsize=8.5)
ax.set_ylabel("% of TOTAL piRNA captured\n(25–32 nt reads, multimapper-weighted)",fontsize=10)
ax.set_ylim(0,max(d.pct_picb.max(),d.pct_trin.max())*1.28)
ax.legend(fontsize=9,frameon=False,loc="upper right")
ax.spines[["top","right"]].set_visible(False)
ax.set_title("DECISIVE — of total piRNA, PICB ≥ Trinity: ~tied at pachytene (P20.5), ~2× ahead at P12.5",fontsize=11.5,fontweight="bold",loc="left")
fig.text(0.5,-0.02,"of total piRNA (25–32 nt reads; samtools view -L + length filter), % mapping inside PICB clusters vs Trinity precursor EXON blocks · italic = PICB/Trinity ratio · sets overlap (not additive) · "
  "pachytene precursors are few + dominant → Trinity ties PICB there; at P12.5 piRNA output is more distributed → PICB captures ~2×.",ha="center",fontsize=6.6,color="#666")
fig.tight_layout()
out=f"{P}/Fig_picb_vs_trinity_readcapture"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
d.to_csv(ROOT+"/analysis/claude_biomni_analysis/source_data/SourceData_Fig_picb_vs_trinity_readcapture.csv",index=False)
print(d.to_string(index=False)); print("wrote",out)
