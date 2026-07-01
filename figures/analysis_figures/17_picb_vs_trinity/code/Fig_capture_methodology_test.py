#!/usr/bin/env python3
"""TEST FIGURE — why the piRNA-capture comparison must use Trinity EXON BLOCKS, not full genomic spans.
40% of Trinity precursor contigs are multi-exon; their minimap2 genomic SPAN includes introns (5-7x the exon
footprint), so counting reads in the span captures intronic reads that are NOT piRNA-precursor reads -> Trinity
capture is falsely inflated to ~=PICB. Using exon blocks (bed12tobed6) removes introns and reveals PICB >= Trinity
(especially at P12.5, where more contigs are multi-exon). PICB is IDENTICAL in both (contiguous clusters, no introns)
-> it is the unaffected reference. Data: read_capture_span_ARTIFACT.csv (span) vs read_capture.csv (exon)."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
P="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity"
sp=pd.read_csv(f"{P}/read_capture_span_ARTIFACT.csv"); ex=pd.read_csv(f"{P}/read_capture.csv")
m=sp.merge(ex,on=["strain","tp"],suffixes=("_span","_exon"))
m["lab"]=m.strain.str.replace("_","/")+"\n"+m.tp
x=np.arange(len(m)); bw=0.26
CT="#CC6677"; CP="#117733"
fig,ax=plt.subplots(figsize=(11,5.6),dpi=300)
ax.bar(x-bw,m.pct_trin_span,bw,color=CT,alpha=0.45,hatch="///",edgecolor="white",linewidth=0.4,label="Trinity — full SPAN (artifact: +introns)",zorder=2)
ax.bar(x,m.pct_trin_exon,bw,color=CT,edgecolor="white",linewidth=0.4,label="Trinity — EXON blocks (correct)",zorder=2)
ax.bar(x+bw,m.pct_picb_exon,bw,color=CP,edgecolor="white",linewidth=0.4,label="PICB clusters (identical both ways)",zorder=2)
for xi,v in zip(x-bw,m.pct_trin_span): ax.text(xi,v+0.4,f"{v:.0f}",ha="center",fontsize=6.5,color=CT)
for xi,v in zip(x,m.pct_trin_exon): ax.text(xi,v+0.4,f"{v:.0f}",ha="center",fontsize=6.8,color=CT,fontweight="bold")
for xi,v in zip(x+bw,m.pct_picb_exon): ax.text(xi,v+0.4,f"{v:.0f}",ha="center",fontsize=6.8,color=CP,fontweight="bold")
# down-arrow span->exon (the correction)
for xi,(_,r) in zip(x,m.iterrows()):
    if r.pct_trin_span-r.pct_trin_exon>0.6:
        ax.annotate("",xy=(xi-bw/2,r.pct_trin_exon+0.2),xytext=(xi-bw/2,r.pct_trin_span-0.2),
                    arrowprops=dict(arrowstyle="->",color="#882255",lw=0.9))
        ax.text(xi-bw/2,(r.pct_trin_span+r.pct_trin_exon)/2,f"−{r.pct_trin_span-r.pct_trin_exon:.1f}pp",ha="right",va="center",fontsize=5.6,color="#882255",rotation=90)
ax.set_xticks(x); ax.set_xticklabels(m.lab,fontsize=8.5)
ax.set_ylabel("% of total piRNA captured",fontsize=10); ax.set_ylim(0,max(m.pct_picb_exon.max(),m.pct_trin_span.max())*1.25)
ax.legend(fontsize=8.3,frameon=False,loc="upper right")
ax.spines[["top","right"]].set_visible(False)
ax.set_title("TEST — Trinity capture must use EXON blocks (not intron-spanning genomic span)",fontsize=12,fontweight="bold",loc="left")
fig.text(0.5,-0.02,"40% of Trinity contigs are multi-exon → full span = 5–7× the exon footprint → captures intronic (non-piRNA) reads → falsely inflates Trinity toward PICB. "
  "Exon blocks reveal PICB ≥ Trinity (largest correction at P12.5). PICB unchanged (contiguous clusters).",ha="center",fontsize=6.8,color="#666")
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{P}/Fig_capture_methodology_test.{e}",bbox_inches="tight")
m.to_csv("/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/source_data/SourceData_Fig_capture_methodology_test.csv",index=False)
print(m[["lab","pct_trin_span","pct_trin_exon","pct_picb_exon"]].to_string(index=False)); print("wrote Fig_capture_methodology_test")
