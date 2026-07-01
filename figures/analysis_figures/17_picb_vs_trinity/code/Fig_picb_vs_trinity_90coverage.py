#!/usr/bin/env python3
"""THEME 17 — do PICB & Trinity agree on the loci that produce 90% of piRNA? Per strain x tp: rank PICB clusters by
all_primary_FPM and Trinity 100/100 precursors by mean rpm; take the top set covering 90% of that method's piRNA;
test genomic overlap of the two 90%-sets. Source: cov90_per_strain_tp.csv (compute_90coverage.py).
Panel A: how FEW loci give 90% (concentration) — PICB collapses onto a handful of dominant clusters at pachytene
(P20.5), Trinity stays spread (fragmentation). Panel B: agreement — at pachytene 98% of PICB's 90%-core clusters are
also Trinity 90%-core precursors -> the DOMINANT loci agree; at pre-pachytene (dispersed) they do not."""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
P="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity"
d=pd.read_csv(f"{P}/cov90_per_strain_tp.csv"); TPO=["E16.5","P12.5","P20.5"]
CP="#117733"; CT="#CC6677"; g=d.groupby("tp")
fig,(axA,axB)=plt.subplots(1,2,figsize=(13,5.4),dpi=300); xb=np.arange(3); bw=0.38
# A: number of loci to reach 90% coverage (log), PICB vs Trinity
mp=[g.get_group(t).n90_picb.mean() for t in TPO]; sp=[g.get_group(t).n90_picb.std() for t in TPO]
mt=[g.get_group(t).n90_trin.mean() for t in TPO]; st=[g.get_group(t).n90_trin.std() for t in TPO]
axA.bar(xb-bw/2,mp,bw,yerr=sp,color=CP,capsize=2,error_kw=dict(elinewidth=0.6),label="PICB clusters")
axA.bar(xb+bw/2,mt,bw,yerr=st,color=CT,capsize=2,error_kw=dict(elinewidth=0.6),label="Trinity precursors")
axA.set_yscale("log")
for xi,v,sd,pcf in zip(xb-bw/2,mp,sp,[g.get_group(t).pct_picb_for90.mean() for t in TPO]): axA.text(xi,(v+sd)*1.15,f"{v:.0f}\n({pcf:.0f}%)",ha="center",va="bottom",fontsize=6.8,color=CP,fontweight="bold")
for xi,v,sd,pcf in zip(xb+bw/2,mt,st,[g.get_group(t).pct_trin_for90.mean() for t in TPO]): axA.text(xi,(v+sd)*1.15,f"{v:.0f}\n({pcf:.0f}%)",ha="center",va="bottom",fontsize=6.8,color=CT,fontweight="bold")
axA.set_xticks(xb); axA.set_xticklabels(TPO,fontsize=10); axA.set_ylabel("loci needed for 90% of piRNA (log; mean±SD)\n(% = fraction of that method's loci)",fontsize=9.5)
axA.set_ylim(top=max(mp+mt)*9); axA.legend(fontsize=8.5,frameon=False,loc="upper right")
axA.set_title("A  PICB concentrates 90% of piRNA into a HANDFUL of clusters at pachytene\n     (73, 2.4%); Trinity stays spread (fragmentation)",fontsize=9.8,fontweight="bold",loc="left")
axA.spines[["top","right"]].set_visible(False)
# B: agreement of the 90%-coverage sets
a1=[g.get_group(t).pct_picb90_on_trin90.mean() for t in TPO]; a1s=[g.get_group(t).pct_picb90_on_trin90.std() for t in TPO]
a2=[g.get_group(t).pct_trin90_on_picb90.mean() for t in TPO]; a2s=[g.get_group(t).pct_trin90_on_picb90.std() for t in TPO]
axB.bar(xb-bw/2,a1,bw,yerr=a1s,color="#0072B2",capsize=2,error_kw=dict(elinewidth=0.6),label="% of PICB 90%-core clusters that are also Trinity 90%-core")
axB.bar(xb+bw/2,a2,bw,yerr=a2s,color="#E69F00",capsize=2,error_kw=dict(elinewidth=0.6),label="% of Trinity 90%-core precursors that are also PICB 90%-core")
for xi,v,sd in zip(xb-bw/2,a1,a1s): axB.text(xi,v+sd+2.5,f"{v:.0f}%",ha="center",fontsize=7.5,color="#0072B2",fontweight="bold")
for xi,v,sd in zip(xb+bw/2,a2,a2s): axB.text(xi,v+sd+2.5,f"{v:.0f}%",ha="center",fontsize=7.5,color="#E69F00",fontweight="bold")
axB.set_xticks(xb); axB.set_xticklabels(TPO,fontsize=10); axB.set_ylabel("% agreement (mean±SD over 16 strains)",fontsize=9.5); axB.set_ylim(0,120)
axB.legend(fontsize=7.2,frameon=False,loc="upper left")
axB.set_title("B  At PACHYTENE the dominant loci AGREE: 98% of PICB's 90%-core\n     clusters are also Trinity 90%-core precursors (asymmetric: Trinity fragments)",fontsize=9.8,fontweight="bold",loc="left")
axB.spines[["top","right"]].set_visible(False)
fig.suptitle("Do PICB & Trinity agree on the loci producing 90% of piRNA?  (16 strains × 3 tp)",fontsize=12.5,fontweight="bold",y=1.0)
fig.text(0.5,-0.02,"PICB clusters ranked by all_primary_FPM, Trinity 100/100 precursors by mean rpm; top set covering 90% of each method's piRNA; genomic overlap of the two 90%-sets · "
  "pachytene piRNA collapses onto few dominant precursors both methods find; pre-pachytene piRNA is dispersed → no small core, low agreement.",ha="center",fontsize=6.6,color="#666")
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{P}/Fig_picb_vs_trinity_90coverage.{e}",bbox_inches="tight")
d.to_csv("/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/source_data/SourceData_Fig_picb_vs_trinity_90coverage.csv",index=False)
print(g[["n90_picb","n90_trin","pct_picb90_on_trin90","pct_trin90_on_picb90"]].mean().round(1).reindex(TPO).to_string()); print("wrote Fig_picb_vs_trinity_90coverage")
