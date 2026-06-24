#!/usr/bin/env python3
"""THEME 18 — effect of WHERE the 27/30 nt length filter is applied vs DESeq2 (data-based test).
Order A = filter BEFORE DESeq2 (production); Order B = DESeq2 on full 25-32 then subset. They differ only in
size factors, dispersion trend, and the BH denominator (27/30 subset vs full 25-32). Result: filter-BEFORE
recovers ~6-7% MORE stage-peak candidates (less multiple-testing burden, length-focused normalization), with
85-93% overlap. Panels: A overall A-vs-B per tp + Jaccard; B per-strain A-vs-B (E16.5).
Data: unique_pirna/deseq16_lenfilt/{tp}.filterorder_perstrain.csv + the printed FO summaries."""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from strain_order import STRAIN_ORDER, WILD, color_wild_labels
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/deseq16_lenfilt"; T=f"{ROOT}/figures/analysis_figures/18_deseq2_stagepeak_unique_piRNA"
TPS=["16.5dpc","12.5dpp","20.5dpp"]; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
# overall A/B/shared/jaccard from the FO run (verbatim from fo_*.log)
summ={"16.5dpc":dict(A=3777,B=3551,shared=3387),"12.5dpp":dict(A=9717,B=9259,shared=9128),"20.5dpp":dict(A=2473,B=2371,shared=2225)}
CA="#0072B2"; CB="#E69F00"   # A=before blue, B=after gold
fig,(axA,axB)=plt.subplots(1,2,figsize=(14,5.6),dpi=300,gridspec_kw=dict(width_ratios=[1,1.7],wspace=0.2))
x=np.arange(len(TPS)); bw=0.38
a=[summ[t]["A"] for t in TPS]; b=[summ[t]["B"] for t in TPS]
axA.bar(x-bw/2,a,bw,color=CA,label="filter BEFORE DESeq2 (Order A, production)")
axA.bar(x+bw/2,b,bw,color=CB,label="filter AFTER DESeq2 (Order B)")
for xi,va,vb in zip(x,a,b):
    axA.text(xi-bw/2,va,f"{va:,}",ha="center",va="bottom",fontsize=8,color=CA,fontweight="bold")
    axA.text(xi+bw/2,vb,f"{vb:,}",ha="center",va="bottom",fontsize=8,color=CB,fontweight="bold")
    j=summ[TPS[xi]]["shared"]/(va+vb-summ[TPS[xi]]["shared"])
    axA.text(xi,max(va,vb)*1.09,f"+{100*(va-vb)/vb:.0f}%\nJ={j:.2f}",ha="center",fontsize=7.5,color="#444")
axA.set_xticks(x); axA.set_xticklabels([TPN[t] for t in TPS]); axA.set_ylabel("stage-peak DESeq2 candidates",fontsize=9.5)
axA.set_ylim(top=max(a)*1.22); axA.legend(fontsize=8,frameon=False,loc="upper left"); axA.spines[["top","right"]].set_visible(False)
axA.set_title("A  Filter BEFORE recovers ~6–7% more candidates\n(less BH burden; length-focused normalization)",fontsize=10.5,fontweight="bold",loc="left")
# per-strain (E16.5)
ps=pd.read_csv(f"{U}/16.5dpc.filterorder_perstrain.csv")
order=[s for s in STRAIN_ORDER if s in set(ps.strain)]; xs=np.arange(len(order)); WPOS=[i for i,s in enumerate(order) if s in WILD]
ps=ps.set_index("strain").reindex(order)
if WPOS: axB.axvspan(min(WPOS)-0.5,max(WPOS)+0.5,color="#C0392B",alpha=0.06,zorder=0)
axB.bar(xs-bw/2,ps.A,bw,color=CA,label="BEFORE"); axB.bar(xs+bw/2,ps.B,bw,color=CB,label="AFTER")
axB.set_xticks(xs); axB.set_xticklabels([s.replace("_","/") for s in order],rotation=45,ha="right",fontsize=8); color_wild_labels(axB,order)
axB.set_ylabel("E16.5 stage-peak candidates (27 nt)",fontsize=9.5); axB.legend(fontsize=8,frameon=False); axB.spines[["top","right"]].set_visible(False)
axB.set_title("B  Per strain (E16.5, 27 nt) — BEFORE ≥ AFTER across strains",fontsize=10.5,fontweight="bold",loc="left")
fig.suptitle("Effect of 27/30 nt filter order relative to DESeq2 — filtering BEFORE is more sensitive for the stage length class",fontsize=12.5,fontweight="bold",y=1.0)
fig.text(0.5,-0.02,"Order A: raw → length-filter → filterByExpr → DESeq2.  Order B: raw → filterByExpr(full 25–32) → DESeq2 → subset to 27/30.  "
  "Both + ≥2-read presence/absence. Difference = size factors, dispersion trend, BH denominator.",ha="center",fontsize=7,color="#666")
out=f"{T}/figures/Fig_filter_order"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
pd.DataFrame([dict(tp=TPN[t],before_A=summ[t]["A"],after_B=summ[t]["B"],shared=summ[t]["shared"]) for t in TPS]).to_csv(f"{T}/data/SourceData_Fig_filter_order.csv",index=False)
print("wrote",out)
