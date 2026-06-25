#!/usr/bin/env python3
"""THEME 19 — EXACT-sequence vs SNP-aware definitions of strain-unique piRNAs.
EXACT: a piRNA is unique iff its EXACT sequence is strain-specific (1-3 SNP variants kept as unique).
SNP-aware (Theme 18): 1-3 mm variants of a same-stage-expressed allele elsewhere are excluded (shared-piRNA alleles).
Exact = SNP-aware + the 4,394 SNP-alleles → 15,118 vs 10,724 (+41 %). Panels: A overall+per-stage counts (exact vs
SNP-aware); B per-strain (canonical, classical/wild) with the SNP-allele addition; C exact-unique composition
(clean-unique vs SNP-allele) per stage. Data: data/exact_stagepeak_classified.csv.gz (make_exact_unique.py)."""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from strain_order import STRAIN_ORDER, WILD, color_wild_labels
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; T=f"{ROOT}/figures/analysis_figures/19_exact_vs_snp_uniqueness"
TPS=["16.5dpc","12.5dpp","20.5dpp"]; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
CBS="unique: conserved-but-silent"; SP="unique: strain-private locus"; SH="unique: stage-shifted (heterochronic)"; GU=[CBS,SP,SH]
CEX="#d95f02"; CSA="#1b9e77"; CADD="#e7298a"   # exact (orange), snp-aware (green), snp-allele add (magenta)
d=pd.read_csv(f"{T}/data/exact_stagepeak_classified.csv.gz")
d["gu_exact"]=d.klass_exact.isin(GU); d["gu_snp"]=d.klass.isin(GU)
fig,(axA,axB,axC)=plt.subplots(1,3,figsize=(16,5.4),dpi=300,gridspec_kw=dict(width_ratios=[1,1.7,1],wspace=0.3))
# A: overall + per-stage
cats=["overall"]+[TPN[t] for t in TPS]
ex=[d.gu_exact.sum()]+[d[d.timepoint==t].gu_exact.sum() for t in TPS]
sa=[d.gu_snp.sum()]+[d[d.timepoint==t].gu_snp.sum() for t in TPS]
x=np.arange(len(cats)); bw=0.38
axA.bar(x-bw/2,ex,bw,color=CEX,label="EXACT-sequence (SNP-variants kept)")
axA.bar(x+bw/2,sa,bw,color=CSA,label="SNP-aware (Theme 18; SNP-variants excluded)")
for xi,e,s in zip(x,ex,sa):
    axA.text(xi-bw/2,e,f"{e:,}",ha="center",va="bottom",fontsize=7.5,color=CEX,fontweight="bold")
    axA.text(xi+bw/2,s,f"{s:,}",ha="center",va="bottom",fontsize=7.5,color=CSA,fontweight="bold")
    axA.annotate(f"+{100*(e-s)/s:.0f}%",(xi,max(e,s)),xytext=(0,13),textcoords="offset points",ha="center",fontsize=7,color=CADD)
axA.set_xticks(x); axA.set_xticklabels(cats,fontsize=8.5); axA.set_ylabel("genuinely-unique piRNAs",fontsize=9.5)
axA.set_ylim(top=max(ex)*1.2); axA.legend(fontsize=7.6,frameon=False); axA.spines[["top","right"]].set_visible(False)
axA.set_title("A  Exact-sequence counts 41% more unique\n(+4,394 SNP-alleles vs SNP-aware)",fontsize=10.5,fontweight="bold",loc="left")
# B: per strain
order=[s for s in STRAIN_ORDER if s in set(d.strain)]; xs=np.arange(len(order)); WPOS=[i for i,s in enumerate(order) if s in WILD]
exs=[d[(d.strain==s)&d.gu_exact].shape[0] for s in order]; sas=[d[(d.strain==s)&d.gu_snp].shape[0] for s in order]
if WPOS: axB.axvspan(min(WPOS)-0.5,max(WPOS)+0.5,color="#C0392B",alpha=0.06,zorder=0)
axB.bar(xs-bw/2,exs,bw,color=CEX,label="EXACT"); axB.bar(xs+bw/2,sas,bw,color=CSA,label="SNP-aware")
axB.set_yscale("log"); axB.set_ylim(bottom=1)
axB.set_xticks(xs); axB.set_xticklabels([s.replace("_","/") for s in order],rotation=45,ha="right",fontsize=8); color_wild_labels(axB,order)
axB.set_ylabel("genuinely-unique (log)",fontsize=9.5); axB.legend(fontsize=8,frameon=False,loc="upper left"); axB.spines[["top","right"]].set_visible(False)
axB.set_title("B  Per strain — SNP-alleles add unique calls across strains (log)",fontsize=10.5,fontweight="bold",loc="left")
# C: exact-unique composition (clean vs SNP-allele) per stage
xb=np.arange(len(TPS))
clean=[d[(d.timepoint==t)&d.gu_snp].shape[0] for t in TPS]
snpadd=[d[(d.timepoint==t)&d.gu_exact&d.was_snp_variant].shape[0] for t in TPS]
axC.bar(xb,clean,0.6,color=CSA,label="clean unique (no 1-3 SNP allele elsewhere)")
axC.bar(xb,snpadd,0.6,bottom=clean,color=CADD,label="SNP-allele (1-3 mm of an expressed allele elsewhere)")
for xi,c,s in zip(xb,clean,snpadd): axC.text(xi,c+s,f"{100*s/(c+s):.0f}%\nSNP",ha="center",va="bottom",fontsize=7.5,color=CADD,fontweight="bold")
axC.set_xticks(xb); axC.set_xticklabels([TPN[t] for t in TPS],fontsize=9); axC.set_ylabel("exact-sequence genuinely-unique",fontsize=9.5)
axC.set_ylim(top=max(np.array(clean)+np.array(snpadd))*1.18); axC.legend(fontsize=7.3,frameon=False,loc="upper right"); axC.spines[["top","right"]].set_visible(False)
axC.set_title("C  What the exact definition adds:\nSNP-alleles of shared piRNAs",fontsize=10.5,fontweight="bold",loc="left")
fig.suptitle("Defining strain-unique piRNAs — EXACT-sequence vs SNP-aware (1-3 SNP variants kept vs excluded)",fontsize=12.5,fontweight="bold",y=1.0)
fig.tight_layout(rect=[0,0,1,0.95])
out=f"{T}/figures/Fig_exact_vs_snp_uniqueness"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
pd.DataFrame({"category":cats,"exact":ex,"snp_aware":sa}).to_csv(f"{T}/data/SourceData_Fig_exact_vs_snp_uniqueness.csv",index=False)
print("wrote",out,"| exact",d.gu_exact.sum(),"snp-aware",d.gu_snp.sum())
