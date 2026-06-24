#!/usr/bin/env python3
"""THEME 18 — DESeq2 stage-peak (27/30 nt) STRICT WITHIN-TIMEPOINT unique piRNAs (BioMNI 3/3-grounded).
Uniqueness judged ONLY within a developmental stage (pre-pachytene E16.5 vs pachytene P20.5 = different
machinery). A piRNA is non-unique only if another strain expresses it AT THE SAME stage. Genuinely-unique =
three biological mechanisms: strain-private locus (insertion), conserved-but-silent (regulatory divergence),
stage-shifted (heterochronic / developmental-timing divergence — exact seq expressed elsewhere only at a
different stage). Panels: A per-strain unique by tp (log; wild dominate); B within-tp class composition;
C the three unique MECHANISMS per stage. Data: deseq16_lenfilt/deseq_stagepeak_classified.csv.gz"""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from strain_order import STRAIN_ORDER, WILD, color_wild_labels, add_classical_wild_companion
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/deseq16_lenfilt"; T=f"{ROOT}/figures/analysis_figures/18_deseq2_stagepeak_unique_piRNA"
d=pd.read_csv(f"{U}/deseq_stagepeak_classified.csv.gz")
TPS=["16.5dpc","12.5dpp","20.5dpp"]; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; TPCOL={"16.5dpc":"#4393C3","12.5dpp":"#E8852B","20.5dpp":"#B2182B"}
CBS="unique: conserved-but-silent"; SP="unique: strain-private locus"; SH="unique: stage-shifted (heterochronic)"
SNP="SNP-variant (1-3mm, same stage)"; LQ="low-quality: no mm0 own-genome locus"; EES="expressed elsewhere (same stage)"
GU=[SP,CBS,SH]   # three unique mechanisms (insertion / regulatory / heterochronic)
KCOL={SP:"#7a3b9a",CBS:"#0072B2",SH:"#009E73",SNP:"#E69F00",LQ:"#cdb892",EES:"#999999"}
KLAB={SP:"strain-private locus (insertion)",CBS:"conserved-but-silent (regulatory divergence)",SH:"stage-shifted (heterochronic)",SNP:"SNP-variant (same-stage allele)",LQ:"low-quality",EES:"expressed elsewhere (same stage)"}
gu=d[d.klass.isin(GU)]
order=[s for s in STRAIN_ORDER if s in set(gu.strain)]; x=np.arange(len(order)); WPOS=[i for i,s in enumerate(order) if s in WILD]
fig=plt.figure(figsize=(15,11),dpi=300); gs=fig.add_gridspec(2,2,height_ratios=[1.15,0.95],hspace=0.5,wspace=0.22,left=0.07,right=0.985,top=0.9,bottom=0.1)
# A: per-strain genuinely-unique (within-tp), grouped by tp, log
axA=fig.add_subplot(gs[0,:]); bw=0.26
piv=gu.groupby(["strain","timepoint"]).size().unstack(fill_value=0).reindex(index=order,columns=TPS,fill_value=0)
if WPOS: axA.axvspan(min(WPOS)-0.5,max(WPOS)+0.5,color="#C0392B",alpha=0.06,zorder=0)
for k,tp in enumerate(TPS): axA.bar(x+(k-1)*bw,piv[tp].values,bw,color=TPCOL[tp],edgecolor="white",linewidth=0.3,label=TPN[tp],zorder=2)
axA.set_yscale("log"); axA.set_ylim(bottom=1)
axA.set_xticks(x); axA.set_xticklabels([s.replace("_","/") for s in order],rotation=45,ha="right",fontsize=8); color_wild_labels(axA,order)
axA.set_ylabel("within-tp genuinely-unique\nstage-peak piRNAs (log)",fontsize=9.5); axA.legend(title="timepoint",fontsize=8.5,title_fontsize=8.5,frameon=False,ncol=3,loc="upper left"); axA.spines[["top","right"]].set_visible(False)
axA.set_title("A  DESeq2 stage-peak (27/30 nt) WITHIN-TIMEPOINT genuinely-unique piRNAs per strain, by stage (log; wild-derived dominate)",fontsize=10.5,fontweight="bold",loc="left")
if WPOS: axA.text(np.mean(WPOS),axA.get_ylim()[1]*0.4,"wild-derived",ha="center",fontsize=8.5,fontweight="bold",color="#C0392B")
add_classical_wild_companion(fig,axA,order,piv.sum(axis=1).values,ylabel="total\n(log)")
# B: within-tp class composition per tp
axB=fig.add_subplot(gs[1,0]); xb=np.arange(len(TPS)); order_k=[SP,CBS,SH,SNP,LQ,EES]
comp=d.groupby(["timepoint","klass"]).size().unstack(fill_value=0).reindex(index=TPS,columns=order_k,fill_value=0)
bottom=np.zeros(len(TPS))
for kk in order_k:
    if comp[kk].sum()==0: continue
    axB.bar(xb,comp[kk].values,bottom=bottom,color=KCOL[kk],edgecolor="white",linewidth=0.4,label=KLAB[kk],zorder=2); bottom+=comp[kk].values
axB.set_xticks(xb); axB.set_xticklabels([TPN[t] for t in TPS],fontsize=9); axB.set_ylabel("candidates",fontsize=9.5)
axB.legend(fontsize=6.6,frameon=False,loc="upper right"); axB.spines[["top","right"]].set_visible(False)
axB.set_title("B  Within-tp class composition (EE-same-stage = 0)",fontsize=10.5,fontweight="bold",loc="left")
# C: the three unique MECHANISMS per stage
axC=fig.add_subplot(gs[1,1]); bottom=np.zeros(len(TPS))
for kk in [SP,CBS,SH]:
    vals=[gu[(gu.timepoint==t)&(gu.klass==kk)].shape[0] for t in TPS]
    axC.bar(xb,vals,bottom=bottom,color=KCOL[kk],edgecolor="white",linewidth=0.4,label=KLAB[kk],zorder=2); bottom+=np.array(vals)
for xi,t in zip(xb,TPS): axC.text(xi,bottom[xi]+max(bottom)*0.01,f"{int(bottom[xi]):,}",ha="center",va="bottom",fontsize=8,fontweight="bold")
axC.set_xticks(xb); axC.set_xticklabels([TPN[t] for t in TPS],fontsize=9); axC.set_ylabel("genuinely-unique piRNAs",fontsize=9.5)
axC.set_ylim(top=max(bottom)*1.15); axC.legend(fontsize=7,frameon=False,loc="upper right"); axC.spines[["top","right"]].set_visible(False)
axC.set_title("C  Three mechanisms of within-tp uniqueness",fontsize=10.5,fontweight="bold",loc="left")
tot=len(d); nu=len(gu); strict=d.klass.isin([CBS,SP]).sum()
fig.suptitle(f"DESeq2 stage-peak (27/30 nt) WITHIN-TIMEPOINT unique piRNAs — {nu:,} of {tot:,} ({100*nu/tot:.0f}%); 16 strains; 3 mechanisms",fontsize=12.5,fontweight="bold",y=0.965)
fig.text(0.5,0.02,"Uniqueness judged ONLY within developmental stage (BioMNI 3/3: stage-specific MILI/MIWI2 vs MIWI). Non-unique = same-stage shared / same-stage SNP-allele / low-quality. "
  f"Unique = insertion (strain-private) + regulatory (conserved-but-silent) + heterochronic (stage-shifted). Strict-sequence subset (excl. heterochronic) = {strict:,}.",ha="center",fontsize=7,color="#666")
out=f"{T}/figures/Fig_deseq_stagepeak_unique"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
piv.assign(total=piv.sum(axis=1)).to_csv(f"{T}/data/SourceData_Fig_deseq_stagepeak_unique_perstrain.csv")
comp.to_csv(f"{T}/data/SourceData_Fig_deseq_stagepeak_unique_byclass.csv")
print("wrote",out,"| within-tp unique:",nu,"of",tot,"| strict(CBS+SP):",strict)
