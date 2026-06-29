#!/usr/bin/env python3
"""THEME 18 — TE-family COVERAGE of the within-tp genuinely-unique piRNAs + the SENSE/ANTISENSE-to-TE proportion
of each piRNA sequence. Companion to Fig_sense_antisense_stageshift (which split by mechanism); here we split by
TE FAMILY. A: how many unique piRNAs derive from each TE family (genomic TE coverage of the unique repertoire).
B: per family, the proportion of piRNAs that are ANTISENSE (silencing-competent, base-pairs the TE transcript)
vs SENSE to the overlapping TE — measured w.r.t. the piRNA sequence. 50% = no strand bias.
Data: deseq16_lenfilt/deseq_stagepeak_classified.csv.gz + sense_antisense/SourceData_sense_antisense16_percand.csv.gz."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; T=f"{ROOT}/figures/analysis_figures/18_deseq2_stagepeak_unique_piRNA"
GU=["unique: conserved-but-silent","unique: strain-private locus","unique: stage-shifted (heterochronic)"]
C_ANTI="#1b9e77"; C_SENSE="#bdbdbd"   # antisense = silencing-competent (green), sense = grey
d=pd.read_csv(f"{U}/deseq16_lenfilt/deseq_stagepeak_classified.csv.gz")
oc=pd.read_csv(f"{U}/sense_antisense/SourceData_sense_antisense16_percand.csv.gz")
n_uniq=int(d.klass.isin(GU).sum())
m=d[d.klass.isin(GU)].merge(oc[["id","family","orientation"]],left_on="cand_id",right_on="id")
n_te=len(m); anti_all=100*(m.orientation=="antisense").mean()
g=m.groupby("family").agg(n=("orientation","size"),anti=("orientation",lambda s:100*(s=="antisense").mean())).sort_values("n",ascending=False)
TOP=g.head(12).iloc[::-1]                              # largest at top after barh
fig,(axA,axB)=plt.subplots(1,2,figsize=(14,6.2),dpi=300,gridspec_kw=dict(width_ratios=[1,1.15],wspace=0.45))
# A: coverage by TE family
y=np.arange(len(TOP)); cols=plt.cm.viridis(np.linspace(0.1,0.9,len(TOP)))
axA.barh(y,TOP.n.values,color=cols,edgecolor="white")
for i,(fam,r) in enumerate(TOP.iterrows()): axA.text(r.n+12,i,f"{int(r.n):,} ({100*r.n/n_te:.0f}%)",va="center",fontsize=8,color="#333")
axA.set_yticks(y); axA.set_yticklabels(TOP.index,fontsize=9); axA.set_xlabel("unique piRNAs derived from this TE family (n)",fontsize=9.5)
axA.set_xlim(0,TOP.n.max()*1.22); axA.spines[["top","right"]].set_visible(False)
axA.set_title(f"A  Genomic TE coverage of the unique piRNAs\n{n_te:,} of {n_uniq:,} ({100*n_te/n_uniq:.0f}%) overlap a TE · top 12 of {g.shape[0]} families",fontsize=10,fontweight="bold",loc="left")
# B: sense / antisense proportion per family
anti=TOP.anti.values; sense=100-anti
axB.barh(y,anti,color=C_ANTI,edgecolor="white",label="antisense to TE (silencing-competent)")
axB.barh(y,sense,left=anti,color=C_SENSE,edgecolor="white",label="sense to TE")
axB.axvline(50,ls="--",color="#333",lw=1.1)
for i,a in enumerate(anti): axB.text(min(a-2,96),i,f"{a:.0f}%",va="center",ha="right",fontsize=7.8,color="white",fontweight="bold")
axB.set_yticks(y); axB.set_yticklabels(TOP.index,fontsize=9); axB.set_xlim(0,100); axB.set_xlabel("% of piRNAs (orientation w.r.t. the overlapping TE)",fontsize=9.5)
axB.legend(fontsize=8,frameon=False,loc="lower center",bbox_to_anchor=(0.5,-0.20),ncol=2); axB.spines[["top","right"]].set_visible(False)
axB.set_title(f"B  Sense/antisense-to-TE proportion per family (w.r.t. piRNA seq)\noverall {anti_all:.0f}% antisense / {100-anti_all:.0f}% sense · dashed = no strand bias (50%)",fontsize=10,fontweight="bold",loc="left")
fig.suptitle("TE-family coverage of the within-tp unique piRNAs and their sense/antisense orientation to the TE",fontsize=12.5,fontweight="bold",y=1.0)
fig.tight_layout(rect=[0,0,1,0.95])
out=f"{T}/figures/Fig_te_coverage_sense_antisense"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
g.assign(sense=100-g.anti).round(1).to_csv(f"{T}/data/source_data/SourceData_Fig_te_coverage_sense_antisense.csv")
print(f"wrote {out} | {n_te} TE-overlapping of {n_uniq} unique; overall antisense {anti_all:.1f}%; top family {TOP.index[-1]} n={int(TOP.n.values[-1])}")
