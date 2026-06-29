#!/usr/bin/env python3
"""THEME 18 — sense/antisense-to-TE orientation of the 3 within-tp unique mechanisms, incl. the new stage-shifted
(heterochronic) class. Orientation is RELATIVE TO THE TE STRAND (theme-08 convention): piRNA OPPOSITE to the TE =
ANTISENSE = silencing-competent (base-pairs the TE transcript); SAME strand = sense; 50% = no strand bias.
Orientation is per production locus (cand_self16 piRNA strand vs STRANDED TE annotation from RM .out) and is
classification-INDEPENDENT, so it is reused from the canonical theme-08 cache
(sense_antisense/SourceData_sense_antisense16_percand.csv.gz); DESeq2 stage-peak candidates are a subset of the
edgeR set by cand_id, joined here. Panels: A antisense% per mechanism (per-strain dots); B antisense% per top
TE family × mechanism.
Data: deseq16_lenfilt/deseq_stagepeak_classified.csv.gz + sense_antisense/SourceData_sense_antisense16_percand.csv.gz"""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
D=f"{U}/deseq16_lenfilt"; T=f"{ROOT}/figures/analysis_figures/18_deseq2_stagepeak_unique_piRNA"; SA=f"{U}/sense_antisense"
SP="unique: strain-private locus"; CBS="unique: conserved-but-silent"; SH="unique: stage-shifted (heterochronic)"
ORDER=[SP,CBS,SH]; LAB={SP:"strain-private\n(insertion)",CBS:"conserved-but-silent\n(regulatory)",SH:"stage-shifted\n(heterochronic)"}; KCOL={SP:"#7a3b9a",CBS:"#0072B2",SH:"#009E73"}
d=pd.read_csv(f"{D}/deseq_stagepeak_classified.csv.gz")
oc=pd.read_csv(f"{SA}/SourceData_sense_antisense16_percand.csv.gz")
m=d[d.klass.isin(ORDER)].merge(oc[["id","family","orientation"]],left_on="cand_id",right_on="id",how="inner")
print("=== coverage (TE-overlapping w/ orientation, of total per mechanism) ===")
for k in ORDER: print(f"  {LAB[k].replace(chr(10),' ')}: {m[m.klass==k].shape[0]} of {d[d.klass==k].shape[0]}")
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig,(axA,axB)=plt.subplots(1,2,figsize=(13.5,5.4),dpi=300,gridspec_kw=dict(width_ratios=[1,1.25],wspace=0.32))
rng=np.random.default_rng(0); x=np.arange(len(ORDER))
src=[]
for i,k in enumerate(ORDER):
    sub=m[m.klass==k]; pooled=100*(sub.orientation=="antisense").mean()
    axA.bar(i,pooled,0.62,color=KCOL[k],edgecolor="white",zorder=2)
    per=sub.groupby("strain").orientation.apply(lambda s:100*(s=="antisense").mean())
    axA.scatter(np.full(len(per),i)+rng.uniform(-0.16,0.16,len(per)),per.values,s=15,c="#333",edgecolor="white",lw=0.3,zorder=3,alpha=0.8)
    axA.text(i,pooled+0.8,f"{pooled:.0f}%\n(n={len(sub):,})",ha="center",va="bottom",fontsize=8,color=KCOL[k],fontweight="bold")
    src.append(dict(mechanism=LAB[k].replace("\n"," "),n_TE=len(sub),pct_antisense=round(pooled,1)))
axA.axhline(50,ls="--",lw=0.8,color="#444"); axA.text(len(ORDER)-0.5,50.6,"no strand bias (50%)",ha="right",va="bottom",fontsize=7,color="#444")
axA.set_xticks(x); axA.set_xticklabels([LAB[k] for k in ORDER],fontsize=8); axA.set_ylim(40,80)
axA.set_ylabel("% antisense to TE (silencing-competent)",fontsize=9.5); axA.spines[["top","right"]].set_visible(False)
axA.set_title("A  Stage-shifted & conserved-but-silent are antisense-biased (silencing-competent);\nstrain-private (new insertions) ≈ 50/50 dual-strand · 16 strains pooled, dots = per-strain",fontsize=8.7,fontweight="bold",loc="left")
# B: antisense% per top family x mechanism
top=m.family.value_counts().head(7).index.tolist()[::-1]; y=np.arange(len(top)); h=0.26
for j,k in enumerate(ORDER):
    vals=[100*(m[(m.klass==k)&(m.family==fm)].orientation=="antisense").mean() if m[(m.klass==k)&(m.family==fm)].shape[0]>=5 else np.nan for fm in top]
    axB.barh(y+(j-1)*h,vals,h,color=KCOL[k],label=LAB[k].replace("\n"," "),zorder=2)
axB.axvline(50,ls="--",lw=0.8,color="#444"); axB.set_yticks(y); axB.set_yticklabels(top,fontsize=8)
axB.set_xlabel("% antisense to TE",fontsize=9.5); axB.set_xlim(0,100); axB.legend(fontsize=7.3,frameon=False,loc="lower right")
axB.spines[["top","right"]].set_visible(False)
axB.set_title("B  Antisense fraction by TE family × mechanism\n(families with ≥5 loci)",fontsize=9.6,fontweight="bold",loc="left")
fig.suptitle("Sense/antisense-to-TE of the within-tp unique mechanisms — is the stage-shifted (heterochronic) class silencing-competent?",fontsize=12,fontweight="bold",y=1.0)
fig.tight_layout(rect=[0,0,1,0.95])
out=f"{T}/figures/Fig_sense_antisense_stageshift"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
pd.DataFrame(src).to_csv(f"{T}/data/source_data/SourceData_Fig_sense_antisense_stageshift.csv",index=False)
print("\nantisense% per mechanism:"); [print(f"  {s['mechanism']}: {s['pct_antisense']}% (n={s['n_TE']})") for s in src]
print("wrote",out)
