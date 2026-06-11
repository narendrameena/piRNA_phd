#!/usr/bin/env python3
"""Figures for the sense/antisense-to-TE analysis (#3). Antisense = silencing-competent (piRNA opposite
to the TE -> can base-pair the TE transcript). (A) antisense fraction by Step-4 class x strain (unique
piRNAs vs common); (B) antisense fraction by top TE family for the genuinely-unique piRNAs. 50% = no
strand bias."""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
SA="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/sense_antisense"
STR=["C57BL_6NJ","CAST_EiJ","SPRET_EiJ"]; SCOL={"C57BL_6NJ":"#0072B2","CAST_EiJ":"#009E73","SPRET_EiJ":"#D55E00"}
CLS=["expressed elsewhere (exact)","SNP-variant of expressed (1-3mm)","unique: conserved-but-silent","unique: strain-private locus"]
CLAB=["expressed\nelsewhere","SNP-variant","unique:\nconserved-silent","unique:\nprivate-locus"]
byc=pd.concat([pd.read_csv(f"{SA}/{X}.antisense_byclass.csv").assign(strain=X) for X in STR])
pivc=byc.pivot(index="klass",columns="strain",values="antisense_pct").reindex(CLS)[STR]
pc=pd.concat([pd.read_csv(f"{SA}/{X}.sense_antisense_percandidate.csv.gz").assign(strain=X) for X in STR])
gu=pc[pc.klass.str.startswith("unique")]
fam_anti=gu.groupby("family").orientation.apply(lambda s:(s=="antisense").mean()*100)
fam_n=gu.groupby("family").size(); top=fam_n.sort_values(ascending=False).head(10).index.tolist()

plt.rcParams.update({"font.family":"Liberation Sans"})
fig,(axA,axB)=plt.subplots(1,2,figsize=(10.4,4.4),dpi=300)
x=np.arange(len(CLS)); w=0.26
for i,X in enumerate(STR):
    axA.bar(x+(i-1)*w,pivc[X].values,w,color=SCOL[X],edgecolor="white",linewidth=0.4,label=X.replace("_","/"),zorder=3)
    for xi,v in zip(x+(i-1)*w,pivc[X].values): axA.text(xi,v+0.4,f"{v:.0f}",ha="center",va="bottom",fontsize=5.6,color=SCOL[X])
axA.axhline(50,ls="--",lw=0.8,color="#444",zorder=2); axA.text(len(CLS)-0.5,50.3,"no strand bias (50%)",ha="right",va="bottom",fontsize=6,color="#444")
axA.set_xticks(x); axA.set_xticklabels(CLAB,fontsize=7); axA.set_ylim(40,62)
axA.set_ylabel("% antisense to TE (silencing-competent)",fontsize=8.5)
axA.set_title("A  Unique piRNAs are more antisense-to-TE than common",fontsize=8.6,fontweight="bold",loc="left")
axA.legend(fontsize=7,frameon=False,loc="upper left"); axA.spines[['top','right']].set_visible(False)
# B: antisense % by TE family (genuinely-unique, pooled)
vals=[fam_anti[f] for f in top]; ns=[fam_n[f] for f in top]
axB.barh(range(len(top))[::-1],vals,color="#C0392B",edgecolor="white",height=0.7,zorder=3)
for i,(v,n) in enumerate(zip(vals,ns)): axB.text(v+0.4,len(top)-1-i,f"{v:.0f}%  (n={n:,})",va="center",fontsize=6.2)
axB.axvline(50,ls="--",lw=0.8,color="#444",zorder=2)
axB.set_yticks(range(len(top))[::-1]); axB.set_yticklabels(top,fontsize=7); axB.set_xlim(40,70)
axB.set_xlabel("% antisense to TE",fontsize=8.5)
axB.set_title("B  Antisense fraction by TE family (genuinely-unique piRNAs)",fontsize=8.6,fontweight="bold",loc="left")
axB.spines[['top','right']].set_visible(False)
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{SA}/Fig_sense_antisense.{e}",bbox_inches="tight")
pivc.to_csv(f"{SA}/SourceData_antisense_byclass.csv")
print(pivc.round(1).to_string()); print("\ntop-family antisense %:"); print(pd.Series(dict(zip(top,np.round(vals,1)))).to_string())
print("wrote Fig_sense_antisense.{png,pdf,svg}")
