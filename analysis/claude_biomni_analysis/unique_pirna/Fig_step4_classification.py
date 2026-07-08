#!/usr/bin/env python3
"""Step 4 result: each strain's DA candidates split by EXPRESSION in the other strains (STAR genome-
anchored, <=3mm). Stacked bars per strain: not-unique (expressed exact / SNP-variant) vs genuinely
unique (conserved-but-silent / strain-private locus). Strain-private locus = TE-driven-novelty set.
"""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD, add_classical_wild_companion
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; S4=f"{U}/pangenome_te"
ORDER=["expressed elsewhere (exact)","SNP-variant (1-3mm)","low-quality: no mm0 own-genome locus","unique: conserved-but-silent","unique: strain-private locus"]
COL={"expressed elsewhere (exact)":"#9e9e9e","SNP-variant (1-3mm)":"#E69F00","low-quality: no mm0 own-genome locus":"#cdb892",
     "unique: conserved-but-silent":"#0072B2","unique: strain-private locus":"#7a3b9a"}
LAB={"expressed elsewhere (exact)":"expressed-elsewhere (not unique)","SNP-variant (1-3mm)":"SNP-variant (allelic — not unique)","low-quality: no mm0 own-genome locus":"low-quality (no mm0 own locus)",
     "unique: conserved-but-silent":"unique: conserved-but-silent","unique: strain-private locus":"unique: strain-private locus (clean)"}
d=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["strain","klass5"])   # 3-strain pilot subset of the ADOPTED 5-class (klass5, ≥2-read)
rows={X: d[d.strain==X].klass5.value_counts().reindex(ORDER,fill_value=0) for X in ["C57BL_6NJ","CAST_EiJ","SPRET_EiJ"]}
tab=pd.DataFrame(rows).T
PILOT=[s for s in STRAIN_ORDER if s in tab.index]; tab=tab.reindex(PILOT)
tab.to_csv(f"{S4}/SourceData_step4_classification.csv")

plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(figsize=(7.0,4.6),dpi=300)
x=np.arange(len(PILOT)); bottom=np.zeros(len(PILOT))
WHITE_TXT={"expressed elsewhere (exact)","SNP-variant (1-3mm)","unique: conserved-but-silent","unique: strain-private locus"}
for k in ORDER:
    vals=tab[k].values
    ax.bar(x,vals,0.6,bottom=bottom,color=COL[k],edgecolor="white",linewidth=0.6,label=LAB[k],zorder=3)
    for xi,v,b in zip(x,vals,bottom):
        if v>0:
            ax.text(xi,b+v/2,f"{int(v):,}",ha="center",va="center",fontsize=6.0,fontweight="bold",
                    color="white" if k in WHITE_TXT else "#222",zorder=4)
    bottom+=vals
for i,X in enumerate(PILOT):
    tot=tab.loc[X].sum(); gu=tab.loc[X,["unique: conserved-but-silent","unique: strain-private locus"]].sum()
    pv=tab.loc[X,"unique: strain-private locus"]
    ax.text(i,tot+4000,f"unique {gu/1000:.0f}k ({100*gu/tot:.0f}%)\nprivate-locus {pv/1000:.0f}k",
            ha="center",va="bottom",fontsize=6.6,fontweight="bold",color="#222")
ax.set_xticks(x); ax.set_xticklabels([])   # strain labels carried by the classical/wild companion below
ax.set_ylabel("strain-specific DA candidates",fontsize=9); ax.set_ylim(0,tab.sum(axis=1).max()*1.20)
ax.set_title("3-strain pilot — unique-piRNA 5-class composition (klass5, ≥2-read; pilot subset of Fig_step4_classification16)",fontsize=8.2,fontweight="bold")
ax.legend(fontsize=6.8,frameon=False,loc="upper left",ncol=1)
ax.spines[['top','right']].set_visible(False)
fig.tight_layout()
# classical(blue)/wild(orange) total-count companion per strain (subspecies colour scheme)
fig.subplots_adjust(bottom=0.34)
_cax=add_classical_wild_companion(fig,ax,PILOT,tab.sum(1).reindex(PILOT).values,gap=0.12,height_frac=0.18,ylabel="total\ncands")
_cax.set_xticks(np.arange(len(PILOT))); _cax.set_xticklabels([s.replace("_","/") for s in PILOT],fontsize=8)
for lab,s in zip(_cax.get_xticklabels(),PILOT): lab.set_color("#C0392B" if s in WILD else "#333")
_cax.set_title("classical (blue) vs wild-derived (orange) — total candidates per strain",fontsize=8,fontweight="bold",loc="left")
for e in ("pdf","svg","png"): fig.savefig(f"{S4}/Fig_step4_classification.{e}",bbox_inches="tight")
print(tab.to_string()); print("\nGenuinely unique:",
   {X:int(tab.loc[X,["unique: conserved-but-silent","unique: strain-private locus"]].sum()) for X in PILOT})
print("wrote Fig_step4_classification + SourceData_step4_classification.csv")
