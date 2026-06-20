#!/usr/bin/env python3
"""SPRET_EiJ unique-piRNA 5-class composition per timepoint (klass5, ≥2-read; pilot-strain view).
Stacked bars per timepoint with the adopted 5-class system: expressed-elsewhere / SNP-variant /
low-quality / unique:conserved-but-silent / unique:strain-private. The two genuinely-unique classes
(blue + purple) are the real strain-private set; grey/amber/tan are NOT truly unique (expressed exactly,
expressed as a SNP-variant, or low-quality with no clean mm0 own-genome locus). This replaces the old
naive min-mismatch (0/1-3/no-hit) split with the pangenome-anchored klass5 used genome-wide."""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; PG=f"{U}/pangenome_te"
ORDER=["expressed elsewhere (exact)","SNP-variant (1-3mm)","low-quality: no mm0 own-genome locus","unique: conserved-but-silent","unique: strain-private locus"]
COL={"expressed elsewhere (exact)":"#9e9e9e","SNP-variant (1-3mm)":"#E69F00","low-quality: no mm0 own-genome locus":"#cdb892",
     "unique: conserved-but-silent":"#0072B2","unique: strain-private locus":"#7a3b9a"}
LAB={"expressed elsewhere (exact)":"expressed-elsewhere (not unique)","SNP-variant (1-3mm)":"SNP-variant (allelic — not unique)","low-quality: no mm0 own-genome locus":"low-quality (no mm0 own locus)",
     "unique: conserved-but-silent":"unique: conserved-but-silent","unique: strain-private locus":"unique: strain-private locus (clean)"}
WHITE={"expressed elsewhere (exact)","SNP-variant (1-3mm)","unique: conserved-but-silent","unique: strain-private locus"}
TPMAP={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; TPO=["E16.5","P12.5","P20.5"]
sp=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["strain","timepoint","klass5"])
sp=sp[sp.strain=="SPRET_EiJ"].copy(); sp["tp"]=sp.timepoint.map(TPMAP)
tab=sp.groupby(["tp","klass5"]).size().unstack(fill_value=0).reindex(TPO)[ORDER]
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(figsize=(6.8,4.8),dpi=300)
x=np.arange(3); bottom=np.zeros(3)
for k in ORDER:
    vals=tab[k].values
    ax.bar(x,vals,width=0.6,bottom=bottom,color=COL[k],edgecolor="white",linewidth=0.6,label=LAB[k],zorder=3)
    for xi,v,b in zip(x,vals,bottom):
        if v>0: ax.text(xi,b+v/2,f"{int(v):,}",ha="center",va="center",fontsize=6.0,fontweight="bold",color="white" if k in WHITE else "#222",zorder=4)
    bottom+=vals
for xi,t in zip(x,TPO):
    gu=tab.loc[t,["unique: conserved-but-silent","unique: strain-private locus"]].sum(); tot=tab.loc[t].sum()
    pv=tab.loc[t,"unique: strain-private locus"]
    ax.text(xi,tot+tab.sum(1).max()*0.02,f"genuinely unique {gu/1000:.0f}k ({100*gu/tot:.0f}%)\nstrain-private locus {pv:,}",ha="center",va="bottom",fontsize=6.6,fontweight="bold",color="#222")
ax.set_xticks(x); ax.set_xticklabels(TPO,fontsize=9); ax.set_ylabel("SPRET_EiJ strain-specific candidates",fontsize=9)
ax.set_ylim(0,tab.sum(1).max()*1.22)
ax.legend(fontsize=6.8,frameon=False,loc="upper left",ncol=1)
ax.set_title("SPRET_EiJ — unique-piRNA 5-class composition per timepoint (klass5, ≥2-read)",fontsize=9.2,fontweight="bold")
fig.text(0.5,-0.04,"Genuinely unique = conserved-but-silent (blue, locus shared but silent elsewhere — divergence) + strain-private locus (purple, locus absent in all other strains — insertion). "
  "Grey/amber/tan are expressed elsewhere exactly / as a SNP-variant / low-quality (no clean mm0 own-genome locus) and are NOT truly unique. Pangenome-anchored klass5 (replaces the old naive ≤3mm split).",
  ha="center",fontsize=5.3,color="#666")
ax.spines[['top','right']].set_visible(False)
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_spret_split.{e}",bbox_inches="tight")
tab.to_csv(f"{PG}/SourceData_spret_split.csv")
print(tab.to_string()); print("\nTOTAL by class:\n", sp.klass5.value_counts().to_string())
print("wrote Fig_spret_split (klass5)")
