#!/usr/bin/env python3
"""TE families that the strain-PRIVATE-LOCUS (genuinely-unique) piRNAs originate from, per strain.
Grouped bars over the top families; annotates TE-derived fraction. CAVEAT: STAR index = main
chromosomes + MT only (no unplaced contigs), so ~20-26% of private piRNAs do not map to the own
genome and cannot be TE-annotated -> TE fraction is a LOWER BOUND (unplaced contigs are TE-rich)."""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
S4="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/step4"
fam=pd.read_csv(f"{S4}/TE_private_families.csv"); summ=pd.read_csv(f"{S4}/TE_private_summary.csv").set_index("strain")
SCOL={"C57BL_6NJ":"#0072B2","CAST_EiJ":"#009E73","SPRET_EiJ":"#D55E00"}
PILOT=[s for s in STRAIN_ORDER if s in set(fam.strain)]
piv=fam.pivot_table(index="classfam",columns="strain",values="n",fill_value=0)
top=piv.sum(axis=1).sort_values(ascending=False).head(12).index.tolist()
piv=piv.reindex(top)[PILOT]

plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(figsize=(8.4,4.4),dpi=300)
x=np.arange(len(top)); w=0.26
ymax=piv.values.max()
for i,s in enumerate(PILOT):
    vals=piv[s].values
    ax.bar(x+(i-1)*w,vals,w,color=SCOL[s],edgecolor="white",linewidth=0.4,
           label=f"{s.replace('_','/')}  (TE-derived {summ.loc[s,'TE_derived']:,} = {summ.loc[s,'TE_frac']}%)",zorder=3)
    for xi,v in zip(x+(i-1)*w,vals):
        if v>0: ax.text(xi,v+ymax*0.012,f"{int(v):,}",ha="center",va="bottom",rotation=90,fontsize=4.8,color=SCOL[s],fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels(top,rotation=40,ha="right",fontsize=7)
ax.set_ylabel("strain-private piRNAs (primary TE family)",fontsize=9); ax.set_ylim(0,ymax*1.22)
ax.set_title("TE families of origin for strain-private (genuinely-unique) piRNAs",fontsize=9.5,fontweight="bold")
ax.legend(fontsize=7,frameon=False,loc="upper right")
ax.spines[['top','right']].set_visible(False)
fig.text(0.5,-0.06,"Primary TE family = largest-overlap RepeatMasker annotation at the piRNA's own-genome locus. "
  "CAVEAT: index = main chr+MT only; ~20–26% of private piRNAs map to unplaced contigs (excluded) → TE fraction is a LOWER BOUND.",
  ha="center",fontsize=5.6,color="#666")
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{S4}/Fig_TE_private_families.{e}",bbox_inches="tight")
print("top families:\n",piv.to_string()); print("\nwrote Fig_TE_private_families.{png,pdf,svg}")
