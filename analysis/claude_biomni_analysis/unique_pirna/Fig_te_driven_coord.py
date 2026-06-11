#!/usr/bin/env python3
"""Coordinate-based pangenome TE-driven result: fold-enrichment of each Step-4 class for having its
PRODUCTION LOCUS inside a strain-private insertion, over the random-locus expectation (private-insertion
bp / genome). The common 'expressed-elsewhere' class is hugely enriched (conserved TE piRNAs map to all
copies of active families, incl. new private insertions) while strain-private-locus piRNAs are only
modestly enriched (C57/CAST) or not (SPRET) -> private TE insertions mostly propagate CONSERVED TE
piRNAs, they do NOT preferentially create strain-private piRNA SEQUENCES. TE-driven candidate counts
(strain-private locus + locus in private insertion + TE-annotated) annotated per strain."""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
PG="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
EXP={"C57BL_6NJ":0.184,"CAST_EiJ":3.585,"SPRET_EiJ":6.961}          # random-locus null (% of genome)
TEDRV={"C57BL_6NJ":42,"CAST_EiJ":1202,"SPRET_EiJ":1259}             # coord TE-driven candidates
SCOL={"C57BL_6NJ":"#0072B2","CAST_EiJ":"#009E73","SPRET_EiJ":"#D55E00"}
PILOT=["C57BL_6NJ","CAST_EiJ","SPRET_EiJ"]
ORD=["expressed elsewhere (exact)","SNP-variant of expressed (1-3mm)","unique: conserved-but-silent","unique: strain-private locus"]
LAB=["expressed elsewhere\n(common)","SNP-variant","unique:\nconserved-silent","unique:\nstrain-private locus"]
df=pd.concat([pd.read_csv(f"{PG}/{x}.coord_byclass.csv") for x in PILOT])
fe=df.copy(); fe["fold"]=fe.apply(lambda r: r["pct"]/EXP[r["strain"]],axis=1)
piv=fe.pivot(index="klass",columns="strain",values="fold").reindex(ORD)[PILOT]
piv.to_csv(f"{PG}/TE_driven_coord_foldenrichment.csv")

plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(figsize=(8.2,4.4),dpi=300)
x=np.arange(len(ORD)); w=0.26
for i,s in enumerate(PILOT):
    ax.bar(x+(i-1)*w,piv[s].values,w,color=SCOL[s],edgecolor="white",linewidth=0.4,
           label=f"{s.replace('_','/')} (TE-driven candidates: {TEDRV[s]:,})",zorder=3)
    for xi,v in zip(x+(i-1)*w,piv[s].values):
        ax.text(xi,v*1.05,f"{v:.1f}×",ha="center",va="bottom",fontsize=5.6,color=SCOL[s],fontweight="bold")
ax.axhline(1.0,ls="--",lw=0.8,color="#444",zorder=2); ax.text(len(ORD)-0.5,1.04,"random-locus expectation (1×)",ha="right",va="bottom",fontsize=6,color="#444")
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(LAB,fontsize=7)
ax.set_ylabel("fold-enrichment: locus inside a strain-private insertion\n(observed ÷ genomic expectation)",fontsize=8)
ax.set_title("Pangenome TE-driven test: only COMMON piRNAs are strongly enriched at strain-private insertions",fontsize=8.4,fontweight="bold")
ax.legend(fontsize=6.6,frameon=False,loc="upper right")
ax.spines[['top','right']].set_visible(False)
fig.text(0.5,-0.04,"Conserved TE piRNAs (expressed in all strains) map to all copies of active families incl. new private insertions -> huge enrichment. "
  "Strain-private piRNA SEQUENCES are only modestly enriched (C57/CAST ~2.5×) or not (SPRET 0.8×): private TE insertions mostly propagate "
  "CONSERVED TE piRNAs, not novel strain-private sequences. ~1.2-1.3k TE-driven candidates/strain (LTR/ERVK, LINE/L1). PROVISIONAL — BioMNI pending.",
  ha="center",fontsize=5.0,color="#666",wrap=True)
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_te_driven_coord.{e}",bbox_inches="tight")
print(piv.round(2).to_string()); print("\nTE-driven candidates:",TEDRV)
print("wrote Fig_te_driven_coord + TE_driven_coord_foldenrichment.csv")
