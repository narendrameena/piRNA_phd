#!/usr/bin/env python3
"""16-strain unique-piRNA classification (pangenome locus-based, final_classified.csv.gz) — the 16-strain
version of the pilot class-count figure. Per strain (CANONICAL thesis order), counts + composition of the three
classes: expressed-elsewhere (not unique), unique conserved-but-silent (expression divergence), unique
strain-private-locus. Wild-derived strains (WSB/CAST/PWK/SPRET) carry far more unique piRNAs — the divergence
signal. Colours match the four-routes figure. Source data (incl. timepoint) saved."""
import warnings; warnings.filterwarnings("ignore")
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; PG=f"{U}/pangenome_te"
d=pd.read_csv(f"{U}/unique16/final_classified.csv.gz")
CANON=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J","NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
KL=["expressed elsewhere (exact)","unique: conserved-but-silent","unique: strain-private locus"]
LAB=["expressed elsewhere (not unique)","unique: conserved-but-silent (expression)","unique: strain-private locus"]
COL=["#9e9e9e","#0072B2","#7a3b9a"]
ct=pd.crosstab(d.strain,d.klass).reindex(CANON)[KL]
pd.crosstab([d.strain,d.timepoint],d.klass).reindex(KL,axis=1).to_csv(f"{PG}/SourceData_unique16_class_breakdown.csv")
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,(ax1,ax2)=plt.subplots(2,1,figsize=(13,9.6),dpi=300); x=np.arange(len(CANON)); bw=0.27
for i,k in enumerate(KL):
    ax1.bar(x+(i-1)*bw,ct[k],bw,color=COL[i],label=LAB[i])
    for xi,v in zip(x+(i-1)*bw,ct[k]): ax1.text(xi,v*1.06,f"{int(v):,}",ha="center",va="bottom",fontsize=4.6,rotation=90,color=COL[i])
ax1.set_yscale("log"); ax1.set_ylim(30,ct.values.max()*2.2); ax1.set_xticks(x); ax1.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=8)
ax1.set_ylabel("piRNA sequences (count, log scale)",fontsize=9); ax1.legend(fontsize=8,frameon=False,ncol=3,loc="upper left")
ax1.set_title("16-strain unique-piRNA classification (pangenome locus-based) — counts per strain, canonical order",fontsize=10.5,fontweight="bold")
ax1.spines[['top','right']].set_visible(False)
prop=ct.div(ct.sum(1),axis=0); bottom=np.zeros(len(CANON))
for i,k in enumerate(KL):
    ax2.bar(x,prop[k],0.74,bottom=bottom,color=COL[i],label=LAB[i]); bottom+=prop[k].values
ax2.set_xticks(x); ax2.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=8)
ax2.set_ylim(0,1); ax2.set_ylabel("class composition (proportion)",fontsize=9); ax2.set_title("class composition per strain",fontsize=10,fontweight="bold")
ax2.spines[['top','right']].set_visible(False)
fig.text(0.5,0.005,"Genuinely unique = conserved-but-silent (locus present elsewhere, expressed only here) + strain-private-locus. Wild-derived strains (WSB/CAST/PWK/SPRET, right) carry far more unique piRNAs "
  "(divergence signal). 16-strain pangenome classification (cand_loci16 -> halLiftover presence). Source: unique16/final_classified.csv.gz; counts summed over E16.5/P12.5/P20.5 (timepoint split in source data).",ha="center",fontsize=6.2,color="#555")
fig.tight_layout(rect=[0,0.02,1,1])
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_unique16_class_breakdown.{e}",bbox_inches="tight")
print("wrote Fig_unique16_class_breakdown.{png,pdf,svg} + source data")
print(ct.assign(total=ct.sum(1)).to_string())
