#!/usr/bin/env python3
"""16-strain Step-4 classification (genome-anchored expression test, 4 classes incl. SNP-variant — the full
pairwise version the pangenome route skipped). Per strain (canonical order): counts + composition of
expressed-elsewhere (exact) / SNP-variant (1-3 mm) = NOT unique, vs unique conserved-but-silent / unique
strain-private-locus = GENUINELY unique by expression. Input = step4_16/{X}.step4_classified16.csv.gz
(map candidates to the 15 other genomes, mm<=3; expressed iff the locus sequence is in that strain's pool)."""
import warnings; warnings.filterwarnings("ignore")
import sys,glob; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; PG=f"{U}/pangenome_te"
d=pd.concat([pd.read_csv(f).assign(strain=f.split("/")[-1].split(".")[0]) for f in glob.glob(f"{U}/step4_16/*.step4_classified16.csv.gz")],ignore_index=True)
CANON=[s for s in STRAIN_ORDER if s in set(d.strain)]
KL=["expressed elsewhere (exact)","SNP-variant of expressed (1-3mm)","unique: conserved-but-silent","unique: strain-private locus"]
LAB=["expressed-elsewhere (not unique)","SNP-variant 1-3mm (not unique)","unique: conserved-but-silent","unique: strain-private locus"]
COL=["#9e9e9e","#E69F00","#0072B2","#7a3b9a"]
ct=pd.crosstab(d.strain,d.klass).reindex(CANON)[KL]; ct.to_csv(f"{PG}/SourceData_step4_classification16.csv")
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
fig,(ax1,ax2)=plt.subplots(2,1,figsize=(13,9.8),dpi=300); x=np.arange(len(CANON)); bw=0.21
for i,(k,lab) in enumerate(zip(KL,LAB)):
    ax1.bar(x+(i-1.5)*bw,ct[k],bw,color=COL[i],label=lab)
    for xi,v in zip(x+(i-1.5)*bw,ct[k]): ax1.text(xi,v*1.06,f"{int(v):,}",ha="center",va="bottom",fontsize=4.0,rotation=90,color=COL[i])
ax1.set_yscale("log"); ax1.set_ylim(20,ct.values.max()*2.4); ax1.set_xticks(x); ax1.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=8)
ax1.set_ylabel("strain-specific piRNA candidates (count, log)",fontsize=9); ax1.legend(fontsize=7.6,frameon=False,ncol=4,loc="upper left")
ax1.set_title("16-strain Step-4 genome-anchored classification (4 classes incl. SNP-variant) — counts per strain, canonical order",fontsize=10.4,fontweight="bold")
ax1.spines[['top','right']].set_visible(False)
prop=ct.div(ct.sum(1),axis=0); bottom=np.zeros(len(CANON))
for i,(k,lab) in enumerate(zip(KL,LAB)):
    ax2.bar(x,prop[k],0.74,bottom=bottom,color=COL[i],label=lab); bottom+=prop[k].values
gu=ct[["unique: conserved-but-silent","unique: strain-private locus"]].sum(1)/ct.sum(1)*100
ax2.set_xticks(x); ax2.set_xticklabels([s.replace("_","/") for s in CANON],rotation=45,ha="right",fontsize=8); ax2.set_ylim(0,1)
ax2.set_ylabel("class composition (proportion)",fontsize=9); ax2.set_title("Class composition per strain (text = % genuinely unique by expression)",fontsize=10,fontweight="bold")
for xi,s in zip(x,CANON): ax2.text(xi,1.01,f"{gu[s]:.0f}%",ha="center",va="bottom",fontsize=5.6,color="#333")
ax2.spines[['top','right']].set_visible(False)
fig.text(0.5,0.005,"Genome-anchored expression test: a candidate is NOT unique if another strain expresses the identical (0 mm) or a ≤3-mm SNP-variant of its locus sequence; GENUINELY unique = locus present elsewhere but silent (conserved-but-silent) or no ≤3-mm homolog (strain-private locus). "
  "All 16 strains x 15 others mapped (the full pairwise Step 4). Source: step4_16/{X}.step4_classified16.csv.gz.",ha="center",fontsize=6.2,color="#555")
fig.tight_layout(rect=[0,0.02,1,1])
for e in ("pdf","svg","png"): fig.savefig(f"{PG}/Fig_step4_classification16.{e}",bbox_inches="tight")
print("wrote Fig_step4_classification16.{png,pdf,svg} + source data"); print(ct.assign(total=ct.sum(1)).to_string())
