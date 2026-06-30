#!/usr/bin/env python3
"""Per-timepoint PCA of piRNA-sequence expression across the 3 pilot strains (thesis Fig 5.21 method:
DESeq2 size-factor normalisation -> PCA). 2 rows x 3 cols: top = all expressed piRNAs (top-500 variable,
Fig 5.21 reproduction); bottom = genuinely-unique (Step-4) piRNAs. Points = samples (3 strains x 3 reps),
coloured by strain (canonical colours); PC1/PC2 with % variance. n_features annotated per panel.
"""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
P=f"{U}/pca"
df=pd.concat([pd.read_csv(f"{P}/{t}.pca.csv") for t in ["16.5dpc","12.5dpp","20.5dpp"]],ignore_index=True)
df.to_csv(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/data/source_data/SourceData_pca_unique.csv",index=False)
TPO=["E16.5","P12.5","P20.5"]; FSO=["all_expressed","unique"]
FTITLE={"all_expressed":"All expressed piRNAs (top-500 variable)","unique":"Genuinely-unique piRNAs (klass5: CBS + strain-private)"}
SCOL={"C57BL_6NJ":"#0072B2","CAST_EiJ":"#009E73","SPRET_EiJ":"#D55E00"}
plt.rcParams.update({"font.family":"Liberation Sans"})
fig,ax=plt.subplots(2,3,figsize=(9.6,6.4),dpi=300)
for r,fs in enumerate(FSO):
    for c,tp in enumerate(TPO):
        a=ax[r,c]; sub=df[(df.feature_set==fs)&(df.timepoint==tp)]
        for s in [x for x in STRAIN_ORDER if x in set(sub.strain)]:
            ss=sub[sub.strain==s]
            a.scatter(ss.PC1,ss.PC2,s=46,color=SCOL[s],edgecolor="white",linewidth=0.6,zorder=3,
                      label=s.replace("_","/") if (r==0 and c==0) else None)
        v1=sub.pc1_var.iloc[0]; v2=sub.pc2_var.iloc[0]; nf=int(sub.n_features.iloc[0])
        a.set_xlabel(f"PC1 ({v1}%)",fontsize=8); a.set_ylabel(f"PC2 ({v2}%)",fontsize=8)
        a.set_title(f"{tp}  (n={nf:,})",fontsize=8.4,fontweight="bold")
        a.tick_params(labelsize=6); a.axhline(0,lw=0.4,color="#ddd",zorder=0); a.axvline(0,lw=0.4,color="#ddd",zorder=0)
    ax[r,0].annotate(FTITLE[fs],xy=(-0.32,0.5),xycoords="axes fraction",rotation=90,va="center",
                     fontsize=9,fontweight="bold")
ax[0,0].legend(fontsize=7,frameon=False,loc="best")
fig.suptitle("PCA of piRNA expression across pilot strains (DESeq2-normalised, per timepoint) — thesis Fig 5.21 method",
             fontsize=9.6,fontweight="bold",y=1.01)
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{U}/pangenome_te/Fig_pca_unique.{e}",bbox_inches="tight")
print(df[["feature_set","timepoint","pc1_var","pc2_var","n_features"]].drop_duplicates().to_string(index=False))
print("wrote Fig_pca_unique.{png,pdf,svg} + SourceData_pca_unique.csv")
