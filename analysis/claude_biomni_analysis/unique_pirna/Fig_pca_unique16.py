#!/usr/bin/env python3
"""16-strain per-timepoint PCA of piRNA-sequence expression (thesis Fig 5.21 / METHODS §7). 2 rows x 3 cols:
top = all expressed piRNAs (top-500 variable); bottom = genuinely-unique (final_classified). Points = samples
(16 strains x up to 3 reps), coloured by strain in canonical order. Input = pca16/{tp}.pca.csv (DESeq2-normalised
prcomp from pca_unique16.R)."""
import sys, os; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]   # 16 strains in the data
PAL=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf','#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5','#9edae5']
SCOL={s:PAL[i] for i,s in enumerate(CANON)}
MARK={s:("*" if s in WILD else "o") for s in CANON}  # wild = star
TPCOL={"E16.5":"#4393C3","P12.5":"#E8852B","P20.5":"#B2182B"}   # timepoint colour (matches the locus figures)
df=pd.concat([pd.read_csv(f"{U}/pca16/{t}.pca.csv") for t in ["16.5dpc","12.5dpp","20.5dpp"]],ignore_index=True)
df.to_csv(f"{U}/SourceData_pca_unique16.csv",index=False)
comb=pd.read_csv(f"{U}/pca16/combined.pca.csv") if os.path.exists(f"{U}/pca16/combined.pca.csv") else None   # NEW: all-timepoints-pooled expressed-piRNA PCA
TPO=["E16.5","P12.5","P20.5"]; FSO=["all_full","all_expressed","unique"]
FTITLE={"all_full":"ALL expressed piRNAs (every feature, no top-500 filter)","all_expressed":"All expressed piRNAs (top-500 variable)","unique":"Genuinely-unique piRNAs (16-strain)"}
plt.rcParams.update({"font.family":"Liberation Sans"})
fig=plt.figure(figsize=(11.5,13.9),dpi=300); gs=fig.add_gridspec(4,3,height_ratios=[1.0,1.0,1.0,1.35],hspace=0.55,wspace=0.34)
ax=np.array([[fig.add_subplot(gs[r,c]) for c in range(3)] for r in range(3)])
for r,fs in enumerate(FSO):
    for c,tp in enumerate(TPO):
        a=ax[r,c]; sub=df[(df.feature_set==fs)&(df.timepoint==tp)]
        for s in [x for x in CANON if x in set(sub.strain)]:
            ss=sub[sub.strain==s]
            a.scatter(ss.PC1,ss.PC2,s=52,color=SCOL[s],marker=MARK[s],edgecolor="black",linewidth=0.4,zorder=3,
                      label=s.replace("_","/") if (r==0 and c==0) else None)
        if len(sub):
            a.set_xlabel(f"PC1 ({sub.pc1_var.iloc[0]}%)",fontsize=8); a.set_ylabel(f"PC2 ({sub.pc2_var.iloc[0]}%)",fontsize=8)
            a.set_title(f"{tp}  (n={int(sub.n_features.iloc[0]):,})",fontsize=8.6,fontweight="bold")
        a.tick_params(labelsize=6); a.axhline(0,lw=0.4,color="#ddd",zorder=0); a.axvline(0,lw=0.4,color="#ddd",zorder=0)
    ax[r,0].annotate(FTITLE[fs],xy=(-0.34,0.5),xycoords="axes fraction",rotation=90,va="center",fontsize=9,fontweight="bold")
ax[0,0].legend(fontsize=5.6,frameon=False,loc="best",ncol=2,handletextpad=0.1,columnspacing=0.6)
# NEW bottom panel: combined all-timepoints expressed-piRNA PCA (colour = timepoint, ★ = wild) — is the dominant axis development or strain?
axb=fig.add_subplot(gs[3,:])
if comb is not None:
    for tp in TPO:
        for s in CANON:
            ss=comb[(comb.tp==tp)&(comb.strain==s)]
            axb.scatter(ss.PC1,ss.PC2,s=60,color=TPCOL[tp],marker=("*" if s in WILD else "o"),edgecolor="black",linewidth=0.4,zorder=3)
    axb.set_xlabel(f"PC1 ({comb.pc1_var.iloc[0]}%)",fontsize=8.5); axb.set_ylabel(f"PC2 ({comb.pc2_var.iloc[0]}%)",fontsize=8.5)
    axb.set_title(f"Combined — all timepoints pooled (ALL shared expressed piRNAs, no top-500 filter; n={int(comb.n_features.iloc[0]):,}; 144 libraries)   ·   colour = timepoint, ★ = wild-derived",fontsize=9.0,fontweight="bold")
    axb.legend(handles=[Line2D([],[],marker='o',color='w',markerfacecolor=TPCOL[t],markeredgecolor='k',label=t,markersize=8) for t in TPO]+[Line2D([],[],marker='*',color='w',markerfacecolor='#aaaaaa',markeredgecolor='k',label='wild-derived',markersize=12)],fontsize=7.2,frameon=False,loc="best",ncol=2)
    axb.tick_params(labelsize=6.5); axb.axhline(0,lw=0.4,color="#ddd",zorder=0); axb.axvline(0,lw=0.4,color="#ddd",zorder=0)
else:
    axb.text(0.5,0.5,"combined.pca.csv not found — run combine_pca16.R",ha="center",va="center",transform=axb.transAxes,color="#999"); axb.axis("off")
fig.suptitle("PCA of piRNA expression across all 16 strains (DESeq2-normalised) — thesis Fig 5.21 / METHODS §7 (wild-derived = star markers); bottom = all-timepoints combined",fontsize=9.8,fontweight="bold",y=1.002)
fig.tight_layout()
for e in ("pdf","svg","png"): fig.savefig(f"{U}/pangenome_te/Fig_pca_unique16.{e}",bbox_inches="tight")
print(df[["feature_set","timepoint","pc1_var","pc2_var","n_features"]].drop_duplicates().to_string(index=False))
print("wrote Fig_pca_unique16.{png,pdf,svg} + SourceData_pca_unique16.csv")
