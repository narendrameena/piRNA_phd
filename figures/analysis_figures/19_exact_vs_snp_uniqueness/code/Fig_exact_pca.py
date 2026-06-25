#!/usr/bin/env python3
"""THEME 19 — PCA of EXACT-sequence genuinely-unique stage-peak piRNA expression, same style as Theme-18
Fig_pca_stagepeak but on the exact set (15,118, incl. SNP-alleles). Per tp: that tp's exact-unique sequences ×
48 libraries (log2 CPM by libsize_window), PCA over libraries (SVD on feature-centred matrix), colour by strain
subspecies (classical/wild). Reads cached per-rep CPM (extract_exact_expression.py). Data: data/exact_cpm_perrep.csv.gz
+ data/exact_stagepeak_classified.csv.gz."""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from strain_order import WILD, WILD_COLOR, CLASSICAL_COLOR
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; T=f"{ROOT}/figures/analysis_figures/19_exact_vs_snp_uniqueness"
TPS=["16.5dpc","12.5dpp","20.5dpp"]; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
GU=["unique: conserved-but-silent","unique: strain-private locus","unique: stage-shifted (heterochronic)"]
g=pd.read_csv(f"{T}/data/exact_stagepeak_classified.csv.gz"); g=g[g.klass_exact.isin(GU)]
E=pd.read_csv(f"{T}/data/exact_cpm_perrep.csv.gz",index_col=0)
res=[]
for tp in TPS:
    cols=[c for c in E.columns if c.split("|")[0]==tp]
    seqs=[s for s in g[g.timepoint==tp].sequence.unique() if s in E.index]
    sub=np.log2(E.loc[seqs,cols].values+1.0); Xc=sub.T - sub.T.mean(0)
    Uu,Sv,_=np.linalg.svd(Xc,full_matrices=False); var=100*Sv**2/np.sum(Sv**2); sc=Uu*Sv
    for i,c in enumerate(cols):
        res.append(dict(sample=c.split("|")[1],strain=c.split("|")[1].rsplit(".",1)[0],tp=TPN[tp],PC1=sc[i,0],PC2=sc[i,1],pc1_var=var[0],pc2_var=var[1],n=len(seqs)))
    print(f"{tp}: {len(seqs)} exact-unique features, PC1 {var[0]:.1f}% PC2 {var[1]:.1f}%",flush=True)
pca=pd.DataFrame(res)
fig,axes=plt.subplots(1,3,figsize=(16,5.4),dpi=300)
for ax,tp in zip(axes,["E16.5","P12.5","P20.5"]):
    s=pca[pca.tp==tp]; yr=(s.PC2.max()-s.PC2.min()) or 1
    for _,r in s.iterrows(): ax.scatter(r.PC1,r.PC2,s=42,color=WILD_COLOR if r.strain in WILD else CLASSICAL_COLOR,edgecolor="white",linewidth=0.4,zorder=3)
    cx,cy=s.PC1.mean(),s.PC2.mean()
    # labels placed in clear space with ARROWS back to the (often coincident) points; collision-prone strains anchored in axes-fraction
    ANCH={"WSB_EiJ":(0.33,0.32),"PWK_PhJ":(0.33,0.19)}
    for st in s.strain.unique():
        if st not in WILD: continue
        c=s[s.strain==st]; px,py=c.PC1.mean(),c.PC2.mean()
        if st in ANCH:
            ax.annotate(st.replace("_","/"),xy=(px,py),xytext=ANCH[st],textcoords="axes fraction",ha="center",va="center",fontsize=6.8,color="#C0392B",fontweight="bold",zorder=4,arrowprops=dict(arrowstyle="->",color="#C0392B",lw=0.5,shrinkB=3))
        else:
            dx,dy=px-cx,py-cy; nrm=(dx*dx+dy*dy)**0.5 or 1.0
            ox,oy=(-8,30) if st=="SPRET_EiJ" else (24*dx/nrm,24*dy/nrm)   # SPRET up (avoid right-edge clip); CAST radial
            ax.annotate(st.replace("_","/"),(px,py),xytext=(ox,oy),textcoords="offset points",ha="center",va="center",fontsize=6.8,color="#C0392B",fontweight="bold",zorder=4,arrowprops=dict(arrowstyle="->",color="#C0392B",lw=0.5))
    cc=s[~s.strain.isin(WILD)]; ax.annotate(f"{cc.strain.nunique()} classical strains\n(clustered)",xy=(cc.PC1.median(),cc.PC2.median()),xytext=(0.54,0.09),textcoords="axes fraction",ha="center",va="center",fontsize=6.5,color="#0072B2",zorder=4,arrowprops=dict(arrowstyle="->",color="#0072B2",lw=0.5,shrinkB=3))
    ax.set_xlabel(f"PC1 ({s.pc1_var.iloc[0]:.0f}%)",fontsize=9.5); ax.set_ylabel(f"PC2 ({s.pc2_var.iloc[0]:.0f}%)",fontsize=9.5)
    ax.set_title(f"{tp}  ({s.n.iloc[0]:,} exact-unique piRNAs)",fontsize=11,fontweight="bold"); ax.spines[["top","right"]].set_visible(False)
axes[0].legend(handles=[Patch(facecolor=CLASSICAL_COLOR,label="classical"),Patch(facecolor=WILD_COLOR,label="wild-derived")],fontsize=8,frameon=False,loc="best")
fig.suptitle("PCA of EXACT-sequence within-tp genuinely-unique stage-peak piRNA expression (48 libraries/tp) — wild-derived strains separate",fontsize=12.5,fontweight="bold",y=1.02)
fig.tight_layout()
out=f"{T}/figures/Fig_exact_pca"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
pca.to_csv(f"{T}/data/SourceData_Fig_exact_pca.csv",index=False); print("wrote",out)
