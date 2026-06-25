#!/usr/bin/env python3
"""THEME 18 — PCA of DESeq2 stage-peak WITHIN-TP genuinely-unique piRNAs. Per tp, extract the unique sequences'
expression across the 48 libraries (stream only those rows from edger16 counts; CPM by libsize_window, log2),
PCA over samples (SVD on feature-centred matrix), colour by strain subspecies (classical/wild). Shows the
strain/subspecies structure of the within-tp unique-piRNA repertoire at each developmental stage.
Data: edger16/{tp}.{counts.tsv.gz,seqs.txt.gz,samples.tsv} + deseq16_lenfilt/deseq_stagepeak_classified.csv.gz"""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import warnings; warnings.filterwarnings("ignore")
import gzip, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from strain_order import STRAIN_ORDER, WILD, WILD_COLOR, CLASSICAL_COLOR
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
ED=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/edger16"; D=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/deseq16_lenfilt"; T=f"{ROOT}/figures/analysis_figures/18_deseq2_stagepeak_unique_piRNA"
TPS=["16.5dpc","12.5dpp","20.5dpp"]; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
GU=["unique: conserved-but-silent","unique: strain-private locus","unique: stage-shifted (heterochronic)"]
cls=pd.read_csv(f"{D}/deseq_stagepeak_classified.csv.gz")
def load_unique(tp, targets):
    counts=f"{ED}/{tp}.counts.tsv.gz"; seqs=f"{ED}/{tp}.seqs.txt.gz"
    with gzip.open(counts,"rt") as cf: header=cf.readline().rstrip("\n").split("\t")
    raw=[]
    with gzip.open(seqs,"rt") as sf, gzip.open(counts,"rt") as cf:
        cf.readline()
        for sl,cl in zip(sf,cf):
            if sl[:-1] in targets: raw.append(cl)
    M=np.array([[int(x) for x in r.rstrip("\n").split("\t")] for r in raw],dtype=float)
    return M, header
import os
_SRC=f"{T}/data/SourceData_Fig_pca_stagepeak.csv"; _DOPCA=not os.path.exists(_SRC); res=[]
for tp in (TPS if _DOPCA else []):
    targets=set(cls[(cls.timepoint==tp)&(cls.klass.isin(GU))].sequence)
    M,header=load_unique(tp,targets)
    s=pd.read_csv(f"{ED}/{tp}.samples.tsv",sep="\t").set_index("sample")
    lib=np.array([s.loc[h,"libsize_window"] for h in header],dtype=float)
    cpm=np.log2(M/lib*1e6+1.0); Xc=cpm.T - cpm.T.mean(0)
    Uu,Sv,_=np.linalg.svd(Xc,full_matrices=False); var=100*Sv**2/np.sum(Sv**2); sc=Uu*Sv
    for i,h in enumerate(header):
        res.append(dict(sample=h,strain=s.loc[h,"strain"],tp=TPN[tp],PC1=sc[i,0],PC2=sc[i,1],pc1_var=var[0],pc2_var=var[1],n_features=M.shape[0]))
    print(f"{tp}: {M.shape[0]} unique features, PC1 {var[0]:.1f}% PC2 {var[1]:.1f}%",flush=True)
pca=pd.DataFrame(res) if _DOPCA else pd.read_csv(_SRC)
fig,axes=plt.subplots(1,3,figsize=(16,5.4),dpi=300)
for ax,tp in zip(axes,["E16.5","P12.5","P20.5"]):
    sub=pca[pca.tp==tp]
    for _,r in sub.iterrows():
        ax.scatter(r.PC1,r.PC2,s=42,color=WILD_COLOR if r.strain in WILD else CLASSICAL_COLOR,
                   edgecolor="white",linewidth=0.4,zorder=3)
    # label only the wild-derived outliers (classical strains overlap at the origin -> annotate collectively)
    yr=(sub.PC2.max()-sub.PC2.min()) or 1
    cx,cy=sub.PC1.mean(),sub.PC2.mean()
    for st in sub.strain.unique():
        if st not in WILD: continue
        c=sub[sub.strain==st]; px,py=c.PC1.mean(),c.PC2.mean()
        dx,dy=px-cx,py-cy; nrm=(dx*dx+dy*dy)**0.5 or 1.0
        ax.annotate(st.replace("_","/"),(px,py),xytext=(20*dx/nrm,20*dy/nrm),textcoords="offset points",ha="center",va="center",fontsize=6.8,color="#C0392B",fontweight="bold",zorder=4,arrowprops=dict(arrowstyle="-",color="#C0392B",lw=0.4))
    cc=sub[~sub.strain.isin(WILD)]
    ax.annotate(f"{cc.strain.nunique()} classical strains\n(clustered)",(cc.PC1.median(),cc.PC2.median()),fontsize=6.5,color="#0072B2",ha="left",va="top",xytext=(16,-20),textcoords="offset points",arrowprops=dict(arrowstyle="-",color="#0072B2",lw=0.5))
    v1,v2=sub.pc1_var.iloc[0],sub.pc2_var.iloc[0]; n=sub.n_features.iloc[0]
    ax.set_xlabel(f"PC1 ({v1:.0f}%)",fontsize=9.5); ax.set_ylabel(f"PC2 ({v2:.0f}%)",fontsize=9.5)
    ax.set_title(f"{tp}  ({n:,} unique piRNAs)",fontsize=11,fontweight="bold"); ax.spines[["top","right"]].set_visible(False)
axes[0].legend(handles=[Patch(facecolor=CLASSICAL_COLOR,label="classical"),Patch(facecolor=WILD_COLOR,label="wild-derived")],fontsize=8,frameon=False,loc="best")
fig.suptitle("PCA of within-tp genuinely-unique stage-peak piRNA expression (48 libraries/tp) — wild-derived strains separate by their divergent unique repertoire",fontsize=12.5,fontweight="bold",y=1.02)
fig.tight_layout()
out=f"{T}/figures/Fig_pca_stagepeak"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
pca.to_csv(f"{T}/data/SourceData_Fig_pca_stagepeak.csv",index=False)
print("wrote",out)
