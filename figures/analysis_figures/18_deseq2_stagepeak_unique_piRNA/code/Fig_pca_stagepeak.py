#!/usr/bin/env python3
"""THEME 18 — PCA of DESeq2 stage-peak WITHIN-TP genuinely-unique piRNAs. Per tp, extract the unique sequences'
expression across the 48 libraries (stream only those rows from edger16 counts; CPM by libsize_window, log2),
PCA over samples (SVD on feature-centred matrix), colour by strain subspecies (classical/wild). A 4th panel pools
ALL timepoints: the union of unique sequences across all 144 libraries (48/tp x 3 tp) in one PCA, markers by
timepoint. Shows the strain/subspecies structure of the unique-piRNA repertoire per stage and pooled.
Data: edger16/{tp}.{counts.tsv.gz,seqs.txt.gz,samples.tsv} + deseq16_lenfilt/deseq_stagepeak_classified.csv.gz"""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import warnings; warnings.filterwarnings("ignore")
import os, gzip, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from strain_order import STRAIN_ORDER, WILD, WILD_COLOR, CLASSICAL_COLOR
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
ED=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/edger16"; D=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/deseq16_lenfilt"; T=f"{ROOT}/figures/analysis_figures/18_deseq2_stagepeak_unique_piRNA"
TPS=["16.5dpc","12.5dpp","20.5dpp"]; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; MK={"E16.5":"o","P12.5":"s","P20.5":"^"}
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
def load_alltp(targets):
    """union targets x 144 libraries (pool all tp); 0 where a seq is absent from a tp's count file."""
    seqlist=sorted(targets); idx={s:i for i,s in enumerate(seqlist)}; blocks=[]; cols=[]; libs=[]
    for tp in TPS:
        with gzip.open(f"{ED}/{tp}.counts.tsv.gz","rt") as cf: header=cf.readline().rstrip("\n").split("\t")
        sub=np.zeros((len(seqlist),len(header)))
        with gzip.open(f"{ED}/{tp}.seqs.txt.gz","rt") as sf, gzip.open(f"{ED}/{tp}.counts.tsv.gz","rt") as cf:
            cf.readline()
            for sl,cl in zip(sf,cf):
                q=sl[:-1]
                if q in idx: sub[idx[q]]=[int(x) for x in cl.rstrip("\n").split("\t")]
        smp=pd.read_csv(f"{ED}/{tp}.samples.tsv",sep="\t").set_index("sample")
        blocks.append(sub)
        for h in header: cols.append((h,smp.loc[h,"strain"],TPN[tp])); libs.append(smp.loc[h,"libsize_window"])
    return np.hstack(blocks), cols, np.array(libs,dtype=float)
def pca_scores(M, lib):
    cpm=np.log2(M/lib*1e6+1.0); Xc=cpm.T-cpm.T.mean(0)
    Uu,Sv,_=np.linalg.svd(Xc,full_matrices=False); var=100*Sv**2/np.sum(Sv**2); return Uu*Sv, var
# ---- per-tp PCA (cached) ----
_SRC=f"{T}/data/source_data/SourceData_Fig_pca_stagepeak.csv"; res=[]
if not os.path.exists(_SRC):
    for tp in TPS:
        targets=set(cls[(cls.timepoint==tp)&(cls.klass.isin(GU))].sequence)
        M,header=load_unique(tp,targets); s=pd.read_csv(f"{ED}/{tp}.samples.tsv",sep="\t").set_index("sample")
        lib=np.array([s.loc[h,"libsize_window"] for h in header],dtype=float); sc,var=pca_scores(M,lib)
        for i,h in enumerate(header):
            res.append(dict(sample=h,strain=s.loc[h,"strain"],tp=TPN[tp],PC1=sc[i,0],PC2=sc[i,1],pc1_var=var[0],pc2_var=var[1],n_features=M.shape[0]))
        print(f"{tp}: {M.shape[0]} unique features, PC1 {var[0]:.1f}% PC2 {var[1]:.1f}%",flush=True)
    pca=pd.DataFrame(res); pca.to_csv(_SRC,index=False)
else: pca=pd.read_csv(_SRC)
# ---- all-timepoints-pooled PCA (cached) ----
_SRC2=f"{T}/data/source_data/SourceData_Fig_pca_stagepeak_alltp.csv"
if not os.path.exists(_SRC2):
    targets_all=set(cls[cls.klass.isin(GU)].sequence); Mall,cols,liball=load_alltp(targets_all); sc,var=pca_scores(Mall,liball)
    allp=pd.DataFrame([dict(sample=h,strain=st,tp=tp,PC1=sc[i,0],PC2=sc[i,1],pc1_var=var[0],pc2_var=var[1],n_features=Mall.shape[0]) for i,(h,st,tp) in enumerate(cols)])
    allp.to_csv(_SRC2,index=False); print(f"ALL: {Mall.shape[0]} unique features x {Mall.shape[1]} libraries, PC1 {var[0]:.1f}% PC2 {var[1]:.1f}%",flush=True)
else: allp=pd.read_csv(_SRC2)
# ---- plot 2x2 ----
def wild_labels(ax,sub,anch=None):
    cx,cy=sub.PC1.mean(),sub.PC2.mean(); ANCH=anch or {"WSB_EiJ":(0.33,0.32),"PWK_PhJ":(0.33,0.19)}
    for st in sub.strain.unique():
        if st not in WILD: continue
        c=sub[sub.strain==st]; px,py=c.PC1.mean(),c.PC2.mean()
        if st in ANCH:
            ax.annotate(st.replace("_","/"),xy=(px,py),xytext=ANCH[st],textcoords="axes fraction",ha="center",va="center",fontsize=6.8,color="#C0392B",fontweight="bold",zorder=4,arrowprops=dict(arrowstyle="->",color="#C0392B",lw=0.5,shrinkB=3))
        else:
            dx,dy=px-cx,py-cy; nrm=(dx*dx+dy*dy)**0.5 or 1.0
            ox,oy=(-8,30) if st=="SPRET_EiJ" else (24*dx/nrm,24*dy/nrm)
            ax.annotate(st.replace("_","/"),(px,py),xytext=(ox,oy),textcoords="offset points",ha="center",va="center",fontsize=6.8,color="#C0392B",fontweight="bold",zorder=4,arrowprops=dict(arrowstyle="->",color="#C0392B",lw=0.5))
    cc=sub[~sub.strain.isin(WILD)]
    ax.annotate(f"{cc.strain.nunique()} classical strains\n(clustered)",xy=(cc.PC1.median(),cc.PC2.median()),xytext=(0.54,0.09),textcoords="axes fraction",ha="center",va="center",fontsize=6.5,color="#0072B2",zorder=4,arrowprops=dict(arrowstyle="->",color="#0072B2",lw=0.5,shrinkB=3))
fig,axes=plt.subplots(2,2,figsize=(13,10.6),dpi=300)
for ax,tp in zip([axes[0,0],axes[0,1],axes[1,0]],["E16.5","P12.5","P20.5"]):
    sub=pca[pca.tp==tp]
    for _,r in sub.iterrows(): ax.scatter(r.PC1,r.PC2,s=42,color=WILD_COLOR if r.strain in WILD else CLASSICAL_COLOR,edgecolor="white",linewidth=0.4,zorder=3)
    wild_labels(ax,sub)
    v1,v2=sub.pc1_var.iloc[0],sub.pc2_var.iloc[0]; n=sub.n_features.iloc[0]
    ax.set_xlabel(f"PC1 ({v1:.0f}%)",fontsize=9.5); ax.set_ylabel(f"PC2 ({v2:.0f}%)",fontsize=9.5)
    ax.set_title(f"{tp}  ({n:,} unique piRNAs, 48 libraries)",fontsize=11,fontweight="bold"); ax.spines[["top","right"]].set_visible(False)
# 4th panel: all timepoints pooled
axP=axes[1,1]
for _,r in allp.iterrows(): axP.scatter(r.PC1,r.PC2,s=34,color=WILD_COLOR if r.strain in WILD else CLASSICAL_COLOR,marker=MK.get(r.tp,"o"),edgecolor="white",linewidth=0.3,zorder=3)
wild_labels(axP,allp,{"WSB_EiJ":(0.33,0.32),"PWK_PhJ":(0.33,0.19),"CAST_EiJ":(0.10,0.82)})  # pooled: CAST centroid sits in the cluster -> anchor it
v1,v2=allp.pc1_var.iloc[0],allp.pc2_var.iloc[0]; n=allp.n_features.iloc[0]
axP.set_xlabel(f"PC1 ({v1:.0f}%)",fontsize=9.5); axP.set_ylabel(f"PC2 ({v2:.0f}%)",fontsize=9.5)
axP.set_title(f"all timepoints pooled  ({n:,} unique piRNAs, 144 libraries)",fontsize=11,fontweight="bold"); axP.spines[["top","right"]].set_visible(False)
axP.legend(handles=[Line2D([0],[0],marker=m,color="#777",ls="",ms=6,label=t) for t,m in MK.items()],fontsize=7.5,frameon=False,loc="lower right",title="timepoint",title_fontsize=7.5)
axes[0,0].legend(handles=[Patch(facecolor=CLASSICAL_COLOR,label="classical"),Patch(facecolor=WILD_COLOR,label="wild-derived")],fontsize=8,frameon=False,loc="best")
fig.suptitle("PCA of genuinely-unique stage-peak piRNA expression — per timepoint (48 libraries) and all timepoints pooled (144) — wild-derived strains separate",fontsize=12.5,fontweight="bold",y=1.0)
fig.tight_layout(rect=[0,0,1,0.97])
out=f"{T}/figures/Fig_pca_stagepeak"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
print("wrote",out)
