#!/usr/bin/env python3
"""THEME 18 — within-tp genuinely-unique stage-peak (27/30 nt) piRNA expression heatmaps SPLIT by subspecies
clade, as SEPARATE plots: CLASSICAL strains (12) and WILD-derived strains (4). Each: rows = piRNAs unique to
that clade's strains (no labels, ordered home stage→strain→class); columns = that clade's strains × E16.5/P12.5/
P20.5 (tp-major); single-hue RED (white = null → dark red = high), SHARED colour scale across both for direct
comparison. Reuses the cached seq×48 matrix from Fig_unique_expression_heatmap (no re-extraction).
NOTE: 10,170 of 10,724 (95%) genuinely-unique piRNAs are wild-derived; only 554 are classical.
Outputs: Fig_unique_expression_heatmap_classical, Fig_unique_expression_heatmap_wild."""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.colors as mc
from matplotlib.patches import Patch
from strain_order import STRAIN_ORDER, WILD, WILD_COLOR
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
D=f"{U}/deseq16_lenfilt"; T=f"{ROOT}/figures/analysis_figures/18_deseq2_stagepeak_unique_piRNA"
TPS=["16.5dpc","12.5dpp","20.5dpp"]; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; TPCOL={"16.5dpc":"#4393C3","12.5dpp":"#E8852B","20.5dpp":"#B2182B"}
ALL=[s for s in STRAIN_ORDER if s!="C57BL_6"]; TR={t:i for i,t in enumerate(TPS)}
SP="unique: strain-private locus"; CBS="unique: conserved-but-silent"; SH="unique: stage-shifted (heterochronic)"
GU=[SP,CBS,SH]; KR={SP:0,CBS:1,SH:2}; KCOL={SP:"#7a3b9a",CBS:"#0072B2",SH:"#009E73"}
g=pd.read_csv(f"{D}/deseq_stagepeak_classified.csv.gz"); g=g[g.klass.isin(GU)].copy(); g["tr"]=g.timepoint.map(TR); g["kr"]=g.klass.map(KR)
E=pd.read_csv(f"{T}/data/SourceData_Fig_unique_expression_heatmap.csv.gz",index_col=0)
vmax=np.percentile(E.values[E.values>0],99.5)   # SHARED scale across both clades
EXPR_CMAP=mc.LinearSegmentedColormap.from_list("onered",["#ffffff","#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d"])
def make(name,gstrains,fig_w,fig_h,lab_fs,bottom):
    sr={s:i for i,s in enumerate(gstrains)}
    sub=g[g.strain.isin(gstrains)].copy(); sub["sr"]=sub.strain.map(sr); sub=sub.sort_values(["tr","sr","kr"]).reset_index(drop=True)
    cols=[f"{tp}|{s}" for tp in TPS for s in gstrains]; nS=len(gstrains)
    Mat=E.reindex(sub.sequence.values)[cols].values
    fig=plt.figure(figsize=(fig_w,fig_h),dpi=300)
    gs=fig.add_gridspec(2,3,width_ratios=[0.05,1,0.03],height_ratios=[0.04,1],wspace=0.012,hspace=0.012,left=0.07,right=0.9,top=0.92,bottom=bottom)
    axH=fig.add_subplot(gs[1,1]); axL=fig.add_subplot(gs[1,0],sharey=axH); axT=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[1,2])
    im=axH.imshow(Mat,aspect="auto",cmap=EXPR_CMAP,vmin=0,vmax=vmax,interpolation="nearest")
    for k in range(1,3): axH.axvline(k*nS-0.5,color="#666",lw=1.2)
    rb=np.cumsum(sub.groupby("tr").size().reindex(range(3),fill_value=0).values)[:-1]
    for r in rb: axH.axhline(r-0.5,color="#666",lw=0.8,alpha=0.8)
    axL.imshow(sub.kr.values.reshape(-1,1),aspect="auto",cmap=mc.ListedColormap([KCOL[SP],KCOL[CBS],KCOL[SH]]),vmin=0,vmax=2,interpolation="nearest"); axL.set_xticks([]); axL.set_yticks([]); axL.set_ylabel(f"{len(sub):,} unique piRNAs (rows; home stage→strain→class; no labels)",fontsize=9)
    axT.imshow(np.array([[TR[c.split('|')[0]] for c in cols]]),aspect="auto",cmap=mc.ListedColormap([TPCOL[t] for t in TPS]),vmin=0,vmax=2,interpolation="nearest"); axT.set_xticks([]); axT.set_yticks([])
    for k,tp in enumerate(TPS): axT.text((k+0.5)*nS-0.5,0,TPN[tp],ha="center",va="center",fontsize=10,fontweight="bold",color="white")
    axH.set_yticks([]); axH.set_xticks(range(len(cols))); axH.set_xticklabels([c.split("|")[1].replace("_","/") for c in cols],rotation=90,fontsize=lab_fs)
    for lab,c in zip(axH.get_xticklabels(),cols): lab.set_color(WILD_COLOR if c.split("|")[1] in WILD else "#333")
    axH.tick_params(axis="x",length=2,pad=1); axH.set_xlabel(f"{nS} {name} strains (canonical order)  ×  E16.5 | P12.5 | P20.5 blocks",fontsize=9,labelpad=6)
    fig.colorbar(im,cax=axC).set_label("expression  log2(CPM+1)",fontsize=8.5); axC.tick_params(labelsize=7)
    axH.legend(handles=[Patch(facecolor=KCOL[SP],label="strain-private (insertion)"),Patch(facecolor=KCOL[CBS],label="conserved-but-silent (regulatory)"),Patch(facecolor=KCOL[SH],label="stage-shifted (heterochronic)")],
        title="unique mechanism (row strip)",fontsize=7.3,title_fontsize=7.6,frameon=False,loc="upper center",bbox_to_anchor=(0.5,-0.10),ncol=3)
    fig.suptitle(f"{name.capitalize()} strains — within-tp genuinely-unique stage-peak (27/30 nt) piRNA expression  ({len(sub):,} piRNAs)",fontsize=12,fontweight="bold",y=0.965)
    out=f"{T}/figures/Fig_unique_expression_heatmap_{name}"
    for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
    print(f"wrote {out} | matrix {Mat.shape}")
make("classical",[s for s in ALL if s not in WILD],12,8.5,7.5,0.22)
make("wild",[s for s in ALL if s in WILD],6.8,14,10.5,0.12)
print(f"shared vmax={vmax:.2f}")
