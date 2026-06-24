#!/usr/bin/env python3
"""THEME 18 — expression heatmap of ALL within-tp genuinely-unique stage-peak (27/30 nt) piRNAs.
Rows = unique piRNA sequences (10,724 entries; NO labels), ordered by home timepoint → home strain (canonical)
→ class. Columns = 48 strain×timepoint combinations (16 strains × E16.5/P12.5/P20.5), tp-major. Colour =
expression log2(CPM+1) (CPM by libsize_window, mean over 3 reps). Left strip = unique mechanism (strain-private/
conserved-but-silent/stage-shifted); top strip = timepoint + classical/wild. The within-tp uniqueness shows as a
strain block-diagonal within each tp; the STAGE-SHIFTED (heterochronic) rows additionally light up in another
strain at a DIFFERENT stage (off-diagonal). Streams the unique-sequence rows from edger16 counts (cached).
Data: deseq16_lenfilt/deseq_stagepeak_classified.csv.gz + edger16/{tp}.{counts,seqs,samples}."""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import warnings; warnings.filterwarnings("ignore")
import os, gzip, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.colors as mc
from matplotlib.patches import Patch
from strain_order import STRAIN_ORDER, WILD, WILD_COLOR, CLASSICAL_COLOR
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
ED=f"{U}/edger16"; D=f"{U}/deseq16_lenfilt"; T=f"{ROOT}/figures/analysis_figures/18_deseq2_stagepeak_unique_piRNA"
TPS=["16.5dpc","12.5dpp","20.5dpp"]; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]; SR={s:i for i,s in enumerate(CANON)}; TR={t:i for i,t in enumerate(TPS)}
SP="unique: strain-private locus"; CBS="unique: conserved-but-silent"; SH="unique: stage-shifted (heterochronic)"
GU=[SP,CBS,SH]; KR={SP:0,CBS:1,SH:2}; KCOL={SP:"#7a3b9a",CBS:"#0072B2",SH:"#009E73"}
TPCOL={"16.5dpc":"#4393C3","12.5dpp":"#E8852B","20.5dpp":"#B2182B"}
g=pd.read_csv(f"{D}/deseq_stagepeak_classified.csv.gz"); g=g[g.klass.isin(GU)].copy()
g["tr"]=g.timepoint.map(TR); g["sr"]=g.strain.map(SR); g["kr"]=g.klass.map(KR)
g=g.sort_values(["tr","sr","kr"]).reset_index(drop=True)
COLS=[f"{tp}|{s}" for tp in TPS for s in CANON]   # 48, tp-major
CACHE=f"{T}/data/SourceData_Fig_unique_expression_heatmap.csv.gz"
if os.path.exists(CACHE):
    E=pd.read_csv(CACHE,index_col=0)
else:
    targets=set(g.sequence); E=pd.DataFrame(0.0,index=sorted(targets),columns=COLS)
    for tp in TPS:
        with gzip.open(f"{ED}/{tp}.counts.tsv.gz","rt") as cf: header=cf.readline().rstrip("\n").split("\t")
        smp=pd.read_csv(f"{ED}/{tp}.samples.tsv",sep="\t").set_index("sample")
        lib=np.array([smp.loc[h,"libsize_window"] for h in header],float); strn=[smp.loc[h,"strain"] for h in header]
        rows={}
        with gzip.open(f"{ED}/{tp}.seqs.txt.gz","rt") as sf, gzip.open(f"{ED}/{tp}.counts.tsv.gz","rt") as cf:
            cf.readline()
            for sl,cl in zip(sf,cf):
                s=sl[:-1]
                if s in targets: rows[s]=np.array(cl.rstrip("\n").split("\t"),float)
        if rows:
            M=np.vstack(list(rows.values())); cpm=M/lib*1e6
            df=pd.DataFrame(cpm,index=list(rows.keys()),columns=strn)
            mean_by_strain=df.T.groupby(level=0).mean().T   # seq x 16 strains (mean over reps)
            for s in CANON:
                if s in mean_by_strain.columns: E.loc[mean_by_strain.index,f"{tp}|{s}"]=mean_by_strain[s].values
        print(f"{tp}: {len(rows)} unique seqs extracted",flush=True)
    E=np.log2(E+1.0); E.to_csv(CACHE)
# ---- build row matrix in the sorted entry order ----
Mat=E.reindex(g.sequence.values)[COLS].values   # 10,724 x 48
vmax=np.percentile(Mat[Mat>0],99.5)
fig=plt.figure(figsize=(15.5,16.5),dpi=300)
gs=fig.add_gridspec(2,3,width_ratios=[0.025,1,0.02],height_ratios=[0.035,1],wspace=0.012,hspace=0.012,left=0.04,right=0.93,top=0.93,bottom=0.16)
axH=fig.add_subplot(gs[1,1]); axL=fig.add_subplot(gs[1,0],sharey=axH); axT=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[1,2])
# main heatmap
EXPR_CMAP=mc.LinearSegmentedColormap.from_list("onered",["#ffffff","#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d"])  # Nature-style single-hue RED: null(white) -> high(dark red)
im=axH.imshow(Mat,aspect="auto",cmap=EXPR_CMAP,vmin=0,vmax=vmax,interpolation="nearest")
nS=len(CANON)
for b in (nS,2*nS): axH.axvline(b-0.5,color="#666",lw=1.2)   # tp block separators (columns)
# row tp-block separators
rb=np.cumsum(g.groupby("tr").size().reindex(range(3),fill_value=0).values)[:-1]
for r in rb: axH.axhline(r-0.5,color="#666",lw=0.8,alpha=0.8)
axH.set_yticks([]); axH.set_xticks(range(48)); axH.set_xticklabels([c.split("|")[1].replace("_","/") for c in COLS],rotation=90,fontsize=8.0)
for lab,c in zip(axH.get_xticklabels(),COLS): lab.set_color(WILD_COLOR if c.split("|")[1] in WILD else "#333")
axH.tick_params(axis="x",length=2,pad=1)
axH.set_xlabel("strain — 16 strains (canonical order) repeated under each timepoint block",fontsize=9,labelpad=8)
# left class strip
cstrip=g.kr.values.reshape(-1,1); cl_cmap=mc.ListedColormap([KCOL[SP],KCOL[CBS],KCOL[SH]])
axL.imshow(cstrip,aspect="auto",cmap=cl_cmap,vmin=0,vmax=2,interpolation="nearest"); axL.set_xticks([]); axL.set_yticks([]); axL.set_ylabel("unique piRNAs (rows, ordered by home stage→strain→class; no labels)",fontsize=9)
# top tp strip
tstrip=np.array([[TR[c.split("|")[0]] for c in COLS]]); tp_cmap=mc.ListedColormap([TPCOL[t] for t in TPS])
axT.imshow(tstrip,aspect="auto",cmap=tp_cmap,vmin=0,vmax=2,interpolation="nearest"); axT.set_xticks([]); axT.set_yticks([])
for k,tp in enumerate(TPS): axT.text((k+0.5)*nS-0.5,0,TPN[tp],ha="center",va="center",fontsize=10,fontweight="bold",color="white")
fig.colorbar(im,cax=axC).set_label("expression  log2(CPM+1)",fontsize=8.5); axC.tick_params(labelsize=7)
# legends
leg1=axH.legend(handles=[Patch(facecolor=KCOL[SP],label="strain-private (insertion)"),Patch(facecolor=KCOL[CBS],label="conserved-but-silent (regulatory)"),Patch(facecolor=KCOL[SH],label="stage-shifted (heterochronic)")],
    title="unique mechanism (row strip)",fontsize=7.5,title_fontsize=7.8,frameon=False,loc="upper left",bbox_to_anchor=(0.0,-0.085),ncol=3)
axH.add_artist(leg1)
axH.legend(handles=[Patch(facecolor=CLASSICAL_COLOR,label="classical"),Patch(facecolor=WILD_COLOR,label="wild-derived")],title="strain x-labels",fontsize=7.5,title_fontsize=7.8,frameon=False,loc="upper right",bbox_to_anchor=(1.0,-0.085),ncol=2)
fig.suptitle(f"Expression of all {len(g):,} within-tp genuinely-unique stage-peak (27/30 nt) piRNAs across 16 strains × 3 stages",fontsize=13,fontweight="bold",y=0.965)
fig.text(0.5,0.025,"rows = unique piRNA sequences (no labels), ordered home stage→strain→class · columns = strain × timepoint (tp-major) · within-tp uniqueness → strain block-diagonal within each stage; "
  "stage-shifted (green strip) rows also light up in another strain at a DIFFERENT stage (heterochronic off-diagonal). CPM by libsize_window, mean of 3 reps.",ha="center",fontsize=7,color="#555")
out=f"{T}/figures/Fig_unique_expression_heatmap"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
print("wrote",out,"| matrix",Mat.shape,"vmax",round(vmax,2))
