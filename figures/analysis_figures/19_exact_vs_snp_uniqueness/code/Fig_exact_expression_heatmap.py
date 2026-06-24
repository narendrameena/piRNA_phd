#!/usr/bin/env python3
"""THEME 19 — expression heatmap of all EXACT-sequence genuinely-unique stage-peak (27/30 nt) piRNAs (15,118),
same style as Theme-18 Fig_unique_expression_heatmap but on the exact set. Rows = unique piRNA sequences (no
labels), ordered home stage → strain → mechanism; left strip distinguishes the 4 categories incl. the SNP-allele
(standing-variation) rows the exact definition adds. Columns = 16 strains (canonical) × E16.5/P12.5/P20.5
(tp-major). Nature-Genetics single-hue RED (white=null → dark red=high). Data: data/exact_cpm_perrep.csv.gz
(extract_exact_expression.py) + data/exact_stagepeak_classified.csv.gz."""
import sys; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import matplotlib.colors as mc
from matplotlib.patches import Patch
from strain_order import STRAIN_ORDER, WILD, WILD_COLOR
plt.rcParams.update({"font.family":"Liberation Sans","pdf.fonttype":42,"svg.fonttype":"none"})
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; T=f"{ROOT}/figures/analysis_figures/19_exact_vs_snp_uniqueness"
TPS=["16.5dpc","12.5dpp","20.5dpp"]; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; TPCOL={"16.5dpc":"#4393C3","12.5dpp":"#E8852B","20.5dpp":"#B2182B"}
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]; TR={t:i for i,t in enumerate(TPS)}; SR={s:i for i,s in enumerate(CANON)}
CBS="unique: conserved-but-silent"; SP="unique: strain-private locus"; SH="unique: stage-shifted (heterochronic)"; GU=[CBS,SP,SH]
# 4 row categories: 0 strain-private, 1 CBS-clean, 2 SNP-allele, 3 stage-shifted
CATCOL=["#7a3b9a","#0072B2","#e7298a","#009E73"]; CATLAB=["strain-private (insertion)","conserved-but-silent (clean)","SNP-allele (standing variation)","stage-shifted (heterochronic)"]
def catof(r): return 0 if r.klass_exact==SP else (3 if r.klass_exact==SH else (2 if r.was_snp_variant else 1))
g=pd.read_csv(f"{T}/data/exact_stagepeak_classified.csv.gz"); g=g[g.klass_exact.isin(GU)].copy()
g["cat"]=[catof(r) for r in g.itertuples()]; g["tr"]=g.timepoint.map(TR); g["sr"]=g.strain.map(SR)
g=g.sort_values(["tr","sr","cat"]).reset_index(drop=True)
E=pd.read_csv(f"{T}/data/exact_cpm_perrep.csv.gz",index_col=0)   # seq x 144 linear CPM, cols tp|strain.rep
strain_of={c:c.split("|")[1].rsplit(".",1)[0] for c in E.columns}; tp_of={c:c.split("|")[0] for c in E.columns}
mean48=pd.DataFrame(index=E.index)
for tp in TPS:
    for s in CANON:
        cc=[c for c in E.columns if tp_of[c]==tp and strain_of[c]==s]
        mean48[f"{tp}|{s}"]=E[cc].mean(axis=1) if cc else 0.0
mean48=np.log2(mean48+1.0)
COLS=[f"{tp}|{s}" for tp in TPS for s in CANON]
Mat=mean48.reindex(g.sequence.values)[COLS].values; vmax=np.percentile(Mat[Mat>0],99.5)
EXPR_CMAP=mc.LinearSegmentedColormap.from_list("onered",["#ffffff","#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d"])
fig=plt.figure(figsize=(15.5,16.5),dpi=300)
gs=fig.add_gridspec(2,3,width_ratios=[0.025,1,0.02],height_ratios=[0.035,1],wspace=0.012,hspace=0.012,left=0.04,right=0.93,top=0.93,bottom=0.16)
axH=fig.add_subplot(gs[1,1]); axL=fig.add_subplot(gs[1,0],sharey=axH); axT=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[1,2])
im=axH.imshow(Mat,aspect="auto",cmap=EXPR_CMAP,vmin=0,vmax=vmax,interpolation="nearest"); nS=len(CANON)
for b in (nS,2*nS): axH.axvline(b-0.5,color="#666",lw=1.2)
rb=np.cumsum(g.groupby("tr").size().reindex(range(3),fill_value=0).values)[:-1]
for r in rb: axH.axhline(r-0.5,color="#666",lw=0.8,alpha=0.8)
axL.imshow(g.cat.values.reshape(-1,1),aspect="auto",cmap=mc.ListedColormap(CATCOL),vmin=0,vmax=3,interpolation="nearest"); axL.set_xticks([]); axL.set_yticks([]); axL.set_ylabel(f"{len(g):,} exact-sequence unique piRNAs (rows; no labels)",fontsize=9)
axT.imshow(np.array([[TR[c.split('|')[0]] for c in COLS]]),aspect="auto",cmap=mc.ListedColormap([TPCOL[t] for t in TPS]),vmin=0,vmax=2,interpolation="nearest"); axT.set_xticks([]); axT.set_yticks([])
for k,tp in enumerate(TPS): axT.text((k+0.5)*nS-0.5,0,TPN[tp],ha="center",va="center",fontsize=10,fontweight="bold",color="white")
axH.set_yticks([]); axH.set_xticks(range(48)); axH.set_xticklabels([c.split("|")[1].replace("_","/") for c in COLS],rotation=90,fontsize=8.0)
for lab,c in zip(axH.get_xticklabels(),COLS): lab.set_color(WILD_COLOR if c.split("|")[1] in WILD else "#333")
axH.tick_params(axis="x",length=2,pad=1); axH.set_xlabel("16 strains (canonical order) repeated under each timepoint block",fontsize=9,labelpad=8)
fig.colorbar(im,cax=axC).set_label("expression  log2(CPM+1)",fontsize=8.5); axC.tick_params(labelsize=7)
axH.legend(handles=[Patch(facecolor=CATCOL[i],label=CATLAB[i]) for i in range(4)],title="row category (left strip)",fontsize=7.3,title_fontsize=7.6,frameon=False,loc="upper center",bbox_to_anchor=(0.5,-0.09),ncol=4)
nsnp=int((g.cat==2).sum())
fig.suptitle(f"EXACT-sequence genuinely-unique stage-peak (27/30 nt) piRNAs — {len(g):,} (incl. {nsnp:,} SNP-alleles) across 16 strains × 3 stages",fontsize=12.5,fontweight="bold",y=0.965)
fig.text(0.5,0.025,"rows = exact-unique piRNAs (no labels), ordered home stage→strain→category · columns = strain × tp (tp-major) · magenta strip = the SNP-alleles (standing variation) the exact definition adds vs SNP-aware. "
  "CPM by libsize_window, mean of 3 reps.",ha="center",fontsize=7,color="#555")
out=f"{T}/figures/Fig_exact_expression_heatmap"
for e in ("pdf","svg","png"): fig.savefig(f"{out}.{e}",bbox_inches="tight")
print("wrote",out,"| matrix",Mat.shape,"| SNP-allele rows",nsnp)
