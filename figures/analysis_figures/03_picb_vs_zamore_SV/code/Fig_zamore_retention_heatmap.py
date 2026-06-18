#!/usr/bin/env python3
"""
Theme 03 (REFRAMED) — per-locus DIVERGENCE heatmap across 16 strains (transposed:
strains = rows, loci = columns). Colour = sequence divergence (100 − % span retained
via cactus pangenome halLiftover). Retained loci stay pale; the few DIVERGENT loci
glow vibrantly and concentrate in the wild-derived strains (bottom rows).

A PowerNorm spreads the small divergences so variation is visible (most loci are
98-100% retained -> a single colour would be uninformative).
"""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from strain_order import STRAIN_ORDER, WILD

plt.rcParams.update({"font.family":"Liberation Sans","font.size":7.5,"axes.linewidth":0.6,
    "pdf.fonttype":42,"svg.fonttype":"none"})

CR="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/all_strains_pangenome/combined_rebuild"
frac=pd.read_csv(f"{CR}/zamore_fraction_lifted.csv")
em=pd.read_csv(f"{CR}/all_strains_expression_matrix.csv")
stage=em.drop_duplicates("locus").set_index("locus")["stage"]
strains=[s for s in STRAIN_ORDER if s in set(frac.strain)]

# divergence (%) = 100 - retention; rows=strains, cols=loci
R=frac.pivot_table(index="strain",columns="locus",values="fraction_lifted").reindex(index=strains)*100
DIV=100-R
STAGE_ORDER=["Prepachytene","Hybrid","Pachytene"]; STAGE_COL={"Prepachytene":"#3BA3EC","Hybrid":"#E84FA8","Pachytene":"#10B981"}
ld=pd.DataFrame({"locus":R.columns,"stage":stage.reindex(R.columns).values,"mr":R.mean(axis=0).values})
ld["stage"]=pd.Categorical(ld["stage"],categories=STAGE_ORDER,ordered=True)
ld=ld.sort_values(["stage","mr"],ascending=[True,True]).reset_index(drop=True)   # divergent first within stage
DIV=DIV[ld.locus.values]

fig=plt.figure(figsize=(11.0,4.2),dpi=300)
gs=fig.add_gridspec(2,1,height_ratios=[0.04,1],hspace=0.02)
axT=fig.add_subplot(gs[0]); axH=fig.add_subplot(gs[1])

# stage top bar
scol=np.array([[*plt.matplotlib.colors.to_rgb(STAGE_COL[ld.stage.iloc[i]])] for i in range(len(ld))])
axT.imshow(scol.reshape(1,-1,3),aspect="auto"); axT.set_xticks([]); axT.set_yticks([])
for sp in axT.spines.values(): sp.set_visible(False)
y0=0
for st in STAGE_ORDER:
    n=int((ld.stage==st).sum())
    if n==0: continue
    axT.text(y0+n/2,-0.7,f"{st} (n={n})",ha="center",va="bottom",fontsize=6.6,color=STAGE_COL[st],fontweight="bold")
    y0+=n
    if y0<len(ld): axH.axvline(y0-0.5,color="white",lw=1.4)

# single-hue: 100% retained = very pale -> divergent = solid (vibrant)
from matplotlib.colors import LinearSegmentedColormap
cmap=LinearSegmentedColormap.from_list("paleSolid",["#fbfcfe","#9ecae1","#3182bd","#08306b"])
cmap.set_bad("#e6e6e6")
im=axH.imshow(DIV.values,aspect="auto",cmap=cmap,norm=PowerNorm(gamma=0.45,vmin=0,vmax=15))
axH.set_yticks(range(len(strains)))
yl=axH.set_yticklabels([s.replace("_","/") for s in strains],fontsize=6.6)
for lab,s in zip(yl,strains): lab.set_color("#C0392B" if s in WILD else "#222")
axH.set_xticks([]); axH.set_xlabel(f"{len(ld)} Zamore piRNA loci (grouped by stage; sorted divergent→retained)",fontsize=7.5)

cb=fig.colorbar(im,ax=[axT,axH],fraction=0.02,pad=0.01,extend="max",ticks=[0,1,2,5,10,15])
cb.set_label("sequence divergence\n(100 − % span retained)",fontsize=7); cb.ax.tick_params(labelsize=6.5)
fig.suptitle("Per-locus divergence of Zamore piRNA loci across 16 strains — divergent loci (solid blue) concentrate in wild strains",
             fontsize=9.6,fontweight="bold",x=0.5,y=1.02)
fig.text(0.5,-0.10,
    "rows = 16 strains (canonical order; red = wild) · columns = 214 loci grouped by Zamore stage (top bar) · "
    "colour = sequence divergence = 100 − % of locus span projecting into the strain via cactus halLiftover (pangenome; PowerNorm spreads small values) · "
    "pale = fully retained; conserved position, divergent sequence (Yu 2021 PMID 33397987)",
    ha="center",fontsize=5.6,color="#666")
base=f"{CR}/Fig_zamore_retention_heatmap"
for ext in ("pdf","svg","png"): fig.savefig(f"{base}.{ext}",bbox_inches="tight")
print("wrote",base)
print("mean divergence by strain (high->low):"); print(DIV.mean(axis=1).sort_values(ascending=False).round(2).head(6).to_string())
