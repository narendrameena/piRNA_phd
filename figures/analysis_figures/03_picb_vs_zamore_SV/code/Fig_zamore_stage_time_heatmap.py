#!/usr/bin/env python3
"""
Theme 03 (REFRAMED) — stage x time x strain expression heatmap (ROTATED 90°):
rows = timepoint blocks (E16.5 | P12.5 | P20.5), each = 16 strains; columns = 214
Zamore loci grouped by stage (top colour bar). Colour = log10(FPM+1), DARK = HIGH.

Developmental WAVE: prepachytene loci fire early (E16.5 rows), pachytene loci fire
late (P20.5 rows), consistently across all 16 strains.
"""
import os, subprocess, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD

plt.rcParams.update({"font.family":"Liberation Sans","font.size":7.5,"axes.linewidth":0.6,
    "pdf.fonttype":42,"svg.fonttype":"none"})

CR="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/all_strains_pangenome/combined_rebuild"
W=f"{CR}/_expr"; TPMAP={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; TPO=["E16.5","P12.5","P20.5"]
# developmental "maturation heat" ramp: gold -> pumpkin -> brick-red (matches Fig_zamore_expression_pangenome)
TP_COL={"E16.5":"#FFC84D","P12.5":"#F07F1A","P20.5":"#C1352A"}
def _pale(hexc,f=0.74):   # blend toward white for the label strip so the dark P20.5 band no longer hides red wild-strain labels
    r,g,b=plt.matplotlib.colors.to_rgb(hexc); return (r+(1-r)*f, g+(1-g)*f, b+(1-b)*f)

em=pd.read_csv(f"{CR}/all_strains_expression_matrix.csv")
stage=em.drop_duplicates("locus").set_index("locus")["stage"]
r=subprocess.run(["bedtools","intersect","-a",f"{W}/loci.bed","-b",f"{W}/clusters.bed","-wa","-wb"],capture_output=True,text=True)
rows=[]
for ln in r.stdout.strip().split("\n"):
    if not ln: continue
    f=ln.split("\t"); locus=f[3]; strain,tp,fpm=f[7].split("|")
    rows.append((locus,strain,TPMAP.get(tp,tp),float(fpm)))
ex=pd.DataFrame(rows,columns=["locus","strain","tp","fpm"]).groupby(["locus","strain","tp"])["fpm"].sum().reset_index()
strains=[s for s in STRAIN_ORDER if s in set(ex.strain)]

# column order: loci grouped by stage, sorted by mean expression
STAGE_ORDER=["Prepachytene","Hybrid","Pachytene"]; STAGE_COL={"Prepachytene":"#3BA3EC","Hybrid":"#E84FA8","Pachytene":"#10B981"}
me=ex.groupby("locus")["fpm"].mean()
ld=pd.DataFrame({"locus":stage.index,"stage":stage.values}); ld["stage"]=pd.Categorical(ld["stage"],categories=STAGE_ORDER,ordered=True)
ld["me"]=ld["locus"].map(me).fillna(0); ld=ld.sort_values(["stage","me"],ascending=[True,False]).reset_index(drop=True)
locus_order=ld.locus.tolist()
piv={t:ex[ex.tp==t].pivot_table(index="locus",columns="strain",values="fpm",aggfunc="sum") for t in TPO}

# matrix rows = (tp, strain); cols = loci
rlab=[(t,s) for t in TPO for s in strains]
M=np.full((len(rlab),len(locus_order)),np.nan)
for i,(t,s) in enumerate(rlab):
    if s in piv[t].columns: M[i,:]=np.log10(piv[t].reindex(locus_order)[s].values.astype(float)+1)

fig=plt.figure(figsize=(12.6,5.6),dpi=300)
gs=fig.add_gridspec(2,3,height_ratios=[0.035,1],width_ratios=[0.05,1,0.016],hspace=0.015,wspace=0.012)
axTop=fig.add_subplot(gs[0,1]); axL=fig.add_subplot(gs[1,0]); axH=fig.add_subplot(gs[1,1]); axCB=fig.add_subplot(gs[1,2])

# stage top bar
scol=np.array([[*plt.matplotlib.colors.to_rgb(STAGE_COL[ld.stage.iloc[i]])] for i in range(len(ld))])
axTop.imshow(scol.reshape(1,-1,3),aspect="auto"); axTop.set_xticks([]); axTop.set_yticks([])
for sp in axTop.spines.values(): sp.set_visible(False)
x0=0
for st in STAGE_ORDER:
    n=int((ld.stage==st).sum())
    if n==0: continue
    axTop.text(x0+n/2,-0.7,f"{st} (n={n})",ha="center",va="bottom",fontsize=7.0,color=STAGE_COL[st],fontweight="bold")
    x0+=n

# timepoint colour strip (left) — PALE tint so the overlaid strain labels (red=wild) stay readable
tcol=np.array([[*_pale(TP_COL[t])] for (t,s) in rlab])
axL.imshow(tcol.reshape(-1,1,3),aspect="auto"); axL.set_xticks([]); axL.set_yticks([])
for sp in axL.spines.values(): sp.set_visible(False)
for k,t in enumerate(TPO):
    axL.text(-1.1, k*len(strains)+len(strains)/2, t, rotation=90, ha="center", va="center", fontsize=8.5, fontweight="bold", color=TP_COL[t])

# main heatmap — single hue: pale = low -> solid = high
from matplotlib.colors import LinearSegmentedColormap
cmap=LinearSegmentedColormap.from_list("paleSolid",["#fbfcfe","#9ecae1","#3182bd","#08306b"]); cmap.set_bad("#d9d9d9")
im=axH.imshow(M,aspect="auto",cmap=cmap,vmin=0,vmax=np.nanpercentile(M,99))
axH.set_yticks(range(len(rlab)))
yl=axH.set_yticklabels([s.replace("_","/") for (t,s) in rlab],fontsize=4.6)
for lab,(t,s) in zip(yl,rlab): lab.set_color("#C0392B" if s in WILD else "#333")
axH.set_xticks([])
for k in (1,2): axH.axhline(k*len(strains)-0.5,color="white",lw=2.4)   # timepoint separators
x0=0
for st in STAGE_ORDER:
    n=int((ld.stage==st).sum()); x0+=n
    if x0<len(ld): axH.axvline(x0-0.5,color="white",lw=1.6)
axH.set_xlabel(f"{len(locus_order)} Zamore piRNA loci (grouped by stage; sorted by expression)",fontsize=8)

cb=fig.colorbar(im,cax=axCB); cb.set_label("log10(FPM+1)",fontsize=7.5); cb.ax.tick_params(labelsize=6.5)
fig.suptitle("Developmental expression wave of Zamore piRNA loci (rows = timepoint × strain; columns = loci by stage; dark = high)",
             fontsize=9.8,fontweight="bold",x=0.5,y=1.0)
fig.text(0.5,-0.04,
    "rows = 3 timepoint blocks × 16 strains (canonical order; red = wild) · columns = 214 loci grouped by Zamore stage (top bar) · "
    "colour = log10(FPM+1) of overlapping combined-run PICB clusters (pangenome) · prepachytene loci fire at E16.5, pachytene at P20.5 — consistently across strains",
    ha="center",fontsize=5.6,color="#666")
base=f"{CR}/Fig_zamore_stage_time_heatmap"
for ext in ("pdf","svg","png"): fig.savefig(f"{base}.{ext}",bbox_inches="tight")
print("wrote",base,"| rows(tp x strain)=",len(rlab),"cols(loci)=",len(locus_order))
