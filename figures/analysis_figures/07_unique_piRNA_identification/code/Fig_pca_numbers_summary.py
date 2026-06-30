#!/usr/bin/env python3
"""STANDALONE summary infographic that explains, intuitively and visually, how all the strain-specific piRNA
candidate numbers fit together (≥2-read adopted) — the numbers behind Fig_pca_classes16 / Fig_pca_unique16.
Four panels, all colour-coded by class:
  A · class composition (donut, n=404,769)
  B · each class total = the sum of its 3 per-timepoint PCA panels (E16.5+P12.5+P20.5)
  C · funnel: pooled candidates -> distinct sequences -> expressed at all 3 timepoints (the 'Combined' PCA set)
  D · how much of each class survives into the 'Combined' PCA (pooled -> all-3-tp, % retained)
Sources: unique16/final_classified_clean_2read.csv.gz  +  pca16/classes_pca.csv (Combined n per class)."""
import sys, os; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from composition_cascade import draw_cascade
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
KL5=["expressed elsewhere (exact)","SNP-variant (1-3mm)","low-quality: no mm0 own-genome locus","unique: conserved-but-silent","unique: strain-private locus"]
KLAB=["expressed-elsewhere","SNP-variant","low-quality","conserved-but-silent","strain-private"]
KCOL=["#9e9e9e","#E69F00","#cdb892","#0072B2","#7a3b9a"]; WHITE={"#9e9e9e","#0072B2","#7a3b9a","#E69F00"}
TP=["16.5dpc","12.5dpp","20.5dpp"]; TPLAB=["E16.5","P12.5","P20.5"]; TPCOL=["#4393C3","#E8852B","#B2182B"]
fc=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["sequence","klass5","timepoint"])
cp=pd.read_csv(f"{U}/pca16/classes_pca.csv")
pooled=[int((fc.klass5==k).sum()) for k in KL5]
unique=[int(fc[fc.klass5==k].sequence.nunique()) for k in KL5]
pertp ={k:[int(((fc.klass5==k)&(fc.timepoint==t)).sum()) for t in TP] for k in KL5}
cbn=cp[cp.view=="combined"].groupby("klass5").n_features.first()
all3=[int(cbn.get(k,0)) for k in KL5]
TOT=sum(pooled); TOT3=sum(all3); GU=all3[3]+all3[4]; GUp=pooled[3]+pooled[4]
pd.DataFrame({"class":KLAB,"pooled":pooled,"E16.5":[pertp[k][0] for k in KL5],"P12.5":[pertp[k][1] for k in KL5],"P20.5":[pertp[k][2] for k in KL5],
  "distinct_sequences":unique,"all3_combined":all3,"pct_of_pooled":[round(100*p/TOT,2) for p in pooled],
  "pct_survive_to_combined":[round(100*a/p,1) for a,p in zip(all3,pooled)]}).to_csv(f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/data/source_data/SourceData_pca_numbers_summary.csv",index=False)
plt.rcParams.update({"font.family":"Liberation Sans"})
fig=plt.figure(figsize=(15,14.5),dpi=300)
fig.text(0.5,0.982,"How the strain-specific piRNA candidate numbers fit together  (≥2-read adopted)",ha="center",va="top",fontsize=14,fontweight="bold")
fig.text(0.5,0.965,f"{TOT:,} candidates (sequence × strain × timepoint) → {sum(unique):,} distinct sequences → {TOT3:,} expressed at all 3 timepoints (enter the pooled 'Combined' PCA)",ha="center",va="top",fontsize=9.5,color="#555")

# ---- A · composition donut ----
axA=fig.add_axes([0.035,0.772,0.255,0.158])
w,_=axA.pie(pooled,colors=KCOL,startangle=90,counterclock=False,wedgeprops=dict(width=0.40,edgecolor="white",linewidth=1.6))
for wi,n in zip(w,pooled):
    if 100*n/TOT>=6:
        ang=np.deg2rad((wi.theta1+wi.theta2)/2); r=0.80
        axA.text(r*np.cos(ang),r*np.sin(ang),f"{100*n/TOT:.0f}%",ha="center",va="center",fontsize=8.5,fontweight="bold",color="white")
axA.text(0,0,f"{TOT:,}\ncandidates",ha="center",va="center",fontsize=11,fontweight="bold")
axA.set_title("A · Class composition (pooled, 100%)",fontsize=10,fontweight="bold")

# ---- B · per-class timepoint split (sum -> class total) ----
axB=fig.add_axes([0.395,0.772,0.575,0.158])
ymax=max(pooled)
for i,k in enumerate(KL5):
    left=0
    for n,tc in zip(pertp[k],TPCOL):
        axB.barh(i,n,left=left,height=0.62,color=tc,edgecolor="white",linewidth=0.6); left+=n
    axB.text(left+ymax*0.008,i,f"{sum(pertp[k]):,}",va="center",ha="left",fontsize=7.6,fontweight="bold",color="#222")
axB.set_yticks(range(5)); axB.set_yticklabels(KLAB,fontsize=8)
for t,c in zip(axB.get_yticklabels(),KCOL): t.set_color(c); t.set_fontweight("bold")
axB.invert_yaxis(); axB.set_xlim(0,ymax*1.13); axB.set_xlabel("number of candidates",fontsize=8); axB.tick_params(labelsize=7)
for s in ("top","right"): axB.spines[s].set_visible(False)
axB.legend(handles=[Line2D([],[],marker="s",color="w",markerfacecolor=c,markersize=9,label=l) for l,c in zip(TPLAB,TPCOL)],
           fontsize=7.5,frameon=False,loc="lower right",title="timepoint",title_fontsize=8,ncol=3)
axB.set_title("B · Each class total = the sum of its 3 per-timepoint PCA panels",fontsize=10,fontweight="bold")

# ---- C · funnel pooled -> unique -> all-3-tp (stacked by class) ----
draw_cascade(fig,[0.065,0.548,0.40,0.165],
    [("pooled over all 3 timepoints",pooled),
     ("distinct sequences (de-duplicated)",unique),
     ("expressed at all 3 timepoints\n→ the pooled 'Combined' PCA",all3)],
    KCOL, title="C · Why the 'Combined' panels are smaller (funnel, by class)",
    note="'all 3 timepoints' is across DEVELOPMENT, not strains")

# ---- D · survival of each class into the Combined PCA (slopegraph, log y) ----
axD=fig.add_axes([0.58,0.553,0.37,0.160])
LNUDGE={0:1.09,2:0.91}   # pooled side: expressed-elsewhere (40,238) up, low-quality (40,011) down — near-equal values
RNUDGE={4:1.13,2:0.88}   # combined side: strain-private (2,091) up, low-quality (1,934) down
for i in range(5):
    axD.plot([0,1],[pooled[i],all3[i]],color=KCOL[i],marker="o",lw=2.2,markersize=6,zorder=3)
    axD.text(-0.04,pooled[i]*LNUDGE.get(i,1.0),f"{pooled[i]:,}",ha="right",va="center",fontsize=7.2,color=KCOL[i],fontweight="bold")
    axD.text(1.04,all3[i]*RNUDGE.get(i,1.0),f"{all3[i]:,}  ({100*all3[i]/pooled[i]:.0f}%)",ha="left",va="center",fontsize=7.2,color=KCOL[i],fontweight="bold")
axD.set_yscale("log"); axD.set_ylim(bottom=1250); axD.set_xlim(-0.5,1.9); axD.set_xticks([0,1]); axD.set_xticklabels(["pooled\n(all candidates)","'Combined' PCA\n(all 3 tp)"],fontsize=8.5)
axD.tick_params(axis="y",labelsize=7); axD.set_ylabel("candidates / sequences (log)",fontsize=8)
for s in ("top","right"): axD.spines[s].set_visible(False)
axD.set_title("D · How much of each class survives into the 'Combined' PCA",fontsize=10,fontweight="bold")
axD.text(0.985,0.965,f"genuinely-unique in 'Combined' =\nconserved-but-silent {all3[3]:,} + strain-private {all3[4]:,} = {GU:,}\n({100*GU/GUp:.0f}% of the {GUp:,} genuinely-unique pooled)",
         transform=axD.transAxes,ha="right",va="top",fontsize=7.0,style="italic",color="#444",bbox=dict(boxstyle="round,pad=0.4",fc="#f3eef7",ec="#7a3b9a",lw=0.8))

# ---- bottom: HIERARCHY tree (UNIQUE = expression test; locus new/shared = subcategory) + glossary cards under the leaves ----
notuq=TOT-GUp
fig.text(0.5,0.520,"How the 5 classes relate  —  UNIQUE is defined by EXPRESSION (1 strain, absent in the rest); the locus then splits it into new vs shared",ha="center",va="bottom",fontsize=10.5,fontweight="bold")
axt=fig.add_axes([0.03,0.345,0.94,0.155]); axt.set_xlim(0.03,0.97); axt.set_ylim(0,1); axt.axis("off"); LC="#999"
axt.text(0.5,0.93,f"{TOT:,} strain-specific candidates",ha="center",va="center",fontsize=8.8,fontweight="bold",color="white",bbox=dict(boxstyle="round,pad=0.35",fc="#444",ec="none"))
axt.plot([0.5,0.5],[0.82,0.76],color=LC,lw=1.2); axt.plot([0.277,0.757],[0.76,0.76],color=LC,lw=1.2)
axt.plot([0.277,0.277],[0.76,0.69],color=LC,lw=1.2)
axt.text(0.277,0.585,f"NOT genuinely unique\n{notuq:,} ({100*notuq/TOT:.0f}%)",ha="center",va="center",fontsize=7.3,fontweight="bold",color="#555",bbox=dict(boxstyle="round,pad=0.3",fc="#ececec",ec="none"))
axt.plot([0.277,0.277],[0.48,0.30],color=LC,lw=1.0); axt.plot([0.085,0.469],[0.30,0.30],color=LC,lw=1.0)
for a in (0.085,0.277,0.469): axt.plot([a,a],[0.30,0.05],color=LC,lw=1.0)
axt.plot([0.757,0.757],[0.76,0.69],color=LC,lw=1.2)
axt.text(0.757,0.585,f"GENUINELY UNIQUE — {GUp:,} ({100*GUp/TOT:.0f}%)\nexpressed in 1 strain, absent in the other 15",ha="center",va="center",fontsize=7.3,fontweight="bold",color="#0b3d66",bbox=dict(boxstyle="round,pad=0.3",fc="#dbeaf6",ec="none"))
axt.plot([0.757,0.757],[0.48,0.40],color=LC,lw=1.0); axt.plot([0.661,0.853],[0.40,0.40],color=LC,lw=1.0)
axt.text(0.661,0.34,"locus SHARED",ha="center",va="center",fontsize=6.6,fontweight="bold",color="#0072B2",bbox=dict(boxstyle="round,pad=0.12",fc="white",ec="none"))
axt.text(0.853,0.34,"locus NEW",ha="center",va="center",fontsize=6.6,fontweight="bold",color="#7a3b9a",bbox=dict(boxstyle="round,pad=0.12",fc="white",ec="none"))
for a in (0.661,0.853): axt.plot([a,a],[0.40,0.05],color=LC,lw=1.0)
GLOSS=["Tagged in one strain, but the IDENTICAL\nsequence is also made in other strains\n→ NOT genuinely strain-unique.",
       "The locus is in every strain; this strain's\ncopy has 1–3 DNA-letter changes (SNPs)\n→ a different-sequence allele, not a new gene.",
       "No perfect (0-mismatch) hit to the strain's\nOWN genome → origin uncertain, so it is\nkept aside as low-confidence.",
       "Locus EXISTS in other strains' genomes, but\nonly THIS strain switches it ON (silent elsewhere)\n→ unique by DIVERGENCE (same locus, new use).",
       "The piRNA locus is ABSENT from every other\nstrain's genome (often a new TE insertion)\n→ unique by a NEW LOCUS + expression."]
TAGS=["","","","   ·  GENUINELY UNIQUE","   ·  GENUINELY UNIQUE"]
gx=[0.035,0.227,0.419,0.611,0.803]
for i,(x0,lab,col,n,desc) in enumerate(zip(gx,KLAB,KCOL,pooled,GLOSS)):
    fig.text(x0,0.340,f" {lab} ",ha="left",va="top",fontsize=8.0,fontweight="bold",color="white" if col in WHITE else "#222",bbox=dict(boxstyle="round,pad=0.3",fc=col,ec="none"))
    fig.text(x0,0.306,f"{n:,} ({100*n/TOT:.1f}%){TAGS[i]}",ha="left",va="top",fontsize=6.6,color=col,fontweight="bold")
    fig.text(x0,0.285,desc,ha="left",va="top",fontsize=6.7,color="#333",linespacing=1.5)
for e in ("pdf","svg","png"): fig.savefig(f"{U}/pangenome_te/Fig_pca_numbers_summary.{e}",bbox_inches="tight")
print("wrote Fig_pca_numbers_summary.{png,pdf,svg}")
print("pooled",pooled,"\nunique",unique,"\nall3 ",all3,"\nTOT",TOT,"TOT3",TOT3,"GU",GU)
