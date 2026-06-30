#!/usr/bin/env python3
"""PCA of EACH piRNA classification class SEPARATELY (the 5 klass5 classes, ≥2-read adopted) — per-timepoint AND
combined. 5 rows (classes, colour-coded) × 4 cols: E16.5 / P12.5 / P20.5 (points coloured by STRAIN, ★=wild) +
Combined (pooled over timepoints, points coloured by TIMEPOINT with soft per-timepoint convex hulls). Companion
to Fig_pca_unique16 (which is organised by feature-set). Source: pca16/classes_pca.csv (combine_pca16_classes.R)."""
import sys, os; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
from composition_cascade import draw_cascade
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]
PAL=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf','#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5','#9edae5']
SCOL={s:PAL[i] for i,s in enumerate(CANON)}; MARK={s:("*" if s in WILD else "o") for s in CANON}
TPCOL={"E16.5":"#4393C3","P12.5":"#E8852B","P20.5":"#B2182B"}; TPO=["E16.5","P12.5","P20.5"]
KL5=["expressed elsewhere (exact)","SNP-variant (1-3mm)","low-quality: no mm0 own-genome locus","unique: conserved-but-silent","unique: strain-private locus"]
KLAB5=["expressed-\nelsewhere","SNP-variant\n(allelic)","low-quality","unique:\nconserved-but-silent","unique:\nstrain-private"]
KCOL=["#9e9e9e","#E69F00","#cdb892","#0072B2","#7a3b9a"]
df=pd.read_csv(f"{U}/pca16/classes_pca.csv")
def _hull(a,sub,tp,col):
    pts=sub[sub.tp==tp][["PC1","PC2"]].values
    if len(pts)>=3:
        try:
            v=ConvexHull(pts).vertices; poly=pts[v]
            a.fill(poly[:,0],poly[:,1],color=col,alpha=0.08,zorder=1,lw=0)
            a.plot(np.append(poly[:,0],poly[0,0]),np.append(poly[:,1],poly[0,1]),color=col,alpha=0.35,lw=0.7,zorder=2)
        except Exception: pass
plt.rcParams.update({"font.family":"Liberation Sans"})
COLS=["E16.5","P12.5","P20.5","Combined"]
fig=plt.figure(figsize=(15,22.0),dpi=300); gs=fig.add_gridspec(5,4,left=0.11,right=0.985,top=0.805,bottom=0.052,hspace=0.5,wspace=0.34)  # explicit extents (deterministic); top<0.81 reserves the header band for the bar + funnel
axcomb0=None
for r,(k,klab,kcol) in enumerate(zip(KL5,KLAB5,KCOL)):
    for c,colname in enumerate(COLS):
        a=fig.add_subplot(gs[r,c])
        if r==0 and colname=="Combined": axcomb0=a
        if colname=="Combined":
            sub=df[(df.klass5==k)&(df.view=="combined")]
            for tp in TPO: _hull(a,sub,tp,TPCOL[tp])
            for tp in TPO:
                for s in CANON:
                    ss=sub[(sub.tp==tp)&(sub.strain==s)]
                    a.scatter(ss.PC1,ss.PC2,s=34,color=TPCOL[tp],marker=("*" if s in WILD else "o"),edgecolor="black",linewidth=0.35,zorder=3)
            a.set_facecolor("#fcfcfd")
        else:
            sub=df[(df.klass5==k)&(df.view=="timepoint")&(df.tp==colname)]
            for s in CANON:
                ss=sub[sub.strain==s]
                a.scatter(ss.PC1,ss.PC2,s=34,color=SCOL[s],marker=MARK[s],edgecolor="black",linewidth=0.35,zorder=3)
        if len(sub):
            a.set_xlabel(f"PC1 ({sub.pc1_var.iloc[0]}%)",fontsize=7); a.set_ylabel(f"PC2 ({sub.pc2_var.iloc[0]}%)",fontsize=7)
            a.set_title(f"{colname}  (n={int(sub.n_features.iloc[0]):,})",fontsize=7.8,fontweight="bold")
        else:
            a.text(0.5,0.5,"(<3 features)",ha="center",va="center",transform=a.transAxes,color="#bbb",fontsize=7); a.set_title(colname,fontsize=7.8,fontweight="bold")
        a.tick_params(labelsize=5.5); a.axhline(0,lw=0.4,color="#ddd",zorder=0); a.axvline(0,lw=0.4,color="#ddd",zorder=0)
        if c==0:
            a.text(-0.46,0.5,klab,transform=a.transAxes,rotation=90,va="center",ha="center",fontsize=10.5,fontweight="bold",color=kcol)
    # subtle colour band behind each row label via the row's leftmost spine colour
fig.suptitle("PCA of each piRNA classification class SEPARATELY (≥2-read adopted) — rows = class · cols = E16.5 / P12.5 / P20.5 (colour=strain, "+r"$\bigstar$"+"=wild) + Combined (colour=timepoint, per-tp hulls)",
             fontsize=10.3,fontweight="bold",y=0.994)
# tp legend in the first Combined panel; strain legend at the bottom
axc0r=axcomb0  # first Combined panel (expressed-elsewhere row)
axc0r.legend(handles=[Line2D([],[],marker='o',color='w',markerfacecolor=TPCOL[t],markeredgecolor='k',label=t,markersize=6) for t in TPO]+[Line2D([],[],marker='*',color='w',markerfacecolor='#aaaaaa',markeredgecolor='k',label='wild',markersize=8)],fontsize=5.6,frameon=False,loc="best",title="timepoint",title_fontsize=6)
strain_handles=[Line2D([],[],marker=MARK[s],color="w",markerfacecolor=SCOL[s],markeredgecolor="black",markeredgewidth=0.4,markersize=6,label=s.replace("_","/")) for s in CANON]
fig.legend(handles=strain_handles,loc="lower center",bbox_to_anchor=(0.5,-0.012),ncol=8,fontsize=7,frameon=False,title="strain (markers in the E16.5 / P12.5 / P20.5 columns)",title_fontsize=7.5)
# ---- TOP header band (above the grid; grid top is fixed at 0.805): composition bar + funnel explaining pooled-vs-'Combined' counts ----
fc2=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["klass5","sequence"])
BARLAB=["expressed-elsewhere","SNP-variant","low-quality","conserved-but-silent","strain-private"]; WHITE_TXT={"#9e9e9e","#0072B2","#7a3b9a"}
cnts=fc2.klass5.value_counts(); tot=int(cnts.sum())
fig.text(0.11,0.972,f"Composition of strain-specific piRNA candidates by class (≥2-read adopted; n={tot:,}; 100% scale)  ·  segment colours = row labels",ha="left",va="bottom",fontsize=8.2,fontweight="bold",color="#222")
axbar=fig.add_axes([0.11,0.902,0.835,0.064]); left=0   # title above the bar; segment labels + GENUINELY-UNIQUE bracket BELOW it (keeps the long suptitle clear)
for k,lab,col in zip(KL5,BARLAB,KCOL):
    n=int(cnts.get(k,0)); frac=100*n/tot; cx=left+n/2
    axbar.barh(0,n,left=left,height=0.45,color=col,edgecolor="white",linewidth=1.2,zorder=3)
    if frac>=18:   # wide segment -> label inside
        axbar.text(cx,0,f"{lab}\n{n:,} · {frac:.0f}%",ha="center",va="center",fontsize=6.5,color="white" if col in WHITE_TXT else "#222",fontweight="bold",zorder=4)
    else:          # narrow segment -> label below with leader
        axbar.plot([cx,cx],[-0.28,-0.52],color=col,lw=0.7,zorder=2)
        axbar.text(cx,-0.58,f"{lab}\n{n:,} ({frac:.1f}%)",ha="center",va="top",fontsize=6.0,color=col,fontweight="bold",zorder=4)
    left+=n
gu=int(cnts.get("unique: conserved-but-silent",0))+int(cnts.get("unique: strain-private locus",0)); gu0=tot-gu
axbar.plot([gu0,tot],[-1.15,-1.15],color="#222",lw=1.2)
for xx in (gu0,tot): axbar.plot([xx,xx],[-1.09,-1.21],color="#222",lw=1.2)
axbar.text((gu0+tot)/2,-1.30,f"GENUINELY UNIQUE — {gu:,} ({100*gu/tot:.0f}%)",ha="center",va="top",fontsize=6.7,fontweight="bold",color="#222")
# complement bracket: the NOT-unique majority = expressed-elsewhere + SNP-variant + low-quality (sequence also present/expressed in the pooled OTHER strains)
axbar.plot([0,gu0],[-1.15,-1.15],color="#8a6d3b",lw=1.2)
for xx in (0,gu0): axbar.plot([xx,xx],[-1.09,-1.21],color="#8a6d3b",lw=1.2)
axbar.text(gu0/2,-1.30,f"EXPRESSED ELSEWHERE / NOT UNIQUE — {gu0:,} ({100*gu0/tot:.0f}%)",ha="center",va="top",fontsize=6.7,fontweight="bold",color="#8a6d3b")
axbar.set_xlim(-tot*0.07,tot*1.07); axbar.set_ylim(-1.85,0.55); axbar.axis("off")
# funnel (left) + plain-language note (right): why the bar > each per-'Combined'-panel n
pooled=[int(cnts.get(k,0)) for k in KL5]
uqn=fc2.groupby("klass5").sequence.nunique(); unique=[int(uqn.get(k,0)) for k in KL5]
cbn=df[df.view=="combined"].groupby("klass5").n_features.first(); all3=[int(cbn.get(k,0)) for k in KL5]
n_pooled=sum(pooled); n_all3=sum(all3)
draw_cascade(fig,[0.085,0.828,0.40,0.056],
    [("strain-specific candidates\npooled over all 3 timepoints",pooled),
     ("distinct sequences\n(de-duplicated across timepoints)",unique),
     ("expressed at all 3 timepoints\n→ each segment = that class's 'Combined' n",all3)],
    KCOL, title="Why bar totals > each 'Combined' panel n  ·  segments coloured by class")
fig.text(0.525,0.884,
         f"Bar = class composition of all {n_pooled:,} strain-specific candidates, pooled across\n"
         "E16.5 + P12.5 + P20.5 (the per-timepoint panels sum to it). A pooled PCA needs a\n"
         "value in all 144 libraries, so 'Combined' keeps only sequences detected at all three\n"
         f"timepoints ({n_all3:,} = {100*n_all3/n_pooled:.0f}% of pooled). Note: 'all 3 timepoints' is across DEVELOPMENT,\n"
         "not strains: a strain-private piRNA stably expressed in its one strain still\n"
         "qualifies — which is why a strain-private 'Combined' panel exists.",
         ha="left",va="top",fontsize=7.3,color="#444",linespacing=1.5)
import os as _os; _SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/07_unique_piRNA_identification/data/source_data"; _os.makedirs(_SD,exist_ok=True)
df.to_csv(f"{_SD}/SourceData_Fig_pca_classes16.csv",index=False)   # per-candidate PCA coordinates (PC1/PC2) by class5 + timepoint view
for e in ("pdf","svg","png"): fig.savefig(f"{U}/pangenome_te/Fig_pca_classes16.{e}",bbox_inches="tight")
print("wrote Fig_pca_classes16.{png,pdf,svg}")
print(df.groupby(["klass5","view"]).size())
