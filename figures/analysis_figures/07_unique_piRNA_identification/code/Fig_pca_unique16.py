#!/usr/bin/env python3
"""16-strain per-timepoint PCA of piRNA-sequence expression (thesis Fig 5.21 / METHODS §7). 2 rows x 3 cols:
top = all expressed piRNAs (top-500 variable); bottom = genuinely-unique (final_classified). Points = samples
(16 strains x up to 3 reps), coloured by strain in canonical order. Input = pca16/{tp}.pca.csv (DESeq2-normalised
prcomp from pca_unique16.R)."""
import sys, os; sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull
def _tphull(a, sub, tp, col):   # soft convex hull around a timepoint's samples — makes the developmental separation pop
    pts=sub[sub.tp==tp][["PC1","PC2"]].values
    if len(pts)>=3:
        try:
            v=ConvexHull(pts).vertices; poly=pts[v]
            a.fill(poly[:,0],poly[:,1],color=col,alpha=0.08,zorder=1,lw=0)
            a.plot(np.append(poly[:,0],poly[0,0]),np.append(poly[:,1],poly[0,1]),color=col,alpha=0.35,lw=0.7,zorder=2)
        except Exception: pass
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
CANON=[s for s in STRAIN_ORDER if s!="C57BL_6"]   # 16 strains in the data
PAL=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf','#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5','#9edae5']
SCOL={s:PAL[i] for i,s in enumerate(CANON)}
MARK={s:("*" if s in WILD else "o") for s in CANON}  # wild = star
TPCOL={"E16.5":"#4393C3","P12.5":"#E8852B","P20.5":"#B2182B"}   # timepoint colour (matches the locus figures)
df=pd.concat([pd.read_csv(f"{U}/pca16/{t}.pca.csv") for t in ["16.5dpc","12.5dpp","20.5dpp"]],ignore_index=True)
df.to_csv(f"{U}/SourceData_pca_unique16.csv",index=False)
comb=pd.read_csv(f"{U}/pca16/combined_byclass.csv") if os.path.exists(f"{U}/pca16/combined_byclass.csv") else None   # all-timepoints-pooled PCA, ONE per feature-set class
TPO=["E16.5","P12.5","P20.5"]; FSO=["all_full","all_expressed","unique"]
FTITLE={"all_full":"ALL expressed piRNAs (every feature, no top-500 filter)","all_expressed":"All expressed piRNAs (top-500 variable)","unique":"Genuinely-unique piRNAs (16-strain)"}
plt.rcParams.update({"font.family":"Liberation Sans"})
fig=plt.figure(figsize=(11.5,21.8),dpi=300); gs=fig.add_gridspec(6,3,height_ratios=[1.0,1.0,1.0,1.0,1.3,0.85],hspace=0.95,wspace=0.34)
ax=np.array([[fig.add_subplot(gs[r,c]) for c in range(3)] for r in range(3)])
for r,fs in enumerate(FSO):
    for c,tp in enumerate(TPO):
        a=ax[r,c]; sub=df[(df.feature_set==fs)&(df.timepoint==tp)]
        for s in [x for x in CANON if x in set(sub.strain)]:
            ss=sub[sub.strain==s]
            a.scatter(ss.PC1,ss.PC2,s=52,color=SCOL[s],marker=MARK[s],edgecolor="black",linewidth=0.4,zorder=3,
                      label=s.replace("_","/") if (r==0 and c==0) else None)
        if len(sub):
            a.set_xlabel(f"PC1 ({sub.pc1_var.iloc[0]}%)",fontsize=8); a.set_ylabel(f"PC2 ({sub.pc2_var.iloc[0]}%)",fontsize=8)
            a.set_title(f"{tp}  (n={int(sub.n_features.iloc[0]):,})",fontsize=8.6,fontweight="bold")
        a.tick_params(labelsize=6); a.axhline(0,lw=0.4,color="#ddd",zorder=0); a.axvline(0,lw=0.4,color="#ddd",zorder=0)
# strain legend handles (placed as a figure legend below the panel grid, after layout — keeps it out of the data)
strain_handles=[Line2D([],[],marker=MARK[s],color="w",markerfacecolor=SCOL[s],markeredgecolor="black",
                       markeredgewidth=0.4,markersize=7,label=s.replace("_","/")) for s in CANON]
# bottom: (row 3) per-class pooled PCA panels + (row 4) the wide all-expressed pooled panel — colour = timepoint, ★ = wild, soft per-tp hulls
CMB=[("all_expressed","top-500 most variable"),("unique","genuinely-unique")]   # all_full pooled dropped here — it is the same data as the wide panel below, so shown only once there
gs3=gs[3,:].subgridspec(1,2,wspace=0.30); axc=[fig.add_subplot(gs3[0,c]) for c in range(2)]; axb=fig.add_subplot(gs[4,:])
tp_handles=[Line2D([],[],marker='o',color='w',markerfacecolor=TPCOL[t],markeredgecolor='k',label=t,markersize=8) for t in TPO]+[Line2D([],[],marker='*',color='w',markerfacecolor='#aaaaaa',markeredgecolor='k',label='wild-derived',markersize=11)]
if comb is not None:
    for c,(fs,lab) in enumerate(CMB):
        a=axc[c]; sub=comb[comb.feature_set==fs]
        for tp in TPO: _tphull(a, sub, tp, TPCOL[tp])
        for tp in TPO:
            for s in CANON:
                ss=sub[(sub.tp==tp)&(sub.strain==s)]
                a.scatter(ss.PC1,ss.PC2,s=40,color=TPCOL[tp],marker=("*" if s in WILD else "o"),edgecolor="black",linewidth=0.4,zorder=3)
        if len(sub):
            a.set_xlabel(f"PC1 ({sub.pc1_var.iloc[0]}%)",fontsize=8); a.set_ylabel(f"PC2 ({sub.pc2_var.iloc[0]}%)",fontsize=8)
            a.set_title(f"{lab}\n(pooled; n={int(sub.n_features.iloc[0]):,})",fontsize=8.3,fontweight="bold")
        a.tick_params(labelsize=6); a.axhline(0,lw=0.4,color="#ddd",zorder=0); a.axvline(0,lw=0.4,color="#ddd",zorder=0)
    # wide RESTORED hero panel: all-shared-expressed pooled (the original detailed view) + per-tp hulls
    aw=comb[comb.feature_set=="all_full"]
    for tp in TPO: _tphull(axb, aw, tp, TPCOL[tp])
    for tp in TPO:
        for s in CANON:
            ss=aw[(aw.tp==tp)&(aw.strain==s)]
            axb.scatter(ss.PC1,ss.PC2,s=62,color=TPCOL[tp],marker=("*" if s in WILD else "o"),edgecolor="black",linewidth=0.4,zorder=3)
    axb.set_xlabel(f"PC1 ({aw.pc1_var.iloc[0]}%)",fontsize=8.5); axb.set_ylabel(f"PC2 ({aw.pc2_var.iloc[0]}%)",fontsize=8.5)
    axb.set_title(f"Combined — all timepoints pooled (ALL shared expressed piRNAs, no top-500 filter; n={int(aw.n_features.iloc[0]):,}; 144 libraries)   ·   colour = timepoint, "+r"$\bigstar$"+" = wild-derived",fontsize=9.0,fontweight="bold")
    axb.legend(handles=tp_handles,fontsize=7.2,frameon=False,loc="best",ncol=2)
    axb.tick_params(labelsize=6.5); axb.axhline(0,lw=0.4,color="#ddd",zorder=0); axb.axvline(0,lw=0.4,color="#ddd",zorder=0)
else:
    for a in axc+[axb]: a.axis("off")
    axb.text(0.5,0.5,"combined_byclass.csv not found — run combine_pca16_byclass.R",ha="center",va="center",transform=axb.transAxes,color="#999")
# 100%-scale class-composition bar (the 5 strain-specific classes, ≥2-read adopted) with actual counts + genuinely-unique bracket
fc2=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["klass5"])
KL5=["expressed elsewhere (exact)","SNP-variant (1-3mm)","low-quality: no mm0 own-genome locus","unique: conserved-but-silent","unique: strain-private locus"]
KLAB5=["expressed-elsewhere","SNP-variant","low-quality","conserved-but-silent","strain-private"]
KCOL=["#9e9e9e","#E69F00","#cdb892","#0072B2","#7a3b9a"]; WHITE_TXT={"#9e9e9e","#0072B2","#7a3b9a"}
cnts=fc2.klass5.value_counts(); tot=int(cnts.sum())
axbar=fig.add_subplot(gs[5,:]); left=0
for k,lab,col in zip(KL5,KLAB5,KCOL):
    n=int(cnts.get(k,0)); frac=100*n/tot; cx=left+n/2
    axbar.barh(0,n,left=left,height=0.5,color=col,edgecolor="white",linewidth=1.2,zorder=3)
    if frac>=18:   # wide segment -> label inside
        axbar.text(cx,0,f"{lab}\n{n:,} · {frac:.0f}%",ha="center",va="center",fontsize=6.6,color="white" if col in WHITE_TXT else "#222",fontweight="bold",zorder=4)
    else:          # narrow segment -> label above with leader (avoids horizontal clipping)
        axbar.plot([cx,cx],[0.25,0.44],color=col,lw=0.7,zorder=2)
        axbar.text(cx,0.48,f"{lab}\n{n:,} ({frac:.1f}%)",ha="center",va="bottom",fontsize=6.0,color=col,fontweight="bold",zorder=4)
    left+=n
gu=int(cnts.get("unique: conserved-but-silent",0))+int(cnts.get("unique: strain-private locus",0)); gu0=tot-gu
axbar.plot([gu0,tot],[-0.4,-0.4],color="#222",lw=1.3)
for xx in (gu0,tot): axbar.plot([xx,xx],[-0.34,-0.46],color="#222",lw=1.3)
axbar.text((gu0+tot)/2,-0.54,f"GENUINELY UNIQUE — {gu:,} ({100*gu/tot:.0f}%)",ha="center",va="top",fontsize=7,fontweight="bold",color="#222")
# complement bracket: the NOT-unique majority = expressed-elsewhere + SNP-variant + low-quality (sequence also present/expressed in the pooled OTHER strains = the 'all expressed pooled' set)
axbar.plot([0,gu0],[-0.4,-0.4],color="#8a6d3b",lw=1.3)
for xx in (0,gu0): axbar.plot([xx,xx],[-0.34,-0.46],color="#8a6d3b",lw=1.3)
axbar.text(gu0/2,-0.54,f"EXPRESSED ELSEWHERE / NOT UNIQUE — {gu0:,} ({100*gu0/tot:.0f}%)",ha="center",va="top",fontsize=7,fontweight="bold",color="#8a6d3b")
axbar.set_xlim(-tot*0.07,tot*1.07); axbar.set_ylim(-0.95,1.25); axbar.axis("off")
axbar.set_title(f"Composition of strain-specific piRNA candidates by class (≥2-read adopted; n={tot:,}; 100% scale)",fontsize=9,fontweight="bold",loc="left")
fig.suptitle("PCA of piRNA expression across all 16 strains (DESeq2-normalised) — thesis Fig 5.21 / METHODS §7 ("+r"$\bigstar$"+"=wild); bottom: per-class pooled PCA + wide all-expressed pooled + class-composition bar",fontsize=9.4,fontweight="bold",y=1.002)
fig.tight_layout()
# centered title above each feature-set row ("panel wise with title"), replacing the rotated side labels that overlapped the y-axis
ROWT={"all_full":"ALL expressed piRNAs (every feature, no top-500 filter)","all_expressed":"All expressed piRNAs (top-500 most variable)","unique":"Genuinely-unique piRNAs (16-strain, pangenome)"}
for r,fs in enumerate(FSO):
    p0=ax[r,0].get_position(); p2=ax[r,2].get_position()
    fig.text((p0.x0+p2.x1)/2, p0.y1+0.033, ROWT[fs], ha="center", va="bottom", fontsize=10.5, fontweight="bold")
# 16-strain legend in the gap between the small-panel grid and the combined panel (out of the data)
ly=(axc[0].get_position().y1+ax[2,0].get_position().y0)/2
fig.legend(handles=strain_handles, loc="center", bbox_to_anchor=(0.5,ly), ncol=8, fontsize=7,
           frameon=False, handletextpad=0.3, columnspacing=1.1, title="strain (markers in the 9 per-timepoint panels above)", title_fontsize=7)
for e in ("pdf","svg","png"): fig.savefig(f"{U}/pangenome_te/Fig_pca_unique16.{e}",bbox_inches="tight")
print(df[["feature_set","timepoint","pc1_var","pc2_var","n_features"]].drop_duplicates().to_string(index=False))
print("wrote Fig_pca_unique16.{png,pdf,svg} + SourceData_pca_unique16.csv")
