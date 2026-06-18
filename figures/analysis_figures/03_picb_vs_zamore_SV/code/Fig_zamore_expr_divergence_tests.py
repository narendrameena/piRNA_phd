#!/usr/bin/env python3
"""
Theme 03 (REFRAMED) — interesting, statistically-tested biology of Zamore pachytene
piRNA loci across 16 strains. Nature-Genetics style; constrained_layout (no overlaps).

A  Expression concentration (Lorenz): a few loci make most piRNAs (Gini~0.90;
   ~25 loci -> 90% of output). Established: Li 2013 Mol Cell PMID 23523368;
   Gainetdinov 2018 (BioMNI-confirmed).
B  Conserved master set: the top loci fire in EVERY strain (cross-strain expression
   Spearman ~0.88) despite rapid sequence divergence (Yu 2021 PMID 33397987).
C  Divergence ~ cross-strain expression VARIABILITY (Spearman rho, p) — divergent
   loci express more variably (level itself is NOT predicted by divergence).
D  Wild vs classical expression concentration (Gini; Mann-Whitney).
"""
import sys
sys.path.insert(0,"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu
from strain_order import STRAIN_ORDER, WILD

plt.rcParams.update({"font.family":"Liberation Sans","font.size":7.5,"axes.linewidth":0.6,
    "axes.spines.top":False,"axes.spines.right":False,"xtick.major.width":0.6,
    "ytick.major.width":0.6,"xtick.major.size":2.5,"ytick.major.size":2.5,
    "pdf.fonttype":42,"svg.fonttype":"none"})
BLUE="#0072B2"; ORANGE="#E69F00"; VERM="#D55E00"; SKY="#56B4E9"; GREEN="#009E73"; PINK="#CC79A7"

CR="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/all_strains_pangenome/combined_rebuild"
M20=pd.read_csv(f"{CR}/zamore_locus_expression_P20.5.csv",index_col=0)
ann=pd.read_csv(f"{CR}/zamore_locus_annotation.csv")
con=pd.read_csv(f"{CR}/SourceData_zamore_concentration.csv")
strains=[s for s in STRAIN_ORDER if s in M20.columns]
M20=M20[strains]

def gini(x):
    x=np.sort(np.asarray(x,float)); x=x[x>=0]; n=len(x)
    if n==0 or x.sum()==0: return np.nan
    return (2*np.sum(np.arange(1,n+1)*x)/(n*x.sum()))-(n+1)/n

fig,axes=plt.subplots(2,2,figsize=(8.4,7.2),dpi=300,constrained_layout=True)
axA,axB,axC,axD=axes.flat

# ── A: Lorenz curves ──────────────────────────────────────────────────────────
for s in strains:
    v=np.sort(M20[s].values)[::-1]; v=v[v>0]
    if len(v)==0: continue
    cx=np.arange(1,len(v)+1)/len(v)*100; cy=np.cumsum(v)/v.sum()*100
    axA.plot(np.concatenate([[0],cx]),np.concatenate([[0],cy]),
             color=(VERM if s in WILD else BLUE),lw=0.7,alpha=0.55,zorder=2)
# median curve over a common loci-fraction grid
grid=np.linspace(0,100,101); curves=[]
for s in strains:
    v=np.sort(M20[s].values)[::-1]; v=v[v>0]
    if len(v)==0: continue
    cx=np.arange(1,len(v)+1)/len(v)*100; cy=np.cumsum(v)/v.sum()*100
    curves.append(np.interp(grid,np.concatenate([[0],cx]),np.concatenate([[0],cy])))
axA.plot(grid,np.median(curves,axis=0),color="#222",lw=2.0,zorder=4,label="median strain")
axA.plot([0,100],[0,100],ls=(0,(3,2)),lw=0.7,color="#999",zorder=1,label="equal output")
medg=np.nanmedian([gini(M20[s].values) for s in strains])
n90=int(np.nanmedian(con.loci_for_90pct)); top10=np.nanmedian(con.top10_pct)
axA.annotate(f"~{n90} loci → 90%\nof all piRNA output",xy=(n90/len(M20)*100,90),
             xytext=(38,55),fontsize=6.6,color="#222",
             arrowprops=dict(arrowstyle="->",lw=0.7,color="#222"))
axA.set_xlabel("loci (%, ranked by output)",fontsize=7.5)
axA.set_ylabel("cumulative piRNA output (%)",fontsize=7.5)
axA.set_title(f"A   A few loci make most piRNAs (Gini={medg:.2f})",fontsize=8.3,fontweight="bold",loc="left")
axA.legend(fontsize=6.3,frameon=False,loc="lower right")
axA.set_xlim(0,100); axA.set_ylim(0,100)

# ── B: conserved master set heatmap (top 25 loci x strains, log FPM) ──────────
top=M20.mean(axis=1).sort_values(ascending=False).head(25).index
H=np.log10(M20.loc[top,strains].values+1)
from matplotlib.colors import LinearSegmentedColormap
_bcm=LinearSegmentedColormap.from_list("paleSolid",["#fbfcfe","#9ecae1","#3182bd","#08306b"])
im=axB.imshow(H,aspect="auto",cmap=_bcm)   # single hue: pale = low -> solid = high
axB.set_xticks(range(len(strains)))
xl=axB.set_xticklabels([s.replace("_","/") for s in strains],rotation=90,fontsize=4.6)
for lab,s in zip(xl,strains): lab.set_color("#C0392B" if s in WILD else "#222")
axB.set_yticks(range(len(top))); axB.set_yticklabels(top,fontsize=4.3)
cc=M20.apply(lambda c: np.log1p(c)).corr(method="spearman").values
mean_rho=cc[np.triu_indices(len(strains),1)].mean()
axB.set_title(f"B   Same master loci fire in every strain (ρ̄={mean_rho:.2f})",fontsize=8.3,fontweight="bold",loc="left")
cb=fig.colorbar(im,ax=axB,fraction=0.046,pad=0.02); cb.set_label("log10(FPM+1)",fontsize=6.3); cb.ax.tick_params(labelsize=5.5)

# ── C: divergence vs expression CV, colour-coded WILD vs CLASSICAL ────────────
# each locus contributes a point per strain-group: divergence & CV computed WITHIN
# that group (so the wild-derived strains' larger divergence is visible).
frac=pd.read_csv(f"{CR}/zamore_fraction_lifted.csv")
wild=[s for s in strains if s in WILD]; classical=[s for s in strains if s not in WILD]
def grp(g,color,label):
    div=(1-frac[frac.strain.isin(g)].groupby("locus")["fraction_lifted"].mean())
    sub=M20[g]; cv=(sub.std(axis=1)/sub.mean(axis=1).replace(0,np.nan))
    d=pd.DataFrame({"dv":div,"cv":cv}).dropna(); d=d[d["dv"]>=0]
    axC.scatter(d["dv"]*100,d["cv"],s=10,color=color,alpha=0.5,linewidths=0,label=label,zorder=3)
    return spearmanr(d["dv"],d["cv"])
rc,pc=grp(classical,BLUE,f"classical (n={len(classical)})")
rw,pw=grp(wild,VERM,f"wild-derived (n={len(wild)})")
axC.set_xlabel("sequence divergence within group (1 − retention, %)",fontsize=7.5)
axC.set_ylabel("within-group expression CV",fontsize=7.5)
axC.set_title(f"C   Divergent loci express more variably (classical ρ={rc:.2f} p={pc:.0e}; wild ρ={rw:.2f})",
              fontsize=7.7,fontweight="bold",loc="left")
axC.legend(fontsize=6.3,frameon=False,loc="upper right")
axC.set_xscale("symlog",linthresh=0.5)

# ── D: wild vs classical concentration (Gini) ────────────────────────────────
gw=con[con.subspecies=="wild-derived"].gini.values; gc=con[con.subspecies=="classical"].gini.values
U,pg=mannwhitneyu(gw,gc)
parts=axD.boxplot([gc,gw],widths=0.5,patch_artist=True,showfliers=False,
                  medianprops=dict(color="#222",lw=1.2))
for patch,c in zip(parts["boxes"],[BLUE,VERM]): patch.set_facecolor(c); patch.set_alpha(0.55)
rng=np.random.default_rng(0)
for i,(g,c) in enumerate([(gc,BLUE),(gw,VERM)],start=1):
    axD.scatter(np.full(len(g),i)+rng.uniform(-0.11,0.11,len(g)),g,s=14,color=c,zorder=3,edgecolors="white",linewidths=0.3)
axD.set_xticks([1,2]); axD.set_xticklabels([f"classical\n(n={len(gc)})",f"wild\n(n={len(gw)})"],fontsize=6.8)
axD.set_ylabel("expression concentration (Gini)",fontsize=7.5)
axD.set_title(f"D   Wild slightly less concentrated (p={pg:.3f})",fontsize=8.3,fontweight="bold",loc="left")

fig.suptitle("Zamore pachytene piRNA loci: a conserved, highly-concentrated program across 16 strains (pangenome)",
             fontsize=9.5,fontweight="bold")
base=f"{CR}/Fig_zamore_expr_divergence_tests"
for ext in ("pdf","svg","png"): fig.savefig(f"{base}.{ext}",bbox_inches="tight")
print("wrote",base)
print(f"Gini median={medg:.3f} | ~{n90} loci=90% | top10={top10:.0f}% | cross-strain rho={mean_rho:.2f} | div-CV classical rho={rc:.2f} wild rho={rw:.2f} | wild-vs-classical Gini p={pg:.3f}")
