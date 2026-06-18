#!/usr/bin/env python3
"""
Theme 03 (REFRAMED, pangenome) — expression of Zamore pachytene piRNA loci across
16 strains x 3 timepoints. Nature-Genetics style (Okabe-Ito / Wong palette,
Liberation Sans, vector + 300 dpi). Cross-strain projection = cactus halLiftover
(pangenome); the ONLY UCSC liftOver upstream is mm10->mm39 of the Zamore loci.

Input : all_strains_pangenome/combined_rebuild/all_strains_expression_matrix.csv
Output: Fig_zamore_expression_pangenome.{pdf,svg,png} + SourceData_*.csv
"""
import sys
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from strain_order import STRAIN_ORDER, WILD

plt.rcParams.update({"font.family":"Liberation Sans","font.size":7.5,"axes.linewidth":0.6,
    "axes.spines.top":False,"axes.spines.right":False,"xtick.major.width":0.6,
    "ytick.major.width":0.6,"xtick.major.size":2.5,"ytick.major.size":2.5,
    "pdf.fonttype":42,"svg.fonttype":"none"})

OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/all_strains_pangenome/combined_rebuild"
df=pd.read_csv(f"{OUT}/all_strains_expression_matrix.csv")
TPO=["E16.5","P12.5","P20.5"]
# developmental "maturation heat" ramp: gold -> pumpkin -> brick-red (E16.5->P12.5->P20.5);
# colourblind-safe YlOrRd-family, vivid, intuitive for a developmental progression
COL={"E16.5":"#FFC84D","P12.5":"#F07F1A","P20.5":"#C1352A"}
order=[s for s in STRAIN_ORDER if s in set(df.strain)]
# per-replicate % expressed (build_replicate_pct_expressed.py) -> +-1 SD across replicates
rep=pd.read_csv(f"{OUT}/replicate_pct_expressed.csv")
sd=(rep.pivot_table(index="strain",columns="timepoint",values="pct_expressed",aggfunc="std")
       .reindex(order)[TPO])

# % of Zamore loci EXPRESSED per strain x timepoint
n_tot=df.locus.nunique()
piv=(df.assign(expr=(df.status=="expressed"))
       .pivot_table(index="strain",columns="timepoint",values="expr",aggfunc="mean")
       .reindex(order)[TPO])*100
# absolute count of expressed loci per strain x timepoint (annotated atop bars)
cnt=(df.assign(expr=(df.status=="expressed"))
       .pivot_table(index="strain",columns="timepoint",values="expr",aggfunc="sum")
       .reindex(order)[TPO]).round().astype(int)

fig,(axA,axB)=plt.subplots(1,2,figsize=(9.0,3.9),dpi=300,gridspec_kw={"width_ratios":[2.5,1]})

# ── Panel A: % expressed, strain x timepoint (grouped bars) ────────────────────
n=len(order); x=np.arange(n); bw=0.27
for j,t in enumerate(TPO):
    xs=x+(j-1)*bw
    axA.bar(xs, piv[t].values, width=bw, color=COL[t], edgecolor="white",
            linewidth=0.3, label=t, zorder=3)
    errs=np.nan_to_num(sd[t].values)                       # +-1 SD across the 3 replicates
    axA.errorbar(xs, piv[t].values, yerr=errs, fmt="none", ecolor="#333",
                 elinewidth=0.6, capsize=1.4, capthick=0.6, zorder=5)
    for xi,v,c,e in zip(xs, piv[t].values, cnt[t].values, errs):   # loci count atop each bar (above the error cap)
        axA.text(xi, v+e+1.2, str(int(c)), rotation=90, ha="center", va="bottom",
                 fontsize=4.2, color="#333", zorder=6)
axA.set_xticks(x)
labs=axA.set_xticklabels(order, rotation=55, ha="right", fontsize=6.6)
for lab,s in zip(labs,order): lab.set_color("#C0392B" if s in WILD else "#222222")
axA.set_ylabel("Zamore pachytene loci expressed (%)", fontsize=8)
axA.set_ylim(0,119); axA.set_yticks([0,20,40,60,80,100]); axA.set_xlim(-0.6,n-0.4)
axA.yaxis.grid(True, lw=0.3, color="#e9e9e9", zorder=0); axA.set_axisbelow(True)
# legend ABOVE the axes (bars fill the panel -> no clean interior space)
axA.legend(title="timepoint", fontsize=7, title_fontsize=7.5, frameon=False,
           loc="lower right", bbox_to_anchor=(1.0, 1.005), ncol=3,
           handlelength=1.1, columnspacing=0.9)
axA.set_title("A   piRNA-precursor expression rises into pachytene (P20.5), all 16 strains",
              fontsize=8.5, fontweight="bold", loc="left", pad=16)

# ── Panel B: developmental rise (mean % expressed per timepoint, dots=strains) ─
rng=np.random.default_rng(0)
for j,t in enumerate(TPO):
    v=piv[t].values
    axB.bar(j, v.mean(), width=0.6, color=COL[t], edgecolor="none", zorder=1)
    axB.scatter(np.full(len(v),j)+rng.uniform(-0.13,0.13,len(v)), v, s=10,
                color="#222", zorder=3, linewidths=0, alpha=0.8)
    axB.text(j, v.max()+3.5, f"{v.mean():.0f}%\n({cnt[t].mean():.0f} loci)", ha="center", va="bottom",
             fontsize=6.3, fontweight="bold", color=COL[t])
axB.set_xticks(range(3)); axB.set_xticklabels(TPO, fontsize=7)
axB.set_ylim(0,113); axB.set_yticks([0,20,40,60,80,100]); axB.set_ylabel("loci expressed (%)", fontsize=8)
axB.set_title("B   developmental\nprogram", fontsize=8.5, fontweight="bold", loc="left", pad=6)

fig.text(0.5,-0.04,
    f"n={n_tot} Zamore pachytene piRNA loci (Li 2013 Mol Cell PMID 23523368; mm10→mm39 UCSC liftOver) · bar = combined-run PICB clusters; number atop each bar = loci expressed (of {n_tot}) · "
    "error bar = ±1 SD across the 3 replicates (per-replicate PICB clusters) · expression = locus overlaps strain PICB clusters · "
    "cross-strain projection = cactus halLiftover (pangenome) · red strain label = wild-derived",
    ha="center", fontsize=5.6, color="#666")
fig.tight_layout()
base=f"{OUT}/Fig_zamore_expression_pangenome"
for ext in ("pdf","svg","png"): fig.savefig(f"{base}.{ext}", bbox_inches="tight")
print("wrote", base+".{pdf,svg,png}")

# ── source data ───────────────────────────────────────────────────────────────
sd=piv.reset_index().rename(columns={t:f"pct_expressed_{t}" for t in TPO})
sd.insert(0,"plot_order",range(1,len(sd)+1))
sd["subspecies"]=np.where(sd.strain.isin(WILD),"wild-derived","classical")
sd.to_csv(f"{OUT}/SourceData_Fig_zamore_expression_pangenome.csv", index=False)
print("mean % expressed by tp:", {t:round(piv[t].mean(),1) for t in TPO})
