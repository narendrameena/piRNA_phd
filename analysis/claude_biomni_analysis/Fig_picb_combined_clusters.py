#!/usr/bin/env python3
"""
Nature-Genetics-style figure: PICB combined-replicate piRNA CLUSTERS per strain.
All 16 strains together, x in canonical figure order (strain_order.py), grouped
bars by timepoint E16.5 -> P12.5 -> P20.5. y = number of PICB clusters (final
`clusters` sheet) from the combined-replicate run.
Input : source_data/SourceData_PICB_cluster_counts.csv  (replicate=='combined')
Writes the plotted source-data subset in the same script.
"""
import sys
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from strain_order import STRAIN_ORDER, WILD, TIMEPOINT_ORDER

plt.rcParams.update({"font.family":"Liberation Sans","font.size":8,"axes.linewidth":0.6,
    "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})

BASE="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis"
df=pd.read_csv(f"{BASE}/source_data/SourceData_PICB_cluster_counts.csv")
df=df[df.replicate=="combined"].copy()

order=[s for s in STRAIN_ORDER if s in set(df.strain)]   # 16 ISV (no C57BL_6 reference)
TPO=TIMEPOINT_ORDER
COL={"E16.5":"#F0C9A0","P12.5":"#E69F00","P20.5":"#B4500A"}
TPMAP={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
df["tp"]=df.timepoint.map(TPMAP)
piv=df.pivot_table(index="strain",columns="tp",values="n_clusters",aggfunc="first").reindex(order)[TPO]

n=len(order); x=np.arange(n); bw=0.27
ymax=np.nanmax(piv.values)
fig,ax=plt.subplots(figsize=(7.8,4.3),dpi=300)
for j,t in enumerate(TPO):
    xs=x+(j-1)*bw; vals=piv[t].values
    ax.bar(xs, vals, width=bw, color=COL[t], edgecolor="white", linewidth=0.3, label=t, zorder=3)
    for xi,v in zip(xs,vals):
        if np.isfinite(v):
            ax.text(xi, v+ymax*0.008, f"{int(round(v))}", rotation=90, ha="center",
                    va="bottom", fontsize=4.2, color="#333", zorder=4)
ax.set_ylim(0, ymax*1.22)
ax.set_xticks(x)
labs=ax.set_xticklabels(order, rotation=55, ha="right", fontsize=6.6)
for lab,s in zip(labs,order): lab.set_color("#C0392B" if s in WILD else "#222222")
ax.set_ylabel("PICB clusters (combined replicates)", fontsize=8.5)
ax.set_xlim(-0.6, n-0.4)
ax.yaxis.grid(True, lw=0.3, color="#e9e9e9", zorder=0); ax.set_axisbelow(True)
ax.legend(title="timepoint", fontsize=7, title_fontsize=7.5, frameon=False,
          loc="upper right", ncol=3, handlelength=1.1, columnspacing=1.0)
ax.set_title("PICB piRNA clusters per strain (combined-replicate run), 16 mouse strains",
             fontsize=9.5, fontweight="bold", pad=8)
fig.text(0.5,-0.10,
    "strains in canonical figure order (thesis Fig 4.4, median P20.5 PC1) · red label = wild-derived (CAST/PWK/SPRET/WSB) · "
    "clusters = PICB final `clusters` sheet, reps pooled before PICB",
    ha="center", fontsize=5.6, color="#666")
fig.tight_layout()
out=f"{BASE}/Fig_picb_combined_clusters"
for ext in ("pdf","svg","png"): fig.savefig(f"{out}.{ext}", bbox_inches="tight")
print("wrote", out+".{pdf,svg,png}")

# source data in plotted order
sd=piv.reset_index().rename(columns={"E16.5":"n_clusters_E16.5","P12.5":"n_clusters_P12.5","P20.5":"n_clusters_P20.5"})
sd.insert(0,"plot_order",range(1,len(sd)+1))
sd["subspecies"]=np.where(sd.strain.isin(WILD),"wild-derived","classical")
sd.to_csv(f"{BASE}/source_data/SourceData_Fig_picb_combined_clusters.csv", index=False)
print("wrote source data; values:\n", piv.round(0).to_string())
