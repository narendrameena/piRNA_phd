#!/usr/bin/env python3
"""
Per-replicate piRNA phasing (+1 nt %): all 16 strains x E16.5/P12.5/P20.5 x replicates.
Horizontal bars; sample order = canonical figure strain order (thesis Fig 4.4,
median P20.5 PC1) -> timepoint (E16.5->P12.5->P20.5) -> replicate. Labels per sample.
Input : phasing_allstrains_1random/ALL_summary.csv
Method: 1 random coord/read (STAR --outSAMmultNmax 1 --outMultimapperOrder Random),
        24-32 nt, GenomicRanges::follow 3'->5' adjacency; +1 nt = phased (Almeida GB2025).
Writes the source-data CSV in the exact plotted order (same script).
"""
import sys
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from strain_order import STRAIN_ORDER, WILD, TIMEPOINT_ORDER, strain_rank

plt.rcParams.update({"font.family":"Liberation Sans","font.size":7,"axes.linewidth":0.5,
    "axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})

BASE="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/phasing_allstrains_1random"
df=pd.read_csv(f"{BASE}/ALL_summary.csv")
def tp(s): return "E16.5" if "16.5dpc" in s else "P12.5" if "12.5dpp" in s else "P20.5" if "20.5dpp" in s else "?"
def strain(s):
    for tk in ("16.5dpc","12.5dpp","20.5dpp"):
        if f"-{tk}." in s: return s.split(f"-{tk}.")[0]
    return s
def rep(s):
    try: return int(s.rsplit(".",1)[-1])
    except ValueError: return 0
df["strain"]=df["sample"].map(strain); df["tp"]=df["sample"].map(tp); df["rep"]=df["sample"].map(rep)
df["pct"]=df["frac_plus1"]*100
TPRANK={t:i for i,t in enumerate(TIMEPOINT_ORDER)}
df["sk"]=df["strain"].map(strain_rank); df["tk"]=df["tp"].map(TPRANK)
df=df.sort_values(["sk","tk","rep"]).reset_index(drop=True)

COL={"E16.5":"#F0C9A0","P12.5":"#E69F00","P20.5":"#B4500A"}
n=len(df); y=np.arange(n)[::-1]      # first sorted row at the top
fig,ax=plt.subplots(figsize=(6.2, max(8.0, n*0.135)), dpi=300)
ax.barh(y, df["pct"].values, color=[COL[t] for t in df["tp"]], height=0.78, edgecolor="none", zorder=2)
for yy, v in zip(y, df["pct"].values):                       # value at the tip of each bar
    ax.text(v + 0.6, yy, f"{v:.0f}", va="center", ha="left", fontsize=4.2, color="#333", zorder=3)
ax.set_yticks(y)
labels=[f"{r.strain}-{r.tp}.{r.rep}" for r in df.itertuples()]
ax.set_yticklabels(labels, fontsize=4.4)
for tick,(_,r) in zip(ax.get_yticklabels(), df.iterrows()):
    tick.set_color("#C0392B" if r.strain in WILD else "#222222")
ax.set_xlabel("+1-nt directly-adjacent piRNA pairs (% of all adjacent pairs)", fontsize=7.5)
ax.set_xlim(0, max(74, df["pct"].max()*1.10)); ax.set_ylim(-0.6, n-0.4)
ax.xaxis.grid(True, lw=0.3, color="#e8e8e8", zorder=0); ax.set_axisbelow(True)
# thin separators + strain name at block centre (left margin)
prev=None; block_start=0
for idx,(yy,(_,r)) in enumerate(zip(y, df.iterrows())):
    if prev is not None and r.strain!=prev:
        ax.axhline(yy+0.5, color="#bbbbbb", lw=0.4, zorder=1)
    prev=r.strain
leg=[Patch(facecolor=COL[t],label=t) for t in TIMEPOINT_ORDER]
ax.legend(handles=leg, title="timepoint", fontsize=6.5, title_fontsize=7, frameon=False, loc="lower right")
ax.set_title("piRNA phasing per replicate — 16 mouse strains × E16.5/P12.5/P20.5",
             fontsize=9.5, fontweight="bold", pad=6)
ax.tick_params(axis="x", labelsize=6.5)
fig.text(0.5,-0.004,
    "sample order = canonical figure strain order (thesis Fig 4.4, median P20.5 PC1) → timepoint → replicate · red label = wild-derived (CAST/PWK/SPRET/WSB)",
    ha="center", fontsize=5.2, color="#666")
fig.tight_layout()
out=f"{BASE}/Fig_phasing_perReplicate"
for ext in ("pdf","svg","png"): fig.savefig(f"{out}.{ext}", bbox_inches="tight")
print("wrote", out+".{pdf,svg,png}  (", n, "samples )")

# ---- source data, exact plotted order ----
SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/source_data"
sd=df[["sample","strain","tp","rep","frac_plus1","pct","zscore_plus1","count_plus1","n_pairs","total_aln"]].copy()
sd.insert(1,"plot_row_top_to_bottom", range(1,len(sd)+1))
sd.insert(4,"subspecies", np.where(sd["strain"].isin(WILD),"wild-derived","classical"))
sd=sd.rename(columns={"tp":"timepoint","rep":"replicate","frac_plus1":"plus1_fraction","pct":"plus1_pct",
    "zscore_plus1":"plus1_zscore","count_plus1":"count_plus1","n_pairs":"n_adjacent_pairs",
    "total_aln":"total_alignments_24_32nt"})
sd.to_csv(f"{SD}/SourceData_Fig_phasing_perReplicate.csv", index=False)
print("wrote source data:", len(sd), "rows ->", f"{SD}/SourceData_Fig_phasing_perReplicate.csv")
