#!/usr/bin/env python3
"""piRNA expression is extremely SKEWED at pachytene: a few 'master' clusters make most piRNAs. Lorenz curves of
PICB-cluster FPM per strain x timepoint (all 16 strains overlaid, coloured by timepoint) + Gini + the % of clusters
that account for 90% of piRNA. PREPACHYTENE (E16.5) = many dispersed TE-silencing loci (near-equal); PACHYTENE
(P20.5) = a handful of A-MYB master clusters dominate (~90% from ~3% of clusters)."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis"); sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
from strain_order import STRAIN_ORDER
import pav_clusters as pc
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]; PG = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"; SD = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/source_data"
TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]; TPC = {"16.5dpc": "#0072B2", "12.5dpp": "#009E73", "20.5dpp": "#D55E00"}; TPL = {"16.5dpc": "E16.5 (prepachytene)", "12.5dpp": "P12.5", "20.5dpp": "P20.5 (pachytene)"}
P = pc._pang()
def gini(f):
    f = np.sort(f); n = len(f); return (2 * np.sum(np.arange(1, n + 1) * f) / (n * f.sum())) - (n + 1) / n
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": 0.8})
fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.6), dpi=300)
ginis = {t: [] for t in TPS}; pct90 = {t: [] for t in TPS}
for X in ORDER:
    for tp in TPS:
        d = P[(P.strain == X) & (P.tp == tp)]
        if len(d) < 20: continue
        f = d.loc[d.groupby(["own_chrom", "own_start"]).all_primary_FPM.idxmax()].all_primary_FPM.values
        fs = np.sort(f); cum = np.concatenate([[0], np.cumsum(fs) / fs.sum()]); x = np.concatenate([[0], np.arange(1, len(fs) + 1) / len(fs)])
        axA.plot(x, cum, color=TPC[tp], lw=0.5, alpha=0.35)
        ginis[tp].append(gini(f)); fd = np.sort(f)[::-1]; cd = np.cumsum(fd) / fd.sum(); pct90[tp].append(100 * (np.searchsorted(cd, 0.9) + 1) / len(fd))
axA.plot([0, 1], [0, 1], "k--", lw=0.9)
for tp in TPS: axA.plot([], [], color=TPC[tp], lw=2.2, label=f"{TPL[tp]} — Gini {np.mean(ginis[tp]):.2f}")
axA.plot([], [], "k--", lw=0.9, label="perfect equality")
axA.set_xlabel("cumulative fraction of piRNA clusters (low→high FPM)", fontsize=8.5); axA.set_ylabel("cumulative fraction of piRNA (FPM)", fontsize=8.5)
axA.legend(fontsize=7.2, frameon=False, loc="upper left"); axA.spines[["top", "right"]].set_visible(False); axA.tick_params(labelsize=8); axA.set_xlim(0, 1); axA.set_ylim(0, 1)
axA.set_title("a   Lorenz curves: piRNA output is concentrated at pachytene\n(each line = one strain × timepoint; 16 strains)", fontsize=9, fontweight="bold", loc="left")
m90 = [np.mean(pct90[t]) for t in TPS]
axB.bar(range(3), m90, color=[TPC[t] for t in TPS], width=0.62, edgecolor="white")
for j, t in enumerate(TPS):
    axB.scatter(np.full(len(pct90[t]), j) + np.random.default_rng(0).uniform(-0.12, 0.12, len(pct90[t])), pct90[t], s=10, color="#333", alpha=0.5, zorder=3)
    axB.text(j, m90[j] + 1.2, f"{m90[j]:.0f}%", ha="center", fontsize=9, fontweight="bold", color=TPC[t])
axB.set_xticks(range(3)); axB.set_xticklabels([TPL[t].replace(" (", "\n(") for t in TPS], fontsize=8); axB.set_ylabel("% of clusters making 90% of piRNA", fontsize=8.5)
axB.spines[["top", "right"]].set_visible(False); axB.tick_params(labelsize=8); axB.set_ylim(0, max(m90) * 1.2)
axB.set_title("b   At pachytene, ~3% of clusters make 90% of piRNA\n(prepachytene is dispersed across many TE loci)", fontsize=9, fontweight="bold", loc="left")
fig.suptitle("A few 'master' clusters dominate the piRNA population at pachytene — the A-MYB-driven pachytene burst (Li 2013); prepachytene TE-silencing piRNAs are dispersed", fontsize=9.6, fontweight="bold", y=1.02)
pd.DataFrame({"timepoint": TPS, "mean_Gini": [np.mean(ginis[t]) for t in TPS], "mean_pct_clusters_for_90pct": m90}).to_csv(f"{SD}/Fig_pirna_lorenz16_sourcedata.csv", index=False)
fig.tight_layout()
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_pirna_lorenz16.{e}", bbox_inches="tight")
print("Gini:", {t: round(np.mean(ginis[t]), 3) for t in TPS}, "| %clusters for 90%:", {t: round(m, 1) for t, m in zip(TPS, m90)})
print("wrote Fig_pirna_lorenz16.{png,pdf,svg} + source data")
