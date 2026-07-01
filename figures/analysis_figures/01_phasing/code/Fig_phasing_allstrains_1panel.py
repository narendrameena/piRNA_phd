#!/usr/bin/env python3
"""ONE-PANEL version of the 16-strain piRNA phasing figure: grouped bars (16 strains × 3 timepoints)
in a SINGLE axis instead of 16 small-multiples. x = strains (canonical order; wild-derived red + shaded
block), 3 bars/strain = E16.5/P12.5/P20.5 (developmental gradient), bar = mean, error = ±SD, dots = replicates.
Input phasing_allstrains_1random/ALL_summary.csv (+1-nt phasing fraction; Almeida GB2025 method)."""
import sys
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from strain_order import STRAIN_ORDER, WILD, add_classical_wild_companion
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
BASE = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/phasing_allstrains_1random"
df = pd.read_csv(f"{BASE}/ALL_summary.csv")
def _tp(s): return "E16.5" if "16.5dpc" in s else "P12.5" if "12.5dpp" in s else "P20.5" if "20.5dpp" in s else "?"
def _strain(s):
    for tk in ("16.5dpc", "12.5dpp", "20.5dpp"):
        if f"-{tk}." in s: return s.split(f"-{tk}.")[0]
    return s
df["tp"] = df["sample"].map(_tp); df["strain"] = df["sample"].map(_strain); df["pct"] = df["frac_plus1"] * 100
TPO = ["E16.5", "P12.5", "P20.5"]; COL = {"E16.5": "#F0C9A0", "P12.5": "#E69F00", "P20.5": "#B4500A"}
TPL = {"E16.5": "E16.5 (prepachytene)", "P12.5": "P12.5", "P20.5": "P20.5 (pachytene)"}
order = [s for s in STRAIN_ORDER if s in set(df["strain"])]
WPOS = [i for i, s in enumerate(order) if s in WILD]
x = np.arange(len(order)); bw = 0.27; rng = np.random.default_rng(0)
fig, ax = plt.subplots(figsize=(15, 5.4), dpi=300)
if WPOS: ax.axvspan(min(WPOS) - 0.5, max(WPOS) + 0.5, color="#C0392B", alpha=0.06, zorder=0)
for j, t in enumerate(TPO):
    xs = x + (j - 1) * bw
    for xi, s in zip(xs, order):
        v = df[(df.strain == s) & (df.tp == t)].pct.values
        m = v.mean() if len(v) else 0
        ax.bar(xi, m, bw, color=COL[t], edgecolor="white", linewidth=0.3, zorder=2)
        if len(v):
            sd = v.std(ddof=1) if len(v) > 1 else 0
            ax.errorbar(xi, m, yerr=sd, fmt="none", ecolor="#333", elinewidth=0.5, capsize=1.5, zorder=3)   # error bar = ±SD
            ax.scatter(np.full(len(v), xi) + rng.uniform(-0.06, 0.06, len(v)), v, s=8, color="#222", zorder=4, linewidths=0)   # dots = replicates
            ax.text(xi, max(v.max(), m + sd) + 1.2, f"{m:.0f}", ha="center", va="bottom", fontsize=4.6, rotation=90, color=COL[t], fontweight="bold", zorder=5)   # mean value on top
ax.set_xticks(x); ax.set_xticklabels([])   # strain labels carried by the classical/wild companion below
ax.set_ylabel("+1 phasing (%)", fontsize=10.5); ax.set_ylim(0, 78)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(handles=[Patch(facecolor=COL[t], label=TPL[t]) for t in TPO], fontsize=9, frameon=False, ncol=3, loc="upper left", title="developmental window", title_fontsize=9)
if WPOS: ax.text(np.mean(WPOS), 70.5, "wild-derived", ha="center", va="top", fontsize=9, fontweight="bold", color="#C0392B")
ax.set_title("piRNA phasing across spermatogenesis — 16 mouse strains × 3 developmental windows", fontsize=12.5, fontweight="bold", loc="left")
fig.text(0.5, -0.17, "bars = mean +1-nt phasing fraction (value on top) · error bar = ±SD · dots = replicates · wild-derived (CAST/PWK/SPRET/WSB) red + shaded · "
  "1 random coordinate/read (STAR) · 25–32 nt · GenomicRanges::follow 3′→5′ adjacency · E16.5→P12.5→P20.5", ha="center", fontsize=7.2, color="#666")
fig.tight_layout()
# classical(blue)/wild(orange) companion: mean +1 phasing per strain (subspecies colour coding)
fig.subplots_adjust(bottom=0.30)
_mean = df.groupby("strain").pct.mean().reindex(order).values
_cax = add_classical_wild_companion(fig, ax, order, _mean, gap=0.12, height_frac=0.22, ylabel="mean\nphasing %", log=False)
_cax.set_xticks(x); _cax.set_xticklabels([s.replace("_", "/") for s in order], rotation=45, ha="right", fontsize=8)
for lab, s in zip(_cax.get_xticklabels(), order): lab.set_color("#C0392B" if s in WILD else "#333")
_cax.set_title("classical (blue) vs wild-derived (orange) — mean +1 phasing per strain", fontsize=8, fontweight="bold", loc="left")
out = f"{BASE}/Fig_phasing_allstrains_1panel"
import os as _os; _SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/01_phasing/data/source_data"; _os.makedirs(_SD,exist_ok=True)
df.to_csv(f"{_SD}/SourceData_Fig_phasing_allstrains_1panel.csv",index=False)
for e in ("pdf", "svg", "png"): fig.savefig(f"{out}.{e}", bbox_inches="tight")
print("wrote", out)
print(df.pivot_table(index="strain", columns="tp", values="pct", aggfunc="mean").reindex(order)[TPO].round(1).to_string())
