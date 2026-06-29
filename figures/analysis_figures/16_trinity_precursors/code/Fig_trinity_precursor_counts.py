#!/usr/bin/env python3
"""THEME 16 — Trinity-based piRNA precursors across 16 strains (thesis Ch.6, Fig 6.3 analog).
A precursor = a Trinity de-novo contig >=500 bp with sRNA coverage 100 RPM AND 100 RPKM (thesis Ch.6, verbatim:
"stringent filters of 100 RPM (piRNA) and 100 RPKM ... to the transcripts"). This equals the saved all_*.csv.gz
(rpm>100 & rpkm>100). NB: the committed piRNA_filter_precursors.sh instead uses length>100 & rpm>100 (filters length,
not rpkm) — a code bug vs the thesis rule. One grouped bar per strain (canonical order) x 3 developmental windows;
bar = mean over the 3 replicates, error = ±SD, dots = replicates; wild-derived strains shaded + red labels.
Developmental pattern: P12.5 peaks, P20.5 fewest (few-but-dominant pachytene precursors — cf. thesis Fig 6.6 /
Li 2013: pachytene piRNA ~95% of piRNA at pachytene). CAVEAT (thesis Fig 6.4/6.5): Trinity FRAGMENTS precursors
vs the Zamore reference, so absolute counts are an OVER-estimate; the per-strain/developmental CONTRAST is the signal.
Input: results/trinity/filter_precursors/all_trinity_filtred_precursors.csv.gz."""
import sys
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from strain_order import STRAIN_ORDER, WILD
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
SD = ROOT + "/analysis/claude_biomni_analysis/source_data"
d = pd.read_csv(ROOT + "/results/trinity/filter_precursors/all_trinity_filtred_precursors.csv.gz",
                header=None, names=["ref", "count", "length", "rpm", "rpkm", "sample"])
d = d[(d["rpm"] > 100) & (d["rpkm"] > 100)].copy()    # thesis Ch.6 coverage rule: 100 RPM & 100 RPKM (== the saved all_ file)
s = d["sample"].str.split("/").str[-1].str.replace(".500.csv", "", regex=False)
def _strain(x):
    for tk in ("16.5dpc", "12.5dpp", "20.5dpp"):
        if f"-{tk}." in x: return x.split(f"-{tk}.")[0]
    return x
def _tp(x): return "E16.5" if "16.5dpc" in x else "P12.5" if "12.5dpp" in x else "P20.5"
def _rep(x): return x.rsplit(".", 1)[-1]
d["strain"], d["tp"], d["rep"] = s.map(_strain), s.map(_tp), s.map(_rep)
pc = d.groupby(["strain", "tp", "rep"]).size().reset_index(name="n")     # precursors per replicate
TPO = ["E16.5", "P12.5", "P20.5"]; COL = {"E16.5": "#F0C9A0", "P12.5": "#E69F00", "P20.5": "#B4500A"}
TPL = {"E16.5": "E16.5 (prepachytene)", "P12.5": "P12.5", "P20.5": "P20.5 (pachytene)"}
order = [x for x in STRAIN_ORDER if x in set(pc.strain)]
WPOS = [i for i, x in enumerate(order) if x in WILD]
x = np.arange(len(order)); bw = 0.27; rng = np.random.default_rng(0)
fig, ax = plt.subplots(figsize=(15, 5.6), dpi=300)
if WPOS: ax.axvspan(min(WPOS) - 0.5, max(WPOS) + 0.5, color="#C0392B", alpha=0.06, zorder=0)
for j, t in enumerate(TPO):
    xs = x + (j - 1) * bw
    for xi, st in zip(xs, order):
        v = pc[(pc.strain == st) & (pc.tp == t)].n.values
        m = v.mean() if len(v) else 0
        ax.bar(xi, m, bw, color=COL[t], edgecolor="white", linewidth=0.3, zorder=2)
        if len(v):
            sd = v.std(ddof=1) if len(v) > 1 else 0
            ax.errorbar(xi, m, yerr=sd, fmt="none", ecolor="#333", elinewidth=0.5, capsize=1.5, zorder=3)
            ax.scatter(np.full(len(v), xi) + rng.uniform(-0.06, 0.06, len(v)), v, s=8, color="#222", zorder=4, linewidths=0)
            ax.text(xi, max(v.max(), m + sd) + 80, f"{m:.0f}", ha="center", va="bottom", fontsize=4.4, rotation=90, color=COL[t], fontweight="bold", zorder=5)
ax.set_xticks(x); ax.set_xticklabels([st.replace("_", "/") for st in order], rotation=45, ha="right", fontsize=8.5)
for lab, st in zip(ax.get_xticklabels(), order): lab.set_color("#C0392B" if st in WILD else "#222")
ax.set_ylabel("Trinity piRNA precursors (≥500 bp, 100 RPM & 100 RPKM — thesis Ch.6)", fontsize=10)
ax.spines[["top", "right"]].set_visible(False)
ax.set_ylim(top=ax.get_ylim()[1] * 1.07)   # headroom so the tall NZO P12.5 ±SD whisker + value labels clear the frame
ax.legend(handles=[Patch(facecolor=COL[t], label=TPL[t]) for t in TPO], fontsize=9, frameon=False, ncol=1, loc="upper left", title="developmental window", title_fontsize=9)
if WPOS: ax.text(np.mean(WPOS), ax.get_ylim()[1] * 0.985, "wild-derived", ha="center", va="top", fontsize=9, fontweight="bold", color="#C0392B")
ax.set_title("Trinity-based piRNA precursors across 16 mouse strains × development  (thesis Ch.6, Fig 6.3 analog)", fontsize=12.5, fontweight="bold", loc="left")
fig.text(0.5, -0.02, "precursor = Trinity de-novo contig ≥500 bp with 100 RPM & 100 RPKM (thesis Ch.6) · bar = mean of 3 replicates · error = ±SD · dots = reps · P12.5 peaks, P20.5 fewest "
  "(few-but-dominant pachytene precursors) · CAVEAT: Trinity fragments precursors vs the Zamore reference → absolute counts over-estimate (thesis Fig 6.4/6.5); the per-strain/stage CONTRAST is the signal.",
  ha="center", fontsize=6.6, color="#666")
fig.tight_layout()
out = ROOT + "/analysis/claude_biomni_analysis/Fig_trinity_precursor_counts"
for e in ("pdf", "svg", "png"): fig.savefig(f"{out}.{e}", bbox_inches="tight")
# source data
src = pc.pivot_table(index="strain", columns=["tp", "rep"], values="n").reindex(order)
src.to_csv(f"{SD}/SourceData_Fig_trinity_precursor_counts.csv")
print(pc.groupby(["strain", "tp"]).n.mean().unstack().reindex(order)[TPO].round(0).to_string())
print("wrote", out)
