#!/usr/bin/env python3
"""
Nature Genetics-standard figure for piRNA phasing distance distribution.

Method: Almeida et al., Genome Biol 2025 (PMID 39844208) — distance from the 3' end
of each piRNA to the 5' end of the next same-strand piRNA; phased (Zucchini-dependent)
biogenesis = enrichment at +1 nt. Computed by workflow/scripts/R/phasing_analysis.R
(25-33 nt mouse window).

Usage:
  python Fig_phasing.py <phasing.csv> <out_prefix> ["Title"]
  (reads <prefix-of-csv>.phasing_summary.csv automatically if present)

Outputs: <out_prefix>.{pdf,svg,png}  (vector + 300 dpi raster)
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Nature Genetics style: Liberation Sans (Arial not installed), small, vector ----
plt.rcParams.update({
    "font.family": "Liberation Sans",
    "font.size": 7,
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "pdf.fonttype": 42,   # editable text in vector output
    "svg.fonttype": "none",
})

# Okabe-Ito colourblind-safe palette
VERMILLION = "#D55E00"   # the +1 phased bar
GREY       = "#BFBFBF"   # background distances

csv      = sys.argv[1]
out_pref = sys.argv[2]
title    = sys.argv[3] if len(sys.argv) > 3 else os.path.basename(out_pref)

df = pd.read_csv(csv)
summ_path = csv.replace(".phasing.csv", ".phasing_summary.csv")
n_pairs = z1 = frac1 = None
if os.path.exists(summ_path):
    s = pd.read_csv(summ_path).iloc[0]
    n_pairs, z1, frac1 = int(s["n_pairs"]), float(s["zscore_plus1"]), float(s["frac_plus1"])

XMAX = 30
d = df[df["distance"] <= XMAX].copy()
d["pct"] = 100.0 * d["count"] / df["count"].sum()
colors = [VERMILLION if x == 1 else GREY for x in d["distance"]]

fig, ax = plt.subplots(figsize=(3.4, 2.5), dpi=300)
ax.bar(d["distance"], d["pct"], color=colors, width=0.82, linewidth=0)

ax.set_xlabel("Distance, piRNA 3′ end → next 5′ end (nt)", fontsize=7.5)
ax.set_ylabel("Adjacent piRNA pairs (%)", fontsize=7.5)
ax.set_title(title, fontsize=8, fontweight="bold", pad=6)
ax.set_xlim(0.3, XMAX + 0.7)
ax.set_xticks([1, 5, 10, 15, 20, 25, 30])

# annotate the +1 phased peak
p1 = d.loc[d["distance"] == 1, "pct"].values[0]
ax.annotate("+1 nt\n(phased)", xy=(1, p1), xytext=(6.5, p1 * 0.86),
            fontsize=7, color=VERMILLION, va="center", ha="left",
            arrowprops=dict(arrowstyle="-", color=VERMILLION, lw=0.7))

# stats box (no overlap: upper-right, framed light)
lines = []
if frac1 is not None: lines.append(f"+1 nt = {frac1*100:.1f}% of pairs")
if z1   is not None: lines.append(f"phasing z(+1) = {z1:.2f}")
if n_pairs is not None: lines.append(f"n = {n_pairs/1e6:.1f}M pairs")
lines.append("25–33 nt · primary")
ax.text(0.97, 0.97, "\n".join(lines), transform=ax.transAxes,
        fontsize=6.3, va="top", ha="right", color="#222222",
        linespacing=1.45,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#CCCCCC", lw=0.5))

fig.tight_layout(pad=0.5)
for ext in ("pdf", "svg", "png"):
    fig.savefig(f"{out_pref}.{ext}", bbox_inches="tight")
print(f"wrote {out_pref}.{{pdf,svg,png}}  (+1={p1:.1f}% of shown, z={z1})")
