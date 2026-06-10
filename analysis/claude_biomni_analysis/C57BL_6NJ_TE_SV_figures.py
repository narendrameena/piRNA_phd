#!/usr/bin/env python3
"""
C57BL_6NJ_TE_SV_figures.py

TE content and SV analysis figures for C57BL_6NJ piRNA clusters.

Figure 4: TE annotation per developmental class (RepeatMasker × PICB clusters)
Figure 5: Timepoint-unique cluster analysis
           - P12.5 exclusive (P12.5_only class)
           - P20.5 exclusive (pachytene class)
           - Shared postnatal
           - TE landscape and SV context

NOTE on SV coordinates:
  PICB cluster BEDs → C57BL_6NJ REL-2205 assembly (UCSC chr1..chrX)
  Pangenome VCF SVs → GRCm39 coordinates
  TE RepeatMasker BED → C57BL_6NJ REL-2205 assembly (same space as clusters)
  The RepeatMasker BED IS in the cluster coordinate space → direct intersection valid.
  The VCF SVs are NOT in the cluster coordinate space → used only genome-wide, not per-cluster.

Outputs: analysis/claude_biomni_analysis/C57BL_6NJ_pangenome/
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────
BASE = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
OUT  = f"{BASE}/analysis/claude_biomni_analysis/C57BL_6NJ_pangenome"

# ── Style ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":        "Liberation Sans",
    "font.size":          7,
    "axes.linewidth":     0.5,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "xtick.major.width":  0.5,
    "ytick.major.width":  0.5,
    "xtick.major.size":   2.5,
    "ytick.major.size":   2.5,
    "xtick.labelsize":    7,
    "ytick.labelsize":    7,
    "figure.dpi":         300,
    "savefig.dpi":        300,
    "pdf.fonttype":       42,
    "svg.fonttype":       "none",
    "figure.facecolor":   "white",
    "axes.facecolor":     "white",
})

WONG = {
    "blue":       "#0072B2",
    "amber":      "#E69F00",
    "green":      "#009E73",
    "pink":       "#CC79A7",
    "sky":        "#56B4E9",
    "yellow":     "#F0E442",
    "vermillion": "#D55E00",
    "black":      "#000000",
}

CLASS_COLORS = {
    "shared_postnatal": WONG["blue"],
    "pachytene":         WONG["amber"],
    "P12.5_only":        WONG["green"],
}

TE_COLORS = {
    "LINE/L1":        WONG["vermillion"],
    "LTR/ERVK":       WONG["pink"],
    "LTR/ERVL-MaLR":  WONG["amber"],
    "LTR/ERV1":       WONG["yellow"],
    "SINE/Alu":       WONG["sky"],
    "SINE/B2":        WONG["blue"],
    "SINE/B4":        WONG["green"],
    "Simple_repeat":  "#BBBBBB",
    "Non-TE":         "#EEEEEE",
}

def save_fig(fig, stem):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(OUT, f"{stem}.{ext}"),
                    bbox_inches="tight", dpi=300)
    print(f"  Saved {stem}.{{pdf,svg,png}}")

# ── Load data ──────────────────────────────────────────────────────────────
print("Loading data...")
te_ann   = pd.read_csv(f"{OUT}/C57BL_6NJ_TE_annotation.csv")
sv_df    = pd.read_csv(f"{OUT}/C57BL_6NJ_TE_sized_SVs_GRCm39.csv")
cs_stats = pd.read_csv(f"{OUT}/C57BL_6NJ_class_stats.csv")

# Load per-class TE breakdown CSVs
te_break = {}
for dc in ("shared_postnatal", "pachytene", "P12.5_only"):
    path = f"{OUT}/C57BL_6NJ_{dc}_TE_breakdown.csv"
    if os.path.exists(path):
        te_break[dc] = pd.read_csv(path)

print("  TE annotation loaded")

# ── FIGURE 4: TE content per developmental class ───────────────────────────
print("\n=== Figure 4: TE content per developmental class ===")

fig, axes = plt.subplots(1, 3, figsize=(7, 2.4), constrained_layout=True)

# Panel A: Total TE% per class (bar chart)
ax = axes[0]
classes  = ["shared_postnatal", "pachytene", "P12.5_only"]
te_pcts  = [float(te_ann[te_ann.dev_class==c]["TE_fraction_pct"].iloc[0]) for c in classes]
n_loci   = [int(te_ann[te_ann.dev_class==c]["n_loci"].iloc[0]) for c in classes]
colors_a = [CLASS_COLORS[c] for c in classes]
bars = ax.bar([0,1,2], te_pcts, color=colors_a, alpha=0.9, width=0.6)
ax.set_xticks([0,1,2])
ax.set_xticklabels(["Shared\npostnatal", "Pachytene", "P12.5\nonly"], fontsize=5.5)
ax.set_ylabel("% cluster bp overlapping TE", fontsize=7)
ax.set_ylim(0, 55)
ax.set_title("A  TE coverage per class", fontsize=7, fontweight="bold", loc="left")
# n= labels inside bars
for i, (p, n) in enumerate(zip(te_pcts, n_loci)):
    ax.text(i, p + 1.2, f"{p:.1f}%", ha="center", va="bottom", fontsize=5.5)
    ax.text(i, 2, f"n={n:,}", ha="center", va="bottom", fontsize=4.5, color="white",
            fontweight="bold")

pd.DataFrame({"dev_class": classes, "TE_fraction_pct": te_pcts, "n_loci": n_loci}).to_csv(
    f"{OUT}/Fig4A_TE_fraction.csv", index=False)

# Panel B: Stacked bar — TE class composition per developmental class
ax = axes[1]
# Key TE families to show (stacked); anything else = "Other TE"
te_families = ["LINE/L1", "LTR/ERVK", "LTR/ERVL-MaLR", "SINE/Alu", "SINE/B2", "SINE/B4"]
te_stack = {}
for dc in classes:
    row = {"dev_class": dc}
    if dc in te_break:
        tb = te_break[dc].set_index("te_class")["pct_of_total_cluster_bp"].to_dict()
    else:
        tb = {}
    for fam in te_families:
        row[fam] = tb.get(fam, 0.0)
    # Non-TE
    te_pct = te_ann[te_ann.dev_class==dc]["TE_fraction_pct"].iloc[0]
    row["Non-TE"] = max(0, 100 - te_pct)
    te_stack[dc] = row

stack_df = pd.DataFrame(list(te_stack.values())).set_index("dev_class")
# Columns: TE families + Non-TE
col_order = te_families + ["Non-TE"]
bottom = np.zeros(3)
x = np.arange(3)
for col in col_order:
    vals = [stack_df.loc[c, col] for c in classes]
    color = TE_COLORS.get(col, "#AAAAAA")
    ax.bar(x, vals, bottom=bottom, color=color,
           label=col, width=0.6, alpha=0.95 if col != "Non-TE" else 0.3)
    bottom += np.array(vals)

ax.set_xticks([0,1,2])
ax.set_xticklabels(["Shared\npostnatal", "Pachytene", "P12.5\nonly"], fontsize=5.5)
ax.set_ylabel("% cluster bp", fontsize=7)
ax.set_title("B  TE composition (stacked)", fontsize=7, fontweight="bold", loc="left")
ax.set_ylim(0, 115)
ax.legend(loc="upper right", fontsize=4.5, frameon=False, ncol=1,
          bbox_to_anchor=(1.0, 1.0))

stack_df.reset_index().to_csv(f"{OUT}/Fig4B_TE_composition.csv", index=False)

# Panel C: TE class enrichment — fold change of TE% vs genomic background
# Use genome-wide TE fraction as background (from RepeatMasker total)
# C57BL_6NJ genome ~2.5 Gb; typical mouse genome is ~42% TE
# Use ratio (class TE%) / (overall cluster TE%) to show relative enrichment
ax = axes[2]
# Compute background: all classified clusters combined
total_bp_all  = sum(te_ann[te_ann.dev_class==c]["total_bp"].iloc[0] for c in classes)
total_te_all  = sum(te_ann[te_ann.dev_class==c]["TE_bp"].iloc[0] for c in classes)
bg_pct = 100 * total_te_all / total_bp_all

enrichments = [(te_pcts[i] / bg_pct) for i in range(3)]
ax.axhline(1.0, color="grey", lw=0.7, ls="--", alpha=0.7)
bars2 = ax.bar([0,1,2], enrichments, color=colors_a, alpha=0.9, width=0.6)
ax.set_xticks([0,1,2])
ax.set_xticklabels(["Shared\npostnatal", "Pachytene", "P12.5\nonly"], fontsize=5.5)
ax.set_ylabel("TE enrichment\n(vs. mean cluster TE%)", fontsize=7)
ax.set_title("C  Relative TE enrichment", fontsize=7, fontweight="bold", loc="left")
ax.set_ylim(0, 1.5)
for i, e in enumerate(enrichments):
    ax.text(i, e + 0.02, f"{e:.2f}×", ha="center", va="bottom", fontsize=5.5)

pd.DataFrame({"dev_class": classes, "TE_fraction_pct": te_pcts,
              "background_pct": [bg_pct]*3,
              "enrichment_fold": enrichments}).to_csv(
    f"{OUT}/Fig4C_TE_enrichment.csv", index=False)

save_fig(fig, "Fig4_TE_content")
plt.close()

# ── FIGURE 5: Timepoint-unique cluster analysis ────────────────────────────
print("\n=== Figure 5: Timepoint-unique cluster analysis ===")

fig, axes = plt.subplots(1, 3, figsize=(7, 2.4), constrained_layout=True)

# Panel A: Cluster counts by timepoint exclusivity — read from class_stats.csv
ax = axes[0]
_cs = cs_stats.set_index("dev_class")
p12_excl = int(_cs.loc["P12.5_only",    "n_loci"])
p12_shar = int(_cs.loc["shared_postnatal", "n_loci"])
p20_excl = int(_cs.loc["pachytene",     "n_loci"])

# Stacked bar showing P12.5 and P20.5 merged clusters by exclusive/shared content
ax.bar([0], [p12_excl], color=WONG["green"],  alpha=0.9, label="P12.5 exclusive", width=0.5)
ax.bar([0], [p12_shar], bottom=[p12_excl], color=WONG["blue"],  alpha=0.9,
       label="Shared postnatal", width=0.5)
ax.bar([1], [p20_excl], color=WONG["amber"],  alpha=0.9, label="P20.5 exclusive", width=0.5)
ax.bar([1], [p12_shar], bottom=[p20_excl], color=WONG["blue"],  alpha=0.9, width=0.5)
ax.set_xticks([0,1])
ax.set_xticklabels(["P12.5\nmerged", "P20.5\nmerged"], fontsize=6)
ax.set_ylabel("N clusters", fontsize=7)
ax.set_title("A  Cluster exclusivity", fontsize=7, fontweight="bold", loc="left")
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
ax.legend(fontsize=4.5, frameon=False, loc="upper right")

pd.DataFrame({
    "timepoint": ["P12.5", "P20.5"],
    "exclusive": [p12_excl, p20_excl],
    "shared":    [p12_shar, p12_shar],
    "total":     [p12_excl+p12_shar, p20_excl+p12_shar]
}).to_csv(f"{OUT}/Fig5A_cluster_exclusivity.csv", index=False)

# Panel B: Cluster size by class (box plot — cleaner than violin for comparison)
ax = axes[1]
from pathlib import Path
classified = pd.read_csv(f"{OUT}/C57BL_6NJ_classified.bed", sep="\t",
                          header=None, names=["chr","start","end","dev_class"])
classified["size_bp"] = classified["end"] - classified["start"]

plot_classes  = ["P12.5_only", "shared_postnatal", "pachytene"]
plot_labels   = ["P12.5\nexclusive", "Shared\npostnatal", "P20.5\nexclusive\n(pachytene)"]
bplot_data    = [np.log10(classified[classified.dev_class==c]["size_bp"].values + 1)
                 for c in plot_classes]
bplot_colors  = [CLASS_COLORS[c] for c in plot_classes]

bp = ax.boxplot(bplot_data, patch_artist=True, notch=False,
                medianprops={"color": "black", "linewidth": 1.0},
                whiskerprops={"linewidth": 0.5},
                capprops={"linewidth": 0.5},
                flierprops={"marker": ".", "markersize": 1.5, "alpha": 0.3},
                boxprops={"linewidth": 0.5})
for patch, c in zip(bp["boxes"], bplot_colors):
    patch.set_facecolor(c); patch.set_alpha(0.8)

ax.set_xticks([1,2,3])
ax.set_xticklabels(plot_labels, fontsize=5)
ax.set_ylabel("log₁₀ cluster size (bp)", fontsize=7)
ax.set_title("B  Cluster size by class", fontsize=7, fontweight="bold", loc="left")

classified[classified.dev_class.isin(plot_classes)][["dev_class","size_bp"]].to_csv(
    f"{OUT}/Fig5B_size_by_class.csv", index=False)

# Panel C: TE superfamily heat-map (row = TE family, col = class)
ax = axes[2]
te_fams = ["LINE/L1", "LTR/ERVK", "LTR/ERVL-MaLR", "SINE/Alu", "SINE/B2", "SINE/B4"]
heat_data = []
for dc in ["P12.5_only", "shared_postnatal", "pachytene"]:
    if dc in te_break:
        tb = te_break[dc].set_index("te_class")["pct_of_total_cluster_bp"].to_dict()
    else:
        tb = {}
    heat_data.append([tb.get(f, 0.0) for f in te_fams])

heat_arr = np.array(heat_data).T   # shape: (n_fams, n_classes)
im = ax.imshow(heat_arr, aspect="auto", cmap="YlOrRd",
               vmin=0, vmax=max(heat_arr.max(), 1))
ax.set_xticks([0,1,2])
ax.set_xticklabels(["P12.5\nexcl.", "Shared", "P20.5\nexcl."], fontsize=5.5)
ax.set_yticks(range(len(te_fams)))
ax.set_yticklabels(te_fams, fontsize=5.5)
plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02, label="% cluster bp")
ax.set_title("C  TE family heatmap", fontsize=7, fontweight="bold", loc="left")
# Annotate cells with values
for i in range(len(te_fams)):
    for j in range(3):
        val = heat_arr[i, j]
        if val > 0.3:
            ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                    fontsize=4, color="black" if val < 7 else "white")

pd.DataFrame(heat_arr, index=te_fams,
             columns=["P12.5_only","shared_postnatal","pachytene"]).to_csv(
    f"{OUT}/Fig5C_TE_heatmap.csv")

save_fig(fig, "Fig5_timepoint_TE_analysis")
plt.close()

print("\n=== ALL DONE ===")
print("Figures: Fig4_TE_content.{pdf,svg,png}")
print("         Fig5_timepoint_TE_analysis.{pdf,svg,png}")
print(f"Output: {OUT}")
