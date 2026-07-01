#!/usr/bin/env python3
"""
C57BL_6NJ_pangenome_figures.py

Nature Genetics quality figures for C57BL_6NJ P12.5/P20.5 pangenome starter analysis.

Outputs: analysis/claude_biomni_analysis/C57BL_6NJ_pangenome/
Each figure: PDF + SVG + PNG (≥300 dpi) + matching input CSV.

Wong 2011 colourblind palette:
  blue   #0072B2  — P12.5 / shared_postnatal
  amber  #E69F00  — P20.5 / pachytene
  green  #009E73  — P12.5_only
  pink   #CC79A7  — SV / INS
  sky    #56B4E9  — DEL / unclassified
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
BASE   = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
OUT    = f"{BASE}/analysis/claude_biomni_analysis/C57BL_6NJ_pangenome"
os.makedirs(OUT, exist_ok=True)

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
    "unclassified":      WONG["sky"],
}

def save_fig(fig, stem):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(OUT, f"{stem}.{ext}"),
                    bbox_inches="tight", dpi=300)
    print(f"  Saved {stem}.{{pdf,svg,png}}")

# ── Load data ──────────────────────────────────────────────────────────────
print("Loading data...")
class_stats = pd.read_csv(f"{OUT}/C57BL_6NJ_class_stats.csv")
classified  = pd.read_csv(f"{OUT}/C57BL_6NJ_classified.bed", sep="\t",
                           header=None, names=["chr","start","end","dev_class"])
classified["size_bp"] = classified["end"] - classified["start"]

p12_merged  = pd.read_csv(f"{OUT}/C57BL_6NJ_P12.5_merged.bed", sep="\t",
                            header=None, names=["chr","start","end","mean_FPM","types"])
p20_merged  = pd.read_csv(f"{OUT}/C57BL_6NJ_P20.5_merged.bed", sep="\t",
                            header=None, names=["chr","start","end","mean_FPM","types"])
sv_df       = pd.read_csv(f"{OUT}/C57BL_6NJ_TE_sized_SVs_GRCm39.csv")

print(f"  Classified loci: {len(classified)}")
print(f"  SVs: {len(sv_df)}")
import os as _os; _SD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/05_C57BL_6NJ_pangenome/data/source_data"; _os.makedirs(_SD,exist_ok=True)  # source-data → theme dir

# ── FIGURE 1: Developmental classification overview ─────────────────────────
print("\n=== Figure 1: Developmental classification ===")

fig, axes = plt.subplots(1, 3, figsize=(7, 2.4), constrained_layout=True)

# Panel A: Merge cascade bar chart — no value labels (crowded), use y-axis instead
ax = axes[0]
stages  = ["P12.5\nall reps", "P12.5\nmerged", "P20.5\nall reps", "P20.5\nmerged", "Master\nunion"]
# all-reps counts come from PICB xlsx (3 reps each, unchanged); merged/union from loaded data
values  = [41448, len(p12_merged), 8052, len(p20_merged), len(classified)]
colors  = [WONG["blue"],      WONG["blue"],      WONG["amber"],     WONG["amber"],    WONG["black"]]
alphas  = [0.4,               1.0,               0.4,               1.0,              0.7]
for i, (v, c, a) in enumerate(zip(values, colors, alphas)):
    ax.bar(i, v, color=c, alpha=a, width=0.7)
ax.set_xticks(range(5))
ax.set_xticklabels(stages, fontsize=5.5)
ax.set_ylabel("N clusters", fontsize=7)
ax.set_title("A  Merge cascade", fontsize=7, fontweight="bold", loc="left")
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
# Legend manually
from matplotlib.patches import Patch as P
ax.legend([P(facecolor=WONG["blue"],alpha=0.4), P(facecolor=WONG["blue"]),
           P(facecolor=WONG["amber"],alpha=0.4), P(facecolor=WONG["amber"]),
           P(facecolor=WONG["black"],alpha=0.7)],
          ["P12.5 reps","P12.5 merged","P20.5 reps","P20.5 merged","Union"],
          fontsize=4.5, frameon=False, ncol=2, loc="upper right")

pd.DataFrame({"stage": stages, "n_clusters": values}).to_csv(
    f"{_SD}/SourceData_Fig1A_merge_cascade.csv", index=False)

# Panel B: Developmental classification donut — clean, no in-plot labels
ax = axes[1]
order    = ["P12.5_only", "shared_postnatal", "pachytene", "unclassified"]
counts   = [classified["dev_class"].value_counts().get(c, 0) for c in order]
colors_b = [CLASS_COLORS[c] for c in order]
ax.pie(counts, colors=colors_b, startangle=90,
       wedgeprops={"linewidth": 0.5, "edgecolor": "white"},
       counterclock=False)
centre_circle = plt.Circle((0,0), 0.58, fc="white")
ax.add_artist(centre_circle)
total_loci = len(classified)
ax.text(0, 0.07, f"{total_loci:,}", ha="center", va="center", fontsize=6.5, fontweight="bold")
ax.text(0, -0.12, "loci", ha="center", va="center", fontsize=6)
ax.set_title("B  Dev. classification", fontsize=7, fontweight="bold", loc="left")
labels_b = [
    f"P12.5 only ({100*counts[0]/total_loci:.1f}%)",
    f"Shared postnatal ({100*counts[1]/total_loci:.1f}%)",
    f"Pachytene ({100*counts[2]/total_loci:.1f}%)",
    f"Unclassified ({100*counts[3]/total_loci:.1f}%)" if counts[3] > 0 else "Unclassified (0%)",
]
ax.legend([P(facecolor=c) for c in colors_b], labels_b,
          loc="lower center", bbox_to_anchor=(0.5, -0.35),
          ncol=1, fontsize=5, frameon=False)

pd.DataFrame({"dev_class": order, "n_loci": counts, "label": labels_b}).to_csv(
    f"{_SD}/SourceData_Fig1B_classification.csv", index=False)

# Panel C: Cluster size violin — medians annotated only if they don't overlap
ax = axes[2]
plot_classes = ["shared_postnatal", "pachytene", "P12.5_only"]
data_c = [np.log10(classified[classified["dev_class"]==c]["size_bp"] + 1).values
          for c in plot_classes]
vp = ax.violinplot(data_c, positions=[0,1,2], showmedians=True, showextrema=False)
for body, c in zip(vp["bodies"], [CLASS_COLORS[cc] for cc in plot_classes]):
    body.set_facecolor(c); body.set_alpha(0.8)
vp["cmedians"].set_color("black"); vp["cmedians"].set_linewidth(1.0)
ax.set_xticks([0,1,2])
ax.set_xticklabels(["Shared\npostnatal", "Pachytene", "P12.5\nonly"], fontsize=5.5)
ax.set_ylabel(r"$\log_{10}$ cluster size (bp)", fontsize=7)   # mathtext subscript (Liberation Sans lacks U+2080-2089)
ax.set_title("C  Cluster size dist.", fontsize=7, fontweight="bold", loc="left")
# Median values as text in the right margin, not overlapping violin
medians = [classified[classified["dev_class"]==c]["size_bp"].median()
           for c in plot_classes]
for i, med in enumerate(medians):
    ax.annotate(f"{int(med):,} bp", xy=(i, np.log10(med+1)),
                xytext=(i + 0.32, np.log10(med+1)),
                fontsize=4.8, va="center",
                arrowprops=dict(arrowstyle="-", color="grey", lw=0.4))

pd.concat([classified[classified["dev_class"]==c][["dev_class","size_bp"]]
           for c in plot_classes]).to_csv(f"{_SD}/SourceData_Fig1C_size_distributions.csv", index=False)

save_fig(fig, "Fig1_developmental_classification")
plt.close()

# ── FIGURE 2: FPM expression characteristics ───────────────────────────────
print("\n=== Figure 2: FPM expression ===")

fig, axes = plt.subplots(1, 3, figsize=(7, 2.4), constrained_layout=True)

# Panel A: P12.5 vs P20.5 FPM violin
ax = axes[0]
log_fpm_p12 = np.log10(p12_merged["mean_FPM"] + 0.01)
log_fpm_p20 = np.log10(p20_merged["mean_FPM"] + 0.01)
vp2 = ax.violinplot([log_fpm_p12, log_fpm_p20], positions=[0,1],
                     showmedians=True, showextrema=False)
for body, c in zip(vp2["bodies"], [WONG["blue"], WONG["amber"]]):
    body.set_facecolor(c); body.set_alpha(0.8)
vp2["cmedians"].set_color("black"); vp2["cmedians"].set_linewidth(1.0)
ax.set_xticks([0,1])
ax.set_xticklabels([f"P12.5\n({len(p12_merged):,})", f"P20.5\n({len(p20_merged):,})"], fontsize=5.5)
ax.set_ylabel(r"$\log_{10}$(mean FPM)", fontsize=7)   # mathtext subscript (Liberation Sans lacks U+2080-2089)
ax.set_title("A  FPM distribution", fontsize=7, fontweight="bold", loc="left")
# Median as side annotation
for i, (data, tp) in enumerate(zip([log_fpm_p12, log_fpm_p20], ["P12.5","P20.5"])):
    med_val = 10 ** float(np.median(data))
    ax.annotate(f"med {med_val:.1f}", xy=(i, np.median(data)),
                xytext=(i + 0.28, np.median(data)),
                fontsize=4.8, va="center",
                arrowprops=dict(arrowstyle="-", color="grey", lw=0.4))

pd.DataFrame({"timepoint": ["P12.5"]*len(log_fpm_p12) + ["P20.5"]*len(log_fpm_p20),
              "log10_FPM": list(log_fpm_p12) + list(log_fpm_p20)}).to_csv(
    f"{_SD}/SourceData_Fig2A_FPM_timepoints.csv", index=False)

# Panel B: Per-class median FPM grouped bar
ax = axes[1]
cs_plot = class_stats[class_stats.dev_class.isin(
    ["shared_postnatal","pachytene","P12.5_only"])].copy()
x = np.arange(len(cs_plot))
w = 0.35
ax.bar(x - w/2, cs_plot["median_FPM_P12.5"].fillna(0), w,
       color=WONG["blue"], label="P12.5", alpha=0.9)
ax.bar(x + w/2, cs_plot["median_FPM_P20.5"].fillna(0), w,
       color=WONG["amber"], label="P20.5", alpha=0.9)
ax.set_xticks(x)
ax.set_xticklabels(["Shared\npostnatal", "Pachytene", "P12.5\nonly"], fontsize=5.5)
ax.set_ylabel("Median FPM", fontsize=7)
ax.set_title("B  Class FPM (median)", fontsize=7, fontweight="bold", loc="left")
ax.legend(fontsize=5.5, frameon=False, loc="upper right")

cs_plot[["dev_class","median_FPM_P12.5","median_FPM_P20.5"]].to_csv(
    f"{_SD}/SourceData_Fig2B_class_FPM.csv", index=False)

# Panel C: N loci per class (bar, no count labels — y-axis is sufficient)
ax = axes[2]
cs_all = class_stats[class_stats.dev_class.isin(
    ["shared_postnatal","pachytene","P12.5_only"])].copy()
colors_c = [CLASS_COLORS[c] for c in cs_all["dev_class"]]
ax.bar(range(len(cs_all)), cs_all["n_loci"], color=colors_c, alpha=0.9, width=0.6)
ax.set_xticks(range(len(cs_all)))
ax.set_xticklabels(["Shared\npostnatal", "Pachytene", "P12.5\nonly"], fontsize=5.5)
ax.set_ylabel("N loci", fontsize=7)
ax.set_title("C  Loci per class", fontsize=7, fontweight="bold", loc="left")
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

cs_all[["dev_class","n_loci"]].to_csv(f"{_SD}/SourceData_Fig2C_loci_per_class.csv", index=False)

save_fig(fig, "Fig2_FPM_expression")
plt.close()

# ── FIGURE 3: Pangenome SV content ─────────────────────────────────────────
print("\n=== Figure 3: Pangenome SV content ===")

fig, axes = plt.subplots(1, 3, figsize=(7, 2.4), constrained_layout=True)

# Panel A: SV size histogram — dual y-axis (INS left, DEL right)
# INS (41,585) overwhelms DEL (74) on shared axis; split y-axes to show both
ax = axes[0]
ins  = sv_df[sv_df.sv_type == "INS"]["sv_size_bp"].values
dels = sv_df[sv_df.sv_type == "DEL"]["sv_size_bp"].values
bins = np.logspace(np.log10(300), np.log10(sv_df["sv_size_bp"].max()), 40)

ax.hist(ins, bins=bins, color=WONG["pink"], alpha=0.8, label=f"INS (n={len(ins):,})")
ax.set_xscale("log")
ax.set_xlabel("SV size (bp)", fontsize=7)
ax.set_ylabel("INS count", fontsize=7, color=WONG["pink"])
ax.tick_params(axis="y", labelcolor=WONG["pink"], labelsize=6)

ax2_sv = ax.twinx()
ax2_sv.hist(dels, bins=bins, color=WONG["sky"], alpha=0.9,
            histtype="step", linewidth=1.2, label=f"DEL (n={len(dels):,})")
ax2_sv.set_ylabel("DEL count", fontsize=7, color=WONG["sky"])
ax2_sv.tick_params(axis="y", labelcolor=WONG["sky"], labelsize=6)
ax2_sv.spines["right"].set_visible(True)
ax2_sv.spines["right"].set_linewidth(0.5)

# Combined legend
lines_a, labs_a = ax.get_legend_handles_labels()
lines_b, labs_b = ax2_sv.get_legend_handles_labels()
ax.legend(lines_a + lines_b, labs_a + labs_b, fontsize=5, frameon=False,
          loc="upper right")
ax.set_title("A  SV size (dual axis)", fontsize=7, fontweight="bold", loc="left")

sv_df[["sv_type","sv_size_bp"]].to_csv(f"{_SD}/SourceData_Fig3A_SV_size.csv", index=False)

# Panel B: SV size categories (bar, no per-bar labels — y-axis suffices)
ax = axes[1]
bins_b = [(300,500,"300–499 bp"), (500,1000,"0.5–1 kb"), (1000,5000,"1–5 kb"),
          (5000,10000,"5–10 kb"), (10000,100000,"10–100 kb"), (100000,10**9,">100 kb")]
labels_b = [l for _, _, l in bins_b]
counts_b = [int(((sv_df["sv_size_bp"] >= lo) & (sv_df["sv_size_bp"] < hi)).sum())
            for lo, hi, _ in bins_b]
te_colors = [WONG["sky"], WONG["blue"], WONG["amber"], WONG["green"],
             WONG["pink"], WONG["vermillion"]]
ax.bar(range(len(labels_b)), counts_b, color=te_colors, alpha=0.9, width=0.7)
ax.set_xticks(range(len(labels_b)))
ax.set_xticklabels(labels_b, fontsize=5.5, rotation=35, ha="right")
ax.set_ylabel("N SVs", fontsize=7)
ax.set_title("B  SV size categories", fontsize=7, fontweight="bold", loc="left")
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

pd.DataFrame({"size_category": labels_b, "n_svs": counts_b}).to_csv(
    f"{_SD}/SourceData_Fig3B_SV_categories.csv", index=False)

# Panel C: Per-chromosome SV counts (bar, every 3rd chr labelled to prevent overlap)
ax = axes[2]
chr_order   = [str(i) for i in range(1,20)] + ["X"]
chr_counts_c = [int(sv_df[sv_df.chr == c].shape[0]) for c in chr_order]
ax.bar(range(len(chr_order)), chr_counts_c, color=WONG["blue"], alpha=0.8, width=0.8)
ax.set_xticks(range(len(chr_order)))
# Label every 3rd chromosome to prevent overlap; always show X
chr_tick_labels = [c if (i % 3 == 0 or c == "X") else "" for i, c in enumerate(chr_order)]
ax.set_xticklabels(chr_tick_labels, fontsize=5, rotation=45, ha="right")
ax.set_xlabel("Chromosome", fontsize=7)
ax.set_ylabel("N SVs (≥300 bp)", fontsize=7)
ax.set_title("C  SVs per chromosome", fontsize=7, fontweight="bold", loc="left")
ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))

pd.DataFrame({"chr": chr_order, "n_svs": chr_counts_c}).to_csv(
    f"{_SD}/SourceData_Fig3C_SV_per_chromosome.csv", index=False)

save_fig(fig, "Fig3_pangenome_SV_content")
plt.close()

print("\n=== ALL FIGURES DONE ===")
print(f"Output directory: {OUT}")
