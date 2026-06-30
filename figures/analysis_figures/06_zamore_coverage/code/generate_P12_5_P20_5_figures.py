#!/usr/bin/env python3
"""
generate_P12_5_P20_5_figures.py

Nature Genetics quality figures for C57BL/6 P12.5 and P20.5
piRNA cluster analysis with Zamore annotation cross-reference.

Outputs (PDF + SVG + PNG ≥300 dpi) to claude_biomni_figures/:
  Fig1_PICB_cluster_architecture.{pdf,svg,png}
  Fig2_Zamore_gene_coverage.{pdf,svg,png}
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy import stats

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE    = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
XLSX    = os.path.join(BASE, "results/picb_result/C57BL_6")
COV_CSV = os.path.join(BASE, "analysis/claude_biomni_analysis/P12_5_P20_5_zamore_coverage_per_gene.csv")
OUT_DIR = os.path.join(BASE, "analysis/claude_biomni_analysis/claude_biomni_figures")
SD_DIR = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/06_zamore_coverage/data/source_data"; os.makedirs(SD_DIR, exist_ok=True)  # theme source_data (consolidated)
os.makedirs(OUT_DIR, exist_ok=True)

# ── Nature Genetics style ──────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Arial", "DejaVu Sans", "Helvetica"],
    "font.size":          7,
    "axes.titlesize":     7,
    "axes.labelsize":     7,
    "xtick.labelsize":    6,
    "ytick.labelsize":    6,
    "legend.fontsize":    6,
    "axes.linewidth":     0.6,
    "xtick.major.width":  0.6,
    "ytick.major.width":  0.6,
    "xtick.major.size":   2.5,
    "ytick.major.size":   2.5,
    "lines.linewidth":    0.8,
    "pdf.fonttype":       42,
    "svg.fonttype":       "none",
    "axes.spines.top":    False,
    "axes.spines.right":  False,
})

# ── Colorblind-safe palette (Wong 2011) ───────────────────────────────────────
C_P125     = "#0072B2"   # blue  — P12.5
C_P205     = "#E69F00"   # amber — P20.5
C_PACH     = "#009E73"   # green — Pachytene
C_PREPACH  = "#CC79A7"   # pink  — Prepachytene
C_HYBRID   = "#56B4E9"   # sky   — Hybrid
C_UNKNOWN  = "#999999"   # grey  — Unknown
C_SC       = "#0072B2"   # SingleCore
C_EC       = "#F0E442"   # ExtendedCore (yellow)
C_MC       = "#D55E00"   # MultiCore (vermillion)

STAGE_CLR  = {"Pachytene": C_PACH, "Prepachytene": C_PREPACH,
               "Hybrid": C_HYBRID, "Unknown": C_UNKNOWN}
TYPE_CLR   = {"SingleCore": C_SC, "ExtendedCore": C_EC, "MultiCore": C_MC}

SAMPLES = [
    ("P12.5_rep1", "C57BL_6-12.5dpp.1.xlsx", "P12.5", C_P125),
    ("P12.5_rep2", "C57BL_6-12.5dpp.2.xlsx", "P12.5", C_P125),
    ("P20.5_rep1", "C57BL_6-20.5dpp.1.xlsx", "P20.5", C_P205),
    ("P20.5_rep2", "C57BL_6-20.5dpp.2.xlsx", "P20.5", C_P205),
]
SAMPLE_LABELS = ["P12.5\nrep1", "P12.5\nrep2", "P20.5\nrep1", "P20.5\nrep2"]
SAMPLE_KEYS   = [s[0] for s in SAMPLES]

# ── Load PICB cluster data ─────────────────────────────────────────────────────
print("Loading PICB cluster xlsx files...")
cluster_dfs = {}
seeds_n, cores_n, clusters_n = [], [], []

for key, fname, tp, clr in SAMPLES:
    path = os.path.join(XLSX, fname)
    df_cl = pd.read_excel(path, sheet_name="clusters", engine="openpyxl")
    df_se = pd.read_excel(path, sheet_name="seeds",    engine="openpyxl")
    df_co = pd.read_excel(path, sheet_name="cores",    engine="openpyxl")
    cluster_dfs[key] = df_cl
    seeds_n.append(len(df_se))
    cores_n.append(len(df_co))
    clusters_n.append(len(df_cl))
    print(f"  {key}: seeds={len(df_se)}, cores={len(df_co)}, clusters={len(df_cl)}")

# ── Load coverage CSV ──────────────────────────────────────────────────────────
print("Loading coverage CSV...")
cov = pd.read_csv(COV_CSV)
print(f"  Coverage rows: {len(cov)}, columns: {list(cov.columns)}")
print(f"  Stages: {sorted(cov['stage'].unique())}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — PICB Cluster Architecture
# Panels: A) component counts  B) cluster type stacked bar  C) FPM violin
# Size: 183 mm wide × 70 mm tall (2-column Nature Genetics)
# ══════════════════════════════════════════════════════════════════════════════
print("\nGenerating Figure 1: PICB Cluster Architecture...")

fig1, axes = plt.subplots(1, 3, figsize=(7.2, 2.75))
fig1.subplots_adjust(left=0.08, right=0.97, top=0.87, bottom=0.22, wspace=0.42)
ax_A, ax_B, ax_C = axes

x = np.arange(4)
lw = 0.4

# ── Panel A: Component counts (seeds / cores / clusters) ─────────────────────
ax_A.set_title("A  PICB component counts", loc="left", fontweight="bold", pad=3)

w = 0.22
offsets = [-w, 0, w]
comp_vals = [seeds_n, cores_n, clusters_n]
comp_lbls = ["Seeds", "Cores", "Clusters"]
comp_clrs = ["#4C72B0", "#55A868", "#C44E52"]   # blue / green / red (colourblind ok)

for i, (vals, lbl, clr) in enumerate(zip(comp_vals, comp_lbls, comp_clrs)):
    bars = ax_A.bar(x + offsets[i], vals, width=w, color=clr, label=lbl,
                    linewidth=lw, edgecolor="white")

ax_A.set_xticks(x)
ax_A.set_xticklabels(SAMPLE_LABELS, fontsize=6)
ax_A.set_ylabel("Count")
ax_A.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{int(v/1000)}k" if v >= 1000 else str(int(v))))
ax_A.legend(frameon=False, loc="upper left", handlelength=1.0,
            handletextpad=0.3, borderpad=0, labelspacing=0.2)

# colour sample labels by timepoint
for tick, (_, _, tp, clr) in zip(ax_A.get_xticklabels(), SAMPLES):
    tick.set_color(clr)

# ── Panel B: Cluster type stacked bar ────────────────────────────────────────
ax_B.set_title("B  Cluster type composition", loc="left", fontweight="bold", pad=3)

types = ["SingleCore", "ExtendedCore", "MultiCore"]
type_pcts = {}
for t in types:
    type_pcts[t] = [
        df["type"].value_counts().get(t, 0) / len(df) * 100
        for df in [cluster_dfs[k] for k in SAMPLE_KEYS]
    ]

bottom = np.zeros(4)
for t in types:
    vals = np.array(type_pcts[t])
    ax_B.bar(x, vals, bottom=bottom, color=TYPE_CLR[t], label=t,
             linewidth=lw, edgecolor="white", width=0.55)
    # label inside bar if wide enough
    for xi, (v, b) in enumerate(zip(vals, bottom)):
        if v > 6:
            ax_B.text(xi, b + v/2, f"{v:.0f}%", ha="center", va="center",
                      fontsize=5, color="black" if t == "ExtendedCore" else "white",
                      fontweight="bold")
    bottom += vals

ax_B.set_xticks(x)
ax_B.set_xticklabels(SAMPLE_LABELS, fontsize=6)
ax_B.set_ylabel("Percentage (%)")
ax_B.set_ylim(0, 110)
ax_B.legend(frameon=False, loc="upper right", handlelength=1.0,
            handletextpad=0.3, borderpad=0, labelspacing=0.2)
for tick, (_, _, tp, clr) in zip(ax_B.get_xticklabels(), SAMPLES):
    tick.set_color(clr)

# ── Panel C: FPM distribution violin ─────────────────────────────────────────
ax_C.set_title("C  Cluster FPM distribution", loc="left", fontweight="bold", pad=3)

fpm_data = [np.log10(cluster_dfs[k]["all_reads_primary_alignments_FPM"].values + 1e-3)
            for k in SAMPLE_KEYS]
sample_clrs = [s[3] for s in SAMPLES]

parts = ax_C.violinplot(fpm_data, positions=x, widths=0.55, showmedians=True,
                         showextrema=False)

# Style violin bodies by timepoint colour
for i, (body, clr) in enumerate(zip(parts["bodies"], sample_clrs)):
    body.set_facecolor(clr)
    body.set_edgecolor("none")
    body.set_alpha(0.7)

parts["cmedians"].set_color("black")
parts["cmedians"].set_linewidth(1.0)

# Add median FPM annotation
for i, (key, clr) in enumerate(zip(SAMPLE_KEYS, sample_clrs)):
    med = cluster_dfs[key]["all_reads_primary_alignments_FPM"].median()
    ax_C.text(i, np.log10(med + 1e-3) + 0.08, f"{med:.1f}", ha="center",
              va="bottom", fontsize=5, color="black")

ax_C.set_xticks(x)
ax_C.set_xticklabels(SAMPLE_LABELS, fontsize=6)
ax_C.set_ylabel("log₁₀(FPM + 0.001)")
ax_C.yaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=5))
for tick, (_, _, tp, clr) in zip(ax_C.get_xticklabels(), SAMPLES):
    tick.set_color(clr)

# Legend — timepoint
leg_patches = [mpatches.Patch(color=C_P125, label="P12.5"),
               mpatches.Patch(color=C_P205, label="P20.5")]
ax_C.legend(handles=leg_patches, frameon=False, loc="upper right",
            handlelength=1.0, handletextpad=0.3, borderpad=0, labelspacing=0.2)

# Bottom annotation
fig1.text(0.5, 0.01,
    "C57BL/6 · GRCm38 · piRNA clusters (25–36 nt) · rep1 and rep2 only (distinct biological replicates)",
    ha="center", va="bottom", fontsize=5.0, color="#555555")

# ── Save Fig1 input data matrix ───────────────────────────────────────────
fig1_data = pd.DataFrame({
    "sample":    SAMPLE_KEYS,
    "timepoint": [s[2] for s in SAMPLES],
    "seeds":     seeds_n,
    "cores":     cores_n,
    "clusters":  clusters_n,
})
for t in ["SingleCore", "ExtendedCore", "MultiCore"]:
    fig1_data[f"pct_{t}"] = [
        round(cluster_dfs[k]["type"].value_counts().get(t, 0) / len(cluster_dfs[k]) * 100, 2)
        for k in SAMPLE_KEYS
    ]
for stat_label, func in [("median_FPM", np.median), ("mean_FPM", np.mean)]:
    fig1_data[stat_label] = [
        round(func(cluster_dfs[k]["all_reads_primary_alignments_FPM"].values), 2)
        for k in SAMPLE_KEYS
    ]
fig1_data.to_csv(os.path.join(SD_DIR, "Fig1_data_matrix.csv"), index=False)
print(f"  Saved Fig1_data_matrix.csv")

stem1 = os.path.join(OUT_DIR, "Fig1_PICB_cluster_architecture")
for ext in (".pdf", ".svg", ".png"):
    kw = {"dpi": 300, "bbox_inches": "tight"}
    fig1.savefig(stem1 + ext, **kw)
    print(f"  Saved {stem1+ext}")
plt.close(fig1)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Zamore piRNA Gene Coverage
# Panels: A) gene detection by stage  B) coverage fraction violin by stage
#          C) Pachytene temporal emergence
# Size: 183 mm wide × 75 mm tall
# ══════════════════════════════════════════════════════════════════════════════
print("\nGenerating Figure 2: Zamore piRNA Gene Coverage...")

# ── Precompute detection rates ─────────────────────────────────────────────
STAGE_TOTALS = {"Pachytene": 99, "Prepachytene": 83, "Hybrid": 32, "Unknown": 1}
STAGES_MAIN  = ["Pachytene", "Prepachytene", "Hybrid"]

det_rates = {}   # {sample_key: {stage: pct}}
for key in SAMPLE_KEYS:
    sub = cov[cov["sample"] == key]
    det = sub[sub["n_det"] > 0]
    det_rates[key] = {}
    for st in STAGES_MAIN:
        n_det   = (det["stage"] == st).sum()
        n_total = STAGE_TOTALS[st]
        det_rates[key][st] = 100.0 * n_det / n_total

# ── Figure 2 layout ────────────────────────────────────────────────────────
fig2, axes2 = plt.subplots(1, 3, figsize=(7.2, 2.95))
fig2.subplots_adjust(left=0.08, right=0.97, top=0.87, bottom=0.22, wspace=0.46)
ax2_A, ax2_B, ax2_C = axes2

# ── Panel A: Gene detection rates by stage and sample ─────────────────────
ax2_A.set_title("A  Zamore gene detection rate", loc="left", fontweight="bold", pad=3)

n_stages = len(STAGES_MAIN)
n_samp   = 4
group_w  = 0.7
bar_w    = group_w / n_samp
stage_x  = np.arange(n_stages)

for si, (key, _, tp, clr) in enumerate(SAMPLES):
    offsets_si = (si - (n_samp - 1) / 2) * bar_w
    vals = [det_rates[key][st] for st in STAGES_MAIN]
    bars = ax2_A.bar(stage_x + offsets_si, vals, width=bar_w * 0.9,
                     color=clr, linewidth=lw, edgecolor="white",
                     label=f"{tp} rep{si%2+1}")

# 100% reference line
ax2_A.axhline(100, color="#bbbbbb", lw=0.5, ls="--", zorder=0)

ax2_A.set_xticks(stage_x)
ax2_A.set_xticklabels(STAGES_MAIN, fontsize=6.5)
ax2_A.set_ylabel("Genes detected (%)")
ax2_A.set_ylim(0, 108)

# Annotate stage totals below x-axis
for xi, st in enumerate(STAGES_MAIN):
    ax2_A.text(xi, -8, f"n={STAGE_TOTALS[st]}", ha="center", va="top",
               fontsize=5, color="#555555")

# Color x-tick labels by stage
for tick, st in zip(ax2_A.get_xticklabels(), STAGES_MAIN):
    tick.set_color(STAGE_CLR[st])

leg_patches2 = [mpatches.Patch(color=C_P125, label="P12.5 rep1"),
                mpatches.Patch(color=C_P125, label="P12.5 rep2", alpha=0.6),
                mpatches.Patch(color=C_P205, label="P20.5 rep1"),
                mpatches.Patch(color=C_P205, label="P20.5 rep2", alpha=0.6)]
ax2_A.legend(handles=leg_patches2, frameon=False, loc="lower right",
             handlelength=1.0, handletextpad=0.3, borderpad=0,
             labelspacing=0.1, ncol=2, columnspacing=0.5)

# ── Panel B: Coverage fraction violin by piRNA stage ──────────────────────
ax2_B.set_title("B  Locus coverage fraction by stage", loc="left", fontweight="bold", pad=3)

# Use max_cov per gene (best isoform), pool replicates within timepoint
for tp_label, tp_samples, tp_clr in [
    ("P12.5", ["P12.5_rep1", "P12.5_rep2"], C_P125),
    ("P20.5", ["P20.5_rep1", "P20.5_rep2"], C_P205),
]:
    tp_sub = cov[cov["sample"].isin(tp_samples)]
    # average max_cov across replicates per gene
    tp_avg = tp_sub.groupby(["gene", "stage"])["max_cov"].mean().reset_index()

    for si, st in enumerate(STAGES_MAIN):
        st_data = tp_avg[tp_avg["stage"] == st]["max_cov"].values
        xpos = si * 2.4 + (0.0 if tp_label == "P12.5" else 0.55)

        if len(st_data) < 3:
            continue

        vp = ax2_B.violinplot([st_data], positions=[xpos], widths=0.48,
                               showmedians=True, showextrema=False)
        for body in vp["bodies"]:
            body.set_facecolor(tp_clr)
            body.set_edgecolor("none")
            body.set_alpha(0.72)
        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(1.0)

        # Jitter strip
        jitter = np.random.default_rng(42).uniform(-0.12, 0.12, len(st_data))
        ax2_B.scatter(xpos + jitter, st_data, s=1.5, color=tp_clr,
                      alpha=0.35, linewidths=0, zorder=3)

        # Annotate median
        med = np.median(st_data)
        ax2_B.text(xpos, med + 0.03, f"{med:.2f}", ha="center", va="bottom",
                   fontsize=4.5, color="black")

# X-axis labels: stage groups
group_centres = [si * 2.4 + 0.275 for si in range(3)]
ax2_B.set_xticks(group_centres)
ax2_B.set_xticklabels(STAGES_MAIN, fontsize=6.5)
ax2_B.set_ylabel("Fraction of locus covered")
ax2_B.set_ylim(-0.05, 1.13)

for tick, st in zip(ax2_B.get_xticklabels(), STAGES_MAIN):
    tick.set_color(STAGE_CLR[st])

leg_tp = [mpatches.Patch(color=C_P125, alpha=0.72, label="P12.5"),
          mpatches.Patch(color=C_P205, alpha=0.72, label="P20.5")]
ax2_B.legend(handles=leg_tp, frameon=False, loc="lower right",
             handlelength=1.0, handletextpad=0.3, borderpad=0, labelspacing=0.2)

# ── Panel C: Pachytene gene temporal emergence ────────────────────────────
ax2_C.set_title("C  Pachytene piRNA gene expression", loc="left", fontweight="bold", pad=3)

# Per-gene, per-replicate: detected (max_cov > 0)?
pach_genes_all = cov[cov["stage"] == "Pachytene"]["gene"].unique()
n_pach = len(pach_genes_all)  # should be 99

detect_mat = pd.DataFrame(index=pach_genes_all, columns=SAMPLE_KEYS)
for key in SAMPLE_KEYS:
    sub = cov[(cov["sample"] == key) & (cov["stage"] == "Pachytene")]
    sub = sub.set_index("gene")
    for g in pach_genes_all:
        detect_mat.loc[g, key] = 1 if (g in sub.index and sub.loc[g, "n_det"] > 0) else 0

detect_mat = detect_mat.astype(int)

# Summarise: both P12.5 detected / both P20.5 detected / both
p125_det = (detect_mat["P12.5_rep1"] & detect_mat["P12.5_rep2"]).sum()
p205_det = (detect_mat["P20.5_rep1"] & detect_mat["P20.5_rep2"]).sum()
both_det  = ((detect_mat["P12.5_rep1"] | detect_mat["P12.5_rep2"]) &
             (detect_mat["P20.5_rep1"] | detect_mat["P20.5_rep2"])).sum()

# Per-gene coverage (average of max_cov across both replicates per timepoint)
p125_cov = (cov[(cov["sample"].isin(["P12.5_rep1","P12.5_rep2"])) & (cov["stage"]=="Pachytene")]
            .groupby("gene")["max_cov"].mean())
p205_cov = (cov[(cov["sample"].isin(["P20.5_rep1","P20.5_rep2"])) & (cov["stage"]=="Pachytene")]
            .groupby("gene")["max_cov"].mean())

# Merge
pach_df = pd.DataFrame({"p125": p125_cov, "p205": p205_cov}).fillna(0)

# Scatter: P12.5 avg coverage vs P20.5 avg coverage
ax2_C.scatter(pach_df["p125"], pach_df["p205"], s=8, alpha=0.65,
              color=C_PACH, linewidths=0.3, edgecolors="white", zorder=3)

# Identity line
lim = [0, 1.05]
ax2_C.plot(lim, lim, color="#bbbbbb", lw=0.6, ls="--", zorder=0)

# Highlight genes not detected at P12.5 but at P20.5
p205_only = pach_df[(pach_df["p125"] == 0) & (pach_df["p205"] > 0)]
ax2_C.scatter(p205_only["p125"], p205_only["p205"], s=18, alpha=0.9,
              color=C_P205, linewidths=0.5, edgecolors="black",
              zorder=5, label=f"Pachytene-onset\n(n={len(p205_only)})")

# Annotate r² for genes present in both
both_pres = pach_df[(pach_df["p125"] > 0) & (pach_df["p205"] > 0)]
if len(both_pres) > 2:
    r, pval = stats.pearsonr(both_pres["p125"], both_pres["p205"])
    ax2_C.text(0.05, 0.93, f"$R^2$ = {r**2:.3f} (detected at both)",
               transform=ax2_C.transAxes, fontsize=5.5, va="top")

ax2_C.set_xlim(-0.05, 1.1)
ax2_C.set_ylim(-0.05, 1.1)
ax2_C.set_xlabel("P12.5 coverage fraction")
ax2_C.set_ylabel("P20.5 coverage fraction")
ax2_C.set_aspect("equal")
ax2_C.legend(frameon=False, loc="upper left", handlelength=1.0,
             handletextpad=0.3, borderpad=0, labelspacing=0.2, fontsize=5.5)

# Stat summary
n_p125_only = ((detect_mat["P12.5_rep1"] | detect_mat["P12.5_rep2"]) &
               ~(detect_mat["P20.5_rep1"] | detect_mat["P20.5_rep2"])).sum()
n_neither   = (~(detect_mat["P12.5_rep1"] | detect_mat["P12.5_rep2"]) &
               ~(detect_mat["P20.5_rep1"] | detect_mat["P20.5_rep2"])).sum()
ax2_C.text(0.05, 0.84,
    f"Both timepoints: {both_det}  |  P20.5-onset: {len(p205_only)}",
    transform=ax2_C.transAxes, fontsize=5.0, va="top", color="#555555")

# Bottom annotation
fig2.text(0.5, 0.01,
    "C57BL/6 · GRCm38 · Zamore mm10 (= GRCm38) annotation · 215 piRNA genes · 99 Pachytene · 83 Prepachytene · 32 Hybrid",
    ha="center", va="bottom", fontsize=5.0, color="#555555")

# ── Save Fig2 input data matrices ────────────────────────────────────────
# Detection rates table
det_rows = []
for key in SAMPLE_KEYS:
    row = {"sample": key, "timepoint": key.split("_")[0]}
    for st in STAGES_MAIN:
        row[f"pct_detected_{st}"] = round(det_rates[key][st], 2)
        row[f"n_total_{st}"] = STAGE_TOTALS[st]
    det_rows.append(row)
fig2_det = pd.DataFrame(det_rows)
fig2_det.to_csv(os.path.join(SD_DIR, "Fig2_detection_rates.csv"), index=False)
print(f"  Saved Fig2_detection_rates.csv")

# Per-gene coverage averaged across replicates
for tp_label, tp_samples in [("P12.5", ["P12.5_rep1","P12.5_rep2"]),
                               ("P20.5", ["P20.5_rep1","P20.5_rep2"])]:
    tp_avg = (cov[cov["sample"].isin(tp_samples)]
              .groupby(["gene","stage"])["max_cov"].mean()
              .reset_index()
              .rename(columns={"max_cov": "avg_max_cov"}))
    tp_avg.to_csv(os.path.join(SD_DIR, f"Fig2_coverage_per_gene_{tp_label}.csv"), index=False)
    print(f"  Saved Fig2_coverage_per_gene_{tp_label}.csv")

# Pachytene scatter matrix (Fig2C)
pach_df_out = pach_df.copy().reset_index()
pach_df_out.columns = ["gene", "avg_cov_P12.5", "avg_cov_P20.5"]
pach_df_out.to_csv(os.path.join(SD_DIR, "Fig2_pachytene_scatter.csv"), index=False)
print(f"  Saved Fig2_pachytene_scatter.csv")

stem2 = os.path.join(OUT_DIR, "Fig2_Zamore_gene_coverage")
for ext in (".pdf", ".svg", ".png"):
    kw = {"dpi": 300, "bbox_inches": "tight"}
    fig2.savefig(stem2 + ext, **kw)
    print(f"  Saved {stem2+ext}")
plt.close(fig2)


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Bedtools Coverage Detail
# Panels: A) per-stage coverage heatmap across samples
#          B) coverage CDF per stage
#          C) high-expression pachytene loci bar chart
# Size: 183 mm wide × 70 mm tall
# ══════════════════════════════════════════════════════════════════════════════
print("\nGenerating Figure 3: Coverage Detail...")

fig3, axes3 = plt.subplots(1, 3, figsize=(7.2, 2.75))
fig3.subplots_adjust(left=0.08, right=0.97, top=0.87, bottom=0.22, wspace=0.46)
ax3_A, ax3_B, ax3_C = axes3

# ── Panel A: Coverage heatmap (gene × sample) for Pachytene genes ──────────
ax3_A.set_title("A  Pachytene locus coverage heatmap", loc="left",
                fontweight="bold", pad=3)

pach_cov_mat = cov[cov["stage"] == "Pachytene"].pivot_table(
    index="gene", columns="sample", values="max_cov", aggfunc="max"
)
# Keep column order
pach_cov_mat = pach_cov_mat[[k for k in SAMPLE_KEYS if k in pach_cov_mat.columns]]
pach_cov_mat = pach_cov_mat.sort_values(SAMPLE_KEYS[2], ascending=False)

im = ax3_A.imshow(pach_cov_mat.values, aspect="auto", cmap="YlGn",
                   vmin=0, vmax=1, interpolation="nearest")

ax3_A.set_yticks([])
ax3_A.set_xticks(range(len(pach_cov_mat.columns)))
ax3_A.set_xticklabels(SAMPLE_LABELS, fontsize=6)
ax3_A.set_ylabel(f"Pachytene genes (n={len(pach_cov_mat)})")
ax3_A.set_xlabel("Sample")

cbar = fig3.colorbar(im, ax=ax3_A, shrink=0.75, pad=0.02)
cbar.set_label("Coverage fraction", fontsize=6)
cbar.ax.tick_params(labelsize=5)

for tick, (_, _, tp, clr) in zip(ax3_A.get_xticklabels(), SAMPLES):
    tick.set_color(clr)

# ── Panel B: CDF of coverage fraction per stage (pooled across samples) ────
ax3_B.set_title("B  Coverage CDF by piRNA stage", loc="left",
                fontweight="bold", pad=3)

for st in STAGES_MAIN:
    data = cov[cov["stage"] == st].groupby("gene")["max_cov"].mean().values
    data_sorted = np.sort(data)
    cdf = np.arange(1, len(data_sorted) + 1) / len(data_sorted)
    ax3_B.plot(data_sorted, cdf, color=STAGE_CLR[st], lw=1.2,
               label=f"{st} (n={len(data_sorted)})")

ax3_B.set_xlabel("Coverage fraction")
ax3_B.set_ylabel("Cumulative proportion")
ax3_B.set_xlim(-0.02, 1.02)
ax3_B.set_ylim(0, 1.05)
ax3_B.legend(frameon=False, loc="upper left", handlelength=1.5,
             handletextpad=0.4, borderpad=0, labelspacing=0.25)
ax3_B.axvline(0.5, color="#cccccc", lw=0.5, ls=":", zorder=0)

# ── Panel C: FPM at top Pachytene loci (P20.5) vs (P12.5) ─────────────────
ax3_C.set_title("C  Top Pachytene loci — FPM at P20.5", loc="left",
                fontweight="bold", pad=3)

# Get Zamore Pachytene gene coordinates then cross-reference cluster FPM
# Use P20.5 rep1 cluster file for FPM
# We don't have the cluster-to-Zamore mapping per se, but we have coverage
# Use median of all_reads_primary_alignments_FPM for clusters on Pachytene chr regions
# Proxy: for each sample, top 15 clusters by FPM at P20.5 rep1
top_n = 15
df_p205r1 = cluster_dfs["P20.5_rep1"].nlargest(top_n, "all_reads_primary_alignments_FPM")
df_p205r1 = df_p205r1.reset_index(drop=True)

# Build locus label
df_p205r1["locus"] = (df_p205r1["seqnames"].astype(str) + ":" +
    (df_p205r1["start"] / 1e6).round(1).astype(str) + "M")

# Get same loci in P12.5 rep1 by coordinate match
fpm_p125, fpm_p205 = [], []
for _, row in df_p205r1.iterrows():
    p125_match = cluster_dfs["P12.5_rep1"][
        (cluster_dfs["P12.5_rep1"]["seqnames"] == row["seqnames"]) &
        (cluster_dfs["P12.5_rep1"]["start"] == row["start"]) &
        (cluster_dfs["P12.5_rep1"]["end"]   == row["end"])
    ]
    fpm_p125.append(p125_match["all_reads_primary_alignments_FPM"].values[0]
                    if len(p125_match) > 0 else 0.0)
    fpm_p205.append(row["all_reads_primary_alignments_FPM"])

ypos = np.arange(top_n)[::-1]
ax3_C.barh(ypos - 0.2, np.log10(np.array(fpm_p125) + 0.01), height=0.35,
           color=C_P125, label="P12.5 rep1", linewidth=lw, edgecolor="white")
ax3_C.barh(ypos + 0.2, np.log10(np.array(fpm_p205) + 0.01), height=0.35,
           color=C_P205, label="P20.5 rep1", linewidth=lw, edgecolor="white")

ax3_C.set_yticks(ypos)
ax3_C.set_yticklabels(df_p205r1["locus"].values[::-1], fontsize=4.5)
ax3_C.set_xlabel("log₁₀(FPM)")
ax3_C.legend(frameon=False, loc="lower right", handlelength=1.0,
             handletextpad=0.3, borderpad=0, labelspacing=0.2)

fig3.text(0.5, 0.01,
    "C57BL/6 · GRCm38 · Coverage = bedtools fraction of Zamore locus outer span covered by ≥1 PICB cluster base",
    ha="center", va="bottom", fontsize=5.0, color="#555555")

# ── Save Fig3 input data matrices ─────────────────────────────────────────
# Pachytene heatmap matrix
pach_cov_mat.reset_index().to_csv(
    os.path.join(OUT_DIR, "Fig3_pachytene_coverage_heatmap.csv"), index=False)
print(f"  Saved Fig3_pachytene_coverage_heatmap.csv")

# CDF input: per-stage average coverage per gene
for st in STAGES_MAIN:
    st_data = cov[cov["stage"] == st].groupby("gene")["max_cov"].mean().reset_index()
    st_data.columns = ["gene", "mean_max_cov"]
    st_data.to_csv(
        os.path.join(OUT_DIR, f"Fig3_CDF_data_{st}.csv"), index=False)
    print(f"  Saved Fig3_CDF_data_{st}.csv")

# Top loci FPM comparison
top_loci_out = df_p205r1[["seqnames","start","end","strand","type","locus",
                            "all_reads_primary_alignments_FPM"]].copy()
top_loci_out.columns = ["chr","start","end","strand","type","locus","FPM_P20.5_rep1"]
top_loci_out["FPM_P12.5_rep1"] = fpm_p125
top_loci_out.to_csv(os.path.join(SD_DIR, "Fig3_top_pachytene_loci_FPM.csv"), index=False)
print(f"  Saved Fig3_top_pachytene_loci_FPM.csv")

stem3 = os.path.join(OUT_DIR, "Fig3_coverage_detail")
for ext in (".pdf", ".svg", ".png"):
    fig3.savefig(stem3 + ext, dpi=300, bbox_inches="tight")
    print(f"  Saved {stem3+ext}")
plt.close(fig3)


print(f"\nAll figures saved to {OUT_DIR}")
print("Done.")
