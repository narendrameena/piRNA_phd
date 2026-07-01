#!/usr/bin/env python3
"""
picb_comparison_figure.py
Compare chr-by-chr vs whole-BAM PICB results.

Usage:
    python3 picb_comparison_figure.py <chrbychr.xlsx> <wholebam.xlsx> <output_dir>

Outputs:
    PICB_chrbychr_vs_wholebam.pdf / .svg / .png
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from scipy import stats

# ── Paths ──────────────────────────────────────────────────────────────────
if len(sys.argv) < 4:
    sys.exit("Usage: picb_comparison_figure.py <chrbychr.xlsx> <wholebam.xlsx> <output_dir>")

CHRBYCHR_XLSX = sys.argv[1]
WHOLEBAM_XLSX = sys.argv[2]
OUT_DIR       = sys.argv[3]
os.makedirs(OUT_DIR, exist_ok=True)

# ── Nature Genetics style ──────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Arial", "DejaVu Sans", "Helvetica"],
    "font.size":         7,
    "axes.titlesize":    7,
    "axes.labelsize":    7,
    "xtick.labelsize":   6,
    "ytick.labelsize":   6,
    "legend.fontsize":   6,
    "axes.linewidth":    0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size":  2.5,
    "ytick.major.size":  2.5,
    "lines.linewidth":   0.8,
    "pdf.fonttype":      42,
    "svg.fonttype":      "none",
})

# Colorblind-safe palette (Wong 2011)
CLR_CHR  = "#0072B2"   # blue  — chr-by-chr
CLR_WBAM = "#E69F00"   # amber — whole-BAM
CLR_DIAG = "#CC79A7"   # pink  — identity diagonal

# ── Load data ──────────────────────────────────────────────────────────────
def load_sheet(path, sheet):
    try:
        df = pd.read_excel(path, sheet_name=sheet, engine="openpyxl")
        if df.shape[0] == 1 and "Message" in df.columns:
            return pd.DataFrame()
        return df
    except Exception:
        return pd.DataFrame()

print(f"Loading {CHRBYCHR_XLSX}")
chr_seeds    = load_sheet(CHRBYCHR_XLSX, "seeds")
chr_cores    = load_sheet(CHRBYCHR_XLSX, "cores")
chr_clusters = load_sheet(CHRBYCHR_XLSX, "clusters")

print(f"Loading {WHOLEBAM_XLSX}")
wb_seeds    = load_sheet(WHOLEBAM_XLSX, "seeds")
wb_cores    = load_sheet(WHOLEBAM_XLSX, "cores")
wb_clusters = load_sheet(WHOLEBAM_XLSX, "clusters")

counts = {
    "seeds":    (len(chr_seeds),    len(wb_seeds)),
    "cores":    (len(chr_cores),    len(wb_cores)),
    "clusters": (len(chr_clusters), len(wb_clusters)),
}
print("Counts chr-by-chr / whole-BAM:")
for k, (a, b) in counts.items():
    print(f"  {k}: {a} / {b}  {'✓ identical' if a == b else '✗ DIFFER'}")

# ── Match clusters by genomic coordinates ──────────────────────────────────
def merge_on_coords(df_a, df_b, suffix_a="chrbychr", suffix_b="wholebam"):
    if df_a.empty or df_b.empty:
        return pd.DataFrame()
    key = ["seqnames", "start", "end", "strand"]
    available = [c for c in key if c in df_a.columns and c in df_b.columns]
    merged = pd.merge(df_a, df_b, on=available, suffixes=(f"_{suffix_a}", f"_{suffix_b}"),
                      how="outer", indicator=True)
    return merged

merged_clusters = merge_on_coords(chr_clusters, wb_clusters)
only_chr = int((merged_clusters["_merge"] == "left_only").sum())  if not merged_clusters.empty else 0
only_wb  = int((merged_clusters["_merge"] == "right_only").sum()) if not merged_clusters.empty else 0
both     = int((merged_clusters["_merge"] == "both").sum())       if not merged_clusters.empty else 0
print(f"\nCluster overlap: both={both}, chr-only={only_chr}, whole-BAM-only={only_wb}")

# ── Per-chromosome counts ──────────────────────────────────────────────────
CHR_ORDER = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"]

def chr_counts(df):
    if df.empty:
        return pd.Series(dtype=int)
    col = "seqnames" if "seqnames" in df.columns else df.columns[0]
    vc = df[col].value_counts()
    return vc.reindex(CHR_ORDER, fill_value=0)

chr_per_chr    = chr_counts(chr_clusters)
chr_wholebam   = chr_counts(wb_clusters)
all_chrs = [c for c in CHR_ORDER if chr_per_chr.get(c, 0) > 0 or chr_wholebam.get(c, 0) > 0]

# ── Figure layout ──────────────────────────────────────────────────────────
# 3 panels: A (counts bar) | B (FPM scatter) | C (per-chr bar)
# Figure width 183 mm (2-col Nature Genetics), height 70 mm

fig = plt.figure(figsize=(7.2, 2.75))   # inches; 7.2" ≈ 183 mm
gs  = fig.add_gridspec(1, 3, wspace=0.42, left=0.07, right=0.97, top=0.88, bottom=0.22)
ax_a = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[1])
ax_c = fig.add_subplot(gs[2])

# ── Panel A: Grouped counts bar ────────────────────────────────────────────
categories = ["Seeds", "Cores", "Clusters"]
vals_chr  = [counts["seeds"][0],    counts["cores"][0],    counts["clusters"][0]]
vals_wb   = [counts["seeds"][1],    counts["cores"][1],    counts["clusters"][1]]

x = np.arange(len(categories))
w = 0.35
bars_chr = ax_a.bar(x - w/2, vals_chr, width=w, color=CLR_CHR,  label="Chr-by-chr",
                    linewidth=0.4, edgecolor="white")
bars_wb  = ax_a.bar(x + w/2, vals_wb,  width=w, color=CLR_WBAM, label="Whole-BAM",
                    linewidth=0.4, edgecolor="white")

# Mark identical bars
for i, (a, b) in enumerate(zip(vals_chr, vals_wb)):
    top_y = max(a, b) * 1.04
    marker = "=" if a == b else "≠"
    ax_a.text(x[i], top_y, marker, ha="center", va="bottom", fontsize=6.5,
              color="#555555")

ax_a.set_xticks(x)
ax_a.set_xticklabels(categories)
ax_a.set_ylabel("Count")
ax_a.set_title("A  PICB component counts", loc="left", fontweight="bold", pad=3)
ax_a.spines["top"].set_visible(False)
ax_a.spines["right"].set_visible(False)
ax_a.yaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=5))
ax_a.legend(frameon=False, loc="upper right", handlelength=1.2, handletextpad=0.4)

# ── Panel B: FPM scatter (chr-by-chr vs whole-BAM) ─────────────────────────
fpm_col = "all_reads_primary_alignments_FPM"
ax_b.set_title("B  Cluster FPM: chr-by-chr vs whole-BAM", loc="left",
               fontweight="bold", pad=3)

if not merged_clusters.empty and f"{fpm_col}_chrbychr" in merged_clusters.columns:
    df_both = merged_clusters[merged_clusters["_merge"] == "both"].copy()
    x_fpm = df_both[f"{fpm_col}_chrbychr"].values
    y_fpm = df_both[f"{fpm_col}_wholebam"].values

    # log10 scale for dynamic range; add tiny pseudocount
    eps = 1e-3
    lx = np.log10(x_fpm + eps)
    ly = np.log10(y_fpm + eps)

    ax_b.scatter(lx, ly, s=4, alpha=0.6, color=CLR_CHR, linewidths=0, rasterized=True)

    # Identity line
    mn = min(lx.min(), ly.min()) - 0.1
    mx = max(lx.max(), ly.max()) + 0.1
    ax_b.plot([mn, mx], [mn, mx], color=CLR_DIAG, lw=0.8, ls="--", zorder=0)

    # R²
    r, _ = stats.pearsonr(lx, ly)
    ax_b.text(0.05, 0.93, f"$R^2$ = {r**2:.6f}", transform=ax_b.transAxes,
              fontsize=6, va="top")
    ax_b.text(0.05, 0.85, f"n = {len(df_both)}", transform=ax_b.transAxes,
              fontsize=6, va="top")

    ax_b.set_xlabel("Chr-by-chr log₁₀(FPM)")
    ax_b.set_ylabel("Whole-BAM log₁₀(FPM)")
else:
    ax_b.text(0.5, 0.5, "No matched clusters\n(run not complete)", ha="center",
              va="center", transform=ax_b.transAxes, fontsize=6, color="grey")

ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)

# ── Panel C: Per-chromosome cluster counts ─────────────────────────────────
ax_c.set_title("C  Clusters per chromosome", loc="left", fontweight="bold", pad=3)

if all_chrs:
    ypos = np.arange(len(all_chrs))
    bw   = 0.35
    ax_c.barh(ypos + bw/2, [chr_per_chr.get(c, 0) for c in all_chrs],
              height=bw, color=CLR_CHR,  linewidth=0.3, edgecolor="white", label="Chr-by-chr")
    ax_c.barh(ypos - bw/2, [chr_wholebam.get(c, 0) for c in all_chrs],
              height=bw, color=CLR_WBAM, linewidth=0.3, edgecolor="white", label="Whole-BAM")
    ax_c.set_yticks(ypos)
    ax_c.set_yticklabels(all_chrs, fontsize=5.5)
    ax_c.set_xlabel("Cluster count")
    ax_c.invert_yaxis()
else:
    ax_c.text(0.5, 0.5, "No cluster data", ha="center", va="center",
              transform=ax_c.transAxes, fontsize=6, color="grey")

ax_c.spines["top"].set_visible(False)
ax_c.spines["right"].set_visible(False)
ax_c.xaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=5))

# ── Summary annotation ─────────────────────────────────────────────────────
overlap_pct = f"{100*both/(both+only_chr+only_wb):.1f}%" if (both+only_chr+only_wb) > 0 else "N/A"
fig.text(
    0.5, 0.01,
    f"E15.5 rep1 · GRCm38 · 25–36 nt · cluster overlap: {overlap_pct}"
    f" ({both}/{both+only_chr+only_wb})",
    ha="center", va="bottom", fontsize=5.5, color="#555555"
)

# ── Source data (plotted values → theme dir) ───────────────────────────────
import os as _os
_SD = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/04_picb_method_qc/data/source_data"
_os.makedirs(_SD, exist_ok=True)
pd.DataFrame({"component": categories, "chrbychr": vals_chr, "wholebam": vals_wb}).to_csv(
    f"{_SD}/SourceData_PICB_chrbychr_vs_wholebam_A_counts.csv", index=False)          # Panel A grouped bars
pd.DataFrame({"chr": all_chrs,
              "chrbychr": [int(chr_per_chr.get(c, 0)) for c in all_chrs],
              "wholebam": [int(chr_wholebam.get(c, 0)) for c in all_chrs]}).to_csv(
    f"{_SD}/SourceData_PICB_chrbychr_vs_wholebam_C_per_chromosome.csv", index=False)  # Panel C per-chr bars
try:  # Panel B scatter only exists when clusters matched
    pd.DataFrame({"FPM_chrbychr": x_fpm, "FPM_wholebam": y_fpm}).to_csv(
        f"{_SD}/SourceData_PICB_chrbychr_vs_wholebam_B_FPM_scatter.csv", index=False)
except NameError:
    pass

# ── Save ───────────────────────────────────────────────────────────────────
stem = os.path.join(OUT_DIR, "PICB_chrbychr_vs_wholebam")
fig.savefig(stem + ".pdf", dpi=300, bbox_inches="tight")
fig.savefig(stem + ".svg", dpi=300, bbox_inches="tight")
fig.savefig(stem + ".png", dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"\nFigure saved:")
print(f"  {stem}.pdf")
print(f"  {stem}.svg")
print(f"  {stem}.png")
print("Done.")
