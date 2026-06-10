#!/usr/bin/env python3
"""
PICB Chr-by-chr vs Whole-BAM Comparison Figure
Compares three PICB processing modes:
  Mode 1 – Correct chr-by-chr  : genome-wide LIBRARY.SIZE passed explicitly
  Mode 2 – Naive chr-by-chr    : LIBRARY.SIZE from single chromosome only
  Mode 3 – Whole-BAM           : genome-wide LIBRARY.SIZE from PICBload
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import matplotlib.gridspec as gridspec

plt.rcParams.update({
    "font.family":       "Liberation Sans",
    "font.size":         7,
    "axes.titlesize":    8,
    "axes.labelsize":    7,
    "xtick.labelsize":   6,
    "ytick.labelsize":   6,
    "lines.linewidth":   1.0,
    "axes.linewidth":    0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size":  2.5,
    "ytick.major.size":  2.5,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "pdf.fonttype":      42,
})

C_CORRECT = "#2166ac"
C_NAIVE   = "#d6604d"
C_GREY    = "#888888"

LS_REP1  = 11_277_390
LS_REP2C = 19_039_720
LS_E165  =  4_961_424

WINDOW_WIDTH = 350
MIN_SEQ_CAP  = 7

def min_aln(lib_size):
    return 2 * (WINDOW_WIDTH / 1000) * (lib_size / 1_000_000)

def min_seq(lib_size):
    return min(min_aln(lib_size), MIN_SEQ_CAP)

CHR_READS_REP1 = {
    "1": 39_510_776, "2": 29_160_252, "3": 29_026_554, "4": 32_045_862,
    "5": 30_235_245, "6": 27_330_548, "7": 34_462_849, "8": 23_099_133,
    "9": 21_383_005, "10": 24_834_389, "11": 19_900_308, "12": 23_714_528,
    "13": 26_398_432, "14": 25_017_294, "15": 16_686_994, "16": 18_479_268,
    "17": 17_792_719, "18": 15_789_088, "19": 9_056_015,
    "X": 41_471_676, "Y": 28_202_642,
}
CHR_READS_E165 = {
    "1": 19_742_015, "2": 14_633_890, "3": 14_970_881, "4": 15_151_553,
    "5": 14_019_698, "6": 13_231_510, "7": 15_606_910, "8": 11_400_786,
    "9": 10_125_868, "10": 12_062_777, "11": 9_457_394, "12": 10_955_634,
    "13": 12_849_855, "14": 12_304_063, "15": 8_068_536, "16": 8_774_282,
    "17": 8_644_161, "18": 7_525_221, "19": 4_190_077,
    "X": 20_844_473, "Y": 10_371_856,
}

CHRS = [str(i) for i in range(1, 20)] + ["X", "Y"]

def naive_ls(reads_dict, genome_ls, c):
    total = sum(reads_dict.values())
    return (reads_dict.get(c, 0) / total) * genome_ls

naive_ls_rep1 = [naive_ls(CHR_READS_REP1, LS_REP1, c) for c in CHRS]
naive_ls_e165 = [naive_ls(CHR_READS_E165, LS_E165, c) for c in CHRS]

thresh_naive_rep1 = [min_aln(ls) for ls in naive_ls_rep1]
thresh_naive_e165 = [min_aln(ls) for ls in naive_ls_e165]
thresh_correct_rep1 = min_aln(LS_REP1)
thresh_correct_e165 = min_aln(LS_E165)

# FPM inflation = genome_LS / chr_LS = total/chr_reads
fpm_infl = [LS_REP1 / ls if ls > 0 else 0 for ls in naive_ls_rep1]

# ── Layout ─────────────────────────────────────────────────────────────────────
# 3 rows, 2 cols: A+B (top), C+D (mid), E+G (bottom)
fig = plt.figure(figsize=(10, 10))
fig.patch.set_facecolor("white")

gs = gridspec.GridSpec(3, 2, figure=fig,
                       top=0.93, bottom=0.06,
                       hspace=0.72, wspace=0.38,
                       left=0.09, right=0.97)

ax_a = fig.add_subplot(gs[0, 0])   # A – schematic
ax_b = fig.add_subplot(gs[0, 1])   # B – library size bars
ax_c = fig.add_subplot(gs[1, 0])   # C – threshold per chr (rep1)
ax_d = fig.add_subplot(gs[1, 1])   # D – FPM inflation
ax_e = fig.add_subplot(gs[2, 0])   # E – threshold comparison rep1 vs E16.5
ax_g = fig.add_subplot(gs[2, 1])   # G – cluster count impact

# ══════════════════════════════════════════════════════════════════════════════
# Panel A – Schematic: three modes
# ══════════════════════════════════════════════════════════════════════════════
ax_a.axis("off")
ax_a.set_title("A   Three PICB processing modes", loc="left",
               fontweight="bold", fontsize=8)

# Use a simple table-style layout: 3 rows (Mode 1, Mode 2, formula)
ax_a.set_xlim(0, 1); ax_a.set_ylim(0, 1)

def hbox(ax, y0, h, color, title, lines, title_color=None):
    rect = FancyBboxPatch((0.02, y0), 0.96, h,
                          boxstyle="round,pad=0.01",
                          facecolor=color, edgecolor="none", zorder=1)
    ax.add_patch(rect)
    tc = title_color if title_color else color
    ax.text(0.05, y0 + h - 0.015, title,
            va="top", ha="left", fontsize=7, fontweight="bold", color=tc, zorder=2)
    for i, line in enumerate(lines):
        ax.text(0.05, y0 + h - 0.015 - (i + 1) * 0.062,
                line, va="top", ha="left", fontsize=6, color="#222222", zorder=2)

hbox(ax_a, 0.65, 0.32, "#dce9f5",
     "Mode 1 & 3  (CORRECT)",
     ["LIBRARY.SIZE = genome-wide  →  uniform threshold across all chromosomes",
      "rep1: LS = 11.28M   MIN.UNIQUE.ALN = 7.89 / 350 nt window",
      "Used by picb_script_chunked.R (Mode 1) and PICBload whole-BAM (Mode 3)"],
     title_color=C_CORRECT)

hbox(ax_a, 0.33, 0.29, "#fce4df",
     "Mode 2  (NAIVE – INCORRECT)",
     ["LIBRARY.SIZE = reads on that chromosome only  →  threshold varies per chr",
      "rep1 chr2: LS ≈ 618K   MIN.UNIQUE.ALN ≈ 0.43  (~18x lower than correct)",
      "Consequence: thousands of spurious clusters pass the threshold"],
     title_color=C_NAIVE)

hbox(ax_a, 0.01, 0.28, "#f5f5f5",
     "Key PICB formula  (PICBbuild.R:62)",
     ["MIN.UNIQUE.ALN = 2 x (350/1000) x (LIBRARY.SIZE / 1e6)",
      "MIN.UNIQUE.SEQ = min(MIN.UNIQUE.ALN,  round(350/50))  =  min(X, 7)",
      "LIBRARY.SIZE must be genome-wide for correct FPKM thresholds"],
     title_color=C_GREY)

# ══════════════════════════════════════════════════════════════════════════════
# Panel B – Genome-wide vs naive LIBRARY.SIZE
# ══════════════════════════════════════════════════════════════════════════════
ax_b.set_title("B   Genome-wide vs naive LIBRARY.SIZE (chr2 example)",
               loc="left", fontweight="bold", fontsize=8)

samples = ["rep1\n(328 cl)", "rep2c\n(8 cl)", "E16.5\n(13,629 cl)"]
ls_gw  = [LS_REP1/1e6, LS_REP2C/1e6, LS_E165/1e6]
ls_chr2 = [naive_ls(CHR_READS_REP1, LS_REP1, "2")/1e6,
           naive_ls({k: CHR_READS_REP1[k]*2 for k in CHR_READS_REP1},
                    LS_REP2C, "2")/1e6,   # use rep1 proportions for rep2c
           naive_ls(CHR_READS_E165, LS_E165, "2")/1e6]

# Recompute chr2 naive for rep2c using actual rep2c data
CHR_READS_REP2C = {
    "1": 78_875_389, "2": 57_314_501, "3": 58_127_389, "4": 63_205_410,
    "5": 60_336_962, "6": 54_791_719, "7": 68_329_595, "8": 45_928_271,
    "9": 42_653_614, "10": 49_338_313, "11": 39_595_353, "12": 48_110_041,
    "13": 53_048_078, "14": 49_412_425, "15": 33_188_880, "16": 36_819_815,
    "17": 35_305_981, "18": 31_583_615, "19": 18_002_744,
    "X": 84_657_764, "Y": 53_418_092,
}
ls_chr2[1] = naive_ls(CHR_READS_REP2C, LS_REP2C, "2") / 1e6

xp = np.arange(3)
bw = 0.32
b1 = ax_b.bar(xp - bw/2, ls_gw,   bw, color=C_CORRECT, label="Genome-wide (correct)", zorder=3)
b2 = ax_b.bar(xp + bw/2, ls_chr2, bw, color=C_NAIVE,   label="Chr2 naive (incorrect)", zorder=3, alpha=0.85)

for bar in b1:
    ax_b.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
              f"{bar.get_height():.1f}M", ha="center", va="bottom", fontsize=5.5)
for bar in b2:
    ax_b.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
              f"{bar.get_height():.2f}M", ha="center", va="bottom",
              fontsize=5.5, color=C_NAIVE)

ax_b.set_xticks(xp); ax_b.set_xticklabels(samples, fontsize=6.5)
ax_b.set_ylabel("LIBRARY.SIZE (x10^6 reads)", fontsize=7)
ax_b.set_ylim(0, 25)
ax_b.legend(fontsize=6, frameon=False, loc="upper left", handlelength=1.2)
ax_b.grid(axis="y", lw=0.3, color="#dddddd", zorder=0)

# ══════════════════════════════════════════════════════════════════════════════
# Panel C – Per-chromosome threshold (rep1)
# ══════════════════════════════════════════════════════════════════════════════
ax_c.set_title("C   Per-chromosome seeding threshold  (rep1, MIN.UNIQUE.ALN)",
               loc="left", fontweight="bold", fontsize=8)

xc = np.arange(len(CHRS))
ax_c.axhline(thresh_correct_rep1, color=C_CORRECT, lw=1.8, ls="-", zorder=4)
ax_c.plot(xc, thresh_naive_rep1, color=C_NAIVE, lw=1.0, ls="--",
          marker="o", ms=2.5, zorder=3)
ax_c.axhline(7, color=C_GREY, lw=0.7, ls=":", zorder=2)
ax_c.fill_between(xc, thresh_naive_rep1, thresh_correct_rep1,
                  color=C_NAIVE, alpha=0.08)

# In-plot text labels instead of legend box
ax_c.text(len(CHRS) - 0.3, thresh_correct_rep1 + 0.15,
          f"Correct: {thresh_correct_rep1:.2f}", ha="right", va="bottom",
          fontsize=6, color=C_CORRECT, fontweight="bold")
ax_c.text(len(CHRS) - 0.3, thresh_naive_rep1[-1] - 0.15,
          "Naive (per chr)", ha="right", va="top",
          fontsize=6, color=C_NAIVE, fontweight="bold")
ax_c.text(len(CHRS) - 0.3, 7.15, "MAX = 7", ha="right", va="bottom",
          fontsize=5.5, color=C_GREY)

# Arrow to chr2
chr2_i = CHRS.index("2")
ax_c.annotate("chr2:\n0.43", xy=(chr2_i, thresh_naive_rep1[chr2_i]),
              xytext=(chr2_i + 3, thresh_naive_rep1[chr2_i] + 2.5),
              fontsize=5.5, color=C_NAIVE, fontweight="bold",
              arrowprops=dict(arrowstyle="->", color=C_NAIVE, lw=0.7),
              ha="left", va="bottom")

ax_c.set_xticks(xc); ax_c.set_xticklabels(CHRS, fontsize=5.5, rotation=45, ha="right")
ax_c.set_ylabel("MIN.UNIQUE.ALN / 350 nt window", fontsize=7)
ax_c.set_ylim(0, 10)
ax_c.grid(axis="y", lw=0.3, color="#dddddd", zorder=0)

# ══════════════════════════════════════════════════════════════════════════════
# Panel D – FPM inflation (cap at 30×, annotate chr19 if clipped)
# ══════════════════════════════════════════════════════════════════════════════
ax_d.set_title("D   FPM/FPKM inflation under naive approach  (rep1)",
               loc="left", fontweight="bold", fontsize=8)

CAP = 30
clipped = [(i, v) for i, v in enumerate(fpm_infl) if v > CAP]
bar_vals = [min(v, CAP) for v in fpm_infl]
bar_colors = [C_NAIVE if v > 1.05 else C_CORRECT for v in fpm_infl]

ax_d.bar(np.arange(len(CHRS)), bar_vals, color=bar_colors, width=0.7, zorder=3,
         edgecolor="none")

# Mark clipped bars with arrow and value
for i, v in clipped:
    ax_d.text(i, CAP + 0.3, f"{v:.0f}x", ha="center", va="bottom",
              fontsize=5, color=C_NAIVE, fontweight="bold")
    ax_d.annotate("", xy=(i, CAP), xytext=(i, CAP - 3),
                  arrowprops=dict(arrowstyle="->", color=C_NAIVE, lw=0.7))

ax_d.axhline(1.0, color=C_CORRECT, lw=1.0, ls="-", zorder=4)
mean_infl = np.mean(fpm_infl)
ax_d.axhline(min(mean_infl, CAP), color=C_GREY, lw=0.8, ls="--", zorder=4)
ax_d.text(len(CHRS) - 0.3, min(mean_infl, CAP) + 0.5,
          f"mean = {mean_infl:.0f}x", ha="right", va="bottom",
          fontsize=5.5, color=C_GREY)

# Highlight chr2 with outline
chr2_i = CHRS.index("2")
ax_d.bar(chr2_i, fpm_infl[chr2_i], color=C_NAIVE, width=0.7, zorder=5,
         edgecolor="#b22222", linewidth=0.8)
ax_d.text(chr2_i, fpm_infl[chr2_i] + 0.5,
          f"{fpm_infl[chr2_i]:.0f}x", ha="center", va="bottom",
          fontsize=5.5, color="#b22222", fontweight="bold")

ax_d.set_xticks(np.arange(len(CHRS)))
ax_d.set_xticklabels(CHRS, fontsize=5.5, rotation=45, ha="right")
ax_d.set_ylabel("FPM inflation factor", fontsize=7)
ax_d.set_ylim(0, CAP + 5)
ax_d.text(0.5, 1.5, "No inflation\n(correct = 1x)", ha="center", va="bottom",
          fontsize=5.5, color=C_CORRECT)
ax_d.grid(axis="y", lw=0.3, color="#dddddd", zorder=0)

# ══════════════════════════════════════════════════════════════════════════════
# Panel E – Threshold comparison: rep1 vs E16.5, correct vs naive
# ══════════════════════════════════════════════════════════════════════════════
ax_e.set_title("E   Threshold per chromosome: rep1 vs E16.5  (correct vs naive)",
               loc="left", fontweight="bold", fontsize=8)

xc = np.arange(len(CHRS))
# rep1
ax_e.axhline(thresh_correct_rep1, color=C_CORRECT, lw=1.8, ls="-", zorder=4)
ax_e.plot(xc, thresh_naive_rep1,  color=C_CORRECT, lw=1.0, ls="--",
          marker="s", ms=2, zorder=3)
# E16.5
ax_e.axhline(thresh_correct_e165, color="#E69F00", lw=1.8, ls="-", zorder=4)
ax_e.plot(xc, thresh_naive_e165,  color="#E69F00", lw=1.0, ls="--",
          marker="o", ms=2, zorder=3)
ax_e.axhline(7, color=C_GREY, lw=0.7, ls=":", zorder=2)

# Right-side labels
ax_e.text(len(CHRS) - 0.2, thresh_correct_rep1 + 0.1,
          f"rep1 correct={thresh_correct_rep1:.2f}", ha="right", va="bottom",
          fontsize=5.5, color=C_CORRECT, fontweight="bold")
ax_e.text(len(CHRS) - 0.2, thresh_correct_e165 + 0.1,
          f"E16.5 correct={thresh_correct_e165:.2f}", ha="right", va="bottom",
          fontsize=5.5, color="#E69F00", fontweight="bold")
ax_e.text(len(CHRS) - 0.2, thresh_naive_rep1[-1] - 0.1,
          "rep1 naive", ha="right", va="top",
          fontsize=5.5, color=C_CORRECT, style="italic")
ax_e.text(len(CHRS) - 0.2, thresh_naive_e165[-1] - 0.1,
          "E16.5 naive", ha="right", va="top",
          fontsize=5.5, color="#E69F00", style="italic")
ax_e.text(len(CHRS) - 0.2, 7.1, "MAX = 7", ha="right", va="bottom",
          fontsize=5.5, color=C_GREY)

ax_e.set_xticks(xc); ax_e.set_xticklabels(CHRS, fontsize=5.5, rotation=45, ha="right")
ax_e.set_ylabel("MIN.UNIQUE.ALN / 350 nt window", fontsize=7)
ax_e.set_ylim(0, 10)
ax_e.grid(axis="y", lw=0.3, color="#dddddd", zorder=0)

# Legend: solid = correct, dashed = naive
handles = [
    mpatches.Patch(color=C_CORRECT,  label="rep1 (11.28M reads)"),
    mpatches.Patch(color="#E69F00",  label="E16.5 (4.96M reads)"),
    plt.Line2D([0], [0], color="black", lw=1.5, ls="-",  label="Correct (genome-wide LS)"),
    plt.Line2D([0], [0], color="black", lw=1.0, ls="--", label="Naive (per-chr LS)"),
]
ax_e.legend(handles=handles, fontsize=5.5, frameon=False, ncol=2,
            loc="upper center", bbox_to_anchor=(0.5, -0.18), handlelength=1.5)

# ══════════════════════════════════════════════════════════════════════════════
# Panel G – Cluster count impact (rep1 chr2 example)
# ══════════════════════════════════════════════════════════════════════════════
ax_g.set_title("G   Threshold impact on cluster detection  (rep1, chr2 example)",
               loc="left", fontweight="bold", fontsize=8)

thresh_naive_chr2 = min_aln(naive_ls(CHR_READS_REP1, LS_REP1, "2"))
xg = np.array([0, 1])
thresh_g = [thresh_correct_rep1, thresh_naive_chr2]

bars_g = ax_g.bar(xg, thresh_g, color=[C_CORRECT, C_NAIVE], width=0.5,
                  zorder=3, edgecolor="none")

# Reference lines
ax_g.axhline(6.4, color="#333333", lw=1.0, ls="-.", zorder=4)
ax_g.axhline(7.0, color=C_GREY,    lw=0.8, ls=":",  zorder=4)

# Compact annotations inside plot
ax_g.text(0, thresh_correct_rep1 / 2,
          f"threshold\n= {thresh_correct_rep1:.2f}",
          ha="center", va="center", fontsize=6, color="white", fontweight="bold")
if thresh_naive_chr2 > 1.5:
    ax_g.text(1, thresh_naive_chr2 / 2,
              f"threshold\n= {thresh_naive_chr2:.2f}",
              ha="center", va="center", fontsize=6, color="white", fontweight="bold")
else:
    ax_g.text(1.27, thresh_naive_chr2 / 2,
              f"threshold={thresh_naive_chr2:.2f}",
              ha="left", va="center", fontsize=6, color=C_NAIVE, fontweight="bold")

# Outcome labels above bars
ax_g.text(0, thresh_correct_rep1 + 0.25,
          "6.4 < 7  -->  FAILS\n(cluster not called)",
          ha="center", va="bottom", fontsize=6, color=C_CORRECT, fontweight="bold")
ax_g.text(1, thresh_naive_chr2 + 0.25,
          f"6.4 > {thresh_naive_chr2:.2f}  -->  PASSES\n(spurious cluster!)",
          ha="center", va="bottom", fontsize=6, color=C_NAIVE, fontweight="bold")

# Horizontal reference labels at right
ax_g.text(1.35, 6.4 + 0.1, "Observed\n(6.4 pos/window)", ha="left",
          va="bottom", fontsize=5.5, color="#333333")
ax_g.text(1.35, 7.1, "MAX.UNIQUE.SEQ = 7", ha="left", va="bottom",
          fontsize=5.5, color=C_GREY)

ax_g.set_xticks(xg)
ax_g.set_xticklabels(["Correct\n(genome-wide LS)", "Naive\n(chr2 only LS)"], fontsize=7)
ax_g.set_ylabel("MIN.UNIQUE.ALN threshold", fontsize=7)
ax_g.set_ylim(0, 13)
ax_g.set_xlim(-0.5, 1.85)
ax_g.grid(axis="y", lw=0.3, color="#dddddd", zorder=0)

# ── Global title ──────────────────────────────────────────────────────────────
fig.suptitle(
    "PICB piRNA cluster pipeline: chr-by-chr vs whole-BAM processing modes\n"
    "E15.5 mPGC C57BL/6  |  GRCm38  |  PICB v0.99.x",
    fontsize=9, fontweight="bold", y=0.975)

# ── Save ──────────────────────────────────────────────────────────────────────
OUT_DIR = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis"
for ext, dpi in [("pdf", 600), ("svg", 150), ("png", 300)]:
    fig.savefig(os.path.join(OUT_DIR, f"PICB_chr_vs_wholebam_comparison.{ext}"),
                bbox_inches="tight", dpi=dpi)
plt.close(fig)
print("Saved PICB_chr_vs_wholebam_comparison.{pdf,svg,png}")
