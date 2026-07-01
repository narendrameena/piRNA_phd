#!/usr/bin/env python3
"""
Figure 7: Pangenome SVs within Zamore piRNA loci and their effect on PICB expression.

Four panels:
  A  Locus schematic — 13-qA5-208: TE insertions at 5′ end → locus absent (chain break)
  B  Locus schematic — 13-qA5-464: large 3′ deletion → locus present, PICB maintained
  C  Genome-wide: PICB detection rate by SV group (case study + all Zamore loci)
  D  50 kb window analysis: PICB class distribution for loci near vs far from SVs
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrow, Rectangle, FancyBboxPatch
from matplotlib.gridspec import GridSpec
import matplotlib.patheffects as pe
from scipy.stats import fisher_exact

BASE = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
OUT  = f"{BASE}/analysis/claude_biomni_analysis/C57BL_6NJ_pangenome"

# ── Colours ───────────────────────────────────────────────────────────────────
C_P12   = "#0072B2"   # blue  — P12.5
C_P20   = "#E69F00"   # amber — P20.5
C_BOTH  = "#009E73"   # green — both
C_NONE  = "#cccccc"   # grey  — no PICB
C_INS   = "#D55E00"   # vermillion — insertions
C_DEL   = "#CC79A7"   # pink  — deletions
C_LOCUS = "#444444"   # dark grey — locus body

plt.rcParams.update({
    "font.family":       "Liberation Sans",
    "font.size":         8,
    "axes.linewidth":    0.6,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        300,
    "savefig.dpi":       300,
    "pdf.fonttype":      42,
})

df = pd.read_csv(f"{OUT}/zamore_SV_PICB_analysis.csv")

fig = plt.figure(figsize=(13, 10))
gs  = GridSpec(2, 2, figure=fig, hspace=0.55, wspace=0.38)

# ═══════════════════════════════════════════════════════════════════════════════
# Panel A — 13-qA5-208 schematic: TE INS → chain break → no PICB
# ═══════════════════════════════════════════════════════════════════════════════
ax_a = fig.add_subplot(gs[0, 0])

locus_start = 50_422_060
locus_end   = 50_471_396
locus_kb    = (locus_end - locus_start) / 1000   # ~49.3 kb

# INS events (GRCm39 positions, sizes from pangenome VCF ALT lengths; TE class from RepeatMasker)
ins_events = [(50_426_422, 18065, "LINE/L1"),
              (50_430_497, 2719,  "LINE/L1"),
              (50_432_524, 1107,  "LTR/ERVK")]
chain_break = 50_426_929   # first chain starts here → gap from locus_start to here

def pos2kb(p):
    return (p - locus_start) / 1000

# Draw locus body
locus_y = 0.55
ax_a.add_patch(Rectangle((0, locus_y - 0.04), locus_kb, 0.08,
                          color=C_LOCUS, zorder=2))
ax_a.text(locus_kb / 2, locus_y + 0.12, "13-qA5-208  (Pachytene, 49 kb)",
          ha="center", va="bottom", fontsize=8, fontweight="bold")

# Chain gap (no coverage) — hatched from locus start to chain_break
gap_end_kb = pos2kb(chain_break)
ax_a.add_patch(Rectangle((0, locus_y - 0.04), gap_end_kb, 0.08,
                          color="white", zorder=3, hatch="////",
                          edgecolor="#aaaaaa", linewidth=0))
# Chain gap annotation below locus to avoid colliding with INS arrows
ax_a.annotate("", xy=(gap_end_kb, locus_y - 0.08), xytext=(0, locus_y - 0.08),
              arrowprops=dict(arrowstyle="<->", color="#aaaaaa", lw=0.8))
ax_a.text(gap_end_kb / 2, locus_y - 0.10, "chain gap\n(no C57BL_6NJ synteny)",
          ha="center", va="top", fontsize=6.5, color="#888888", style="italic")

# INS arrows above locus — stagger heights so labels don't overlap
ins_tops = [locus_y + 0.24, locus_y + 0.36, locus_y + 0.24]
te_colors = {"LINE/L1": "#D55E00", "LTR/ERVK": "#CC79A7", "LTR/ERVL-MaLR": "#E69F00"}
for i, (pos, size, te_class) in enumerate(ins_events):
    kb = pos2kb(pos)
    top_y = ins_tops[i]
    ax_a.annotate("", xy=(kb, locus_y + 0.04), xytext=(kb, top_y - 0.03),
                  arrowprops=dict(arrowstyle="-|>", color=C_INS,
                                  lw=1.5, mutation_scale=10))
    if size >= 10000:
        size_str = f"+{round(size/1000)}kb"
    elif size >= 1000:
        size_str = f"+{size/1000:.1f}kb"
    else:
        size_str = f"+{size}bp"
    ax_a.text(kb, top_y, size_str,
              ha="center", va="bottom", fontsize=6.5, color=C_INS, fontweight="bold")
    ax_a.text(kb, top_y + 0.035, te_class,
              ha="center", va="bottom", fontsize=5.5,
              color=te_colors.get(te_class, "#888888"), style="italic")

# PICB signal area (below locus) — empty (no signal)
picb_y = 0.2
ax_a.add_patch(Rectangle((0, picb_y - 0.06), locus_kb, 0.06,
                          facecolor="#f0f0f0", edgecolor="#cccccc", lw=0.5, zorder=1))
ax_a.text(-1.5, picb_y - 0.03, "PICB\nP12.5", va="center", ha="right",
          fontsize=6.5, color=C_P12)
ax_a.add_patch(Rectangle((0, picb_y - 0.13), locus_kb, 0.06,
                          facecolor="#f0f0f0", edgecolor="#cccccc", lw=0.5, zorder=1))
ax_a.text(-1.5, picb_y - 0.1, "PICB\nP20.5", va="center", ha="right",
          fontsize=6.5, color=C_P20)
ax_a.text(locus_kb / 2, picb_y - 0.18,
          "Locus ABSENT from C57BL_6NJ — no PICB signal",
          ha="center", va="top", fontsize=7, color=C_INS, fontweight="bold")

ax_a.set_xlim(-3, locus_kb + 3)
ax_a.set_ylim(0.02, 1.12)
ax_a.set_xlabel("Position within locus (kb)", fontsize=8)
ax_a.set_xticks(np.arange(0, locus_kb + 1, 10))
ax_a.set_yticks([])
ax_a.spines["left"].set_visible(False)
ax_a.set_title("A   TE insertions at 5′ end → chain break → piRNA loss",
               fontsize=8.5, fontweight="bold", loc="left")

ins_patch = mpatches.Patch(color=C_INS, label="TE insertion (INS)")
gap_patch = mpatches.Patch(facecolor="white", edgecolor="#aaaaaa",
                           hatch="////", label="No syntenic chain")
ax_a.legend(handles=[ins_patch, gap_patch], fontsize=7, frameon=False,
            loc="upper right", bbox_to_anchor=(1.0, 0.98))

# ═══════════════════════════════════════════════════════════════════════════════
# Panel B — 13-qA5-464 schematic: large 3′ DEL → chain intact → PICB maintained
# ═══════════════════════════════════════════════════════════════════════════════
ax_b = fig.add_subplot(gs[0, 1])

locus_start2 = 50_655_589
locus_end2   = 50_749_523
locus_kb2    = (locus_end2 - locus_start2) / 1000  # ~93.9 kb

del_start    = 50_723_521
del_end_clip = 50_749_523  # DEL extends past locus; clip to locus end
del_kb_s     = (del_start   - locus_start2) / 1000
del_kb_e     = (del_end_clip - locus_start2) / 1000

# Draw locus body (intact 5′ part)
ax_b.add_patch(Rectangle((0, locus_y - 0.04), locus_kb2, 0.08,
                          color=C_LOCUS, zorder=2))
# Overlay DEL region (pink)
ax_b.add_patch(Rectangle((del_kb_s, locus_y - 0.04), del_kb_e - del_kb_s, 0.08,
                          facecolor=C_DEL, alpha=0.7, zorder=3))
ax_b.annotate("", xy=(del_kb_e, locus_y - 0.08), xytext=(del_kb_s, locus_y - 0.08),
              arrowprops=dict(arrowstyle="<->", color=C_DEL, lw=1.2))
ax_b.text((del_kb_s + del_kb_e) / 2, locus_y - 0.12,
          "−71.7 kb DEL  (LINE/L1, LTR/ERVK)\n(3′ portion removed)", ha="center", va="top",
          fontsize=6.5, color=C_DEL, fontweight="bold")
ax_b.text(locus_kb2 / 2, locus_y + 0.12, "13-qA5-464  (Pachytene, 94 kb)",
          ha="center", va="bottom", fontsize=8, fontweight="bold")

# PICB signal (below locus) — both P12.5 and P20.5 present
picb_coverage_start = 0
picb_coverage_end   = del_kb_s  # PICB only over the conserved 5′ part
for i, (yoff, col, tp) in enumerate([(-0.06, C_P12, "P12.5"), (-0.13, C_P20, "P20.5")]):
    ax_b.add_patch(Rectangle((0, picb_y + yoff), del_kb_s, 0.06,
                              facecolor=col, alpha=0.6, edgecolor=col, lw=0.5, zorder=2))
    ax_b.add_patch(Rectangle((del_kb_s, picb_y + yoff), locus_kb2 - del_kb_s, 0.06,
                              facecolor="#f0f0f0", edgecolor="#cccccc", lw=0.5, zorder=1))
    ax_b.text(-2, picb_y + yoff + 0.03, f"PICB\n{tp}", va="center", ha="right",
              fontsize=6.5, color=col)
ax_b.text(del_kb_s / 2, picb_y - 0.17,
          "piRNA production maintained (P12.5 & P20.5)",
          ha="center", va="top", fontsize=7, color=C_BOTH, fontweight="bold")

ax_b.set_xlim(-4, locus_kb2 + 3)
ax_b.set_ylim(0.02, 1.0)
ax_b.set_xlabel("Position within locus (kb)", fontsize=8)
ax_b.set_xticks(np.arange(0, locus_kb2 + 1, 20))
ax_b.set_yticks([])
ax_b.spines["left"].set_visible(False)
ax_b.set_title("B   3′ deletion → promoter intact → piRNA maintained",
               fontsize=8.5, fontweight="bold", loc="left")

del_patch = mpatches.Patch(color=C_DEL, alpha=0.7, label="DEL in C57BL_6NJ")
p12_patch = mpatches.Patch(color=C_P12, alpha=0.6, label="PICB P12.5")
p20_patch = mpatches.Patch(color=C_P20, alpha=0.6, label="PICB P20.5")
ax_b.legend(handles=[del_patch, p12_patch, p20_patch], fontsize=7, frameon=False,
            loc="upper left", bbox_to_anchor=(0.0, 0.98))

# ═══════════════════════════════════════════════════════════════════════════════
# Panel C — PICB detection rate by SV group
# ═══════════════════════════════════════════════════════════════════════════════
ax_c = fig.add_subplot(gs[1, 0])

# Groups:
#   1. Liftover failed (n=4): all have 0 PICB
#   2. Direct INS (13-qA5-208, in not_lifted, n=1): 0 PICB
#   3. Direct DEL (13-qA5-464): both timepoints
#   4. No direct SV, lifted (n~208): varies

not_lifted = df[df['picb_class'] == 'not_lifted']        # 4 loci
ins_direct = df[(df['SV_direct']) & (df['INS_direct'] > 0)]  # 1: 13-qA5-208
del_direct = df[(df['SV_direct']) & (df['DEL_direct'] > 0)]  # 1: 13-qA5-464
no_sv_lift = df[(~df['SV_direct']) & (df['picb_class'] != 'not_lifted')]  # 208

groups = [
    ("Liftover\nfailed\n(no large SV)", not_lifted, "#888888"),
    ("Direct INS\n13-qA5-208\n(chain break)", ins_direct, C_INS),
    ("Direct DEL\n13-qA5-464\n(3′ truncation)", del_direct, C_DEL),
    ("No direct SV\n(lifted, n=209)", no_sv_lift, C_LOCUS),
]

picb_order = ['P12.5_only', 'P20.5_only', 'both', 'none']
picb_cols  = [C_P12, C_P20, C_BOTH, C_NONE]
picb_labels = ['P12.5 only', 'P20.5 only', 'Both', 'No PICB']

x = np.arange(len(groups))
bottoms = np.zeros(len(groups))
for cls, col, label in zip(picb_order, picb_cols, picb_labels):
    vals = []
    for _, grp, _ in groups:
        n_cls = (grp['picb_class'] == cls).sum()
        pct   = 100.0 * n_cls / max(len(grp), 1)
        vals.append(pct)
    ax_c.bar(x, vals, bottom=bottoms, color=col, label=label, width=0.55,
             edgecolor='white', lw=0.4)
    bottoms += np.array(vals)

ax_c.set_xticks(x)
ax_c.set_xticklabels([g[0] for g in groups], fontsize=7)
ax_c.set_ylabel("% of Zamore piRNA loci", fontsize=8)
ax_c.set_ylim(0, 112)
ax_c.set_title("C   PICB expression status by SV category", fontsize=8.5,
               fontweight="bold", loc="left")
for i, (_, grp, _) in enumerate(groups):
    ax_c.text(i, bottoms[i] + 1.5, f"n={len(grp)}", ha="center", va="bottom",
              fontsize=7)
ax_c.legend(title="PICB class", fontsize=7, title_fontsize=7,
            loc="upper center", bbox_to_anchor=(0.5, -0.18),
            ncol=2, frameon=False,
            handles=[mpatches.Patch(color=c, label=l)
                     for c, l in zip(picb_cols, picb_labels)])

# ═══════════════════════════════════════════════════════════════════════════════
# Panel D — 50 kb window: PICB detection rate vs no nearby SV
# ═══════════════════════════════════════════════════════════════════════════════
ax_d = fig.add_subplot(gs[1, 1])

# Only lifted loci
lifted = df[df['picb_class'] != 'not_lifted'].copy()
lifted['expressed'] = lifted['picb_class'] != 'none'

sv_grps  = lifted[lifted['SV_50kb']]
nosv_grps = lifted[~lifted['SV_50kb']]

# Per-stage detection rate
stages = ['Prepachytene', 'Pachytene', 'Hybrid']
x2 = np.arange(len(stages))
w  = 0.35

for offset, grp, label, col, n_yoff in [
        (-w/2, nosv_grps, 'No SV ≤50 kb',    "#888888", 1),
        (+w/2, sv_grps,   'SV within 50 kb',  C_INS,     8)]:
    rates = []
    ns    = []
    for st in stages:
        sub = grp[grp['stage'] == st]
        n   = len(sub)
        k   = sub['expressed'].sum()
        rates.append(100.0 * k / n if n else 0)
        ns.append(n)
    bars = ax_d.bar(x2 + offset, rates, width=w, color=col,
                    label=label, alpha=0.85, edgecolor='white', lw=0.3)
    for b, n in zip(bars, ns):
        ax_d.text(b.get_x() + b.get_width()/2, b.get_height() + n_yoff,
                  f"n={n}", ha="center", va="bottom", fontsize=6.5)

# Fisher exact p-value for Pachytene (most data)
for i, st in enumerate(stages):
    a = sv_grps[sv_grps['stage']==st]
    b = nosv_grps[nosv_grps['stage']==st]
    if len(a) > 0 and len(b) > 0:
        tbl = [[a['expressed'].sum(), len(a)-a['expressed'].sum()],
               [b['expressed'].sum(), len(b)-b['expressed'].sum()]]
        try:
            _, p = fisher_exact(tbl)
            if p < 0.1:
                ymax = max(100*a['expressed'].mean(), 100*b['expressed'].mean())
                ax_d.text(i, ymax + 8, f"p={p:.2f}", ha="center",
                          fontsize=6, color="#555555")
        except Exception:
            pass

ax_d.set_xticks(x2)
ax_d.set_xticklabels(stages, fontsize=8)
ax_d.set_ylabel("% loci with any PICB coverage", fontsize=8)
ax_d.set_ylim(0, 115)
ax_d.set_title("D   PICB detection: SVs within 50 kb vs none\n(lifted loci only)",
               fontsize=8.5, fontweight="bold", loc="left")
ax_d.legend(fontsize=7, frameon=False, loc="lower right")

fig.suptitle("Pangenome SVs within Zamore piRNA loci and effect on PICB expression\n"
             "C57BL_6NJ REL-2205 vs GRCm39 reference", fontsize=9.5, y=1.01)

for ext in ("pdf", "svg", "png"):
    fig.savefig(os.path.join(OUT, f"Fig7_SV_expression.{ext}"),
                dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved Fig7_SV_expression.{pdf,svg,png}")
