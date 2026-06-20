"""
Canonical strain order for EVERY figure in this project.

Source: thesis Figure 4.4 (RNA-seq PCA) strain legend; strains ordered by median
P20.5 PC1 position from the RNA-seq PCA (thesis p.136). Same order used in the
piRNA chapter (Fig 5.7/5.9). This is NOT the Snakemake/project order and NOT
alphabetical. Import this everywhere instead of redefining order per script.

  from strain_order import STRAIN_ORDER, WILD, TIMEPOINT_ORDER, strain_rank
"""

STRAIN_ORDER = [
    'C57BL_6', 'C57BL_6NJ', 'BALB_cJ', 'A_J', 'FVB_NJ', 'C3H_HeJ', 'LP_J',
    '129S1_SvImJ', 'DBA_2J', 'AKR_J', 'CBA_J', 'NZO_HlLtJ', 'NOD_ShiLtJ',
    'WSB_EiJ', 'CAST_EiJ', 'PWK_PhJ', 'SPRET_EiJ',
]

# C57BL_6 (pos 1) = reference C57BL/6J (T2T), distinct from C57BL_6NJ.

WILD = {'CAST_EiJ', 'PWK_PhJ', 'SPRET_EiJ', 'WSB_EiJ'}

TIMEPOINT_ORDER = ['E16.5', 'P12.5', 'P20.5']  # developmental, never alphabetical

_RANK = {s: i for i, s in enumerate(STRAIN_ORDER)}

def strain_rank(s):
    """Sort key: position in STRAIN_ORDER; unknown strains sort to the end."""
    return _RANK.get(s, len(STRAIN_ORDER))

# ── Canonical classical / wild-derived colour scheme (Okabe-Ito, colourblind-safe) ──
# Use for per-strain bar plots so colour = subspecies (classical vs wild-derived).
CLASSICAL_COLOR = "#0072B2"   # blue
WILD_COLOR      = "#D55E00"   # vermillion/orange
WILD_LABEL_COLOR = "#C0392B"  # red x-axis labels for wild strains

def strain_bar_colors(strains):
    """List of bar fill colours for the given strains (wild=orange, classical=blue)."""
    return [WILD_COLOR if s in WILD else CLASSICAL_COLOR for s in strains]

def color_wild_labels(ax, strains, classical="#222222"):
    """Colour the x-tick labels red for wild strains (call after set_xticklabels)."""
    for lab, s in zip(ax.get_xticklabels(), strains):
        lab.set_color(WILD_LABEL_COLOR if s in WILD else classical)

def classical_wild_legend_handles():
    """Two Patch handles labelled classical / wild-derived for a legend."""
    from matplotlib.patches import Patch
    return [Patch(facecolor=CLASSICAL_COLOR, label="classical"),
            Patch(facecolor=WILD_COLOR, label="wild-derived")]

def add_classical_wild_companion(fig, main_ax, strains, totals, height_frac=0.12, gap=0.045,
                                 ylabel="total\n(log)", log=True, fontsize=6.0):
    """Add a slim companion bar panel BELOW `main_ax`, x-aligned to it, showing `totals`
    per strain coloured classical/wild (subspecies colour scheme alongside a klass5/composition panel).
    Non-destructive: the main panel keeps its own colours. Returns the companion axes."""
    import numpy as np
    pos = main_ax.get_position()
    cax = fig.add_axes([pos.x0, pos.y0 - gap - height_frac * pos.height,
                        pos.width, height_frac * pos.height])
    x = np.arange(len(strains))
    cax.bar(x, totals, width=0.8, color=strain_bar_colors(strains), edgecolor="white", linewidth=0.3, zorder=3)
    if log: cax.set_yscale("log")
    cax.set_xlim(*main_ax.get_xlim())
    cax.set_xticks([])
    cax.set_ylabel(ylabel, fontsize=fontsize)
    cax.tick_params(labelsize=fontsize - 0.5)
    cax.spines[["top", "right"]].set_visible(False)
    return cax
