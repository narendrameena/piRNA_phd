"""Shared 'HOW TO READ' corner key for the 16-strain piRNA circos figures.
Each figure calls draw_readkey() with `rows` describing its concentric tracks (outer->inner). Every row is
(label, [colour,...], description); the helper draws the colour swatch(es) + text + an outer->inner arrow, so each
circos is self-explanatory. Colour always = a categorical attribute (family/strand/identity/# strains); HEIGHT
always = a log amount — the key states both for every track."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch


def draw_readkey(fig, rows, note, title="HOW TO READ — one strain block", loc=(0.008, 0.735, 0.27, 0.23)):
    gax = fig.add_axes(list(loc)); gax.axis("off"); gax.set_xlim(0, 1); gax.set_ylim(0, 1)
    gax.add_patch(plt.Rectangle((0, 0), 1, 1, facecolor="#f7f7f7", edgecolor="#bbb", lw=1.2, zorder=0))
    gax.text(0.5, 0.955, title, ha="center", va="top", fontsize=9.5, fontweight="bold")
    n = len(rows); top = 0.85; bot = 0.30; dy = (top - bot) / max(n, 1)
    for i, (label, colors, desc) in enumerate(rows):
        yy = top - (i + 0.5) * dy
        x = 0.055; seg = 0.20 / len(colors)
        for c in colors:
            gax.add_patch(plt.Rectangle((x, yy - 0.030), seg, 0.060, color=c, ec="none")); x += seg
        gax.text(0.285, yy, f"{label} — {desc}", fontsize=6.6, va="center")
    gax.annotate("", xy=(0.045, top - 0.02), xytext=(0.045, bot + 0.02), arrowprops=dict(arrowstyle="<->", color="#777", lw=1.0))
    gax.text(0.017, (top + bot) / 2, "outer→inner", rotation=90, ha="center", va="center", fontsize=6.0, color="#777")
    gax.text(0.5, bot - 0.02, note, ha="center", va="top", fontsize=6.3, style="italic")


def zoom_inset_6nj(fig, ax, rings, hl_theta, r_lo, r_hi, loc=(0.625, 0.805, 0.18, 0.15)):
    """Small magnified ZOOM INSET of the C57BL/6NJ block, placed next to the START of its tracks (chr1), for a
    many-track block (e.g. integrated16) where in-place names would cram together. Highlights 6NJ's chr1-start
    wedge (red box, r_lo..r_hi over hl_theta) and connects it to a small inset that redraws the track stack
    enlarged + labelled. rings = [(label, [colours]), ...] outer->inner."""
    tt = np.linspace(hl_theta[0], hl_theta[1], 12)
    for r in (r_lo, r_hi):
        ax.plot(tt, [r] * 12, color="#C0392B", lw=1.0, zorder=13, solid_capstyle="butt")
    for th in hl_theta:
        ax.plot([th, th], [r_lo, r_hi], color="#C0392B", lw=1.0, zorder=13)
    gax = fig.add_axes(list(loc)); gax.axis("off"); gax.set_xlim(0, 1); gax.set_ylim(0, 1)
    gax.add_patch(Rectangle((0, 0), 1, 1, facecolor="white", edgecolor="#C0392B", lw=1.3, zorder=0))
    n = len(rings); h = 0.90 / n; y0 = 0.95
    for i, (label, colours) in enumerate(rings):
        y = y0 - (i + 0.5) * h; x = 0.05; seg = 0.24 / len(colours)
        for c in colours:
            gax.add_patch(Rectangle((x, y - h * 0.34), seg, h * 0.64, color=c, ec="none")); x += seg
        gax.text(0.34, y, label, fontsize=5.9, va="center", ha="left", color="#222")
    fig.add_artist(ConnectionPatch(xyA=(0.0, 0.5), coordsA=gax.transAxes,
                                   xyB=((hl_theta[0] + hl_theta[1]) / 2, r_hi), coordsB=ax.transData,
                                   color="#C0392B", lw=0.7, ls="--", zorder=12))


def zoom_6nj(ax, rings, theta_c, fs=4.8, spread=False):
    """Name 6NJ's tracks in place: write each track's short NAME in the empty label gap (at angle theta_c, on the
    open chr1 side of the 6NJ strain-name spoke) at that track's own radius — aligned with 6NJ's real tracks. No
    colour, no header (the strain name already sits at the spoke). `rings` = [(label, r_mid), ...] outer->inner.
    If spread=True (blocks with many thin sub-tracks, e.g. integrated16): fan the names over a wider radial range
    and draw a thin arrow from each name to its real track, so the names don't overlap — readability only."""
    if not spread:
        for label, r in rings:
            ax.text(theta_c, r, label, fontsize=fs, ha="right", va="center", color="#222", zorder=12)
        return
    rs = [r for _, r in rings]; rc = sum(rs) / len(rs); span = (max(rs) - min(rs)) * 3.0
    labr = [rc + span / 2 - i * span / (len(rings) - 1) for i in range(len(rings))]
    tl = theta_c - 0.085
    for (label, r), lr in zip(rings, labr):
        ax.annotate(label, xy=(theta_c, r), xytext=(tl, lr), fontsize=fs, ha="right", va="center",
                    color="#222", zorder=12, arrowprops=dict(arrowstyle="->", color="#888", lw=0.6,
                    shrinkA=2, shrinkB=1, connectionstyle="arc3,rad=0.25"))
