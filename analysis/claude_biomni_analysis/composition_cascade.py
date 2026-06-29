"""Shared figure helper: a compact 3-step *funnel* explaining why the 100%-scale class-composition bar totals are
LARGER than the per-'Combined'/pooled-PCA panel n. Each level is a horizontal bar STACKED BY CLASS (segment colours
match the composition bar / row labels), so the reader sees per-class shrinkage:

  (1) strain-specific candidates pooled over all 3 timepoints   (= the composition bar; the 3 per-timepoint panels sum to it)
  (2) distinct sequences after de-duplicating timepoints        (a sequence flagged at >1 tp is counted once)
  (3) sequences expressed at all 3 timepoints                   (the intersection a single pooled PCA can use → each bottom
                                                                  segment = that class's 'Combined' panel n)

Used by Fig_pca_classes16 (top) and Fig_pca_unique16 (bottom). Numbers are passed in (computed in the figure script
from the same source tables), never hard-coded here."""

def draw_cascade(fig, rect, levels, colors, title="Why bar totals > each 'Combined'/pooled panel n", note=None):
    """`levels` = [(label, [per-class counts]), ...] top→bottom; `colors` = per-class colour list (len == each counts
    list). Bars are left-aligned and stacked by class on a shared x-scale, so widths shrink down the funnel. `note`
    (optional) is a small italic line under the last level."""
    ax = fig.add_axes(rect); ax.axis("off")
    top = sum(levels[0][1]); ny = len(levels)
    for i, (lab, arr) in enumerate(levels):
        y = ny - 1 - i; ti = sum(arr); L = 0
        for n, c in zip(arr, colors):
            ax.barh(y, n, left=L, height=0.56, color=c, edgecolor="white", linewidth=0.5, zorder=3); L += n
        ax.text(-top*0.013, y, lab, ha="right", va="center", fontsize=6.2, color="#333")
        pct = "" if i == 0 else f"   ({100*ti/top:.0f}% of pooled)"
        ax.text(ti + top*0.013, y, f"{ti:,}{pct}", ha="left", va="center", fontsize=6.9,
                fontweight="bold", color="#222")
        if i < ny-1:   # light funnel edges (left aligned at 0; right edge slopes inward)
            n2 = sum(levels[i+1][1])
            ax.plot([ti, n2], [y-0.28, (y-1)+0.28], color="#bbb", lw=0.8, zorder=2)
            ax.plot([0, 0],   [y-0.28, (y-1)+0.28], color="#ddd", lw=0.8, zorder=2)
    if note:
        ax.text(0, -0.88, note, ha="left", va="center", fontsize=6.0, style="italic", color="#555")
    ax.set_xlim(-top*0.66, top*1.34); ax.set_ylim(-1.18 if note else -0.62, ny-0.30)
    if title:
        ax.set_title(title, fontsize=8.0, fontweight="bold", loc="left")
    return ax
