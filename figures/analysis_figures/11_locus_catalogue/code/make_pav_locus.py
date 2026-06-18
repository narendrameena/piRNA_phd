#!/usr/bin/env python3
"""Locus plot for a PICB piRNA cluster (top-expressing strain detail). DATA SOURCES ONLY: pangenome cluster
projection (picb_pangenome_clusters.tsv = PICB-COMBINED clusters with own+GRCm39 coords), sRNA BAMs (PRIMARY
reads inside the cluster), RepeatMasker TE annotation, pangenome genome-PAV. NO per-figure liftover; NO spillover.
  (A) PANGENOME cross-strain × timepoint PICB FPM (all 16 strains) + genome presence (● in genome / ○ absent);
  (B) the cluster in the TOP strain: PRIMARY-read coverage INSIDE the cluster (plus ↑ green / minus ↓ purple) +
      TE track + dominant TE; right margin = per-timepoint ON (strand split) / OFF (no cluster);
  (C) base resolution in the top strain.
TERMINOLOGY: genomic +/- = architecture ≠ sense/antisense (relative to TE; antisense-to-TE = silencing).
Usage: make_pav_locus.py <g39chrom> <g39start> <g39end> <gene_label> <unused> <outbase>."""
import sys, os, numpy as np
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
from strain_order import STRAIN_ORDER, WILD
from collections import Counter
import pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch, Patch
import pav_clusters as pc
from pav_clusters import TPS, TPLAB, TPCOL, clusters_at, present_strains, fpm_by_tp, cluster_extent, fetch_primary, te_at, dom_te_family, genome_pav
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]
NT = {"A": "#33a02c", "C": "#1f78b4", "G": "#ff7f00", "T": "#e31a1c", "N": "#999"}; comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
TECOL = {"LINE/L1": "#E69F00", "LTR/ERVK": "#6a3d9a", "LTR/ERVL-MaLR": "#b15928", "SINE/Alu": "#33a02c", "SINE/B2": "#1f78b4", "LTR/ERVL": "#a6cee3", "LTR/ERV1": "#cab2d6", "Simple_repeat": "#bbbbbb"}
G39C, G39S, G39E, GENE, OUT = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[6]
nb = 200
sub = clusters_at(G39C, G39S, G39E)
FPM = sub.groupby(["strain", "tp"])["all_primary_FPM"].max().unstack(fill_value=0.0).reindex(index=ORDER, columns=TPS).fillna(0.0)
present = present_strains(sub, ORDER)
TOP = max(present, key=lambda X: FPM.loc[X].max()) if present else ORDER[0]
PAV = {X: ("present" if FPM.loc[X].max() > 0 else genome_pav(G39C, G39S, G39E, X)) for X in ORDER}   # expression PROVES presence: a PICB cluster with reads => locus IS there (halLiftover g39->strain can false-negative)
nabs = sum(1 for X in ORDER if PAV[X] == "absent"); nsil = sum(1 for X in ORDER if PAV[X] == "present" and FPM.loc[X].max() == 0)
print(f"{GENE}: PICB present strains={present}; top={TOP} (FPM={FPM.loc[TOP].max():.1f}); PAV absent={nabs} silent={nsil}")
SF = {}                                                          # (strain, tp) -> +strand fraction, for the Panel A strand height-split (solid + / pale −)
for _X in present:
    _ex = cluster_extent(sub, _X); _fX = fpm_by_tp(sub, _X)
    if not _ex: continue
    for _tp in TPS:
        if _tp in _fX:
            _dd = fetch_primary(_X, _ex[0], _ex[1], _ex[2], _tp, nb)
            if _dd and _dd["ntot"] > 0: SF[(_X, _tp)] = _dd["nplus"] / _dd["ntot"]
fps = fpm_by_tp(sub, TOP); ext = cluster_extent(sub, TOP); CH, ps, pe = ext[0], ext[1], ext[2]; BAMC = f"{TOP}#1#chr{CH}"; N = max(1, pe - ps)
per = {}; plus = np.zeros(nb); minus = np.zeros(nb); reads = []; nat = nte = n1u = npl = nmi = 0; tes = None
for tp in TPS:
    if tp not in fps: per[tp] = None; continue
    d = fetch_primary(TOP, CH, ps, pe, tp, nb); per[tp] = d
    plus += d["plus"]; minus += d["minus"]; reads += d["reads"]; nat += d["nat"]; nte += d["nte"]; n1u += d["n1u"]; npl += d["nplus"]; nmi += d["nminus"]; tes = d["tes"]
if tes is None: tes = te_at(TOP, CH, ps, pe)
ntot = npl + nmi; pct_minus = 100 * nmi / max(1, ntot); arch = "dual-strand" if min(nmi, ntot - nmi) / max(1, ntot) > 0.2 else "uni-strand"; pct_at = 100 * nat / max(1, nte); domTE = dom_te_family(tes, ps, pe)
def testr(p):
    for ts, te, st, fam in tes:
        if ts <= p < te: return st
    return None
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(13, 11.5), dpi=300); gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 2.0, 1.4], hspace=0.55); fig.subplots_adjust(top=0.9, bottom=0.07)
axA, axB, axC = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
x = np.arange(len(ORDER)); bw = 0.26
for j, tp in enumerate(TPS):
    h = FPM[tp].values; axA.bar(x + (j - 1) * bw, np.maximum(h, 1e-3), width=bw, color=TPCOL[tp], edgecolor="white", linewidth=0.2, label=TPLAB[tp])
    for xi, X in enumerate(ORDER):       # split each FPM bar by HEIGHT → solid + (bottom) / pale − (top)
        if h[xi] > 0 and (X, tp) in SF: pc.panelA_strand(axA, xi + (j - 1) * bw, h[xi], bw, SF[(X, tp)], TPCOL[tp])
    for xi in range(len(ORDER)):
        if h[xi] > 0: axA.text(xi + (j - 1) * bw, h[xi] * 1.18, f"{h[xi]:.0f}" if h[xi] >= 1 else f"{h[xi]:.1f}", ha="center", va="bottom", fontsize=4.8, rotation=90, color=TPCOL[tp], fontweight="bold")
axA.set_yscale("log"); axA.set_ylim(0.1, max(FPM.values.max(), 1) * 9); axA.set_xticks(x); axA.set_xticklabels([s.replace("_", "/") for s in ORDER], rotation=90, fontsize=7.5)
for t, X in zip(axA.get_xticklabels(), ORDER):
    if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
for xi, X in enumerate(ORDER):
    pv = PAV[X]; fc = ("#333" if pv == "present" else "white"); ec = ("#333" if pv == "present" else ("#C0392B" if pv == "absent" else "#ccc"))
    axA.plot(xi, 0.135, marker="o", markersize=4.2, markerfacecolor=fc, markeredgecolor=ec, markeredgewidth=0.9, clip_on=False)
axA.set_ylabel("PICB cluster\nFPM (log)", fontsize=8); axA.spines[["top", "right"]].set_visible(False)
axA.legend([TPLAB[t] for t in TPS], fontsize=7.5, frameon=False, ncol=1, loc="upper left", bbox_to_anchor=(1.005, 1.0), title="timepoint", title_fontsize=7, handlelength=1.2)
axA.plot([], [], "o", markersize=4.2, markerfacecolor="#333", markeredgecolor="#333", label="locus in genome")
axA.plot([], [], "o", markersize=4.2, markerfacecolor="white", markeredgecolor="#C0392B", label="locus absent (loss)")
axA.legend(fontsize=6.0, frameon=False, ncol=1, loc="upper left", bbox_to_anchor=(1.005, 0.42), handlelength=1.0, handletextpad=0.3)
pc.pbadge(axA, "A", "Pangenome × timepoint — PICB cluster expression across 16 strains (FPM, log)   ·   ● locus present in genome   ○ lost", fs=8.4, y=1.10)
axA.text(0.012, 1.02, "each FPM bar split by HEIGHT (proportion):  solid = + strand  ·  pale = − strand", transform=axA.transAxes, fontsize=6.0, color="#777", style="italic", ha="left", va="center")
# B: TOP strain — per-timepoint sRNA coverage block (canonical multi-figure style, shared helper)
N = max(1, pe - ps)
axB.set_xlim(0, 1); axB.set_ylim(-3.5, 2.1); axB.axis("off"); fig.canvas.draw()   # 0–1 fractional x (like the multi figure); lims fixed before pc.rtext()
_f5 = Counter(r[0] for r in reads if not r[2]); _r5 = Counter(r[1] - 1 for r in reads if r[2]); _rn = lambda i: sum(_r5[j] for j in range(i - 2, i + 13))   # PRIORITISE Panel C zoom on a BOTH-strands region (sense+antisense), so both red and grey show
_cd = [(min(_f5[i], _rn(i)), _f5[i] + _rn(i), i + 5) for i in _f5 if _rn(i) > 0]   # rank by min(fwd,rev) → both strands strongest, then by total
pk = max(_cd)[2] if _cd else ((Counter((r[1] - 1 if r[2] else r[0]) for r in reads).most_common(1)[0][0]) if reads else ps + N // 2); z0, z1 = pk - 30, pk + 50
ztpL = " + ".join(TPLAB[t] for t in TPS if per[t] is not None and any(r[0] < z1 and r[1] > z0 for r in per[t]["reads"])) or "—"   # timepoints pooled into Panel C
_ntot, _domTE, _arch, _pas, _p1u = pc.draw_strain_block(axB, per, tes, CH, ps, pe, fps, 0.0, z0, z1, name=TOP, is_top=True, wild=(TOP in WILD), TECOL=TECOL)
_fams = list(dict.fromkeys((f.split("|")[-1] if "|" in f else f) for _, _, _, f in tes))[:5]
_lh = [Patch(facecolor=pc.PLUS_COL[t], label=TPLAB[t]) for t in TPS] + [Patch(facecolor="#6a3d9a", label="solid = + strand"), Patch(facecolor=pc.pale("#6a3d9a", 0.55), label="pale = − strand"), Patch(facecolor="#efefef", label="non-TE piRNA"), Patch(facecolor="#C0392B", label="antisense-to-TE = silencing (red outline); bars = TE family"), Patch(facecolor="#cfcfcf", label="sense-to-TE"), Patch(facecolor="#c9d6ea", edgecolor=pc.C_GENE, label="gene model")] + [Patch(facecolor=pc.famcol(f), label=f) for f in _fams]
pc.color_legend(axB.legend(handles=_lh, fontsize=6.0, frameon=False, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.06), columnspacing=1.0, handlelength=1.0), _lh)   # right-margin element: TOP = TE families ranked × strand + non-TE; BOTTOM = antisense/sense-to-TE (silencing) of TE-mapped piRNA
pc.pbadge(axB, "B", f"Top strain {TOP.replace('_','/')} ({BAMC}:{ps:,}–{pe:,}) — per-timepoint sRNA coverage (height = expression, colour = timepoint) above TE + gene tracks · {_ntot:,} primary piRNAs · {_p1u:.0f}% 1U", fs=7.4, y=1.07)
axB.text(0.012, 1.02, "‘primary reads’ = each sRNA read (24–32 nt) counted once at its STAR primary locus (multimappers kept, not double-counted) · architecture (genomic strand) ≠ sense/antisense (relative to TE)", transform=axB.transAxes, fontsize=5.2, color="#8a8a8a", style="italic", ha="left", va="center")
def anti_te(rs, re, isrev):
    s = testr((rs + re) // 2); return None if s is None else ((s == "-") != isrev)
zr = [r for r in reads if r[0] < z1 and r[1] > z0]; pr, mr = Counter(), Counter()
for rs, re, isrev, seq in zr: (mr if isrev else pr)[(rs, re, seq, isrev)] += 1
def draw(items, y0, dirn):
    y = y0
    for (rs, re, seq, isrev), cnt in items:
        for k, ch in enumerate(seq):
            xx = rs + k
            if z0 - 2 <= xx <= z1 + 2: axC.text(xx, y, ch, fontsize=4.6, ha="center", va="center", family="monospace", color="white", bbox=dict(boxstyle="square,pad=0.02", fc=NT.get(ch, "#999"), ec="none"))
        a = anti_te(rs, re, isrev); acol = "#C0392B" if a else ("#888" if a is not None else "#ccc"); fp = re - 1 if isrev else rs
        axC.annotate("", xy=(fp + (0.5 if not isrev else -0.5), y), xytext=(fp + (-2.6 if not isrev else 2.6), y), arrowprops=dict(arrowstyle="-|>", color=acol, lw=1.0))
        axC.text(z1 + 4, y, f"×{cnt}", fontsize=4.6, va="center", color="#666"); y += dirn
    return y
ytop = draw(pr.most_common(7), 1, 1); axC.axhline(0, color="#333", lw=0.7); ybot = draw(mr.most_common(7), -1, -1)
axC.set_xlim(z0 - 1, z1 + 8); axC.set_ylim(ybot - 1.6, ytop + 0.6)
for sp in ("top", "left", "right"): axC.spines[sp].set_visible(False)
axC.set_yticks([]); axC.spines["bottom"].set_position(("data", ybot - 1.0))
tk = np.linspace(z0, z1, 5).astype(int); axC.set_xticks(tk); axC.set_xticklabels([f"{t:,}" for t in tk], fontsize=6.5); axC.tick_params(axis="x", length=3)
axC.set_xlabel(f"{TOP.replace('_','/')} chr{CH} position (bp) — every base at its true genomic coordinate", fontsize=7)
axC.text(z0 - 0.6, ytop, "+ strand", fontsize=6.5, color="#33a02c", fontweight="bold", ha="right", va="center"); axC.text(z0 - 0.6, ybot, "− strand", fontsize=6.5, color="#6a3d9a", fontweight="bold", ha="right", va="center")
fig.add_artist(ConnectionPatch(xyA=((z0 - ps) / N, -1.30), coordsA=axB.transData, xyB=(z0, ytop + 0.6), coordsB=axC.transData, color="#E8A33D", lw=0.8, ls=(0, (3, 2))))
fig.add_artist(ConnectionPatch(xyA=((z1 - ps) / N, -1.30), coordsA=axB.transData, xyB=(z1, ytop + 0.6), coordsB=axC.transData, color="#E8A33D", lw=0.8, ls=(0, (3, 2))))
pc.pbadge(axC, "C", f"Base resolution, {TOP.replace('_','/')} — pooled {ztpL}   ·   5′-U = 1U   ·   5′ arrow RED = antisense-to-TE (silencing), grey = sense-to-TE", fs=7.6)
fig.suptitle(f"{GENE} — PICB piRNA cluster: pangenome cross-strain × timepoint + nucleotide resolution", fontsize=12, fontweight="bold", y=0.965)
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/{OUT}.{e}", bbox_inches="tight")
print(f"   wrote {OUT}.png")
