#!/usr/bin/env python3
"""SINGLE-strain locus plot for a PICB piRNA cluster, resolved by TIMEPOINT (the developmental dimension).
DATA SOURCES ONLY: pangenome cluster projection (picb_pangenome_clusters.tsv = PICB-COMBINED clusters with
own+GRCm39 coords), sRNA BAMs (PRIMARY reads inside the cluster), RepeatMasker TE annotation, pangenome PAV.
NO per-figure liftover; NO multimapping spillover.
  (A) this strain's PICB-combined cluster FPM across E16.5/P12.5/P20.5 — WHEN the cluster is ON (★ = nucleotide tp);
  (B) for EACH timepoint: PRIMARY-read coverage INSIDE the PICB cluster (plus ↑ green / minus ↓ purple) when the
      cluster is ON, or 'no PICB cluster (OFF)' when off — the developmental 'why'; + RepeatMasker TE track;
  (C) base resolution in the chosen (top-FPM) timepoint.
TERMINOLOGY: genomic +/- = architecture, ≠ sense/antisense (relative to TE, RM .out strand; antisense-to-TE = silencing).
Usage: make_pav_locus_single.py <g39chrom> <g39start> <g39end> <gene> <outbase> [strain]."""
import sys, os, numpy as np
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
from strain_order import STRAIN_ORDER, WILD
from collections import Counter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch, Patch
import pav_clusters as pc
from pav_clusters import TPS, TPLAB, TPCOL, clusters_at, present_strains, fpm_by_tp, cluster_extent, fetch_primary, te_at, dom_te_family, genome_pav
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]
NT = {"A": "#33a02c", "C": "#1f78b4", "G": "#ff7f00", "T": "#e31a1c", "N": "#999"}; comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
TECOL = {"LINE/L1": "#E69F00", "LTR/ERVK": "#6a3d9a", "LTR/ERVL-MaLR": "#b15928", "SINE/Alu": "#33a02c", "SINE/B2": "#1f78b4", "LTR/ERVL": "#a6cee3", "LTR/ERV1": "#cab2d6", "Simple_repeat": "#bbbbbb"}
G39C, G39S, G39E, GENE, OUT = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
WANT = sys.argv[6] if len(sys.argv) > 6 else None
nb = 200
sub = clusters_at(G39C, G39S, G39E)
present = present_strains(sub, ORDER)
import pandas as pd
FPM = sub.groupby(["strain", "tp"])["all_primary_FPM"].max().unstack(fill_value=0.0).reindex(index=ORDER, columns=TPS).fillna(0.0)
STRAIN = WANT if (WANT and WANT in present) else (max(present, key=lambda X: FPM.loc[X].max()) if present else (WANT or ORDER[0]))
fps = fpm_by_tp(sub, STRAIN); fpm_tp = FPM.loc[STRAIN]
CHOSEN = max(fps, key=fps.get) if fps else TPS[0]
ext = cluster_extent(sub, STRAIN)
if ext is None:
    print(f"{GENE}: {STRAIN} has NO PICB cluster at this locus — nothing to plot"); sys.exit(0)
CH, ps, pe = ext[0], ext[1], ext[2]; BAMC = f"{STRAIN}#1#chr{CH}"; N = max(1, pe - ps)
PAVst = "present" if (fps and max(fps.values()) > 0) else genome_pav(G39C, G39S, G39E, STRAIN)   # expression PROVES presence (halLiftover g39->strain can false-negative)
DAT = {tp: (fetch_primary(STRAIN, CH, ps, pe, tp, nb) if tp in fps else None) for tp in TPS}
tes = next((DAT[tp]["tes"] for tp in TPS if DAT[tp]), te_at(STRAIN, CH, ps, pe))
domTE = dom_te_family(tes, ps, pe)
CHOSEN = max((tp for tp in TPS if DAT[tp] and DAT[tp]["ntot"] > 0), key=lambda tp: DAT[tp]["ntot"], default=CHOSEN)  # nucleotide tp must have data (some tp BAMs missing, e.g. C57BL_6NJ E16.5)
print(f"{GENE}: single-strain={STRAIN}; FPM by tp={dict(fpm_tp.round(1))}; cluster {BAMC}:{ps:,}-{pe:,}; domTE={domTE}; genome={PAVst}")
def testr(p):
    for ts, te, st, fam in tes:
        if ts <= p < te: return st
    return None
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(12.5, 11.5), dpi=300); gs = fig.add_gridspec(3, 1, height_ratios=[0.85, 2.3, 1.4], hspace=0.6); fig.subplots_adjust(top=0.9, bottom=0.06)
axA, axB, axC = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
xb = np.arange(3)
axA.bar(xb, [max(fpm_tp[t], 1e-2) for t in TPS], width=0.6, color=[TPCOL[t] for t in TPS], edgecolor="white")
for _i, _t in enumerate(TPS):       # split each FPM bar by HEIGHT: solid + strand (bottom) / pale − strand (top)
    if DAT[_t] and DAT[_t]["ntot"] > 0: pc.panelA_strand(axA, _i, fpm_tp[_t], 0.6, DAT[_t]["nplus"] / DAT[_t]["ntot"], TPCOL[_t], ylo=0.05)
axA.set_yscale("log"); axA.set_ylim(0.05, max(fpm_tp.max() * 6, 1)); axA.set_xticks(xb); axA.set_xticklabels([TPLAB[t] for t in TPS], fontsize=9)
for i, t in enumerate(TPS):
    if fpm_tp[t] > 0: axA.text(i, fpm_tp[t] * 1.15, f"{fpm_tp[t]:.1f}" + ("  ●" if t == CHOSEN else ""), ha="center", va="bottom", fontsize=8, fontweight="bold", color=TPCOL[t])
    else: axA.text(i, 0.06, "OFF\n(no cluster)", ha="center", va="bottom", fontsize=6.5, color="#bbb", style="italic")
axA.set_ylabel("PICB FPM (log)", fontsize=8.5); axA.spines[["top", "right"]].set_visible(False)
pc.pbadge(axA, "A", f"{STRAIN.replace('_','/')} — PICB-cluster FPM across development (log); ON/OFF per timepoint; ● = nucleotide tp (C); locus in genome: {PAVst}", fs=8.2, y=1.10)
axB.set_xlim(0, 1); axB.set_ylim(-3.5, 2.1); axB.axis("off"); fig.canvas.draw()   # 0–1 fractional x (like the multi figure)
_dc = DAT[CHOSEN]; reads = _dc["reads"] if _dc else []
_f5 = Counter(r[0] for r in reads if not r[2]); _r5 = Counter(r[1] - 1 for r in reads if r[2]); _rn = lambda i: sum(_r5[j] for j in range(i - 2, i + 13))   # PRIORITISE Panel C zoom on a BOTH-strands region (sense+antisense), so both red and grey show
_cd = [(min(_f5[i], _rn(i)), _f5[i] + _rn(i), i + 5) for i in _f5 if _rn(i) > 0]   # rank by min(fwd,rev) → both strands strongest, then by total
pk = max(_cd)[2] if _cd else ((Counter((r[1] - 1 if r[2] else r[0]) for r in reads).most_common(1)[0][0]) if reads else ps + N // 2); z0, z1 = pk - 30, pk + 50
_ntot, _domTE, _arch, _pas, _p1u = pc.draw_strain_block(axB, DAT, tes, CH, ps, pe, fpm_tp, 0.0, z0, z1, name=STRAIN, is_top=True, wild=(STRAIN in WILD), TECOL=TECOL, only_tp=CHOSEN)
_fams = list(dict.fromkeys((f.split("|")[-1] if "|" in f else f) for _, _, _, f in tes))[:5]
_lh = [Patch(facecolor=pc.PLUS_COL[t], label=TPLAB[t]) for t in TPS] + [Patch(facecolor="#6a3d9a", label="solid = + strand"), Patch(facecolor=pc.pale("#6a3d9a", 0.55), label="pale = − strand"), Patch(facecolor="#efefef", label="non-TE piRNA"), Patch(facecolor="#C0392B", label="antisense-to-TE = silencing (red outline); bars = TE family"), Patch(facecolor="#cfcfcf", label="sense-to-TE"), Patch(facecolor="#c9d6ea", edgecolor=pc.C_GENE, label="gene model")] + [Patch(facecolor=pc.famcol(f), label=f) for f in _fams]
pc.color_legend(axB.legend(handles=_lh, fontsize=6.0, frameon=False, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.06), columnspacing=1.0, handlelength=1.0), _lh)
pc.pbadge(axB, "B", f"{STRAIN.replace('_','/')} {BAMC}:{ps:,}–{pe:,} — per-timepoint sRNA coverage (height = expression, colour = timepoint) above TE + gene tracks · {_ntot:,} primary piRNAs", fs=7.4, y=1.07)
axB.text(0.012, 1.02, "‘primary reads’ = each sRNA read (24–32 nt) counted once at its STAR primary locus (multimappers kept, not double-counted) · architecture (genomic strand) ≠ sense/antisense (relative to TE)", transform=axB.transAxes, fontsize=5.2, color="#8a8a8a", style="italic", ha="left", va="center")
def anti_te(rs, re, isrev):
    st = testr((rs + re) // 2); return None if st is None else ((st == "-") != isrev)
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
axC.set_xlabel(f"{STRAIN.replace('_','/')} chr{CH} position (bp) — every base at its true genomic coordinate", fontsize=7)
axC.text(z0 - 0.6, ytop, "+ strand", fontsize=6.5, color="#33a02c", fontweight="bold", ha="right", va="center"); axC.text(z0 - 0.6, ybot, "− strand", fontsize=6.5, color="#6a3d9a", fontweight="bold", ha="right", va="center")
for x in (z0, z1):
    fig.add_artist(ConnectionPatch(xyA=((x - ps) / N, -1.30), coordsA=axB.transData, xyB=(x, ytop + 0.6), coordsB=axC.transData, color="#E8A33D", lw=0.8, ls=(0, (3, 2))))
pc.pbadge(axC, "C", f"Base resolution, {STRAIN.replace('_','/')} at {TPLAB[CHOSEN]} (top-FPM tp)   ·   5′-U = 1U   ·   5′ arrow RED = antisense-to-TE (silencing), grey = sense-to-TE", fs=7.6)
fig.suptitle(f"{GENE} — SINGLE-strain PICB cluster across development ({STRAIN.replace('_','/')})", fontsize=12, fontweight="bold", y=0.965)
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/{OUT}.{e}", bbox_inches="tight")
print(f"   wrote {OUT}.png")
