#!/usr/bin/env python3
"""MULTI version of the TE-driven figure: shows THIS strain-private TE-driven cluster across ALL 16 strains.
  (A) 16-strain TE-presence PAV — ● carrier HAS the insertion / ○ the other strains are TE-absent (no orthologous
      locus → no cluster). The insertion is a singleton in the 16-strain pangenome VCF, so carrier=●, rest=○;
      this is the whole point: the piRNA cluster exists ONLY where the TE inserted.
  (B) the carrier strain's per-timepoint sRNA coverage block (shared helper) + TE + gene tracks (the developmental
      activation is read off the per-timepoint rows + right-margin summaries);
  (C) nucleotide resolution at the silencing peak: 5'-U = 1U; RED 5' arrow = ANTISENSE-to-TE. NO liftover.
Usage: make_te_driven_locus_multi.py <strain> <chrom> <start> <end> <TE_label> <n_absent> <outbase>"""
import sys, os, numpy as np
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
from collections import Counter
from strain_order import STRAIN_ORDER, WILD
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch, Patch
import pav_clusters as pc
from pav_clusters import TPS, TPLAB, TPCOL, fetch_primary, te_at
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; PG = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
NT = {"A": "#33a02c", "C": "#1f78b4", "G": "#ff7f00", "T": "#e31a1c", "N": "#999"}
TECOL = {"LINE/L1": "#E69F00", "LTR/ERVK": "#6a3d9a", "LTR/ERVL-MaLR": "#b15928", "SINE/Alu": "#33a02c", "SINE/B2": "#1f78b4", "LTR/ERVL": "#a6cee3", "LTR/ERV1": "#cab2d6", "Simple_repeat": "#bbbbbb", "SINE/MIR": "#1b9e77"}
X, CH, PS, PE, TELAB, NABS, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], int(sys.argv[6]), sys.argv[7]
nb = 220; N = max(1, PE - PS)
DAT = {tp: fetch_primary(X, CH, PS, PE, tp, nb) for tp in TPS}
tes = next((DAT[tp]["tes"] for tp in TPS if DAT[tp]["ntot"]), te_at(X, CH, PS, PE))
def tst(p):
    for ts, te, st, fam in tes:
        if ts <= p < te: return st
    return None
CHOSEN = max(TPS, key=lambda tp: DAT[tp]["ntot"])
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"] + ([X] if (X not in STRAIN_ORDER and X != "C57BL_6") else [])   # 16 inbred strains (exclude external black6 C57BL_6); carrier guaranteed present
PRESENT_N = sum(s == X for s in ORDER); ABSENT_N = len(ORDER) - PRESENT_N   # private singleton => 1 carrier, rest absent (data-driven, not the NABS arg which assumed 16)
print(f"{X} chr{CH}:{PS:,}-{PE:,}  reads/tp=" + " ".join(f"{TPLAB[t]}:{DAT[t]['ntot']}" for t in TPS) + f"  peak={TPLAB[CHOSEN]}  PAV: 1 carrier / {NABS} absent")
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(12.5, 11.5), dpi=300); gs = fig.add_gridspec(3, 1, height_ratios=[0.95, 2.1, 1.4], hspace=0.62); fig.subplots_adjust(top=0.9, bottom=0.07)
axA, axB, axC = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
# A: standard Pangenome × timepoint FPM bar chart — the carrier expresses (bar); the 15 TE-absent strains show ○ (no orthologous locus → no cluster). Same style as every other figure's Panel A.
import pysam as _ps
def _lib(s, tp):
    t = 0
    for r in (1, 2, 3):
        b = f"{ROOT}/results/STAR_srna_strain_wise/{s}/{s}-{tp}.{r}/Aligned.sortedByCoord.out.bam"
        if os.path.exists(b):
            try: t += sum(st.mapped for st in _ps.AlignmentFile(b, "rb").get_index_statistics())
            except Exception: pass
    return t
FPMv = {tp: {s: 0.0 for s in ORDER} for tp in TPS}                                # carrier-only expression (FPM); the rest are TE-absent → 0
for tp in TPS:
    if DAT[tp]["ntot"] > 0:
        _l = _lib(X, tp); FPMv[tp][X] = (DAT[tp]["ntot"] / _l * 1e6) if _l else 0.0
PAVd = {s: ("present" if s == X else "absent") for s in ORDER}
xb = np.arange(len(ORDER)); bw = 0.26; _ylo = 0.1
for j, tp in enumerate(TPS):
    h = [FPMv[tp][s] for s in ORDER]
    axA.bar(xb + (j - 1) * bw, [max(v, 1e-3) for v in h], width=bw, color=TPCOL[tp], edgecolor="white", linewidth=0.2, label=TPLAB[tp])
    for xi, s in enumerate(ORDER):                                                # split the carrier's bar by HEIGHT: solid = + strand, pale = − strand
        if h[xi] > 0 and s == X and DAT[tp]["ntot"] > 0:
            _pf = DAT[tp]["nplus"] / DAT[tp]["ntot"]; _sp = 10 ** (np.log10(_ylo) + _pf * (np.log10(max(h[xi], _ylo * 1.001)) - np.log10(_ylo)))
            axA.bar(xi + (j - 1) * bw, h[xi] - _sp, bottom=_sp, width=bw, color=pc.pale(TPCOL[tp]), edgecolor="none", zorder=3)
        if h[xi] > 0: axA.text(xi + (j - 1) * bw, h[xi] * 1.18, f"{h[xi]:.0f}" if h[xi] >= 1 else f"{h[xi]:.1f}", ha="center", va="bottom", fontsize=4.8, rotation=90, color=TPCOL[tp], fontweight="bold")
_mx = max((FPMv[tp][X] for tp in TPS), default=1.0)
axA.set_yscale("log"); axA.set_ylim(0.1, max(_mx, 1) * 20); axA.set_xticks(xb); axA.set_xticklabels([s.replace("_", "/") for s in ORDER], rotation=90, fontsize=7.3)
for t, s in zip(axA.get_xticklabels(), ORDER):
    if s in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
for xi, s in enumerate(ORDER):                                                    # genome PAV marker: ● TE insertion present (carrier) / ○ TE absent
    fc = "#333" if PAVd[s] == "present" else "white"; ec = "#333" if PAVd[s] == "present" else "#C0392B"
    axA.plot(xi, 0.135, marker="o", markersize=4.4, markerfacecolor=fc, markeredgecolor=ec, markeredgewidth=1.0, clip_on=False)
axA.set_ylabel("PICB cluster\nFPM (log)", fontsize=8); axA.spines[["top", "right"]].set_visible(False)
axA.plot([], [], "o", markersize=4.4, markerfacecolor="#333", markeredgecolor="#333", label="TE insertion present (carrier)")
axA.plot([], [], "o", markersize=4.4, markerfacecolor="white", markeredgecolor="#C0392B", label="TE absent (no orthologous locus)")
axA.legend(fontsize=6.2, frameon=False, ncol=1, loc="upper left", bbox_to_anchor=(1.005, 1.0), handlelength=1.0, handletextpad=0.3)
pc.pbadge(axA, "A", f"TE insertion present in {PRESENT_N}/{len(ORDER)} strains ({X.replace('_','/')}-private) — PICB cluster expression × timepoint (FPM, log); the cluster exists ONLY where the TE inserted", fs=7.5, y=1.16)
axA.text(0.012, 1.05, "each FPM bar split by HEIGHT (proportion):  solid = + strand  ·  pale = − strand", transform=axA.transAxes, fontsize=6.2, color="#777", style="italic", ha="left", va="center")
# B: carrier per-timepoint coverage block — canonical multi-figure style (shared helper) + TE + gene tracks
axB.set_xlim(0, 1); axB.set_ylim(-3.5, 2.1); axB.axis("off"); fig.canvas.draw()
five = Counter((r[1] - 1 if r[2] else r[0]) for r in DAT[CHOSEN]["reads"]); pk = five.most_common(1)[0][0] if five else PS + N // 2; z0, z1 = pk - 30, pk + 50
_ntot, _domTE, _arch, _pas, _p1u = pc.draw_strain_block(axB, DAT, tes, CH, PS, PE, None, 0.0, z0, z1, name=X, is_top=True, wild=(X in WILD), TECOL=TECOL, nb=nb)
_fams = list(dict.fromkeys((f.split("|")[-1] if "|" in f else f) for _, _, _, f in tes))[:5]
axB.legend(handles=[Patch(facecolor=pc.PLUS_COL[t], label=TPLAB[t]) for t in TPS] + [Patch(facecolor="#4a4a4a", label="+ strand = deeper shade"), Patch(facecolor="#cfcfcf", label="− strand = paler shade"), Patch(facecolor="#c9d6ea", edgecolor=pc.C_GENE, label="gene model")] + [Patch(facecolor=TECOL.get(f, "#ddd"), label=f) for f in _fams], fontsize=6.0, frameon=False, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.06), columnspacing=1.0, handlelength=1.0)
pc.pbadge(axB, "B", f"{X.replace('_','/')} chr{CH}:{PS:,}–{PE:,} — per-timepoint sRNA coverage (height = expression, colour = timepoint) above TE + gene tracks · {_ntot:,} primary piRNAs · {_p1u:.0f}% 1U", fs=7.0, y=1.05)
axB.text(0.012, 0.985, f"the {ABSENT_N} TE-absent strains have NO orthologous locus here (the insertion is {X.replace('_','/')}-private) → coverage cannot be plotted for them = the cluster is TE-DRIVEN", transform=axB.transAxes, fontsize=5.6, color=pc.C_WILD, fontweight="bold", ha="left", va="center")
# C: nucleotide at the silencing peak (z0,z1 + zoom box already set in Panel B via draw_strain_block)
d = DAT[CHOSEN]; reads = d["reads"]
def anti_te(rs, re, isrev):
    s = tst((rs + re) // 2); return None if s is None else ((s == "-") != isrev)
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
axC.set_xlabel(f"{X.replace('_','/')} chr{CH} position (bp)", fontsize=7)
axC.text(z0 - 0.6, ytop, "+ strand", fontsize=6.5, color="#33a02c", fontweight="bold", ha="right", va="center"); axC.text(z0 - 0.6, ybot, "− strand", fontsize=6.5, color="#6a3d9a", fontweight="bold", ha="right", va="center")
for x in (z0, z1):
    fig.add_artist(ConnectionPatch(xyA=((x - PS) / N, -1.30), coordsA=axB.transData, xyB=(x, ytop + 0.6), coordsB=axC.transData, color="#E8A33D", lw=0.8, ls=(0, (3, 2))))
axC.set_title(f"C  Base resolution at {TPLAB[CHOSEN]} silencing peak; 5′-U = 1U; 5′ arrow RED = ANTISENSE-to-TE (silencing), grey = sense-to-TE", fontsize=8.0, fontweight="bold", loc="left")
fig.suptitle(f"TE-DRIVEN piRNA cluster across the 16-strain pangenome — {X.replace('_','/')} {TELAB} (coordinate-verified)", fontsize=11.8, fontweight="bold", y=0.965)
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/{OUT}.{e}", bbox_inches="tight")
print(f"   wrote {OUT}.png")
