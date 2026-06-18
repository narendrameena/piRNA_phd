#!/usr/bin/env python3
"""Strain-private piRNA SOURCE LOCUS (individual piRNAs, NOT a PICB cluster) rendered in the PANGENOME-MULTI
LAYOUT (matches make_pav_locus_multi/single so the two figure families share one look while staying distinct
in content):
  (A) Pangenome × timepoint — source-locus RPM across the 16 strains; for a strain-PRIVATE insertion the carrier
      is present (●) and the other 15 are genetically absent (○) → 'present in 1/16';
  (B) per-timepoint sRNA coverage in the carrier (height = expression, deep ↑+ / pale ↓−) over TE + gene tracks;
  (C) base resolution at true genomic coordinates (5'-U = 1U; 5' arrow RED = antisense-to-TE = silencing).
Driven DIRECTLY from the carrier's strain-wise BAMs (these loci are BELOW PICB's cluster threshold, so they are
not in the PICB pangenome table), own coordinates, NO per-figure liftover.
Usage: make_source_pav.py <STRAIN> <CHROM_own e.g. CAST_EiJ#1#chr15> <start> <end> <te_label> <te_strand> <outbase>"""
import sys, os, numpy as np, pysam
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
from strain_order import STRAIN_ORDER, WILD
from collections import Counter
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch, Patch
import pav_clusters as pc
from pav_clusters import TPS, TPLAB, TPCOL, fetch_primary, te_at, dom_te_family
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"; RB = f"{ROOT}/results/STAR_srna_strain_wise"
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]
NT = {"A": "#33a02c", "C": "#1f78b4", "G": "#ff7f00", "T": "#e31a1c", "N": "#999"}; comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
TECOL = {"LINE/L1": "#E69F00", "LTR/ERVK": "#6a3d9a", "LTR/ERVL-MaLR": "#b15928", "SINE/Alu": "#33a02c", "SINE/B2": "#1f78b4", "LTR/ERVL": "#a6cee3", "LTR/ERV1": "#cab2d6", "Simple_repeat": "#bbbbbb"}
STRAIN, CHROM, S, E, TELAB, TEST, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6], sys.argv[7]
CH = CHROM.split("chr")[-1]; ps, pe = S, E; N = max(1, pe - ps); BAMC = CHROM; nb = 200
# per-timepoint PRIMARY coverage from the carrier BAMs (own coords) — these loci are not PICB clusters
DAT = {tp: fetch_primary(STRAIN, CH, ps, pe, tp, nb) for tp in TPS}
DAT = {tp: (d if d and d["ntot"] > 0 else None) for tp, d in DAT.items()}
tes = next((DAT[tp]["tes"] for tp in TPS if DAT[tp]), te_at(STRAIN, CH, ps, pe))
domTE = dom_te_family(tes, ps, pe)
def total_mapped(strain, tp):   # genome-mapped reads (pooled reps) for RPM normalisation
    tot = 0
    for r in (1, 2, 3):
        b = f"{RB}/{strain}/{strain}-{tp}.{r}/Aligned.sortedByCoord.out.bam"
        if os.path.exists(b):
            try: tot += sum(s.mapped for s in pysam.AlignmentFile(b, "rb").get_index_statistics())
            except Exception: pass
    return tot
rpm = {tp: ((DAT[tp]["ntot"] / total_mapped(STRAIN, tp) * 1e6) if (DAT[tp] and total_mapped(STRAIN, tp)) else 0.0) for tp in TPS}
CHOSEN = max((tp for tp in TPS if DAT[tp]), key=lambda tp: DAT[tp]["ntot"], default=TPS[0])
_rc = {t: (DAT[t]["ntot"] if DAT[t] else 0) for t in TPS}; _rp = {t: round(rpm[t], 1) for t in TPS}
print(f"{OUT}: {STRAIN} {BAMC}:{ps:,}-{pe:,}; reads/tp={_rc}; RPM={_rp}; domTE={domTE}")
def testr(p):
    for ts, te, st, fam in tes:
        if ts <= p < te: return st
    return None
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(13, 11.5), dpi=300); gs = fig.add_gridspec(3, 1, height_ratios=[1.05, 2.3, 1.4], hspace=0.62); fig.subplots_adjust(top=0.9, bottom=0.07)
axA, axB, axC = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
# ---------- Panel A: pangenome x timepoint (carrier present, others absent for a private insertion) ----------
xb = np.arange(len(ORDER)); w = 0.26
for j, tp in enumerate(TPS):
    vals = [rpm[tp] if s == STRAIN else 0.0 for s in ORDER]
    axA.bar(xb + (j - 1) * w, [max(v, 1e-2) for v in vals], width=w, color=TPCOL[tp], edgecolor="white", lw=0.3, label=TPLAB[tp])
    if DAT[tp] and rpm[tp] > 0: pc.panelA_strand(axA, ORDER.index(STRAIN) + (j - 1) * w, rpm[tp], w, DAT[tp]["nplus"] / DAT[tp]["ntot"], TPCOL[tp], ylo=0.05)   # split carrier bar by HEIGHT: solid + / pale −
axA.set_yscale("log"); axA.set_ylim(0.05, max(max(rpm.values()), 1) * 8)
for i, s in enumerate(ORDER):
    axA.plot(i, 0.066, "o", ms=4.2, mfc=("#111" if s == STRAIN else "white"), mec=("#111" if s == STRAIN else "#bbb"), mew=0.7)
axA.set_xticks(xb); axA.set_xticklabels([s.replace("_", "/") for s in ORDER], rotation=45, ha="right", fontsize=6.4)
for t, s in zip(axA.get_xticklabels(), ORDER):
    t.set_color(pc.C_WILD if s in WILD else "#333"); t.set_fontweight("bold" if s == STRAIN else "normal")
axA.set_ylabel("source-locus RPM (log)", fontsize=8.5); axA.spines[["top", "right"]].set_visible(False)
axA.legend(fontsize=6.2, frameon=False, loc="upper left", ncol=3, columnspacing=1.0, handlelength=1.0)
pc.pbadge(axA, "A", f"Pangenome × timepoint — strain-private piRNA SOURCE LOCUS (individual piRNAs) · present in 1/16 ({STRAIN.replace('_','/')}, private insertion) · ● present  ○ genetically absent", fs=8.0, y=1.13)
# ---------- Panel B: per-timepoint coverage in the carrier (shared draw_strain_block) ----------
axB.set_xlim(0, 1); axB.set_ylim(-3.5, 2.1); axB.axis("off"); fig.canvas.draw()
_dc = DAT[CHOSEN]; reads = _dc["reads"] if _dc else []
_f5 = Counter(r[0] for r in reads if not r[2]); _r5 = Counter(r[1] - 1 for r in reads if r[2]); _rn = lambda i: sum(_r5[j] for j in range(i - 2, i + 13))   # zoom on a BOTH-strands region so red+grey both show in C
_cd = [(min(_f5[i], _rn(i)), _f5[i] + _rn(i), i + 5) for i in _f5 if _rn(i) > 0]
pk = max(_cd)[2] if _cd else ((Counter((r[1] - 1 if r[2] else r[0]) for r in reads).most_common(1)[0][0]) if reads else ps + N // 2); z0, z1 = pk - 30, pk + 50
_ntot, _domTE, _arch, _pas, _p1u = pc.draw_strain_block(axB, DAT, tes, CH, ps, pe, None, 0.0, z0, z1, name=STRAIN, is_top=True, wild=(STRAIN in WILD), TECOL=TECOL, only_tp=CHOSEN)
_fams = list(dict.fromkeys((f.split("|")[-1] if "|" in f else f) for _, _, _, f in tes))[:5]
_lh = [Patch(facecolor=pc.PLUS_COL[t], label=TPLAB[t]) for t in TPS] + [Patch(facecolor="#6a3d9a", label="solid = + strand"), Patch(facecolor=pc.pale("#6a3d9a", 0.55), label="pale = − strand"), Patch(facecolor="#efefef", label="non-TE piRNA"), Patch(facecolor="#C0392B", label="antisense-to-TE = silencing (red outline); bars = TE family"), Patch(facecolor="#cfcfcf", label="sense-to-TE"), Patch(facecolor="#c9d6ea", edgecolor=pc.C_GENE, label="gene model")] + [Patch(facecolor=pc.famcol(f), label=f) for f in _fams]
pc.color_legend(axB.legend(handles=_lh, fontsize=6.0, frameon=False, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.06), columnspacing=1.0, handlelength=1.0), _lh)
pc.pbadge(axB, "B", f"{STRAIN.replace('_','/')} {BAMC}:{ps:,}–{pe:,} — per-timepoint sRNA coverage (height = expression, colour = timepoint) above TE + gene tracks · {_ntot:,} primary piRNAs · {_p1u:.0f}% 1U", fs=7.4, y=1.07)
axB.text(0.012, 1.02, "‘primary reads’ = each sRNA read (24–32 nt) counted once at its STAR primary locus (multimappers kept, not double-counted) · architecture (genomic strand) ≠ sense/antisense (relative to TE)", transform=axB.transAxes, fontsize=5.2, color="#8a8a8a", style="italic", ha="left", va="center")
# ---------- Panel C: base resolution at true genomic coordinates ----------
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
pc.pbadge(axC, "C", f"Base resolution, {STRAIN.replace('_','/')} at {TPLAB[CHOSEN]} (top tp)   ·   5′-U = 1U   ·   5′ arrow RED = antisense-to-TE (silencing), grey = sense-to-TE", fs=7.6)
fig.suptitle(f"{TELAB} → strain-private piRNA SOURCE LOCUS (individual piRNAs, NOT a PICB cluster) — pangenome layout ({STRAIN.replace('_','/')} chr{CH}:{ps:,}-{pe:,})", fontsize=11.5, fontweight="bold", y=0.965)
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/{OUT}.{e}", bbox_inches="tight")
print(f"   wrote {OUT}.png")
