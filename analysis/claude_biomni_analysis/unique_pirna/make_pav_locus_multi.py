#!/usr/bin/env python3
"""MULTI-strain locus plot for a PICB piRNA cluster: pangenome comparison + per-strain coverage + nucleotide detail.
  (A) PANGENOME cross-strain × timepoint comparison — PICB-combined cluster FPM at this locus for every strain ×
      timepoint, from picb_pangenome_fpm.tsv (all strains' PICB clusters projected once into the GRCm39 frame via
      the cactus HAL — NOT pairwise locus liftover). Grouped bars, FPM atop each, wild red, data-driven label.
  (B) coverage in EVERY PRESENT (PICB-expressing) strain: genomic-strand read coverage (plus ↑ green / minus ↓
      purple) on a shared relative axis; each track labelled with PICB FPM, architecture, and % ANTISENSE-to-TE.
  (C) base resolution in the TOP-expressing strain — piRNAs as coloured bases on a true-coordinate ruler.
TERMINOLOGY: genomic +/- strand != sense/antisense; SENSE/ANTISENSE is relative to the TE (RM .out strand),
antisense-to-TE = silencing-competent. PICB COMBINED run; reads from own-genome STAR_srna_strain_wise BAMs; each
present strain's locus = GRCm39 locus back-lifted via the HAL (cached). Does NOT use the old coverage PAV.
Usage: make_pav_locus.py <g39chrom> <g39start> <g39end> <gene_label> <unused> <outbase>."""
import sys, os, glob, subprocess, numpy as np, pysam
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
from collections import Counter
import pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch, Patch
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"; CP = f"{U}/cluster_pav"
SIF = f"{ROOT}/cactus_v2.9.3.sif"; HAL = f"{ROOT}/results/pangenome/output/mouse_17strain_pangenome.full.hal"
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]; TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]
TPLAB = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}; TPCOL = {"16.5dpc": "#4393C3", "12.5dpp": "#FDB863", "20.5dpp": "#B2182B"}
NT = {"A": "#33a02c", "C": "#1f78b4", "G": "#ff7f00", "T": "#e31a1c", "N": "#999"}; comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
TECOL = {"LINE/L1": "#E69F00", "LTR/ERVK": "#6a3d9a", "LTR/ERVL-MaLR": "#b15928", "SINE/Alu": "#33a02c", "SINE/B2": "#1f78b4", "LTR/ERVL": "#a6cee3", "LTR/ERV1": "#cab2d6", "Simple_repeat": "#bbbbbb"}
G39C, G39S, G39E, GENE, OUT = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[6]
nb = 180
# ---- (1) pangenome FPM (cross-strain x timepoint) ----
P = pd.read_csv(f"{CP}/picb_pangenome_fpm.tsv", sep="\t", dtype={"g39_chrom": str})
sub = P[(P.g39_chrom == G39C) & (P.start < G39E) & (P.end > G39S)]
FPM = sub.groupby(["strain", "tp"])["all_primary_FPM"].max().unstack(fill_value=0.0).reindex(index=ORDER, columns=TPS).fillna(0.0)
present = [X for X in ORDER if FPM.loc[X].max() > 0]   # CANONICAL order (same as Panel A)
TOP = present[-1] if present else max(ORDER, key=lambda X: FPM.loc[X].max())   # last/bottom strain = Panel C example
pattern_auto = f"PICB cluster present in {len(present)}/16 strains: " + ", ".join(s.replace("_", "/") for s in present)
print(f"{GENE}: present={present}; top={TOP}")
# ---- (2) per-present-strain coverage (own genome; back-lift cached) ----
def bams_for(X): return [b for tp in TPS for r in (1, 2, 3) for b in [f"{ROOT}/results/STAR_srna_strain_wise/{X}/{X}-{tp}.{r}/Aligned.sortedByCoord.out.bam"] if os.path.exists(b)]
def lift(X):
    lc = f"{PG}/.lift_{G39C}_{G39S}_{X}.bed"
    if not os.path.exists(lc):
        open(f"/tmp/_g39_{OUT}_{X}.bed", "w").write(f"{G39C}\t{G39S}\t{G39E}\tx\n")
        subprocess.run(f"singularity exec --bind /mnt {SIF} halLiftover {HAL} GRCm39 /tmp/_g39_{OUT}_{X}.bed {X} {lc}", shell=True, capture_output=True)
    rows = [l.split("\t") for l in open(lc) if l.strip()]
    if not rows: return None
    ch = Counter(r[0] for r in rows).most_common(1)[0][0]
    return ch, min(int(r[1]) for r in rows if r[0] == ch), max(int(r[2]) for r in rows if r[0] == ch)
def te_at(X, BAMC, ps, pe):
    outf = glob.glob(f"{ROOT}/resources/repeatMasker/{X}_*.out"); tes = []
    if outf:
        aw = subprocess.run("awk -v c=%s -v s=%d -v e=%d 'NR>3 && $5==c && $7>s && $6<e{st=($9==\"C\")?\"-\":\"+\"; print $6\"\\t\"$7\"\\t\"st\"\\t\"$11}' %s" % (BAMC, ps, pe, outf[0]), shell=True, capture_output=True, text=True).stdout
        for ln in aw.splitlines():
            f = ln.split("\t"); tes.append((int(f[0]), int(f[1]), f[2], f[3]))
    return sorted(tes)
def collect(X):
    lf = lift(X)
    if lf is None: return None
    ch, ps, pe = lf; BAMC = f"{X}#1#chr{ch}"; N = max(1, pe - ps); tes = te_at(X, BAMC, ps, pe)
    def tst(p):
        for ts, te, s, fam in tes:
            if ts <= p < te: return s
        return None
    plus = np.zeros(nb); minus = np.zeros(nb); reads = []; nat = nte = n1u = 0
    for b in bams_for(X):
        try:
            bam = pysam.AlignmentFile(b, "rb")
            for a in bam.fetch(BAMC, ps, pe):
                if a.is_unmapped or not a.query_sequence: continue
                L = a.reference_end - a.reference_start
                if not (24 <= L <= 32): continue
                b0 = max(0, min(nb - 1, int((a.reference_start - ps) / N * nb))); b1 = max(0, min(nb, int((a.reference_end - ps) / N * nb)))
                seq = a.query_sequence.upper(); (minus if a.is_reverse else plus)[b0:b1] += 1
                reads.append((a.reference_start, a.reference_end, a.is_reverse, seq)); n1u += ((comp.get(seq[-1], "N") if a.is_reverse else seq[0]) == "T")
                s = tst((a.reference_start + a.reference_end) // 2)
                if s is not None: nte += 1; nat += ((s == "-") != a.is_reverse)
            bam.close()
        except (OSError, ValueError):
            try: bam.close()
            except Exception: pass
    ntot = len(reads); nminus = sum(1 for r in reads if r[2])
    return dict(ch=ch, ps=ps, pe=pe, BAMC=BAMC, N=N, tes=tes, tst=tst, plus=plus, minus=minus, reads=reads,
                ntot=ntot, nminus=nminus, n1u=n1u, nat=nat, nte=nte, pct_at=100 * nat / max(1, nte),
                arch="dual-strand" if min(nminus, ntot - nminus) / max(1, ntot) > 0.2 else "uni-strand")
COV = {X: collect(X) for X in present}; COV = {X: d for X, d in COV.items() if d}
present = [X for X in present if X in COV]; TOP = present[-1] if present else TOP
# ---- (3) figure ----
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
nP = max(1, len(present)); fig = plt.figure(figsize=(14, 7.2 + 0.95 * nP), dpi=300)
gs = fig.add_gridspec(3, 1, height_ratios=[0.95, 0.9 * nP + 1.0, 2.1], hspace=0.78); fig.subplots_adjust(top=0.85, bottom=0.075)
axA, axB, axC = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
x = np.arange(len(ORDER)); bw = 0.26
for j, tp in enumerate(TPS):
    h = FPM[tp].values; axA.bar(x + (j - 1) * bw, np.maximum(h, 1e-3), width=bw, color=TPCOL[tp], edgecolor="white", linewidth=0.2, label=TPLAB[tp])
    for xi in range(len(ORDER)):
        if h[xi] > 0: axA.text(xi + (j - 1) * bw, h[xi] * 1.18, f"{h[xi]:.0f}" if h[xi] >= 1 else f"{h[xi]:.1f}", ha="center", va="bottom", fontsize=4.8, rotation=90, color=TPCOL[tp], fontweight="bold")
axA.set_yscale("log"); axA.set_ylim(0.1, max(FPM.values.max(), 1) * 4.5); axA.set_xticks(x); axA.set_xticklabels([s.replace("_", "/") for s in ORDER], rotation=90, fontsize=7.5)
for t, X in zip(axA.get_xticklabels(), ORDER):
    if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
axA.set_ylabel("PICB cluster\nFPM (log)", fontsize=8); axA.spines[["top", "right"]].set_visible(False)
axA.legend([TPLAB[t] for t in TPS], fontsize=7.5, frameon=False, ncol=1, loc="upper left", bbox_to_anchor=(1.005, 1.0), title="timepoint", title_fontsize=7, handlelength=1.2)
axA.set_title("A  Pangenome cross-strain × timepoint comparison — PICB-combined cluster FPM (projected to GRCm39)", fontsize=8.8, fontweight="bold", loc="left")
# zoom window from TOP strain (for nucleotide panel C + the zoom bar in B)
dT = COV[TOP]; psT = dT["ps"]; NT_ = dT["N"]
fiveT = Counter((r[1] - 1 if r[2] else r[0]) for r in dT["reads"]); pk = fiveT.most_common(1)[0][0] if fiveT else psT + NT_ // 2; z0, z1 = pk - 30, pk + 50
# B: per-present-strain coverage + per-strain TE track; EXAMPLE (zoomed) strain at the BOTTOM (next to C)
SP = 3.0; off_top = -(nP - 1) * SP   # canonical order top->bottom; last (TOP) at the bottom, next to C
for i, X in enumerate(present):
    d = COV[X]; off = -i * SP; xr = np.linspace(0, 1, nb); pm = max(d["plus"].max(), d["minus"].max(), 1)
    if X == TOP:                                          # zoom HIGHLIGHTER bar on the example (bottom) track
        zr0 = (z0 - d["ps"]) / d["N"]; zr1 = (z1 - d["ps"]) / d["N"]
        axB.add_patch(Rectangle((zr0, off - 1.65), zr1 - zr0, 2.75, facecolor="#FDB863", alpha=0.32, edgecolor="#E8A33D", lw=0.5, zorder=0))
    axB.fill_between(xr, off, off + d["plus"] / pm, color="#33a02c", alpha=0.85, step="mid"); axB.fill_between(xr, off, off - d["minus"] / pm, color="#6a3d9a", alpha=0.85, step="mid")
    axB.axhline(off, color="#888", lw=0.4)
    for (ts, te, st, fam) in d["tes"]:                    # this strain's TE annotation track + STRAND arrow
        x0 = max(0.0, (ts - d["ps"]) / d["N"]); x1 = min(1.0, (te - d["ps"]) / d["N"])
        if x1 > x0:
            axB.add_patch(Rectangle((x0, off - 1.6), x1 - x0, 0.32, facecolor=TECOL.get(fam, "#dddddd"), edgecolor="none"))
            if x1 - x0 > 0.012:                            # TE strand: arrow points 5'->3' of the TE
                af, at = (x0, x1) if st == "+" else (x1, x0)
                axB.annotate("", xy=(at, off - 1.44), xytext=(af, off - 1.44), arrowprops=dict(arrowstyle="-|>", color="#111", lw=0.4, mutation_scale=4.5))
    axB.text(-0.015, off + 0.35, X.replace("_", "/") + ("  ▼zoom" if X == TOP else ""), fontsize=7.0, ha="right", va="center", fontweight="bold", color="#C0392B" if X in WILD else "#222")
    axB.text(-0.015, off - 0.78, f"chr{d['ch']}:{d['ps']:,}", fontsize=4.6, ha="right", va="center", color="#999")
    axB.text(1.015, off + 0.25, f"FPM {FPM.loc[X].max():.1f} · {d['arch']}", fontsize=6.0, ha="left", va="center", color="#444")
    axB.text(1.015, off - 0.45, f"{d['pct_at']:.0f}% AS→TE · {100*d['n1u']/max(1,d['ntot']):.0f}% 1U · TE↓", fontsize=5.4, ha="left", va="center", color="#777")
axB.set_xlim(0, 1); axB.set_ylim(off_top - 1.8, 1.5); axB.axis("off")
axB.text(0.5, 1.04, "B  Per-PRESENT-strain genomic-strand coverage (plus ↑ green / minus ↓ purple, POOLED over all timepoints — per-tp FPM in A) + each strain's TE track; example strain at bottom → zoom into C", fontsize=7.9, fontweight="bold", ha="center", transform=axB.transAxes)
famset = list(dict.fromkeys(f for X in present for (_, _, _, f) in COV[X]["tes"]))[:6]
fig.legend(handles=[Patch(facecolor="#33a02c", label="+ strand"), Patch(facecolor="#6a3d9a", label="− strand")] + [Patch(facecolor=TECOL.get(f, "#ddd"), label=f) for f in famset], fontsize=6.4, frameon=False, loc="lower center", ncol=8, bbox_to_anchor=(0.5, 0.012), title="genomic strand + TE families", title_fontsize=6.5)
# C: nucleotide in TOP strain
d = COV[TOP]; ps = d["ps"]; reads = d["reads"]; tst = d["tst"]; N = d["N"]; CH = d["ch"]
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
axC.set_xlabel(f"{TOP.replace('_','/')} chr{CH} position (bp) — every base at its true genomic coordinate", fontsize=7)
axC.text(z0 - 0.6, ytop, "+ strand", fontsize=6.5, color="#33a02c", fontweight="bold", ha="right", va="center"); axC.text(z0 - 0.6, ybot, "− strand", fontsize=6.5, color="#6a3d9a", fontweight="bold", ha="right", va="center")
axC.set_title(f"C  Base resolution in example strain {TOP.replace('_','/')} (zoom of the highlighted bar in B); 5′-U = 1U; 5′ arrow RED = ANTISENSE-to-TE (silencing), grey = sense-to-TE", fontsize=7.7, fontweight="bold", loc="left")
for (xbr, xc) in [((z0 - psT) / NT_, z0), ((z1 - psT) / NT_, z1)]:   # zoom-out callout: bottom (example) track -> Panel C
    fig.add_artist(ConnectionPatch(xyA=(xbr, off_top - 1.65), coordsA=axB.transData, xyB=(xc, ytop + 0.6), coordsB=axC.transData, color="#E8A33D", lw=0.9, ls=(0, (3, 2))))
fig.suptitle(f"{GENE} — multi-strain PICB piRNA cluster", fontsize=12.5, fontweight="bold", y=0.985)
fig.text(0.5, 0.95, pattern_auto if len(pattern_auto) < 120 else pattern_auto[:117] + "…", ha="center", fontsize=7.2, color="#555")
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/{OUT}.{e}", bbox_inches="tight")
print(f"   wrote {OUT}.png ({len(present)} present strains)")
