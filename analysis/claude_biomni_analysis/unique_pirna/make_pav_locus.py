#!/usr/bin/env python3
"""Nucleotide-resolution figure for a strain-variable piRNA CLUSTER from the pangenome PAV. Three linked panels:
  (A) cross-strain presence/absence — per-strain cluster coverage at the GRCm39 locus (the PAV pattern), wild red;
  (B) the REAL piRNA cluster in a PRESENT strain (own genome): coverage by GENOMIC STRAND (plus up / minus down,
      24-32 nt, reps+timepoints pooled) — this describes cluster ARCHITECTURE (uni- vs dual-strand) — with the
      strain RepeatMasker TE track + legend beneath;
  (C) base resolution at the dominant-piRNA peak: individual piRNAs as coloured nucleotide letters on a genomic
      coordinate ruler (every base at its true coordinate), 5'-U (1U) = primary signature.
IMPORTANT terminology: genomic + / - strand is NOT sense/antisense. SENSE/ANTISENSE is defined RELATIVE TO THE TE
(from the RepeatMasker .out strand): a piRNA is ANTISENSE-to-TE (silencing-competent) if it lies on the opposite
strand to the TE it overlaps. Panels show genomic strand (architecture); the antisense-to-TE fraction (silencing)
is reported separately and the 5' arrows in C are coloured by it. GRCm39 locus halLiftover'd (cactus HAL) to the
present strain; reads from own-genome STAR_srna_strain_wise BAMs. Usage: make_pav_locus.py <g39chrom> <g39start>
<g39end> <present_strain> <gene_label> <pattern_label> <outbase>. Biology BioMNI triple-verified (VERIFICATION_QUEUE)."""
import sys, os, glob, subprocess, numpy as np, pysam
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
from collections import Counter
import pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch, Patch
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"
SIF = f"{ROOT}/cactus_v2.9.3.sif"; HAL = f"{ROOT}/results/pangenome/output/mouse_17strain_pangenome.full.hal"
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]; TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]
NT = {"A": "#33a02c", "C": "#1f78b4", "G": "#ff7f00", "T": "#e31a1c", "N": "#999"}
comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
TECOL = {"LINE/L1": "#E69F00", "LTR/ERVK": "#6a3d9a", "LTR/ERVL-MaLR": "#b15928", "SINE/Alu": "#33a02c", "SINE/B2": "#1f78b4", "LTR/ERVL": "#a6cee3", "LTR/ERV1": "#cab2d6", "Simple_repeat": "#bbbbbb", "LTR/ERVL?": "#a6cee3"}
G39C, G39S, G39E, STR, GENE, PATT, OUT = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]

# ---- (1) halLiftover GRCm39 locus -> present strain (cached) ----
lc = f"{PG}/.lift_{OUT}_{STR}.bed"
if not os.path.exists(lc):
    open(f"/tmp/_g39_{OUT}.bed", "w").write(f"{G39C}\t{G39S}\t{G39E}\t{OUT}\n")
    subprocess.run(f"singularity exec --bind /mnt {SIF} halLiftover {HAL} GRCm39 /tmp/_g39_{OUT}.bed {STR} {lc}", shell=True, capture_output=True)
lift = [l.split("\t") for l in open(lc) if l.strip()]
CH = Counter(x[0] for x in lift).most_common(1)[0][0]
ps = min(int(x[1]) for x in lift if x[0] == CH); pe = max(int(x[2]) for x in lift if x[0] == CH)
BAMC = f"{STR}#1#chr{CH}"; N = pe - ps; nb = 200
print(f"{GENE}: GRCm39 {G39C}:{G39S}-{G39E} -> {STR} {BAMC}:{ps}-{pe} ({N}bp)")

# ---- (2) RepeatMasker TE at locus WITH STRAND (from .out) -> for orientation-relative-to-TE ----
outf = glob.glob(f"{ROOT}/resources/repeatMasker/{STR}_*.out")
tes = []  # (start, end, te_strand, family)
if outf:
    aw = subprocess.run("awk -v c=%s -v s=%d -v e=%d 'NR>3 && $5==c && $7>s && $6<e{st=($9==\"C\")?\"-\":\"+\"; print $6\"\\t\"$7\"\\t\"st\"\\t\"$11}' %s" % (BAMC, ps, pe, outf[0]), shell=True, capture_output=True, text=True).stdout
    for ln in aw.splitlines():
        f = ln.split("\t"); tes.append((int(f[0]), int(f[1]), f[2], f[3]))
tes.sort()
def te_strand_at(pos):
    for ts, te, st, fam in tes:
        if ts <= pos < te: return st
    return None

# ---- (3) own-genome reads: genomic-strand coverage + 5' + seq + orientation-relative-to-TE ----
plus = np.zeros(nb); minus = np.zeros(nb); reads = []
for tp in TPS:
    for rep in (1, 2, 3):
        b = f"{ROOT}/results/STAR_srna_strain_wise/{STR}/{STR}-{tp}.{rep}/Aligned.sortedByCoord.out.bam"
        if not os.path.exists(b): continue
        bam = pysam.AlignmentFile(b, "rb")
        try: it = bam.fetch(BAMC, ps, pe)
        except Exception: continue
        for a in it:
            if a.is_unmapped or not a.query_sequence: continue
            L = a.reference_end - a.reference_start
            if not (24 <= L <= 32): continue
            b0 = max(0, min(nb - 1, int((a.reference_start - ps) / N * nb))); b1 = max(0, min(nb, int((a.reference_end - ps) / N * nb)))
            seq = a.query_sequence.upper()
            if a.is_reverse: minus[b0:b1] += 1; b5 = comp.get(seq[-1], "N")
            else: plus[b0:b1] += 1; b5 = seq[0]
            reads.append((a.reference_start, a.reference_end, a.is_reverse, seq, b5))
        bam.close()
ntot = len(reads); nminus = sum(1 for r in reads if r[2]); n1u = sum(1 for r in reads if r[4] == "T")
pct_minus = 100 * nminus / max(1, ntot); minor = min(ntot - nminus, nminus) / max(1, ntot)
arch = "dual-strand (bidirectional)" if minor > 0.2 else "uni-strand"
# CORRECT sense/antisense = relative to TE: piRNA antisense-to-TE if on opposite strand to its TE
nte = nat = 0
def anti_te(rs, re, isrev):
    st = te_strand_at((rs + re) // 2)
    return None if st is None else ((st == "-") != isrev)
for rs, re, isrev, seq, b5 in reads:
    a = anti_te(rs, re, isrev)
    if a is not None: nte += 1; nat += a
pct_at = 100 * nat / max(1, nte)
print(f"   reads={ntot} plus/minus={100-pct_minus:.0f}/{pct_minus:.0f} ({arch}); 1U={100*n1u/max(1,ntot):.0f}%; antisense-to-TE={pct_at:.0f}% (of {nte})")
# zoom centred on dominant piRNA 5' end
five = Counter((r[1] - 1 if r[2] else r[0]) for r in reads)
pk = five.most_common(1)[0][0] if five else ps + N // 2
z0, z1 = pk - 30, pk + 50

# ---- (4) PAV presence row ----
M = pd.read_csv(f"{U}/cluster_pav/pan_cluster_coverage_matrix.tsv.gz", sep="\t", low_memory=False)
rowm = M[(M.chrom.astype(str) == G39C) & (M.start <= G39E) & (M.end >= G39S)]
pav = rowm[ORDER].max().reindex(ORDER).fillna(0) if len(rowm) else pd.Series(0.0, index=ORDER)

# ---- (5) figure ----
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(13.5, 10.2), dpi=300)
gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.5, 1.6], hspace=0.62); fig.subplots_adjust(top=0.91, bottom=0.06)
axA, axB, axC = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
# A: presence/absence
xs = np.arange(len(ORDER))
axA.bar(xs, pav.values, width=0.82, color=["#C0392B" if pav.values[i] >= 0.5 else "#bbb" for i in range(len(ORDER))], edgecolor="white", linewidth=0.4)
axA.axhline(0.5, color="#888", lw=0.7, ls="--"); axA.set_ylim(0, 1.05); axA.set_xticks(xs)
axA.set_xticklabels([s.replace("_", "/") for s in ORDER], rotation=90, fontsize=7.5)
for t, X in zip(axA.get_xticklabels(), ORDER):
    if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
axA.set_ylabel("cluster\ncoverage", fontsize=8); axA.spines[["top", "right"]].set_visible(False)
axA.set_title(f"A  Cross-strain presence/absence (pangenome PAV, GRCm39) — {PATT}", fontsize=8.8, fontweight="bold", loc="left")
# B: genomic-strand coverage + TE track
xg = ps + (np.arange(nb) + 0.5) / nb * N; ymax = max(plus.max(), minus.max(), 1)
axB.fill_between(xg, 0, plus, color="#33a02c", alpha=0.8, step="mid", label="plus (+) strand")
axB.fill_between(xg, 0, -minus, color="#6a3d9a", alpha=0.8, step="mid", label="minus (−) strand")
axB.axhline(0, color="#333", lw=0.8); axB.set_xlim(ps, pe)
axB.set_ylabel("piRNA reads\n(genomic strand)", fontsize=8)
for (ts, te, st, fam) in tes:
    axB.add_patch(Rectangle((ts, -ymax * 1.2), te - ts, ymax * 0.12, facecolor=TECOL.get(fam, "#dddddd"), edgecolor="none"))
    axB.annotate("", xy=(te if st == "+" else ts, -ymax * 1.14), xytext=(ts if st == "+" else te, -ymax * 1.14), arrowprops=dict(arrowstyle="-|>", color="#444", lw=0.5))
axB.set_ylim(-ymax * 1.28, ymax * 1.12)
famset = sorted({f for _, _, _, f in tes}, key=lambda f: -sum(te - ts for ts, te, _, ff in tes if ff == f))[:6]
axB.legend(handles=[Patch(facecolor="#33a02c", label="+ strand"), Patch(facecolor="#6a3d9a", label="− strand")] + [Patch(facecolor=TECOL.get(f, "#ddd"), label=f) for f in famset],
           fontsize=6.4, frameon=False, loc="lower left", ncol=4, bbox_to_anchor=(0, -0.32), columnspacing=1.2, handlelength=1.1)
axB.axvspan(z0, z1, color="#FDB863", alpha=0.35, zorder=0)
axB.set_title(f"B  piRNA cluster in {STR.replace('_','/')} ({BAMC}:{ps:,}-{pe:,}) — {ntot:,} piRNAs · {arch} (+{100-pct_minus:.0f}%/−{pct_minus:.0f}%) · {pct_at:.0f}% ANTISENSE-to-TE (silencing) · {100*n1u/max(1,ntot):.0f}% 1U",
             fontsize=8.6, fontweight="bold", loc="left")
axB.set_xlabel(f"{STR.replace('_','/')} chr{CH} position (bp) · TE track below (arrow = TE orientation)", fontsize=7.5)
# C: base resolution
zreads = [r for r in reads if r[0] < z1 and r[1] > z0]
plus_r, minus_r = Counter(), Counter()
for rs, re, isrev, seq, b5 in zreads:
    (minus_r if isrev else plus_r)[(rs, re, seq, isrev)] += 1
def draw(items, y0, dirn):
    y = y0
    for (rs, re, seq, isrev), cnt in items:
        for k, ch in enumerate(seq):
            x = rs + k
            if z0 - 2 <= x <= z1 + 2:
                axC.text(x, y, ch, fontsize=4.6, ha="center", va="center", family="monospace", color="white", bbox=dict(boxstyle="square,pad=0.02", fc=NT.get(ch, "#999"), ec="none"))
        at = anti_te(rs, re, isrev); acol = "#C0392B" if at else ("#888" if at is not None else "#ccc")
        fp = re - 1 if isrev else rs
        axC.annotate("", xy=(fp + (0.5 if not isrev else -0.5), y), xytext=(fp + (-2.6 if not isrev else 2.6), y), arrowprops=dict(arrowstyle="-|>", color=acol, lw=1.0))
        axC.text(z1 + 4, y, f"×{cnt}", fontsize=4.6, va="center", color="#666")
        y += dirn
    return y
ytop = draw(plus_r.most_common(8), 1, 1); axC.axhline(0, color="#333", lw=0.7); ybot = draw(minus_r.most_common(8), -1, -1)
axC.set_xlim(z0 - 1, z1 + 8); axC.set_ylim(ybot - 1.6, ytop + 0.6)
for sp in ("top", "left", "right"): axC.spines[sp].set_visible(False)
axC.set_yticks([]); axC.spines["bottom"].set_position(("data", ybot - 1.0))
tk = np.linspace(z0, z1, 5).astype(int); axC.set_xticks(tk); axC.set_xticklabels([f"{t:,}" for t in tk], fontsize=6.5); axC.tick_params(axis="x", length=3)
axC.set_xlabel(f"{STR.replace('_','/')} chr{CH} position (bp) — every base at its true genomic coordinate", fontsize=7)
axC.text(z0 - 0.6, ytop, "+ strand", fontsize=6.5, color="#33a02c", fontweight="bold", ha="right", va="center")
axC.text(z0 - 0.6, ybot, "− strand", fontsize=6.5, color="#6a3d9a", fontweight="bold", ha="right", va="center")
for xb in (z0, z1):
    fig.add_artist(ConnectionPatch(xyA=(xb, -ymax * 1.28), coordsA=axB.transData, xyB=(xb, ytop + 0.6), coordsB=axC.transData, color="#E8A33D", lw=0.8, ls=(0, (3, 2))))
axC.set_title(f"C  Base resolution (zoom of band in B); bases at true coordinates; 5′-U = 1U signature; 5′ arrow RED = ANTISENSE-to-TE (silencing), grey = sense-to-TE", fontsize=8.2, fontweight="bold", loc="left")
fig.suptitle(f"{GENE} — strain-variable piRNA cluster at nucleotide resolution", fontsize=12.5, fontweight="bold", y=0.965)
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/{OUT}.{e}", bbox_inches="tight")
print(f"   wrote {OUT}.png")
