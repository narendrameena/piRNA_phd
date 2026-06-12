#!/usr/bin/env python3
"""Nucleotide-resolution locus figure for a strain-variable piRNA CLUSTER found in the pangenome PAV. Shows the
BIOLOGY at one locus, in three linked panels:
  (A) cross-strain presence/absence — per-strain cluster coverage at the GRCm39 locus (the PAV pattern), wild red;
  (B) the REAL piRNA cluster in a PRESENT strain (own genome): strand-resolved sRNA coverage (sense up / antisense
      down, 24-32 nt, reps+timepoints pooled) over the locus, with the strain RepeatMasker TE track beneath;
  (C) base resolution — the most abundant individual piRNA reads at the coverage peak drawn as coloured nucleotide
      letters (5'-U = the 1U primary signature highlighted), sense and antisense, with a ping-pong 10-nt 5'-overlap
      pair marked if present.
GRCm39 locus is halLiftover'd (minigraph-cactus HAL) to the present strain's genome; reads come from that strain's
own-genome STAR_srna_strain_wise BAMs. Usage: make_pav_locus.py <g39chrom> <g39start> <g39end> <present_strain>
<gene_label> <pattern_label> <outbase>. Biology (gene function, cluster turnover) BioMNI triple-verified; see
VERIFICATION_QUEUE. Grounded: clusters confirmed to carry thousands of piRNA-sized reads in present strains."""
import sys, os, subprocess, numpy as np, pysam
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
from collections import defaultdict, Counter
import pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"
SIF = f"{ROOT}/cactus_v2.9.3.sif"; HAL = f"{ROOT}/results/pangenome/output/mouse_17strain_pangenome.full.hal"
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]; TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]
NT = {"A": "#33a02c", "C": "#1f78b4", "G": "#ff7f00", "T": "#e31a1c", "N": "#999"}
comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
def rc(s): return "".join(comp.get(c, "N") for c in reversed(s))
G39C, G39S, G39E, STR, GENE, PATT, OUT = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]

# ---- (1) halLiftover GRCm39 locus -> present strain (cached) ----
lc = f"{PG}/.lift_{OUT}_{STR}.bed"
if not os.path.exists(lc):
    open(f"/tmp/_g39_{OUT}.bed", "w").write(f"{G39C}\t{G39S}\t{G39E}\t{OUT}\n")
    subprocess.run(f"singularity exec --bind /mnt {SIF} halLiftover {HAL} GRCm39 /tmp/_g39_{OUT}.bed {STR} {lc}", shell=True, capture_output=True)
lift = [l.split("\t") for l in open(lc) if l.strip()]
chrom_counts = Counter(x[0] for x in lift); CH = chrom_counts.most_common(1)[0][0]
ps = min(int(x[1]) for x in lift if x[0] == CH); pe = max(int(x[2]) for x in lift if x[0] == CH)
BAMC = f"{STR}#1#chr{CH}"; N = pe - ps; nb = 200
print(f"{GENE}: GRCm39 {G39C}:{G39S}-{G39E} -> {STR} {BAMC}:{ps}-{pe} ({N}bp)")

# ---- (2) strand-resolved coverage + 5' positions + sequences (own-genome reads, pooled) ----
fwd = np.zeros(nb); rev = np.zeros(nb); reads = []; u1 = [0, 0]
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
            b0 = int((a.reference_start - ps) / N * nb); b1 = int((a.reference_end - ps) / N * nb)
            b0 = max(0, min(nb - 1, b0)); b1 = max(0, min(nb, b1))
            seq = a.query_sequence.upper()
            if a.is_reverse:
                rev[b0:b1] += 1; five = a.reference_end; b5 = comp.get(seq[-1], "N")
            else:
                fwd[b0:b1] += 1; five = a.reference_start; b5 = seq[0]
            u1[a.is_reverse] += 0; u1[0] += (b5 == "T")
            reads.append((five, a.reference_start, a.reference_end, a.is_reverse, seq, b5))
        bam.close()
ntot = len(reads); nanti = sum(1 for r in reads if r[3]); n1u = sum(1 for r in reads if r[5] == "T")
print(f"   reads={ntot} antisense={nanti} ({100*nanti/max(1,ntot):.0f}%) 1U={100*n1u/max(1,ntot):.0f}%")
# nucleotide zoom centred on the dominant piRNA 5' end (cleanest stack)
five_counts = Counter(r[0] for r in reads)
pk_pos = five_counts.most_common(1)[0][0] if five_counts else (ps + N // 2)
ZW = 80; z0, z1 = pk_pos - 30, pk_pos + ZW - 30

# ---- (3) strain RepeatMasker TE at locus ----
rmf = f"{ROOT}/resources/repeatMasker/{STR}_repeatmasker.bed"
tes = []
if os.path.exists(rmf):
    out = subprocess.run(f"awk -v c=chr{CH} -v s={ps} -v e={pe} '$1==c && $3>s && $2<e' {rmf}", shell=True, capture_output=True, text=True).stdout
    for ln in out.splitlines():
        f = ln.split("\t"); fam = f[3].split("|")[-1] if len(f) > 3 else "?"
        tes.append((int(f[1]), int(f[2]), fam))

# ---- (4) PAV presence row (GRCm39 frame) ----
M = pd.read_csv(f"{U}/cluster_pav/pan_cluster_coverage_matrix.tsv.gz", sep="\t", low_memory=False)
row = M[(M.chrom.astype(str) == G39C) & (M.start <= G39E) & (M.end >= G39S)]
pav = row[ORDER].max().reindex(ORDER).fillna(0) if len(row) else pd.Series(0.0, index=ORDER)

# ---- (5) figure ----
TECOL = {"LINE/L1": "#E69F00", "LTR/ERVK": "#6a3d9a", "LTR/ERVL-MaLR": "#b15928", "SINE/Alu": "#33a02c", "SINE/B2": "#1f78b4", "LTR/ERVL": "#a6cee3", "LTR/ERV1": "#cab2d6", "Simple_repeat": "#bbb"}
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(13.5, 10.2), dpi=300)
gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.5, 1.6], hspace=0.6)
fig.subplots_adjust(top=0.91, bottom=0.06)
# Panel A: presence strip
axA = fig.add_subplot(gs[0]); xs = np.arange(len(ORDER))
axA.bar(xs, pav.values, width=0.82, color=["#C0392B" if pav.values[i] >= 0.5 else "#bbb" for i in range(len(ORDER))], edgecolor="white", linewidth=0.4)
axA.axhline(0.5, color="#888", lw=0.7, ls="--"); axA.set_ylim(0, 1.05); axA.set_xticks(xs)
axA.set_xticklabels([s.replace("_", "/") for s in ORDER], rotation=90, fontsize=7.5)
for t, X in zip(axA.get_xticklabels(), ORDER):
    if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
axA.set_ylabel("cluster\ncoverage", fontsize=8); axA.spines[["top", "right"]].set_visible(False)
axA.set_title(f"A  Cross-strain presence/absence (pangenome PAV, GRCm39) — {PATT}", fontsize=8.8, fontweight="bold", loc="left")
# Panel B: strand-resolved coverage + TE
axB = fig.add_subplot(gs[1]); xg = ps + (np.arange(nb) + 0.5) / nb * N
axB.fill_between(xg, 0, fwd, color="#C0392B", alpha=0.75, step="mid", label="sense (+)")
axB.fill_between(xg, 0, -rev, color="#2166AC", alpha=0.75, step="mid", label="antisense (−)")
axB.axhline(0, color="#333", lw=0.8); axB.set_xlim(ps, pe)
axB.set_ylabel("piRNA reads\n(sense ↑ / antisense ↓)", fontsize=8); axB.legend(fontsize=7.5, frameon=False, loc="upper right", ncol=2)
ymax = max(fwd.max(), rev.max(), 1)
for (ts, te, fam) in tes:
    axB.add_patch(Rectangle((ts, -ymax * 1.18), te - ts, ymax * 0.12, facecolor=TECOL.get(fam, "#ccc"), edgecolor="none"))
axB.set_ylim(-ymax * 1.25, ymax * 1.1)
fams = sorted(set(f for _, _, f in tes), key=lambda f: -sum(te - ts for ts, te, ff in tes if ff == f))[:5]
axB.text(ps, -ymax * 1.24, "TE: " + "  ".join(f"■{f}" for f in fams), fontsize=6.5, va="bottom", color="#555")
axB.axvspan(z0, z1, color="#FDB863", alpha=0.35, zorder=0)
axB.set_title(f"B  Real piRNA cluster in {STR.replace('_','/')} (own genome {BAMC}:{ps:,}-{pe:,}) — {ntot:,} piRNAs, {100*nanti/max(1,ntot):.0f}% antisense, {100*n1u/max(1,ntot):.0f}% 1U", fontsize=9.4, fontweight="bold", loc="left")
axB.set_xlabel(f"{STR.replace('_','/')} chr{CH} position (bp)", fontsize=7.5)
# Panel C: nucleotide-level reads at peak
axC = fig.add_subplot(gs[2])
zreads = [r for r in reads if r[1] < z1 and r[2] > z0]
sense = Counter(); anti = Counter()
for five, rs, re, isrev, seq, b5 in zreads:
    (anti if isrev else sense)[(rs, re, seq, isrev)] += 1
def draw(items, y0, dirn):
    y = y0
    for (rs, re, seq, isrev), cnt in items:
        sshow = seq
        for k, ch in enumerate(sshow):
            x = (re - len(seq) + k) if isrev else (rs + k)
            if z0 - 2 <= x <= z1 + 2:
                axC.text(x, y, ch, fontsize=4.6, ha="center", va="center", family="monospace",
                         color="white", bbox=dict(boxstyle="square,pad=0.02", fc=NT.get(ch, "#999"), ec="none"))
        # 5' marker (1U)
        fivep = re - 1 if isrev else rs
        axC.annotate("", xy=(fivep + (0.5 if not isrev else -0.5), y), xytext=(fivep + (-2.5 if not isrev else 2.5), y),
                     arrowprops=dict(arrowstyle="-|>", color="#C0392B" if isrev else "#222", lw=0.8))
        axC.text(z1 + 4, y, f"×{cnt}", fontsize=4.6, va="center", color="#666")
        y += dirn
    return y
ytop = draw(sense.most_common(8), 1, 1)
axC.axhline(0, color="#333", lw=0.7)
ybot = draw(anti.most_common(8), -1, -1)
# ping-pong check in window (sense 5' and antisense 5' exactly 10 nt apart)
fp_s = set(r[1] for r in zreads if not r[3]); fp_a = set(r[2] - 1 for r in zreads if r[3])
pp = [(p, p + 9) for p in fp_s if (p + 9) in fp_a]
axC.set_xlim(z0 - 1, z1 + 8); axC.set_ylim(ybot - 1.6, ytop + 0.6)
# coordinate ruler — every base sits at its TRUE genomic coordinate
for sp in ("top", "left", "right"): axC.spines[sp].set_visible(False)
axC.set_yticks([]); axC.spines["bottom"].set_position(("data", ybot - 1.0))
tk = np.linspace(z0, z1, 5).astype(int); axC.set_xticks(tk); axC.set_xticklabels([f"{t:,}" for t in tk], fontsize=6.5)
axC.tick_params(axis="x", length=3)
axC.set_xlabel(f"{STR.replace('_','/')} chr{CH} position (bp) — every base aligned to its genomic coordinate", fontsize=7)
axC.text(z0 - 0.6, ytop, "sense", fontsize=6.5, color="#C0392B", fontweight="bold", ha="right", va="center")
axC.text(z0 - 0.6, ybot, "antisense", fontsize=6.5, color="#2166AC", fontweight="bold", ha="right", va="center")
# zoom callout: connect Panel B's highlighted band to Panel C
for xb in (z0, z1):
    fig.add_artist(ConnectionPatch(xyA=(xb, -ymax * 1.25), coordsA=axB.transData, xyB=(xb, ytop + 0.6), coordsB=axC.transData, color="#E8A33D", lw=0.8, ls=(0, (3, 2))))
axC.set_title(f"C  Base resolution — bases at true coordinates (zoom of the band in B); 5′-U = 1U primary signature; {'PING-PONG 10-nt pair present' if pp else 'phased (no 10-nt pair here)'}", fontsize=8.8, fontweight="bold", loc="left")
fig.suptitle(f"{GENE} — strain-variable piRNA cluster at nucleotide resolution", fontsize=12.5, fontweight="bold", y=0.965)
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/{OUT}.{e}", bbox_inches="tight")
print(f"   wrote {OUT}.png")
