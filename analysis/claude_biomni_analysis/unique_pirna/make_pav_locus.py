#!/usr/bin/env python3
"""Nucleotide-resolution figure for a strain-variable piRNA CLUSTER (pangenome PAV), showing ALL strains.
  (A) per-strain piRNA EXPRESSION at the locus = RPM (piRNA-sized reads at the syntenic locus / library size x1e6,
      library-normalised; ~ PICB FPM) — expression, NOT mere coverage; all 16 strains, wild red, log scale.
  (B) the cluster in EVERY PRESENT strain: genomic-strand coverage (plus up / minus down) on a shared relative-
      position axis — shows the biology is shared across strains; each track labelled with RPM and % ANTISENSE-to-TE.
  (C) base resolution in the top-expressing present strain: piRNAs as coloured bases on a true-coordinate ruler.
TERMINOLOGY (correct): genomic + / - strand != sense/antisense. SENSE/ANTISENSE is RELATIVE TO THE TE (RepeatMasker
.out strand): a piRNA is ANTISENSE-to-TE (silencing-competent) if on the opposite strand to its overlapping TE.
Panels show genomic strand (architecture) AND report sense/antisense-to-TE %. Each strain's locus = GRCm39 locus
halLiftover'd to that strain (cactus HAL); reads from own-genome STAR_srna_strain_wise BAMs.
Usage: make_pav_locus.py <g39chrom> <g39start> <g39end> <gene_label> <pattern_label> <outbase>."""
import sys, os, glob, json, subprocess, numpy as np, pysam
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
from collections import Counter
import pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch, Patch
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"
SIF = f"{ROOT}/cactus_v2.9.3.sif"; HAL = f"{ROOT}/results/pangenome/output/mouse_17strain_pangenome.full.hal"
SAMT = sorted(glob.glob("/mnt/home3/miska/nm667/miniconda3/envs/*/bin/samtools"))[0]
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]; TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]
NT = {"A": "#33a02c", "C": "#1f78b4", "G": "#ff7f00", "T": "#e31a1c", "N": "#999"}
comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
TECOL = {"LINE/L1": "#E69F00", "LTR/ERVK": "#6a3d9a", "LTR/ERVL-MaLR": "#b15928", "SINE/Alu": "#33a02c", "SINE/B2": "#1f78b4", "LTR/ERVL": "#a6cee3", "LTR/ERV1": "#cab2d6", "Simple_repeat": "#bbbbbb"}
G39C, G39S, G39E, GENE, PATT, OUT = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6]

def bams_for(X): return [f"{ROOT}/results/STAR_srna_strain_wise/{X}/{X}-{tp}.{r}/Aligned.sortedByCoord.out.bam" for tp in TPS for r in (1, 2, 3) if os.path.exists(f"{ROOT}/results/STAR_srna_strain_wise/{X}/{X}-{tp}.{r}/Aligned.sortedByCoord.out.bam")]
# ---- library size = total mapped READS (each read once, from STAR Log.final.out), cached ----
LIBF = f"{PG}/.libsize_srna_reads.json"; LIB = json.load(open(LIBF)) if os.path.exists(LIBF) else {}
def libsize(X):
    if X in LIB: return LIB[X]
    tot = 0
    for b in bams_for(X):
        log = b.replace("Aligned.sortedByCoord.out.bam", "Log.final.out")
        if os.path.exists(log):
            for l in open(log):
                if "Uniquely mapped reads number" in l or "Number of reads mapped to multiple loci" in l:
                    tot += int(l.strip().split("\t")[-1])
    LIB[X] = tot; json.dump(LIB, open(LIBF, "w")); return tot
def lift(X):  # GRCm39 locus -> strain X coords (cached); returns (chrom, ps, pe) or None
    lc = f"{PG}/.lift_{OUT}_{X}.bed"
    if not os.path.exists(lc):
        open(f"/tmp/_g39_{OUT}.bed", "w").write(f"{G39C}\t{G39S}\t{G39E}\t{OUT}\n")
        subprocess.run(f"singularity exec --bind /mnt {SIF} halLiftover {HAL} GRCm39 /tmp/_g39_{OUT}.bed {X} {lc}", shell=True, capture_output=True)
    rows = [l.split("\t") for l in open(lc) if l.strip()]
    if not rows: return None
    ch = Counter(r[0] for r in rows).most_common(1)[0][0]
    return ch, min(int(r[1]) for r in rows if r[0] == ch), max(int(r[2]) for r in rows if r[0] == ch)
def te_at(X, BAMC, ps, pe):  # stranded TEs from RM .out
    outf = glob.glob(f"{ROOT}/resources/repeatMasker/{X}_*.out"); tes = []
    if outf:
        aw = subprocess.run("awk -v c=%s -v s=%d -v e=%d 'NR>3 && $5==c && $7>s && $6<e{st=($9==\"C\")?\"-\":\"+\"; print $6\"\\t\"$7\"\\t\"st\"\\t\"$11}' %s" % (BAMC, ps, pe, outf[0]), shell=True, capture_output=True, text=True).stdout
        for ln in aw.splitlines():
            f = ln.split("\t"); tes.append((int(f[0]), int(f[1]), f[2], f[3]))
    return sorted(tes)
def collect(X, nb=160):  # per-strain: RPM, strand coverage, seqs, antisense-to-TE
    lf = lift(X)
    if lf is None: return None
    ch, ps, pe = lf; BAMC = f"{X}#1#chr{ch}"; N = max(1, pe - ps)
    tes = te_at(X, BAMC, ps, pe)
    def testr(pos):
        for ts, te, st, fam in tes:
            if ts <= pos < te: return st
        return None
    plus = np.zeros(nb); minus = np.zeros(nb); reads = []; nat = nte = 0; fcount = 0.0
    for b in bams_for(X):
        try:
            bam = pysam.AlignmentFile(b, "rb")
            for a in bam.fetch(BAMC, ps, pe):
                if a.is_unmapped or not a.query_sequence: continue
                L = a.reference_end - a.reference_start
                if not (24 <= L <= 32): continue
                fcount += 1.0 / (a.get_tag("NH") if a.has_tag("NH") else 1)   # multimapping-fair (1/NH)
                b0 = max(0, min(nb - 1, int((a.reference_start - ps) / N * nb))); b1 = max(0, min(nb, int((a.reference_end - ps) / N * nb)))
                (minus if a.is_reverse else plus)[b0:b1] += 1
                reads.append((a.reference_start, a.reference_end, a.is_reverse, a.query_sequence.upper()))
                st = testr((a.reference_start + a.reference_end) // 2)
                if st is not None: nte += 1; nat += ((st == "-") != a.is_reverse)
            bam.close()
        except (OSError, ValueError):
            try: bam.close()
            except Exception: pass
            continue
    ntot = len(reads); rpm = ntot / max(1, libsize(X)) * 1e6
    return dict(ch=ch, ps=ps, pe=pe, BAMC=BAMC, N=N, plus=plus, minus=minus, reads=reads, tes=tes,
                ntot=ntot, rpm=rpm, nminus=sum(1 for r in reads if r[2]), n1u=sum(1 for r in reads if (comp.get(r[3][-1], "N") if r[2] else r[3][0]) == "T"),
                pct_at=100 * nat / max(1, nte), nte=nte, testr=testr)

print(f"{GENE}: collecting all 16 strains at GRCm39 {G39C}:{G39S}-{G39E} ...")
D = {X: collect(X) for X in ORDER}
rpm = {X: (D[X]["rpm"] if D[X] else 0.0) for X in ORDER}
present = [X for X in ORDER if rpm[X] >= max(2.0, 0.05 * max(rpm.values()))]  # expressed strains
present = sorted(present, key=lambda X: -rpm[X])
print("  present (expressing) strains:", [(X, round(rpm[X], 1)) for X in present])
TOP = present[0] if present else max(ORDER, key=lambda X: rpm[X])

# ---- figure ----
plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
nB = max(1, len(present)); fig = plt.figure(figsize=(13.5, 4.6 + 0.5 * nB), dpi=300)
gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.42 * nB + 0.5, 1.5], hspace=0.6); fig.subplots_adjust(top=0.92, bottom=0.06)
axA, axB, axC = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
# A: per-strain RPM (expression)
xs = np.arange(len(ORDER)); vals = np.array([max(rpm[X], 1e-3) for X in ORDER])
axA.bar(xs, vals, width=0.82, color=["#C0392B" if X in present else "#bbb" for X in ORDER], edgecolor="white", linewidth=0.4)
axA.set_yscale("log"); axA.set_ylim(max(0.1, min(vals) * 0.7), max(vals) * 2)
axA.set_xticks(xs); axA.set_xticklabels([s.replace("_", "/") for s in ORDER], rotation=90, fontsize=7.5)
for t, X in zip(axA.get_xticklabels(), ORDER):
    if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
axA.set_ylabel("piRNA RPM\n(log)", fontsize=8); axA.spines[["top", "right"]].set_visible(False)
axA.set_title(f"A  Per-strain piRNA EXPRESSION at the locus (RPM = piRNA reads/library ×1e6; ≈ PICB FPM) — {PATT}", fontsize=8.7, fontweight="bold", loc="left")
# B: every present strain, strand coverage on shared relative axis
axB.axhline(0, color="#ccc", lw=0.4)
for i, X in enumerate(present):
    dd = D[X]; off = -i * 2.4; xr = np.linspace(0, 1, len(dd["plus"]))
    pmax = max(dd["plus"].max(), dd["minus"].max(), 1)
    axB.fill_between(xr, off, off + dd["plus"] / pmax, color="#33a02c", alpha=0.8, step="mid")
    axB.fill_between(xr, off, off - dd["minus"] / pmax, color="#6a3d9a", alpha=0.8, step="mid")
    axB.axhline(off, color="#888", lw=0.4)
    arch = "dual" if min(dd["nminus"], dd["ntot"] - dd["nminus"]) / max(1, dd["ntot"]) > 0.2 else "uni"
    axB.text(-0.01, off, f"{X.replace('_','/')}", fontsize=7, ha="right", va="center", fontweight="bold", color="#C0392B" if X in WILD else "#222")
    axB.text(1.012, off, f"RPM {dd['rpm']:.0f} · {arch}-strand · {dd['pct_at']:.0f}% AS→TE", fontsize=6, ha="left", va="center", color="#555")
axB.set_xlim(0, 1); axB.set_ylim(-2.4 * nB, 1.3); axB.axis("off")
axB.text(0.5, 1.15, f"B  piRNA cluster in every PRESENT strain — genomic strand (plus ↑ green / minus ↓ purple); RPM + architecture + % ANTISENSE-to-TE (silencing) per strain",
         fontsize=8.7, fontweight="bold", ha="center", transform=axB.get_xaxis_transform() if False else axB.transAxes)
# C: base resolution in top strain
dd = D[TOP]; ps = dd["ps"]; reads = dd["reads"]; tes = dd["tes"]; testr = dd["testr"]
five = Counter((r[1] - 1 if r[2] else r[0]) for r in reads); pk = five.most_common(1)[0][0] if five else ps + dd["N"] // 2
z0, z1 = pk - 30, pk + 50
def anti_te(rs, re, isrev):
    st = testr((rs + re) // 2); return None if st is None else ((st == "-") != isrev)
zreads = [r for r in reads if r[0] < z1 and r[1] > z0]; pr, mr = Counter(), Counter()
for rs, re, isrev, seq in zreads: (mr if isrev else pr)[(rs, re, seq, isrev)] += 1
def draw(items, y0, dirn):
    y = y0
    for (rs, re, seq, isrev), cnt in items:
        for k, ch in enumerate(seq):
            x = rs + k
            if z0 - 2 <= x <= z1 + 2: axC.text(x, y, ch, fontsize=4.6, ha="center", va="center", family="monospace", color="white", bbox=dict(boxstyle="square,pad=0.02", fc=NT.get(ch, "#999"), ec="none"))
        at = anti_te(rs, re, isrev); acol = "#C0392B" if at else ("#888" if at is not None else "#ccc")
        fp = re - 1 if isrev else rs
        axC.annotate("", xy=(fp + (0.5 if not isrev else -0.5), y), xytext=(fp + (-2.6 if not isrev else 2.6), y), arrowprops=dict(arrowstyle="-|>", color=acol, lw=1.0))
        axC.text(z1 + 4, y, f"×{cnt}", fontsize=4.6, va="center", color="#666"); y += dirn
    return y
ytop = draw(pr.most_common(7), 1, 1); axC.axhline(0, color="#333", lw=0.7); ybot = draw(mr.most_common(7), -1, -1)
axC.set_xlim(z0 - 1, z1 + 8); axC.set_ylim(ybot - 1.6, ytop + 0.6)
for sp in ("top", "left", "right"): axC.spines[sp].set_visible(False)
axC.set_yticks([]); axC.spines["bottom"].set_position(("data", ybot - 1.0))
tk = np.linspace(z0, z1, 5).astype(int); axC.set_xticks(tk); axC.set_xticklabels([f"{t:,}" for t in tk], fontsize=6.5); axC.tick_params(axis="x", length=3)
axC.set_xlabel(f"{TOP.replace('_','/')} chr{dd['ch']} position (bp) — every base at its true genomic coordinate", fontsize=7)
axC.text(z0 - 0.6, ytop, "+ strand", fontsize=6.5, color="#33a02c", fontweight="bold", ha="right", va="center")
axC.text(z0 - 0.6, ybot, "− strand", fontsize=6.5, color="#6a3d9a", fontweight="bold", ha="right", va="center")
axC.set_title(f"C  Base resolution in top strain {TOP.replace('_','/')}; bases at true coordinates; 5′-U = 1U signature; 5′ arrow RED = ANTISENSE-to-TE (silencing), grey = sense-to-TE", fontsize=8.0, fontweight="bold", loc="left")
fig.suptitle(f"{GENE} — strain-variable piRNA cluster: expression across all strains + nucleotide resolution", fontsize=12.2, fontweight="bold", y=0.975)
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/{OUT}.{e}", bbox_inches="tight")
print(f"   wrote {OUT}.png ({len(present)} present strains; top={TOP})")
