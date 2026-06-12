#!/usr/bin/env python3
"""SINGLE-strain locus plot for a PICB piRNA cluster, resolved by TIMEPOINT (the developmental dimension).
One strain (top-expressing by PICB FPM, or passed as arg), three panels:
  (A) this strain's PICB-combined cluster FPM across E16.5 / P12.5 / P20.5 — shows WHEN the cluster is expressed
      and which timepoint is chosen for the nucleotide panel (starred);
  (B) genomic-strand read coverage at the locus for EACH timepoint (3 stacked tracks, reps pooled) + RepeatMasker
      TE track; each track labelled with its PICB FPM and % ANTISENSE-to-TE (silencing);
  (C) base resolution in the CHOSEN (top-FPM) timepoint — piRNAs as coloured bases on a true-coordinate ruler.
TERMINOLOGY: genomic +/- strand != sense/antisense; SENSE/ANTISENSE is relative to the TE (RM .out strand),
antisense-to-TE = silencing-competent. PICB COMBINED run; reads from own-genome STAR_srna_strain_wise BAMs; locus
back-lifted GRCm39 -> strain via cactus HAL. Usage: make_pav_locus_single.py <g39chrom> <g39start> <g39end> <gene>
<outbase> [strain]."""
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
G39C, G39S, G39E, GENE, OUT = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
WANT = sys.argv[6] if len(sys.argv) > 6 else None

P = pd.read_csv(f"{CP}/picb_pangenome_fpm.tsv", sep="\t", dtype={"g39_chrom": str})
sub = P[(P.g39_chrom == G39C) & (P.start < G39E) & (P.end > G39S)]
FPM = sub.groupby(["strain", "tp"])["all_primary_FPM"].max().unstack(fill_value=0.0).reindex(index=ORDER, columns=TPS).fillna(0.0)
STRAIN = WANT if WANT else max(ORDER, key=lambda X: FPM.loc[X].max())
fpm_tp = FPM.loc[STRAIN]; CHOSEN = fpm_tp.idxmax()
print(f"{GENE}: single-strain={STRAIN}; FPM by tp={dict(fpm_tp.round(1))}; nucleotide tp={TPLAB[CHOSEN]}")

lc = f"{PG}/.lift_{G39C}_{G39S}_{STRAIN}.bed"
if not os.path.exists(lc):
    open(f"/tmp/_g39s_{OUT}.bed", "w").write(f"{G39C}\t{G39S}\t{G39E}\tx\n")
    subprocess.run(f"singularity exec --bind /mnt {SIF} halLiftover {HAL} GRCm39 /tmp/_g39s_{OUT}.bed {STRAIN} {lc}", shell=True, capture_output=True)
lift = [l.split("\t") for l in open(lc) if l.strip()]
CH = Counter(r[0] for r in lift).most_common(1)[0][0]; ps = min(int(r[1]) for r in lift if r[0] == CH); pe = max(int(r[2]) for r in lift if r[0] == CH)
BAMC = f"{STRAIN}#1#chr{CH}"; N = max(1, pe - ps); nb = 200
outf = glob.glob(f"{ROOT}/resources/repeatMasker/{STRAIN}_*.out"); tes = []
if outf:
    aw = subprocess.run("awk -v c=%s -v s=%d -v e=%d 'NR>3 && $5==c && $7>s && $6<e{st=($9==\"C\")?\"-\":\"+\"; print $6\"\\t\"$7\"\\t\"st\"\\t\"$11}' %s" % (BAMC, ps, pe, outf[0]), shell=True, capture_output=True, text=True).stdout
    for ln in aw.splitlines():
        f = ln.split("\t"); tes.append((int(f[0]), int(f[1]), f[2], f[3]))
tes.sort()
def testr(pos):
    for ts, te, st, fam in tes:
        if ts <= pos < te: return st
    return None
def fetch(tp):
    plus = np.zeros(nb); minus = np.zeros(nb); reads = []; nat = nte = n1u = 0
    for r in (1, 2, 3):
        b = f"{ROOT}/results/STAR_srna_strain_wise/{STRAIN}/{STRAIN}-{tp}.{r}/Aligned.sortedByCoord.out.bam"
        if not os.path.exists(b): continue
        try:
            bam = pysam.AlignmentFile(b, "rb")
            for a in bam.fetch(BAMC, ps, pe):
                if a.is_unmapped or not a.query_sequence: continue
                L = a.reference_end - a.reference_start
                if not (24 <= L <= 32): continue
                b0 = max(0, min(nb - 1, int((a.reference_start - ps) / N * nb))); b1 = max(0, min(nb, int((a.reference_end - ps) / N * nb)))
                seq = a.query_sequence.upper(); (minus if a.is_reverse else plus)[b0:b1] += 1
                reads.append((a.reference_start, a.reference_end, a.is_reverse, seq)); n1u += ((comp.get(seq[-1], "N") if a.is_reverse else seq[0]) == "T")
                st = testr((a.reference_start + a.reference_end) // 2)
                if st is not None: nte += 1; nat += ((st == "-") != a.is_reverse)
            bam.close()
        except (OSError, ValueError):
            try: bam.close()
            except Exception: pass
    return dict(plus=plus, minus=minus, reads=reads, ntot=len(reads), nat=nat, nte=nte, n1u=n1u)
DAT = {tp: fetch(tp) for tp in TPS}

plt.rcParams.update({"font.family": "Liberation Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(12.5, 10), dpi=300); gs = fig.add_gridspec(3, 1, height_ratios=[0.8, 1.7, 1.5], hspace=0.6); fig.subplots_adjust(top=0.9, bottom=0.06)
axA, axB, axC = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
xb = np.arange(3)
axA.bar(xb, [max(fpm_tp[t], 1e-2) for t in TPS], width=0.6, color=[TPCOL[t] for t in TPS], edgecolor="white")
axA.set_yscale("log"); axA.set_ylim(0.05, max(fpm_tp.max() * 3, 1)); axA.set_xticks(xb); axA.set_xticklabels([TPLAB[t] for t in TPS], fontsize=9)
for i, t in enumerate(TPS):
    if fpm_tp[t] > 0: axA.text(i, fpm_tp[t] * 1.15, f"{fpm_tp[t]:.1f}" + ("  ★" if t == CHOSEN else ""), ha="center", va="bottom", fontsize=8, fontweight="bold", color=TPCOL[t])
axA.set_ylabel("PICB FPM (log)", fontsize=8.5); axA.spines[["top", "right"]].set_visible(False)
axA.set_title(f"A  {STRAIN.replace('_','/')} PICB-cluster FPM across development — ★ = timepoint shown at nucleotide resolution (C)", fontsize=9.2, fontweight="bold", loc="left")
xg = np.linspace(ps, pe, nb); GM = max(max(DAT[t]["plus"].max(), DAT[t]["minus"].max()) for t in TPS) or 1
for i, tp in enumerate(TPS):
    d = DAT[tp]; off = -i * 2.6
    axB.fill_between(xg, off, off + d["plus"] / GM, color="#33a02c", alpha=0.8, step="mid"); axB.fill_between(xg, off, off - d["minus"] / GM, color="#6a3d9a", alpha=0.8, step="mid")
    axB.axhline(off, color="#888", lw=0.4)
    at = 100 * d["nat"] / max(1, d["nte"])
    axB.text(ps - 0.012 * N, off, TPLAB[tp], fontsize=8, ha="right", va="center", fontweight="bold", color=TPCOL[tp])
    axB.text(pe + 0.004 * N, off, f"FPM {fpm_tp[tp]:.1f} · {d['ntot']:,} reads · {at:.0f}% AS→TE · {100*d['n1u']/max(1,d['ntot']):.0f}% 1U", fontsize=6.2, ha="left", va="center", color="#555")
tey = -2.6 * 3 + 0.4
for (ts, te, st, fam) in tes:
    axB.add_patch(Rectangle((ts, tey), te - ts, 0.3, facecolor=TECOL.get(fam, "#ddd"), edgecolor="none"))
    axB.annotate("", xy=(te if st == "+" else ts, tey + 0.5), xytext=(ts if st == "+" else te, tey + 0.5), arrowprops=dict(arrowstyle="-|>", color="#444", lw=0.5))
axB.set_xlim(ps - 0.13 * N, pe + 0.16 * N); axB.set_ylim(tey - 0.5, 1.4); axB.axis("off")
famset = sorted({f for _, _, _, f in tes}, key=lambda f: -sum(te - ts for ts, te, _, ff in tes if ff == f))[:5]
axB.legend(handles=[Patch(facecolor="#33a02c", label="+ strand"), Patch(facecolor="#6a3d9a", label="− strand")] + [Patch(facecolor=TECOL.get(f, "#ddd"), label=f) for f in famset], fontsize=6.2, frameon=False, loc="lower center", ncol=7, bbox_to_anchor=(0.5, -0.07))
axB.set_title(f"B  {STRAIN.replace('_','/')} {BAMC}:{ps:,}-{pe:,} — genomic-strand coverage per timepoint (plus ↑ green / minus ↓ purple); TE track at bottom", fontsize=9.0, fontweight="bold", loc="left")
d = DAT[CHOSEN]; reads = d["reads"]
five = Counter((r[1] - 1 if r[2] else r[0]) for r in reads); pk = five.most_common(1)[0][0] if five else ps + N // 2; z0, z1 = pk - 30, pk + 50
axB.add_patch(Rectangle((z0, tey - 0.15), z1 - z0, 1.55 - (tey - 0.15), facecolor="#FDB863", alpha=0.22, edgecolor="#E8A33D", lw=0.5, zorder=0))   # zoom highlighter bar across Panel B
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
    fig.add_artist(ConnectionPatch(xyA=(x, tey), coordsA=axB.transData, xyB=(x, ytop + 0.6), coordsB=axC.transData, color="#E8A33D", lw=0.8, ls=(0, (3, 2))))
axC.set_title(f"C  Base resolution in {STRAIN.replace('_','/')} at {TPLAB[CHOSEN]} (top-FPM timepoint); 5′-U = 1U; 5′ arrow RED = ANTISENSE-to-TE (silencing), grey = sense-to-TE", fontsize=8.0, fontweight="bold", loc="left")
fig.suptitle(f"{GENE} — SINGLE-strain PICB cluster across development ({STRAIN.replace('_','/')})", fontsize=12, fontweight="bold", y=0.965)
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/{OUT}.{e}", bbox_inches="tight")
print(f"   wrote {OUT}.png")
