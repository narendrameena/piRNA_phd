#!/usr/bin/env python3
"""Sense/antisense-to-TE of strain-specific piRNAs across 16 strains (klass5 ≥2-read). ORIENTATION IS RELATIVE TO
THE TE FEATURE STRAND, NOT genomic +/-: piRNA on the SAME strand as the TE = sense (from the TE's own transcript);
OPPOSITE = antisense (silencing-competent — can base-pair the TE mRNA). 50% = no strand bias.
Per strain: piRNA locus+strand (cand_self16 BAM) vs stranded TE annotation (RM .out -> chr-named stranded BED).
(A) antisense fraction by klass5 class, pooled over 16 strains (dots = per-strain). (B) antisense fraction by top TE
family, SPLIT by the two genuinely-unique subcategories (conserved-but-silent vs strain-private). Cached per-candidate."""
import warnings; warnings.filterwarnings("ignore")
import sys, os, subprocess, tempfile
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
from collections import defaultdict
import pandas as pd, numpy as np, pysam
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
SA = f"{U}/sense_antisense"; PG = f"{U}/pangenome_te"; BT = "/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"
CANON = [s for s in STRAIN_ORDER if s != "C57BL_6"]
CLS = ["expressed elsewhere (exact)", "SNP-variant (1-3mm)", "unique: conserved-but-silent", "unique: strain-private locus"]
CLAB = ["expressed-\nelsewhere", "SNP-variant\n(allelic)", "unique:\nconserved-silent", "unique:\nprivate-locus"]
CCOL = ["#9e9e9e", "#E69F00", "#0072B2", "#7a3b9a"]
cache = f"/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/08_unique_piRNA_TE_origin_and_strand/data/source_data/SourceData_sense_antisense16_percand.csv.gz"
if os.path.exists(cache):
    pc = pd.read_csv(cache); print("loaded cached per-candidate orientation table")
else:
    d = pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz", usecols=["strain", "sequence", "timepoint", "klass5"])
    recs = []
    for X in CANON:
        outf = f"{ROOT}/resources/repeatMasker/{X}_chromosomes_MT_unplaced.fasta.out"
        if not os.path.exists(outf): print("no .out", X); continue
        teb = f"{SA}/{X}.TE_stranded16.bed"   # chr-named stranded TE BED (matches cand_self16 ref split)
        subprocess.run(f"""awk 'NR>3 && NF>=11 {{n=$5; sub(/.*#/,"",n); st=($9=="C")?"-":"+"; print n"\\t"($6-1)"\\t"$7"\\t"$10"|"$11"\\t0\\t"st}}' "{outf}" | sort -k1,1 -k2,2n > "{teb}" """, shell=True)   # main chroms in .out already 'chrN' -> do NOT prepend (matches cand_self16 ref split)
        g = d[d.strain == X]; kl = {X + "|" + r.timepoint + "|" + r.sequence: r.klass5 for r in g.itertuples()}
        bam = pysam.AlignmentFile(f"{U}/cand_self16/{X}.cand_self16.bam", "rb")
        pb = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False, dir=SA); seen = set()
        for a in bam.fetch(until_eof=True):
            if a.is_unmapped or a.query_name not in kl or a.query_name in seen: continue
            seen.add(a.query_name); pb.write(f"{a.reference_name.split('#')[-1]}\t{a.reference_start}\t{a.reference_end}\t{a.query_name}\t0\t{'-' if a.is_reverse else '+'}\n")
        pb.close(); bam.close()
        out = subprocess.run(f"sort -k1,1 -k2,2n {pb.name} | {BT} intersect -a - -b {teb} -wo", shell=True, capture_output=True, text=True).stdout
        os.unlink(pb.name)
        best = {}   # id -> (overlap, family, orientation)  [orientation = piRNA strand vs TE strand]
        for ln in out.splitlines():
            f = ln.split("\t")
            if len(f) < 12: continue
            cid = f[3]; p_str = f[5]; fam = f[9].split("|")[-1]; t_str = f[11]; ov = int(f[-1])
            orient = "sense" if p_str == t_str else "antisense"
            if cid not in best or ov > best[cid][0]: best[cid] = (ov, fam, orient)
        for cid, (ov, fam, orient) in best.items(): recs.append((X, cid, kl[cid], fam, orient))
        print(f"{X}: TE-overlapping candidates = {len(best):,}")
    pc = pd.DataFrame(recs, columns=["strain", "id", "klass5", "family", "orientation"]); pc.to_csv(cache, index=False)
# ---- plot ----
plt.rcParams.update({"font.family": "Liberation Sans"})
fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.2), dpi=300)
# A: antisense % by class, pooled over 16 strains; dots = per-strain
x = np.arange(len(CLS)); rng = np.random.default_rng(0)
for i, (k, col) in enumerate(zip(CLS, CCOL)):
    sub = pc[pc.klass5 == k]; pooled = (sub.orientation == "antisense").mean() * 100 if len(sub) else np.nan
    axA.bar(i, pooled, 0.62, color=col, edgecolor="white", zorder=2)
    per = sub.groupby("strain").orientation.apply(lambda s: (s == "antisense").mean() * 100)
    axA.scatter(np.full(len(per), i) + rng.uniform(-0.16, 0.16, len(per)), per.values, s=14, c="#333", edgecolor="white", lw=0.3, zorder=3, alpha=0.8)
    if pd.notna(pooled): axA.text(i, pooled + 0.6, f"{pooled:.0f}%", ha="center", va="bottom", fontsize=8, color=col, fontweight="bold")
axA.axhline(50, ls="--", lw=0.8, color="#444"); axA.text(len(CLS) - 0.5, 50.4, "no strand bias (50%)", ha="right", va="bottom", fontsize=7, color="#444")
axA.set_xticks(x); axA.set_xticklabels(CLAB, fontsize=7.5); axA.set_ylim(35, 70)
axA.set_ylabel("% antisense to TE (silencing-competent)", fontsize=9); axA.spines[["top", "right"]].set_visible(False)
axA.set_title("A  Unique piRNAs are more antisense-to-TE than common\n(16 strains pooled; dots = per-strain)", fontsize=9.4, fontweight="bold", loc="left")
# B: antisense % by top TE family, split by the 2 genuinely-unique subcategories
SUB = [("unique: strain-private locus", "strain-private", "#7a3b9a"), ("unique: conserved-but-silent", "conserved-but-silent", "#0072B2")]
gu = pc[pc.klass5.isin([k for k, _, _ in SUB])]
top = gu.family.value_counts().head(9).index.tolist()[::-1]
y = np.arange(len(top)); h = 0.38
for j, (k, lab, col) in enumerate(SUB):
    sk = gu[gu.klass5 == k]
    vals = [(sk[sk.family == f].orientation == "antisense").mean() * 100 if (sk.family == f).any() else np.nan for f in top]
    ns = [int((sk.family == f).sum()) for f in top]
    axB.barh(y + (0.5 - j) * h, vals, h, color=col, edgecolor="white", label=lab, zorder=3)
    for yi, v, nn in zip(y + (0.5 - j) * h, vals, ns):
        if pd.notna(v): axB.text(v + 0.4, yi, f"{v:.0f}% (n={nn:,})", va="center", fontsize=5.6, color=col)
axB.axvline(50, ls="--", lw=0.8, color="#444"); axB.set_yticks(y); axB.set_yticklabels(top, fontsize=7.5); axB.set_xlim(35, 75)
axB.set_xlabel("% antisense to TE", fontsize=9); axB.legend(fontsize=7.5, frameon=False, loc="lower right")
axB.set_title("B  Antisense fraction by TE family — split by subcategory\n(strain-private vs conserved-but-silent)", fontsize=9.4, fontweight="bold", loc="left")
axB.spines[["top", "right"]].set_visible(False)
fig.text(0.5, 0.005, "Orientation is relative to the TE feature strand (sense = same strand as the TE; antisense = opposite = silencing-competent), NOT genomic +/-. klass5 ≥2-read; 16 strains.", ha="center", fontsize=6.5, color="#555")
fig.tight_layout(rect=[0, 0.02, 1, 1])
for e in ("pdf", "svg", "png"): fig.savefig(f"{SA}/Fig_sense_antisense.{e}", bbox_inches="tight")
print("wrote Fig_sense_antisense (16-strain, by subcategory)")
print(pc.groupby("klass5").orientation.apply(lambda s: round((s == "antisense").mean() * 100, 1)))
