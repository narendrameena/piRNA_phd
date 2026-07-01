#!/usr/bin/env python3
"""Build per-2Mb-bin, per-family active-TE EXPRESSION from the RNA-seq featureCounts (results/TE_expression_rna/
{strain}.TE.featureCounts). Each row = one TE copy; family is read from the Geneid (L1Md->L1, IAP->IAP,
ERVK/ERVB->ERVK), chrom from Chr (#chrN main chroms only), counts summed across that strain's 9 RNA libraries.
Output active_te_expression_byfamily.tsv (strain, chrom, bin, family, RNA_counts) — same layout as
active_te_byfamily.tsv but transcription instead of genomic bp, so te16/integrated16 can read it identically."""
import glob, os
from collections import defaultdict
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
CHROMS = set([str(i) for i in range(1, 20)] + ["X"])
out = defaultdict(float); nstrain = 0
for f in sorted(glob.glob(f"{ROOT}/results/TE_expression_rna/*.TE.featureCounts")):
    strain = os.path.basename(f).split(".")[0]; nstrain += 1
    for line in open(f):
        if line[0] == "#" or line.startswith("Geneid"): continue
        cols = line.rstrip("\n").split("\t")
        g = cols[0]
        if "L1Md" in g: fam = "L1"
        elif "IAP" in g: fam = "IAP"
        elif "ERVK" in g or "ERVB" in g: fam = "ERVK"
        else: continue
        c = cols[1].split("#")[-1]
        if not c.startswith("chr"): continue
        c = c[3:]
        if c not in CHROMS: continue
        try: b = int(cols[2]) // 2_000_000
        except ValueError: continue
        out[(strain, c, b, fam)] += sum(float(x) for x in cols[6:])
with open(f"{U}/active_te_expression_byfamily.tsv", "w") as o:
    for (s, c, b, fam), v in sorted(out.items()):
        o.write(f"{s}\t{c}\t{b}\t{fam}\t{v:.2f}\n")
fam_tot = defaultdict(float)
for (s, c, b, fam), v in out.items(): fam_tot[fam] += v
print(f"strains={nstrain}; wrote active_te_expression_byfamily.tsv: {len(out)} (strain,chrom,bin,family) entries")
print("family RNA-expression totals:", {k: f"{v:.0f}" for k, v in sorted(fam_tot.items())})
