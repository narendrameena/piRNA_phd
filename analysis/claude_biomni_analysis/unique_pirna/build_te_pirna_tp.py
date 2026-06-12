#!/usr/bin/env python3
"""Per-TIMEPOINT, per-family active-TE EXPRESSION (RNA-seq) and piRNA-on-TE (sRNA), re-parsed from the existing
featureCounts (no new counting). RNA: each strain's featureCounts has 9 sample columns (3 tp x 3 reps); the column
header carries the timepoint, so sum reps within a tp. sRNA: all_featureCounts_TE.tab has one row per copy per
sample, sample (col8) carries the tp. Output two TSVs (strain, tp, chrom, bin, family, value)."""
import glob, os
from collections import defaultdict
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
CHROMS = set([str(i) for i in range(1, 20)] + ["X"]); TPMAP = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}
def fam(g):
    if "L1Md" in g: return "L1"
    if "IAP" in g: return "IAP"
    if "ERVK" in g or "ERVB" in g: return "ERVK"
    return None
def chrom(x):
    c = x.split("#")[-1]; c = c[3:] if c.startswith("chr") else c
    return c if c in CHROMS else None
def tpof(s):
    for k, v in TPMAP.items():
        if k in s: return v
    return None
# --- RNA TE expression per tp ---
rna = defaultdict(float)
for f in sorted(glob.glob(f"{ROOT}/results/TE_expression_rna/*.TE.featureCounts")):
    strain = os.path.basename(f).split(".")[0]; lines = open(f); next(lines)
    hdr = next(lines).rstrip("\n").split("\t"); coltp = [tpof(h) for h in hdr]
    for line in lines:
        cols = line.rstrip("\n").split("\t"); fm = fam(cols[0]); c = chrom(cols[1])
        if not fm or not c: continue
        b = int(cols[2]) // 2_000_000
        for i in range(6, len(cols)):
            if coltp[i] and cols[i] != "0":
                rna[(strain, coltp[i], c, b, fm)] += float(cols[i])
# --- sRNA piRNA-on-TE per tp ---
pi = defaultdict(float)
for line in open(f"{ROOT}/results/TEtranscriptCount/all_featureCounts_TE.tab"):
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 8: continue
    fm = fam(cols[0]); c = chrom(cols[1]); tp = tpof(cols[7]); strain = cols[1].split("#")[0]
    if not fm or not c or not tp: continue
    try: pi[(strain, tp, c, int(cols[2]) // 2_000_000, fm)] += float(cols[6])
    except ValueError: continue
for name, D in [("active_te_expression_byfamily_tp.tsv", rna), ("active_pirna_on_te_byfamily_tp.tsv", pi)]:
    with open(f"{U}/{name}", "w") as o:
        for (s, tp, c, b, fm), val in sorted(D.items()): o.write(f"{s}\t{tp}\t{c}\t{b}\t{fm}\t{val:.2f}\n")
    print(f"wrote {name}: {len(D)} entries, {len(set(k[0] for k in D))} strains, tps={sorted(set(k[1] for k in D))}")
