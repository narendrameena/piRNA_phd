#!/usr/bin/env python3
"""Build per-2Mb-bin, per-family piRNA-on-TE abundance from the sRNA featureCounts (results/TEtranscriptCount/
all_featureCounts_TE.tab = small-RNA/piRNA reads on each TE copy; VERIFIED provenance: STAR_srna_strain_wise
BAMs, piRNA params). Columns: Geneid(with family) | Chr(strain#1#chrN) | Start | End | Strand | Length | count |
sample. Family from Geneid (L1Md->L1, IAP->IAP, ERVK/ERVB->ERVK), counts summed per (strain, chrom, bin, family)
over all that strain's libraries. Output active_pirna_on_te_byfamily.tsv — the DIRECT piRNA response to each TE
family per region (the defence side of the TE arms-race), distinct from PICB cluster loci."""
import os
from collections import defaultdict
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
CHROMS = set([str(i) for i in range(1, 20)] + ["X"])
F = f"{ROOT}/results/TEtranscriptCount/all_featureCounts_TE.tab"
out = defaultdict(float); n = 0
for line in open(F):
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 8: continue
    g = cols[0]
    if "L1Md" in g: fam = "L1"
    elif "IAP" in g: fam = "IAP"
    elif "ERVK" in g or "ERVB" in g: fam = "ERVK"
    else: continue
    parts = cols[1].split("#")
    if len(parts) < 3: continue
    strain = parts[0]; c = parts[-1]
    c = c[3:] if c.startswith("chr") else c
    if c not in CHROMS: continue
    try: b = int(cols[2]) // 2_000_000; cnt = float(cols[6])
    except ValueError: continue
    out[(strain, c, b, fam)] += cnt; n += 1
with open(f"{U}/active_pirna_on_te_byfamily.tsv", "w") as o:
    for (s, c, b, fam), v in sorted(out.items()):
        o.write(f"{s}\t{c}\t{b}\t{fam}\t{v:.2f}\n")
fam_tot = defaultdict(float)
for (s, c, b, fam), v in out.items(): fam_tot[fam] += v
strains = sorted(set(k[0] for k in out))
print(f"rows parsed={n}; strains={len(strains)}; wrote active_pirna_on_te_byfamily.tsv ({len(out)} entries)")
print("piRNA-on-TE family totals:", {k: f"{v:.0f}" for k, v in sorted(fam_tot.items())})
