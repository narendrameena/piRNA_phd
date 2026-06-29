#!/usr/bin/env python3
"""Phase 4 selection: top INSERTION-DRIVEN source loci per strain (all 16) from TE_driven_COORDINATE16_{X}.csv,
centred on a ~5 kb window over the driver TE, for rendering with make_source_pav_multi.py. Driver TE comes from
the strain's own RepeatMasker .out (te_at) so no per-strain TE_stranded.bed is needed."""
import sys, pandas as pd
from collections import Counter
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
from strain_order import STRAIN_ORDER
import pav_clusters as pc
PG = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]
DRIVERS = {"LINE/L1", "LTR/ERVK", "LTR/ERV1", "LTR/ERVL-MaLR", "LTR/ERVL", "SINE/B2"}
N_PER_STRAIN = 3
rows = []
for X in ORDER:
    try: d = pd.read_csv(f"{PG}/TE_driven_COORDINATE16_{X}.csv")
    except FileNotFoundError: continue
    d = d.drop_duplicates(["chrom", "start"]).sort_values("FPM", ascending=False)
    n = 0; used = []
    for _, r in d.iterrows():
        if n >= N_PER_STRAIN: break
        chrom = str(r.chrom); S, E = int(r.start), int(r.end)
        mid = (S + E) // 2; S2, E2 = (mid - 2500, mid + 2500) if (E - S) > 6000 else (S, E)
        if any(c == chrom and abs(m - mid) < 30000 for c, m in used): continue   # dedup nearby
        tes = pc.te_at(X, chrom, S2, E2)
        dt = [(s, e, st, f.split("|")[-1]) for s, e, st, f in tes if f.split("|")[-1] in DRIVERS]
        if not dt: continue
        cov = Counter()
        for s, e, st, f in dt: cov[f] += min(e, E2) - max(s, S2)
        domfam = cov.most_common(1)[0][0]; strand = [st for s, e, st, f in dt if f == domfam][0]
        rows.append((X, f"{X}#1#chr{chrom}", S2, E2, domfam, strand, round(float(r.FPM), 1), "insertion")); used.append((chrom, mid)); n += 1
out = pd.DataFrame(rows, columns=["carrier", "chrom_own", "start", "end", "te", "strand", "FPM", "class"])
out.to_csv(f"{PG}/source_loci_master_insertion.tsv", sep="\t", index=False)
print(f"selected {len(out)} insertion-driven loci across {out.carrier.nunique()}/16 strains")
print(out.groupby("carrier").size().to_string())
