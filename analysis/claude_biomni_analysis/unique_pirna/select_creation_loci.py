#!/usr/bin/env python3
"""Confound-fix step 3a: re-select source-locus figures to CREATION only (true new strain-private loci, cluster
breadth<=3 by clusters_at) — excludes PROPAGATION (private insertion landed in a conserved cluster, e.g. chr17
pachytene master). Top-FPM creation loci per strain, over a driver TE. -> source_loci_master_creation.tsv."""
import sys; sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis"); sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
import pav_clusters as pc, pandas as pd
from collections import Counter
from strain_order import STRAIN_ORDER
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]; PG = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
DRIVERS = {"LINE/L1", "LTR/ERVK", "LTR/ERV1", "LTR/ERVL-MaLR", "LTR/ERVL", "SINE/B2"}; N_PER = 3; MAXSCAN = 150
P = pc._pang(); rows = []
for X in ORDER:
    d = pd.read_csv(f"{PG}/TE_driven_COORDINATE16_{X}.csv").drop_duplicates(["chrom", "start"]).sort_values("FPM", ascending=False)
    n = 0; used = []
    for _, r in d.head(MAXSCAN).iterrows():
        if n >= N_PER: break
        chrom = str(r.chrom); S, E = int(r.start), int(r.end)
        row = P[(P.strain == X) & (P.own_chrom == chrom) & (P.own_start >= S - 3000) & (P.own_start <= E + 3000)]
        if len(row) == 0: continue
        if pc.clusters_at(row.iloc[0].g39_chrom, int(row.start.min()), int(row.end.max())).strain.nunique() > 3: continue   # CREATION only (breadth<=3)
        mid = (S + E) // 2; S2, E2 = (mid - 2500, mid + 2500) if (E - S) > 6000 else (S, E)
        if any(c == chrom and abs(m - mid) < 30000 for c, m in used): continue
        tes = pc.te_at(X, chrom, S2, E2); dt = [(s, e, st, f.split("|")[-1]) for s, e, st, f in tes if f.split("|")[-1] in DRIVERS]
        if not dt: continue
        cov = Counter()
        for s, e, st, f in dt: cov[f] += min(e, E2) - max(s, S2)
        domfam = cov.most_common(1)[0][0]; strand = [st for s, e, st, f in dt if f == domfam][0]
        rows.append((X, f"{X}#1#chr{chrom}", S2, E2, domfam, strand, round(float(r.FPM), 1), "creation")); used.append((chrom, mid)); n += 1
out = pd.DataFrame(rows, columns=["carrier", "chrom_own", "start", "end", "te", "strand", "FPM", "class"])
out.to_csv(f"{PG}/source_loci_master_creation.tsv", sep="\t", index=False)
print(f"selected {len(out)} CREATION loci across {out.carrier.nunique()}/16 strains")
print(out.groupby("carrier").size().to_string())
