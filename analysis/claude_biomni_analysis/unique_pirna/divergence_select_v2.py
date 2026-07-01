#!/usr/bin/env python3
"""Phase 3 (v2): select DIVERGENCE-driven loci = genome-CONSERVED but STRAIN-RESTRICTED CLUSTER expression.
Expression breadth is verified with clusters_at (authoritative; the interval-merge only UNDER-counts, so it is a
safe pre-filter superset). Genome conservation is added later by HAL projection. Output g39 loci -> render via
make_pav_locus_multi (cluster-based expression, correct for divergence)."""
import sys, pandas as pd
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
import pav_clusters as pc
PG = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
P = pc._pang()
P2 = P.sort_values(["g39_chrom", "start"]).reset_index(drop=True)
_cmax = P2.groupby("g39_chrom")["end"].cummax()
P2 = P2.assign(locid=((P2.g39_chrom != P2.g39_chrom.shift()) | (P2.start > _cmax.shift() + 5000)).cumsum())
g = P2.groupby("locid").agg(g39_chrom=("g39_chrom", "first"), start=("start", "min"), end=("end", "max"),
                            mergebreadth=("strain", "nunique"), maxFPM=("all_primary_FPM", "max")).reset_index()
cand = g[(g.mergebreadth <= 3) & (g.maxFPM >= 10)].sort_values("maxFPM", ascending=False).head(600)   # superset (merge under-counts); top-FPM for clear examples
print(f"merge pre-filter candidates (mergebreadth<=3, FPM>=10): {len(g[(g.mergebreadth<=3)&(g.maxFPM>=10)])}; verifying top {len(cand)} with clusters_at")
rows = []
for _, r in cand.iterrows():
    sub = pc.clusters_at(str(r.g39_chrom), int(r.start), int(r.end))
    b = sub.strain.nunique()
    if b == 0 or b > 3: continue                       # clusters_at AUTHORITATIVE breadth: strain-restricted only
    carrier = sub.loc[sub.all_primary_FPM.idxmax(), "strain"]
    rows.append((str(r.g39_chrom), int(r.start), int(r.end), b, carrier, round(float(sub.all_primary_FPM.max()), 1), ",".join(sorted(sub.strain.unique()))))
out = pd.DataFrame(rows, columns=["g39_chrom", "g39_start", "g39_end", "breadth", "carrier", "maxFPM", "expressing"]).sort_values("maxFPM", ascending=False)
out.to_csv(f"{PG}/divergence_candidates.tsv", sep="\t", index=False)
print(f"clusters_at-verified strain-restricted (breadth<=3): {len(out)}")
print(out.head(25).to_string(index=False))
