#!/usr/bin/env python3
"""Phase 3: select DIVERGENCE-driven piRNA loci — genome-CONSERVED (present in many strains) but STRAIN-RESTRICTED
in expression (piRNAs in only a few strains) = the 'present-but-silent' pattern, complementing the insertion-driven
(locus-absent-elsewhere) set. Output g39 loci -> rendered by the existing make_pav_locus_multi.py.
NOTE: divergence loci are WEAKLY expressed (the chr1wildtrio archetype is maxFPM~27), so the FPM floor is LOW;
loci are merged by INTERVAL overlap (not fixed bins) so expression breadth is counted correctly (== clusters_at)."""
import sys, pandas as pd
U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"; CP = f"{U}/cluster_pav"; PG = f"{U}/pangenome_te"
P = pd.read_csv(f"{CP}/picb_pangenome_clusters.tsv", sep="\t", dtype={"g39_chrom": str, "own_chrom": str}, low_memory=False)
P = P.sort_values(["g39_chrom", "start"]).reset_index(drop=True)   # interval-merge clusters into g39 loci
_cmax = P.groupby("g39_chrom")["end"].cummax()
P["locid"] = ((P.g39_chrom != P.g39_chrom.shift()) | (P.start > _cmax.shift() + 5000)).cumsum()
g = P.groupby("locid").agg(g39_chrom=("g39_chrom", "first"), start=("start", "min"), end=("end", "max"), n_expr=("strain", "nunique"),
                           maxFPM=("all_primary_FPM", "max"), strains=("strain", lambda s: ",".join(sorted(set(s))))).reset_index()
cand = g[(g.n_expr <= 3) & (g.maxFPM >= 10)].sort_values("maxFPM", ascending=False).copy()   # strain-restricted + visible (LOW FPM floor: divergence loci are weak)
print(f"candidates (n_expr<=3, maxFPM>=10): {len(cand)}")
cand = cand.head(800)   # genome_n is the slow per-locus catalogue lookup → only the top-FPM candidates
C = pd.read_csv(f"{CP}/cluster_PAV_catalogue.csv.gz", dtype={"chrom": str})
def genome_n(chrom, s, e):
    seg = C[(C.chrom == chrom) & (C.start < e) & (C.end > s)]
    return int(seg.n_strains.max()) if len(seg) else 0
cand["genome_n"] = [genome_n(r.g39_chrom, r.start, r.end) for _, r in cand.iterrows()]
def top_strain(chrom, s, e):
    sub = P[(P.g39_chrom == chrom) & (P.start < e) & (P.end > s)]
    return sub.loc[sub.all_primary_FPM.idxmax(), "strain"] if len(sub) else "NA"
div = cand[cand.genome_n >= 12].copy()                  # CONSERVED genome (>=12/16) but expressed in <=3 strains = present-but-silent
div["carrier"] = [top_strain(r.g39_chrom, r.start, r.end) for _, r in div.iterrows()]
div = div.sort_values(["genome_n", "maxFPM"], ascending=False)
sel = []; seen = {}
for _, r in div.iterrows():                              # diversify: <=2 per carrier strain, top 12
    if seen.get(r.carrier, 0) >= 2: continue
    seen[r.carrier] = seen.get(r.carrier, 0) + 1; sel.append(r)
    if len(sel) >= 12: break
print(f"conserved divergence (genome_n>=12): {len(div)} | selected: {len(sel)}")
if not sel:
    print("NO divergence loci at these thresholds (genome-conserved + strain-restricted is rare)"); sys.exit(0)
out = pd.DataFrame(sel)[["g39_chrom", "start", "end", "carrier", "n_expr", "genome_n", "maxFPM", "strains"]]
out.to_csv(f"{PG}/divergence_loci_selected.tsv", sep="\t", index=False)
print(out.to_string(index=False))
