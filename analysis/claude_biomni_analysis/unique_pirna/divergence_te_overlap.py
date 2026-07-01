#!/usr/bin/env python3
"""Is the DIVERGENCE class TE-associated? Fraction of divergence loci (genome-conserved, strain-restricted) that
overlap a driver-family TE in the carrier's own genome, vs a baseline of random PICB clusters. TE overlap =
correlation (TE-derived locus), NOT proof of causation; the private-INSERTION mechanism is excluded by construction
(divergence loci are genome-conserved)."""
import sys, random; sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis"); sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
import pav_clusters as pc, pandas as pd
from collections import Counter
PG = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"; CP = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"
DRIVERS = {"LINE/L1", "LTR/ERVK", "LTR/ERV1", "LTR/ERVL-MaLR", "LTR/ERVL", "SINE/B2", "LINE/L1?", "LTR/ERVK?"}
cand = pd.read_csv(f"{PG}/divergence_candidates.tsv", sep="\t", dtype={"g39_chrom": str})
sup = pd.read_csv(f"{CP}/locus_genome_pav_divergence.tsv", sep="\t", dtype={"g39_chrom": str})
gn = sup[sup.present].groupby(["g39_chrom", "g39_start"]).size().rename("genome_n").reset_index()
div = cand.merge(gn, on=["g39_chrom", "g39_start"], how="left")
div = div[div.genome_n >= 10]
def te_profile(rows_iter, get_owncoords):
    anyte = hit = n = 0; fams = Counter()
    for r in rows_iter:
        oc, s, e, carrier = get_owncoords(r)
        if oc is None: continue
        tes = pc.te_at(carrier, oc, s, e)
        if not tes: n += 1; continue
        n += 1; anyte += 1
        dr = [f.split("|")[-1] for _, _, _, f in tes if f.split("|")[-1] in DRIVERS]
        if dr: hit += 1; fams[Counter(dr).most_common(1)[0][0]] += 1
    return n, anyte, hit, fams
# divergence sample
random.seed(1); ds = div.sample(min(200, len(div)), random_state=1)
def div_coords(r):
    sub = pc.clusters_at(r.g39_chrom, int(r.g39_start), int(r.g39_end)); cc = sub[sub.strain == r.carrier]
    if len(cc) == 0: return None, 0, 0, None
    return str(cc.iloc[0].own_chrom), int(cc.own_start.min()), int(cc.own_end.max()), r.carrier
nd, ad, hd, fd = te_profile((r for _, r in ds.iterrows()), div_coords)
# baseline: random PICB clusters
P = pc._pang(); bs = P.sample(200, random_state=2)
def base_coords(r): return str(r.own_chrom), int(r.own_start), int(r.own_end), r.strain
nb, ab, hb, fb = te_profile((r for _, r in bs.iterrows()), base_coords)
print(f"DIVERGENCE loci (n={nd}): any-TE {100*ad/nd:.0f}% | DRIVER-TE {100*hd/nd:.0f}% | top families {fd.most_common(5)}")
print(f"BASELINE random clusters (n={nb}): any-TE {100*ab/nb:.0f}% | DRIVER-TE {100*hb/nb:.0f}% | top families {fb.most_common(5)}")
from scipy.stats import fisher_exact
odds, p = fisher_exact([[hd, nd - hd], [hb, nb - hb]])
print(f"DRIVER-TE enrichment divergence vs baseline: OR={odds:.2f}, Fisher P={p:.3g}")
