#!/usr/bin/env python3
"""VERIFY the propagation/creation confound with a TIGHT read threshold (no background). A strain counts as
expressing the orthologous locus only if its pooled primary-piRNA reads there are >=10% of the carrier's reads
(a comparable cluster, not background; absolute floor 100). Balanced sample: clusters_at-CREATION (breadth<=3)
+ clusters_at-PROPAGATION (breadth>=10), carrier FPM>=50. Compares clusters_at breadth vs tight READ breadth ->
does the confound hold under read evidence?"""
import sys, subprocess, tempfile, os
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis"); sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
import pav_clusters as pc, pandas as pd
from strain_order import STRAIN_ORDER
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; SIF = f"{ROOT}/cactus_v2.9.3.sif"; HAL = f"{ROOT}/results/pangenome/output/mouse_17strain_pangenome.full.hal"; PG = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]; TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]
THR_REL = 0.10; THR_ABS = 100; FPM_MIN = 50
P = pc._pang()
def reads_at(X, ch, s, e): return sum((pc.fetch_primary(X, ch, s, e, tp, 50) or {}).get("ntot", 0) for tp in TPS)
samp = []
for X in ["CAST_EiJ", "SPRET_EiJ", "PWK_PhJ", "WSB_EiJ", "NOD_ShiLtJ", "BALB_cJ", "AKR_J", "C57BL_6NJ"]:
    d = pd.read_csv(f"{PG}/TE_driven_COORDINATE16_{X}.csv").drop_duplicates(["chrom", "start"]).sort_values("FPM", ascending=False)
    cre = []; prop = []
    for r in d.head(100).itertuples():
        if r.FPM < FPM_MIN or (len(cre) >= 2 and len(prop) >= 2): break
        chrom = str(r.chrom); S, E = int(r.start), int(r.end)
        row = P[(P.strain == X) & (P.own_chrom == chrom) & (P.own_start >= S - 3000) & (P.own_start <= E + 3000)]
        if len(row) == 0: continue
        ca = pc.clusters_at(row.iloc[0].g39_chrom, int(row.start.min()), int(row.end.max())).strain.nunique()
        if ca <= 3 and len(cre) < 2: cre.append((X, chrom, S, E, r.FPM, ca, "creation"))
        elif ca >= 10 and len(prop) < 2: prop.append((X, chrom, S, E, r.FPM, ca, "propagation"))
    samp += cre + prop
rows = []
for X, chrom, S, E, FPM, ca, cls in samp:
    cr = reads_at(X, chrom, S, E)
    if cr < THR_ABS: continue
    thr = max(THR_ABS, THR_REL * cr)
    bed = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False); bed.write(f"{chrom}\t{S}\t{E}\tl\n"); bed.close()
    rb = 1
    for Y in ORDER:
        if Y == X: continue
        out = subprocess.run(f"singularity exec --bind /mnt {SIF} halLiftover {HAL} {X} {bed.name} {Y} /dev/stdout", shell=True, capture_output=True, text=True).stdout
        frs = [(f[0], int(f[1]), int(f[2])) for f in (l.split("\t") for l in out.splitlines()) if len(f) >= 3]
        if not frs: continue
        ch2 = frs[0][0]; s2 = min(f[1] for f in frs); e2 = max(f[2] for f in frs)
        if reads_at(Y, ch2, s2, e2) >= thr: rb += 1
    os.unlink(bed.name)
    rows.append((X, chrom, S, round(FPM), ca, cls, rb, int(cr)))
    print(f"{X} chr{chrom}:{S} FPM={FPM:.0f} clustersAt={ca}({cls}) carrierReads={cr} -> TIGHT read_breadth={rb}", flush=True)
res = pd.DataFrame(rows, columns=["carrier", "chrom", "start", "FPM", "clusters_at_breadth", "clusters_at_class", "read_breadth_tight", "carrier_reads"])
res.to_csv(f"{PG}/confound_verification_tight.csv", index=False)
cre = res[res.clusters_at_class == "creation"].read_breadth_tight; prop = res[res.clusters_at_class == "propagation"].read_breadth_tight
print(f"\nn={len(res)} | TIGHT (>=10% of carrier reads):")
print(f"clusters_at-CREATION (breadth<=3): tight read_breadth median {cre.median():.0f} (range {cre.min()}-{cre.max()}), n={len(cre)}")
print(f"clusters_at-PROPAGATION (breadth>=10): tight read_breadth median {prop.median():.0f} (range {prop.min()}-{prop.max()}), n={len(prop)}")
print(f"AGREEMENT: creation→read≤3: {(cre<=3).mean()*100:.0f}% | propagation→read≥10: {(prop>=10).mean()*100:.0f}% → confound {'REAL' if (cre.median()<=3 and prop.median()>=10) else 'NEEDS REVIEW'}")
