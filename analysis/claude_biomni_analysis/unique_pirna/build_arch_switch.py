#!/usr/bin/env python3
"""Catalogue PICB piRNA CLUSTERS whose genomic-strand ARCHITECTURE (uni- vs dual-strand) SWITCHES across
developmental timepoints. For every strain, take the PICB-combined clusters (clusters_fpm.bed, own genome),
find clusters present in >=2 timepoints (coordinate overlap) with FPM>=15 in at least one, and for each present
timepoint compute the minus-strand fraction of own-genome sRNA reads (24-32 nt) -> architecture (dual if
minor-strand fraction >20%, else uni). Records clusters where the call CHANGES across timepoints (>=50 reads/tp).
Output: cluster_pav/arch_switch_catalogue.tsv. (Architecture = genomic strand, NOT sense/antisense.)"""
import os, glob, subprocess, pysam, numpy as np, pandas as pd
from collections import defaultdict
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; CP = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"
BT = "/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/bedtools"
STRAINS = ["WSB_EiJ", "CAST_EiJ", "BALB_cJ", "C57BL_6NJ", "A_J", "DBA_2J", "NOD_ShiLtJ", "SPRET_EiJ", "AKR_J", "C3H_HeJ", "CBA_J", "PWK_PhJ", "NZO_HlLtJ", "FVB_NJ", "129S1_SvImJ", "LP_J"]
TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]; TPLAB = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}
def arch_at(X, chrom, ps, pe, tp):
    p = m = 0; bamc = f"{X}#1#chr{chrom}"
    for r in (1, 2, 3):
        b = f"{ROOT}/results/STAR_srna_strain_wise/{X}/{X}-{tp}.{r}/Aligned.sortedByCoord.out.bam"
        if not os.path.exists(b): continue
        try:
            bam = pysam.AlignmentFile(b, "rb")
            for a in bam.fetch(bamc, ps, pe):
                if a.is_unmapped or not a.query_sequence: continue
                if 24 <= a.reference_end - a.reference_start <= 32: m += a.is_reverse; p += (not a.is_reverse)
            bam.close()
        except (OSError, ValueError):
            try: bam.close()
            except Exception: pass
    tot = p + m
    if tot < 50: return None
    mf = m / tot; return ("dual" if min(mf, 1 - mf) > 0.2 else "uni", round(mf, 3), tot)
rows = []
for X in STRAINS:
    fb = f"{CP}/{X}.clusters_fpm.bed"
    if not os.path.exists(fb): continue
    d = pd.read_csv(fb, sep="\t", header=None, names=["chrom", "start", "end", "allFPM", "uniqFPM", "strand", "tp"])
    d = d[d.chrom.astype(str).str.match(r"^(\d+|X)$")].copy()
    # union intervals across all tp (merge), then which tp present + max FPM in each
    raw = f"/tmp/.arch_{X}.bed"; d[["chrom", "start", "end", "tp", "allFPM"]].sort_values(["chrom", "start"]).to_csv(raw, sep="\t", header=False, index=False)
    uni = subprocess.run(f"sort -k1,1 -k2,2n {raw} | {BT} merge -d 200 -c 4,5 -o distinct,max", shell=True, capture_output=True, text=True).stdout
    ncand = nsw = 0
    for ln in uni.splitlines():
        f = ln.split("\t")
        tps = set(f[3].split(",")); maxfpm = float(f[4])
        if len(tps) < 2 or maxfpm < 15: continue
        ncand += 1; chrom = f[0]; ps = int(f[1]); pe = int(f[2])
        if pe - ps > 50000: continue
        archs = {tp: arch_at(X, chrom, ps, pe, tp) for tp in TPS if tp in tps}
        archs = {k: v for k, v in archs.items() if v}
        cats = {v[0] for v in archs.values()}
        if len(cats) > 1:
            nsw += 1
            rows.append(dict(strain=X, chrom=chrom, start=ps, end=pe, maxFPM=round(maxfpm, 1),
                             **{f"{TPLAB[tp]}_arch": archs[tp][0] if tp in archs else "" for tp in TPS},
                             **{f"{TPLAB[tp]}_minusfrac": archs[tp][1] if tp in archs else "" for tp in TPS}))
    print(f"[{X}] multi-tp candidates(FPM>=15)={ncand}; architecture SWITCHES={nsw}", flush=True)
T = pd.DataFrame(rows); T.to_csv(f"{CP}/arch_switch_catalogue.tsv", sep="\t", index=False)
print(f"\nTotal architecture-switch clusters: {len(T)} across {T.strain.nunique() if len(T) else 0} strains")
if len(T):
    pat = T.apply(lambda r: "→".join(r[f"{TPLAB[tp]}_arch"][0].upper() if r[f"{TPLAB[tp]}_arch"] else "·" for tp in TPS), axis=1)
    print("switch patterns (E16.5→P12.5→P20.5):"); print(pat.value_counts().to_string())
