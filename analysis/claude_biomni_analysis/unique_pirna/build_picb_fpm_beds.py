#!/usr/bin/env python3
"""PICB-identified piRNA clusters WITH FPM, per strain, in HAL/BAM coordinates (chrN -> N) — the foundation for
using PICB clusters + their pipeline FPM (multimapping-aware) in every figure and for the pangenome cross-strain
comparison. Reads each strain's combined PICB xlsx 'clusters' sheet for all 3 timepoints (same source as
extract_clusters.py / the cluster PAV), keeps all_reads_primary_alignments_FPM + uniq_reads_FPM + strand + timepoint.
Output: cluster_pav/{X}.clusters_fpm.bed (chrom, start0, end, all_primary_FPM, uniq_FPM, strand, tp). Coordinates are
identical to the BAM/HAL system (verified: CAST chr1:124,197,401 PICB cluster == back-lift of the GRCm39 wild-trio locus)."""
import os, glob
import pandas as pd
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
RES = f"{ROOT}/results/picb_result_combined"; OUT = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"
STRAINS = ["WSB_EiJ", "CAST_EiJ", "BALB_cJ", "C57BL_6NJ", "A_J", "DBA_2J", "NOD_ShiLtJ", "SPRET_EiJ", "AKR_J",
           "C3H_HeJ", "CBA_J", "PWK_PhJ", "NZO_HlLtJ", "FVB_NJ", "129S1_SvImJ", "LP_J"]
TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]
for X in STRAINS:
    rows = []
    for tp in TPS:
        f = f"{RES}/{X}/{X}-{tp}.combined.xlsx"
        if not os.path.exists(f): continue
        d = pd.read_excel(f, sheet_name="clusters")
        for _, r in d.iterrows():
            c = str(r["seqnames"]); c = c[3:] if c.lower().startswith("chr") else c
            rows.append((c, int(r["start"]) - 1, int(r["end"]),
                         round(float(r["all_reads_primary_alignments_FPM"]), 3), round(float(r["uniq_reads_FPM"]), 3),
                         str(r["strand"]), tp))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "all_primary_FPM", "uniq_FPM", "strand", "tp"]).sort_values(["chrom", "start"])
    df.to_csv(f"{OUT}/{X}.clusters_fpm.bed", sep="\t", header=False, index=False)
    print(f"[{X}] {len(df)} PICB clusters (3 tp) -> {X}.clusters_fpm.bed  (max all_primary_FPM={df.all_primary_FPM.max():.0f})")
print("done: PICB clusters+FPM BEDs in HAL coords")
