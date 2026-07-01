#!/usr/bin/env python3
"""Project every strain's PICB-combined piRNA clusters (with FPM + timepoint) into the COMMON GRCm39 pangenome
frame, via the minigraph-cactus HAL (halLiftover strain -> GRCm39) — the SAME pangenome alignment that built the
cluster PAV. This is the pangenome-based cross-strain (AND cross-timepoint) comparison: NOT pairwise liftover of a
locus to each strain, but one projection of all clusters into the shared reference, where strains are then compared
directly. FPM/timepoint/strand are carried through halLiftover in the BED name field. Output:
cluster_pav/picb_pangenome_fpm.tsv = (g39_chrom, start, end, strain, tp, all_primary_FPM, uniq_FPM, strand)."""
import os, glob, subprocess
import pandas as pd
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; CP = f"{U}/cluster_pav"
SIF = f"{ROOT}/cactus_v2.9.3.sif"; HAL = f"{ROOT}/results/pangenome/output/mouse_17strain_pangenome.full.hal"
STRAINS = ["WSB_EiJ", "CAST_EiJ", "BALB_cJ", "C57BL_6NJ", "A_J", "DBA_2J", "NOD_ShiLtJ", "SPRET_EiJ", "AKR_J",
           "C3H_HeJ", "CBA_J", "PWK_PhJ", "NZO_HlLtJ", "FVB_NJ", "129S1_SvImJ", "LP_J"]
allrows = []
for X in STRAINS:
    fb = f"{CP}/{X}.clusters_fpm.bed"
    if not os.path.exists(fb): print(f"[skip] {X} no fpm bed"); continue
    d = pd.read_csv(fb, sep="\t", header=None, names=["chrom", "start", "end", "allFPM", "uniqFPM", "strand", "tp"])
    inb = f"{CP}/.{X}.fpm_liftin.bed"
    with open(inb, "w") as o:                     # name = strain;tp;allFPM;uniqFPM;strand;OWNchrom;OWNstart;OWNend  (carried through halLiftover)
        for _, r in d.iterrows():
            o.write(f"{r.chrom}\t{r.start}\t{r.end}\t{X};{r.tp};{r.allFPM};{r.uniqFPM};{r.strand};{r.chrom};{r.start};{r.end}\n")
    outb = f"{CP}/.{X}.fpm_g39.bed"
    rc = subprocess.run(f"singularity exec --bind /mnt {SIF} halLiftover {HAL} {X} {inb} GRCm39 {outb}", shell=True, capture_output=True, text=True)
    n = 0
    if os.path.exists(outb):
        for ln in open(outb):
            f = ln.rstrip("\n").split("\t")
            if len(f) < 4: continue
            nm = f[3].split(";")
            if len(nm) < 8: continue
            allrows.append((f[0], int(f[1]), int(f[2]), nm[0], nm[1], float(nm[2]), float(nm[3]), nm[4], nm[5], int(nm[6]), int(nm[7]))); n += 1
    os.remove(inb)
    print(f"[{X}] {len(d)} clusters -> {n} GRCm39-projected fragments  (rc={rc.returncode})")
T = pd.DataFrame(allrows, columns=["g39_chrom", "start", "end", "strain", "tp", "all_primary_FPM", "uniq_FPM", "strand", "own_chrom", "own_start", "own_end"])
T = T[T.g39_chrom.astype(str).str.match(r"^(\d+|X|MT)$")].sort_values(["g39_chrom", "start"])
T.to_csv(f"{CP}/picb_pangenome_clusters.tsv", sep="\t", index=False)   # NEW superset (own+g39 coords); old picb_pangenome_fpm.tsv preserved
print(f"wrote picb_pangenome_clusters.tsv: {len(T)} rows; strains={T.strain.nunique()}; tps={sorted(T.tp.unique())}")
