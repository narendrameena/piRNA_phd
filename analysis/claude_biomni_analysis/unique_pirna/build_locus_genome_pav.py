#!/usr/bin/env python3
"""Genome PRESENCE/ABSENCE (PAV) of each cluster-PAV example locus in each strain's genome, via the pangenome HAL
(halLiftover GRCm39 -> strain). present = the GRCm39 locus sequence aligns into the strain's genome (>=20% of the
window covered). Distinguishes genetic LOSS (absent) from regulatory SILENCING (present, but no PICB cluster).
One-time pangenome pass (not per-figure liftover). Output cluster_pav/locus_genome_pav.tsv."""
import os, csv, subprocess
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; CP = f"{U}/cluster_pav"
SIF = f"{ROOT}/cactus_v2.9.3.sif"; HAL = f"{ROOT}/results/pangenome/output/mouse_17strain_pangenome.full.hal"
STRAINS = ["C57BL_6NJ", "BALB_cJ", "A_J", "FVB_NJ", "C3H_HeJ", "LP_J", "129S1_SvImJ", "DBA_2J", "AKR_J", "CBA_J",
           "NZO_HlLtJ", "NOD_ShiLtJ", "WSB_EiJ", "CAST_EiJ", "PWK_PhJ", "SPRET_EiJ"]
loci = []  # (g39c, g39s, g39e, label)
for r in csv.DictReader(open(f"{U}/pangenome_te/SourceData_cluster_pav_examples.csv")):
    loci.append((r["chrom"], int(r["start"]), int(r["end"]), r["label"]))
inb = f"{CP}/.pav_loci_g39.bed"
with open(inb, "w") as o:
    for i, (c, s, e, lab) in enumerate(loci):
        o.write(f"{c}\t{s}\t{e}\tL{i}\n")
rows = [("g39_chrom", "g39_start", "g39_end", "strain", "present", "frac_covered", "label")]
for X in STRAINS:
    outb = f"{CP}/.pav_{X}_own.bed"
    subprocess.run(f"singularity exec --bind /mnt {SIF} halLiftover {HAL} GRCm39 {inb} {X} {outb}", shell=True, capture_output=True, text=True)
    cov = {i: 0 for i in range(len(loci))}
    if os.path.exists(outb):
        for ln in open(outb):
            f = ln.rstrip("\n").split("\t")
            if len(f) < 4: continue
            i = int(f[3][1:]); cov[i] += int(f[2]) - int(f[1])
        os.remove(outb)
    for i, (c, s, e, lab) in enumerate(loci):
        frac = cov[i] / max(1, e - s); rows.append((c, s, e, X, int(frac >= 0.2), round(frac, 3), lab))
    print(f"[{X}] present {sum(1 for i in cov if cov[i]/max(1,loci[i][2]-loci[i][1])>=0.2)}/{len(loci)} loci", flush=True)
os.remove(inb)
with open(f"{CP}/locus_genome_pav.tsv", "w") as o:
    for r in rows: o.write("\t".join(map(str, r)) + "\n")
print(f"wrote locus_genome_pav.tsv: {len(rows)-1} rows ({len(loci)} loci x {len(STRAINS)} strains)")
