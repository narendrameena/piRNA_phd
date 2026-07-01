#!/usr/bin/env python3
"""Genome PAV for divergence candidates via pangenome HAL (GRCm39 g39 window -> each strain). Writes a supplement
cluster_pav/locus_genome_pav_divergence.tsv (read by pav_clusters.genome_pav) and selects the genome-CONSERVED
ones (present in >=10 strains) -> divergence_loci_final.tsv for make_pav_locus_multi. One bed, 15 liftovers."""
import sys, os, subprocess, tempfile, pandas as pd
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"; CP = f"{U}/cluster_pav"
SIF = f"{ROOT}/cactus_v2.9.3.sif"; HAL = f"{ROOT}/results/pangenome/output/mouse_17strain_pangenome.full.hal"
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]
cand = pd.read_csv(f"{PG}/divergence_candidates.tsv", sep="\t", dtype={"g39_chrom": str})
bed = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False); locinfo = {}
for _, r in cand.iterrows():
    lid = f"{r.g39_chrom}:{int(r.g39_start)}:{int(r.g39_end)}"; bed.write(f"{r.g39_chrom}\t{int(r.g39_start)}\t{int(r.g39_end)}\t{lid}\n"); locinfo[lid] = (r.g39_chrom, int(r.g39_start), int(r.g39_end))
bed.close()
present = {lid: set() for lid in locinfo}
for X in ORDER:
    out = subprocess.run(f"singularity exec --bind /mnt {SIF} halLiftover {HAL} GRCm39 {bed.name} {X} /dev/stdout", shell=True, capture_output=True, text=True).stdout
    agg = {}
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) >= 4: agg[f[3]] = agg.get(f[3], 0) + (int(f[2]) - int(f[1]))
    for lid, (c, s, e) in locinfo.items():
        if agg.get(lid, 0) >= 0.30 * (e - s): present[lid].add(X)
    print(f"[{X}] projected", flush=True)
os.unlink(bed.name)
rows = []
for _, r in cand.iterrows():
    lid = f"{r.g39_chrom}:{int(r.g39_start)}:{int(r.g39_end)}"; expr = set(str(r.expressing).split(","))
    for X in ORDER:
        rows.append((r.g39_chrom, int(r.g39_start), int(r.g39_end), X, (X in present[lid]) or (X in expr)))   # expression proves presence
sup = pd.DataFrame(rows, columns=["g39_chrom", "g39_start", "g39_end", "strain", "present"])
sup.to_csv(f"{CP}/locus_genome_pav_divergence.tsv", sep="\t", index=False)
gn = sup[sup.present].groupby(["g39_chrom", "g39_start"]).size()
cand["genome_n"] = [int(gn.get((r.g39_chrom, int(r.g39_start)), 0)) for _, r in cand.iterrows()]
div = cand[cand.genome_n >= 10].sort_values("maxFPM", ascending=False)
sel = []; seen = {}
for _, r in div.iterrows():
    if seen.get(r.carrier, 0) >= 2: continue
    seen[r.carrier] = seen.get(r.carrier, 0) + 1; sel.append(r)
    if len(sel) >= 12: break
fin = pd.DataFrame(sel)
fin.to_csv(f"{PG}/divergence_loci_final.tsv", sep="\t", index=False)
print(f"\ncandidates={len(cand)} | genome-conserved (>=10): {len(div)} | selected: {len(fin)}")
if len(fin): print(fin[["g39_chrom", "g39_start", "g39_end", "breadth", "genome_n", "carrier", "maxFPM", "expressing"]].to_string(index=False))
