#!/usr/bin/env python3
"""THEME 23 — batch-generate per-locus figures for the top INJECTED strain-PRIVATE TE-insertion piRNA clusters
(theme-22 data/private_te_loci.csv): clusters present in ONE strain only, absent from every other strain AND GRCm39 —
the young, strain-specific TE insertions that found brand-new piRNA clusters (the leading edge). Same figure as the
shared-subset atlas (make_nonref_locus.py); panel A reads the private clusters' graph PAV (appended to
graph_check_pav3.tsv). Run with the biomni_e1 python (pysam). The private signature: carrier coverage ≈ 1.0, every
other strain + GRCm39 ≈ 0."""
import pandas as pd, subprocess, sys
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
T22 = f"{ROOT}/figures/analysis_figures/22_odgi_inject_cluster_pav"; T23 = f"{ROOT}/figures/analysis_figures/23_nonref_locus_atlas"
sel = pd.read_csv(f"{T22}/data/private_te_loci.csv")
for _, r in sel.iterrows():
    name = f"Fig_private_locus_{r.strain}_chr{r.chrom}_{int(r.start)//1000}kb"
    title = (f"{r.strain.replace('_','/')} chr{r.chrom} — a {r.kb}-kb strain-PRIVATE piRNA cluster ({r.fpm:.0f} FPM): "
             f"a young TE insertion present in ONLY this strain, absent from all others and GRCm39 (odgi inject+pav)")
    print(f">>> {name}  ({r.strain} chr{r.chrom}:{int(r.start):,} {r.fpm:.0f} FPM)", flush=True)
    subprocess.run([sys.executable, f"{T23}/code/make_nonref_locus.py", r.strain, str(r.chrom), str(int(r.start)), str(int(r.end)), r.cid, title, name])
print("=== PRIVATE BATCH DONE ===", flush=True)
