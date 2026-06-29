#!/usr/bin/env python3
"""THEME 23/22 — regenerate the 3 FEATURED non-reference locus figures (skipped by gen_nonref_loci_batch.py via the
FEATURED set) and the CONSERVED WSB chr7 liftover-break figure (theme 22), with their curated output names. This
documents the previously ad-hoc command-line calls so every locus figure AND its per-figure SourceData table is
reproducible from a committed driver. Coords for the featured 3 come from shared_subset_loci.csv (rep_cid); the
conserved WSB chr7 cluster is chr7:3,216,400-3,298,600 (31,188 FPM, ~82 kb; minimap2 -> GRCm39 chr7:6.57 Mb).
Run with the biomni_e1 python (pysam). Usage: gen_featured_loci.py"""
import subprocess, sys, pandas as pd
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; PY = sys.executable
T22 = f"{ROOT}/figures/analysis_figures/22_odgi_inject_cluster_pav"
T23 = f"{ROOT}/figures/analysis_figures/23_nonref_locus_atlas"
T24 = f"{ROOT}/figures/analysis_figures/24_strain_private_te_loci"
d = pd.read_csv(f"{T22}/data/shared_subset_loci.csv").set_index("rep_cid")
FEAT = [("WSB_EiJ|10683", "Fig_nonref_locus_WSB_chr4_L1"),
        ("FVB_NJ|1226", "Fig_nonref_locus_FVB_chr1_ERVK"),
        ("SPRET_EiJ|14463", "Fig_nonref_locus_SPRET_chr6_SINEB4")]
for cid, nm in FEAT:
    r = d.loc[cid]
    title = f"{r.strain.replace('_','/')} chr{r.chrom} — a {r.te_family} piRNA cluster genetically ABSENT from GRCm39 ({r.n_strains}/16 strains · odgi inject+pav)"
    print(f">>> {nm}  ({r.strain} chr{r.chrom}:{int(r.start):,})", flush=True)
    subprocess.run([PY, f"{T23}/code/make_nonref_locus.py", r.strain, str(r.chrom), str(int(r.start)), str(int(r.end)), cid, title, nm])
# CONSERVED WSB chr7 (make_conserved_locus.py, theme 22) — a major conserved cluster mis-called by halLiftover
ctitle = "WSB/EiJ chr7 — a CONSERVED, highly-expressed pachytene cluster mis-called 'non-reference' by halLiftover (present in all 16 strains + GRCm39 by minimap2, at strain-specific coordinates)"
print(">>> Fig_conserved_locus_liftover_break_WSB_chr7", flush=True)
subprocess.run([PY, f"{T22}/code/make_conserved_locus.py", "WSB_EiJ", "7", "3216400", "3298600",
                f"{T24}/data/wsb_conserved.tsv", ctitle, "Fig_conserved_locus_liftover_break_WSB_chr7"])
print("=== FEATURED + CONSERVED DONE ===", flush=True)
