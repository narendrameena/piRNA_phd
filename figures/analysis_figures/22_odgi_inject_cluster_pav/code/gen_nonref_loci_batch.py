#!/usr/bin/env python3
"""THEME 22 — batch-generate the per-locus figure (make_nonref_locus.py) for EVERY mapped non-reference shared-subset
locus (those present in graph_check_pav3.tsv, so the graph-PAV panel A works). Each figure: pangenome x timepoint FPM
(strand-split) + graph-coverage strip (how partial) + carrier sRNA/TE/gene tracks + base resolution.
Usage: gen_nonref_loci_batch.py   (run with the biomni_e1 python that has pysam)"""
import pandas as pd, subprocess, re, sys
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; T22 = f"{ROOT}/figures/analysis_figures/22_odgi_inject_cluster_pav"
PY = sys.executable
d = pd.read_csv(f"{T22}/data/shared_subset_loci.csv"); gp = set(pd.read_csv(f"{T22}/data/graph_check_pav3.tsv", sep="\t").name)
FEATURED = {"WSB_EiJ|10683", "FVB_NJ|1226", "SPRET_EiJ|14463"}   # already made with curated names — skip
sub = d[d.rep_cid.isin(gp)].sort_values("uniqFPM", ascending=False)
def teshort(t): return re.sub(r"[^A-Za-z0-9]", "", str(t)) if pd.notna(t) else "nonTE"
seen = set()
for _, r in sub.iterrows():
    if r.rep_cid in FEATURED: continue
    base = f"{r.strain}_chr{r.chrom}_{teshort(r.te_family)}"
    if base in seen: base = f"{r.strain}_chr{r.chrom}_{int(r.start)//1000}kb_{teshort(r.te_family)}"
    seen.add(base); name = f"Fig_nonref_locus_{base}"
    title = f"{r.strain.replace('_','/')} chr{r.chrom} — a {r.te_family} piRNA cluster genetically ABSENT from GRCm39 ({r.n_strains}/16 strains · odgi inject+pav)"
    print(f">>> {name}  ({r.strain} chr{r.chrom}:{int(r.start):,} {r.uniqFPM:.1f} FPM)", flush=True)
    subprocess.run([PY, f"{T22}/code/make_nonref_locus.py", r.strain, str(r.chrom), str(int(r.start)), str(int(r.end)), r.rep_cid, title, name])
print("=== BATCH DONE ===", flush=True)
