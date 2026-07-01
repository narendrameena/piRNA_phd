#!/usr/bin/env python3
"""Render ALL cluster-PAV example loci (from SourceData_cluster_pav_examples.csv) as the THREE locus-figure
types (original / multi-strain / single-strain), each with the per-timepoint READS-BY-STRAND Panel-B update.
Coordinates + labels are read straight from the catalogue CSV (data-driven, not hard-coded). Slugs are explicit
so filenames are stable and match the 3 already done (Col13a1/Dgcr8/Fthl17)."""
import csv, subprocess, sys, os
U = "/mnt/beegfs/scratch/miska/nm667/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
PY = "/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python"
CSV = f"{U}/pangenome_te/SourceData_cluster_pav_examples.csv"
SLUG = {  # catalogue label -> stable filesystem slug (chr1wildtrio/Col13a1/Dgcr8/Fthl17 overwrite the existing figures)
    "near Dars (intergenic)": "chr1wildtrio", "Col13a1": "Col13a1", "Dgcr8 (+NOD)": "Dgcr8",
    "Ccr7": "Ccr7", "Gm20275": "Gm20275", "Fthl17 family (X, germline)": "Fthl17", "Prmt9": "Prmt9",
    "chr15 70kb (Gm53029)": "Gm53029_chr15_70kb", "Hsbp1": "Hsbp1", "Snx29": "Snx29",
    "LP_J 136kb": "LPJ_chr12_136kb", "C57BL/6NJ 109kb": "C57BL6NJ_chr12_109kb",
    "WSB_EiJ 95kb (D430041D05Rik)": "D430041D05Rik_WSB_95kb", "NZO 93kb": "NZO_chr12_93kb",
}
LABEL_OVERRIDE = {  # richer, BioMNI-verified display titles for the loci that already had them (cosmetic only)
    "near Dars (intergenic)": "chr1:128Mb (intergenic, near Dars)", "Col13a1": "Col13a1 (collagen XIII)",
    "Dgcr8 (+NOD)": "Dgcr8 (Microprocessor / miRNA biogenesis)",
    "Fthl17 family (X, germline)": "Fthl17 family (X-linked germline multicopy)",
}
ONLY = set(sys.argv[1:])  # optional: restrict to given slugs (else all)
rows = list(csv.DictReader(open(CSV)))
for r in rows:
    label, c, s, e = r["label"], r["chrom"], r["start"], r["end"]
    slug = SLUG.get(label)
    if slug is None:
        print(f"!! no slug for {label!r} — SKIP", flush=True); continue
    if ONLY and slug not in ONLY: continue
    disp = LABEL_OVERRIDE.get(label, label)   # richer title where available; SLUG lookup still uses raw catalogue label
    print(f"=== {slug}  ({c}:{s}-{e})  [{r['cat']}] {disp} ===", flush=True)
    for script, extra in (("make_pav_locus.py", f"Fig_pav_locus_{slug}"),
                          ("make_pav_locus_multi.py", f"Fig_pav_locus_{slug}_multi")):
        # original/multi: argv = c s e label PATT(dummy) outbase
        p = subprocess.run([PY, f"{U}/{script}", c, s, e, disp, "_", extra], capture_output=True, text=True, cwd=U)
        out = "\n".join(l for l in (p.stdout + p.stderr).splitlines() if ("wrote" in l or "Error" in l or "Traceback" in l or "present strains" in l))
        print(f"  [{script}] {out or 'NO OUTPUT (rc=%d)' % p.returncode}", flush=True)
    # single: argv = c s e label outbase [strain]
    p = subprocess.run([PY, f"{U}/make_pav_locus_single.py", c, s, e, disp, f"Fig_pav_locus_{slug}_single"], capture_output=True, text=True, cwd=U)
    out = "\n".join(l for l in (p.stdout + p.stderr).splitlines() if ("wrote" in l or "Error" in l or "Traceback" in l))
    print(f"  [single] {out or 'NO OUTPUT (rc=%d)' % p.returncode}", flush=True)
print("RENDER_ALL_DONE", flush=True)
