#!/usr/bin/env python3
"""CANONICAL driver for the strain-private TE-source-locus catalogue (theme 09: Fig_locus_full_*).
Renders each locus in the tracked loci list `09_TE_driven_evolution/data/te_source_loci_catalogue.tsv`
(cols: slug, strain, chrom, start, end, te_label, te_strand) via make_locus_full.py (reconciled 25-32 coverage
caption + fixes). The loci were the delivered curated example set; because find_example_locus's live candidate set
has since drifted (0/19 reproducible), the EXACT coords were recovered from the delivered figures' own titles and
frozen into that TSV (see render_delivered_locus_full.py for the recovery). This driver is now the reproducible
source of truth. Usage: render_te_source_catalogue.py [slug ...]  (restrict to given slugs)."""
import subprocess, sys, os
import pandas as pd
U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
PY = "/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python"
TSV = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/09_TE_driven_evolution/data/te_source_loci_catalogue.tsv"
ONLY = set(sys.argv[1:])
d = pd.read_csv(TSV, sep="\t")
if ONLY: d = d[d.slug.isin(ONLY)]
print(f"rendering {len(d)} TE-source loci from {os.path.basename(TSV)}:", flush=True)
ok = 0
for _, r in d.iterrows():
    cp = f"{r.strain}#1#{r.chrom}"
    print(f"=== {r.slug}  {cp}:{int(r.start):,}-{int(r.end):,}  TE={r.te_label}({r.te_strand})  {r.strain} ===", flush=True)
    p = subprocess.run([PY, f"{U}/make_locus_full.py", r.slug, cp, str(int(r.start)), str(int(r.end)),
                        str(r.te_label), str(r.te_strand), r.strain], capture_output=True, text=True, cwd=U)
    msg = "\n".join(l for l in (p.stdout + p.stderr).splitlines() if any(k in l for k in ("wrote", "Error", "Traceback")))
    print(f"  {msg or 'NO OUTPUT (rc=%d)' % p.returncode}", flush=True); ok += (p.returncode == 0)
print(f"TE_SOURCE_LOCI_DONE — {ok}/{len(d)}", flush=True)
