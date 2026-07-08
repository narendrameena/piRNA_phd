#!/usr/bin/env python3
"""Re-render the DELIVERED strain-private TE-source-locus catalogue (theme 09: Fig_locus_full_*) with the reconciled
code (25-32 coverage caption + fixes), PRESERVING the exact delivered loci. The exact coords are recovered from each
delivered figure's OWN title (chrom:start-end, rendered by make_locus_full) — NOT from find_example_locus, whose
candidate set has since changed. TELAB is read from the SVG title; the TE-strand (load-bearing: sets antisense-to-TE
colouring) is re-derived from {strain}.TE_stranded.bed as the largest-overlap TE at the locus (as find_example_locus did).
Strain = CAST_EiJ for 'CAST_'-prefixed slugs, else SPRET_EiJ (make_locus_full's default). Usage: [slug ...] to restrict."""
import re, glob, os, subprocess, sys
import pandas as pd
U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
SA = f"{U}/sense_antisense"; PY = "/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/python"
ARCH = f"{U}/pangenome_te/_archive_pilot_locus_full"
_TE = {}


def te_strand_and_fam(strain, chrom_pansn, S, E):
    """largest-overlap TE at the locus -> (strand, 'class/family') from {strain}.TE_stranded.bed."""
    if strain not in _TE:
        f = f"{SA}/{strain}.TE_stranded.bed"
        _TE[strain] = pd.read_csv(f, sep="\t", header=None, names=["c", "s", "e", "fam", "sc", "st"]) if os.path.exists(f) else None
    te = _TE[strain]
    if te is None: return "+", None
    sub = te[(te.c == chrom_pansn) & (te.s < E) & (te.e > S)]
    if not len(sub): return "+", None
    row = sub.iloc[(sub.e.clip(upper=E) - sub.s.clip(lower=S)).values.argmax()]
    return row.st, str(row.fam).split("|")[-1]


ONLY = set(sys.argv[1:])
todo = []
for svg in sorted(glob.glob(f"{ARCH}/Fig_locus_full_*.svg")):
    slug = os.path.basename(svg).replace("Fig_locus_full_", "").replace(".svg", "")
    if ONLY and slug not in ONLY: continue
    txt = re.sub(r"<[^>]+>", "", open(svg, errors="ignore").read())
    m = re.search(r"(chr[0-9XY]+):([0-9,]{5,})[-–]([0-9,]{5,})", txt)
    if not m: print(f"  SKIP {slug}: no coord in SVG", flush=True); continue
    chrom = m.group(1); S = int(m.group(2).replace(",", "")); E = int(m.group(3).replace(",", ""))
    mt = re.search(r"Strain-private\s+(\S+)\s", txt); telab = mt.group(1) if mt else None
    strain = "CAST_EiJ" if slug.startswith("CAST_") else "SPRET_EiJ"
    cp = f"{strain}#1#{chrom}"
    test, fam = te_strand_and_fam(strain, cp, S, E)
    if telab is None: telab = fam or "TE"                      # fall back to the annotated family if the title had no TELAB
    todo.append((slug, cp, S, E, telab, test, strain))

print(f"re-rendering {len(todo)} delivered loci (exact coords recovered from the figures):", flush=True)
ok = 0
for slug, cp, S, E, telab, test, strain in todo:
    print(f"=== {slug}  {cp}:{S:,}-{E:,}  TE={telab}({test})  {strain} ===", flush=True)
    p = subprocess.run([PY, f"{U}/make_locus_full.py", slug, cp, str(S), str(E), telab, test, strain],
                       capture_output=True, text=True, cwd=U)
    msg = "\n".join(l for l in (p.stdout + p.stderr).splitlines() if any(k in l for k in ("wrote", "Error", "Traceback")))
    print(f"  {msg or 'NO OUTPUT (rc=%d)' % p.returncode}", flush=True); ok += (p.returncode == 0)
print(f"DELIVERED_RERENDER_DONE — {ok}/{len(todo)}", flush=True)
