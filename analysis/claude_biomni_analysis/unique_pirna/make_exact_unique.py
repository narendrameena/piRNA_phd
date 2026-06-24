#!/usr/bin/env python3
"""THEME 19 — EXACT-SEQUENCE uniqueness (user-directed 2026-06-23).
A piRNA is unique to strain X at stage T iff its EXACT sequence is not expressed in another strain at T; 1-3 SNP
variants are treated as DIFFERENT sequences and are KEPT as unique (NO SNP-refinement). This is the Theme-18
within-tp scheme WITHOUT the D4 SNP step: candidates classed "SNP-variant" (a 1-3 mm allele of a same-stage-
expressed sequence elsewhere) revert to conserved-but-silent (their EXACT sequence is still strain-specific).
Genuinely-unique (exact) = conserved-but-silent + strain-private + stage-shifted. Compared against the SNP-aware
Theme-18 set (SNP-variants excluded). Output: theme-19 data/exact_stagepeak_classified.csv.gz with
klass_exact + was_snp_variant flag."""
import pandas as pd
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/19_exact_vs_snp_uniqueness/data"
SNP="SNP-variant (1-3mm, same stage)"; CBS="unique: conserved-but-silent"; SP="unique: strain-private locus"; SH="unique: stage-shifted (heterochronic)"
GU=[CBS,SP,SH]
d=pd.read_csv(f"{U}/deseq16_lenfilt/deseq_stagepeak_classified.csv.gz")
d["was_snp_variant"]=d.klass==SNP
d["klass_exact"]=d.klass.where(d.klass!=SNP, CBS)        # exact-only: drop D4 -> SNP-variant becomes conserved-but-silent
d.to_csv(f"{OUT}/exact_stagepeak_classified.csv.gz",index=False)
na=d.klass.isin(GU).sum(); ex=d.klass_exact.isin(GU).sum()
print(f"SNP-aware genuinely-unique (Theme 18): {na:,}")
print(f"EXACT-sequence genuinely-unique (Theme 19): {ex:,}  (+{d.was_snp_variant.sum():,} SNP-alleles kept as unique)")
print(f"  => exact inflates unique by {100*(ex-na)/na:.0f}%")
print("\nexact klass per timepoint:")
print(d.groupby(["timepoint","klass_exact"]).size().to_string())
