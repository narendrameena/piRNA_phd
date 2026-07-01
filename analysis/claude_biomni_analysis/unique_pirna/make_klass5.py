#!/usr/bin/env python3
"""STEP 5 (merge) — build the ADOPTED 5-class unique-piRNA set from committed Step-4 determinants.

Closes the reproducibility gap flagged in DATA_INTEGRITY_AUDIT_2026-06-19.md: previously the klass4/klass5
columns and the >=2-read `final_classified_clean*.csv.gz` files had no committed producer. This script
reconstructs them deterministically from already-committed inputs. The logic was reverse-verified on
2026-06-19 to reproduce the existing files at 100% (klass4 0/451,972 mismatch; klass5 low-quality split
60,857/60,857; >=2-read filter 451,972/451,972).

Inputs (all committed):
  unique16/final_classified.csv.gz            <- classify_unique16_locus.py  (3-class `klass`)
  unique16/snp_variant_refinement.csv         <- classify_step416.py / step4_16 (1-3mm SNP-variant set, by cand_id)
  cand_self16/{X}.cand_self16.bam (16)        <- cand_self16.sh  (candidate -> OWN genome, --outFilterMismatchNmax 0,
                                                 so ANY alignment == a clean mm0 own-genome locus)
  edger16/{tp}.strain_specific_DA_2read.csv.gz <- the reproducible >=2-read edgeR DA candidate set

Refinement:
  klass4 = klass, but conserved-but-silent candidates that are a 1-3mm SNP-variant of an expressed allele
           at the orthologous locus elsewhere (in snp_variant_refinement) -> "SNP-variant (1-3mm)".
  klass5 = klass4, but strain-private candidates with NO mm0 own-genome locus (not mapped in cand_self16)
           -> "low-quality: no mm0 own-genome locus".

Outputs:
  unique16/final_classified_clean.csv.gz        (ALL candidates, 5-class klass5)
  unique16/final_classified_clean_2read.csv.gz  (ADOPTED: rows whose (strain,timepoint,sequence) is in the
                                                 >=2-read DA candidate set)

Usage: make_klass5.py [OUTDIR]   (OUTDIR default = unique16/; pass a temp dir to verify without overwriting)
"""
import sys, glob, os
import pandas as pd, pysam
U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
IN = f"{U}/unique16"
OUTDIR = sys.argv[1] if len(sys.argv) > 1 else IN
os.makedirs(OUTDIR, exist_ok=True)
CBS = "unique: conserved-but-silent"; SP = "unique: strain-private locus"
KL4_SNP = "SNP-variant (1-3mm)"; KL5_LQ = "low-quality: no mm0 own-genome locus"
COLS = ["sequence", "strain", "timepoint", "length", "logFC", "FDR", "n_other_exact",
        "other_exact_strains", "expr_class", "klass", "homolog_strains", "klass4", "klass5"]

d = pd.read_csv(f"{IN}/final_classified.csv.gz")
d["cand_id"] = d.strain + "|" + d.timepoint + "|" + d.sequence

# klass4: SNP-variant refinement of conserved-but-silent
snp = set(pd.read_csv(f"{IN}/snp_variant_refinement.csv", usecols=["cand_id"]).cand_id)
d["klass4"] = d.klass.where(~((d.klass == CBS) & (d.cand_id.isin(snp))), KL4_SNP)

# klass5: low-quality split of strain-private (no mm0 own-genome locus = unmapped in cand_self16)
mapped = set()
for bam in sorted(glob.glob(f"{U}/cand_self16/*.cand_self16.bam")):
    b = pysam.AlignmentFile(bam, "rb")
    for a in b.fetch(until_eof=True):
        if not a.is_unmapped: mapped.add(a.query_name)
    b.close()
d["klass5"] = d.klass4.where(~((d.klass4 == SP) & (~d.cand_id.isin(mapped))), KL5_LQ)

d[COLS].to_csv(f"{OUTDIR}/final_classified_clean.csv.gz", index=False)

# ADOPTED >=2-read subset: membership in the reproducible >=2-read DA candidate set
da = set()
for tp in ["16.5dpc", "12.5dpp", "20.5dpp"]:
    t = pd.read_csv(f"{U}/edger16/{tp}.strain_specific_DA_2read.csv.gz", usecols=["strain", "sequence"])
    da |= set(t.strain + "|" + tp + "|" + t.sequence)
d2 = d[d.cand_id.isin(da)]
d2[COLS].to_csv(f"{OUTDIR}/final_classified_clean_2read.csv.gz", index=False)

print(f"wrote {OUTDIR}/final_classified_clean.csv.gz ({len(d):,} rows) and final_classified_clean_2read.csv.gz ({len(d2):,} rows)")
print("klass5 (clean):", d.klass5.value_counts().to_dict())
print("klass5 (2read):", d2.klass5.value_counts().to_dict())
