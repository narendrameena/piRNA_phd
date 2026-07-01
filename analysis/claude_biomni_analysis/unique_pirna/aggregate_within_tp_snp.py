#!/usr/bin/env python3
"""Aggregate the 16 per-strain within-tp SNP-variant record files (step4_16/{X}.snp_withintp.csv.gz) into
unique16/snp_variant_refinement_withintp.csv, matching the cross-tp file's columns
(cand_id,home,variant_strain,home_seq,Y_allele,mm) so within_tp_rederive.py and the SNP figures consume it
unchanged. Reports the within-tp vs cross-tp SNP-variant cand_id counts (within-tp is stricter -> fewer)."""
import glob, pandas as pd
U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
fs = sorted(glob.glob(f"{U}/step4_16/*.snp_withintp.csv.gz"))
assert len(fs) == 16, f"expected 16 per-strain files, found {len(fs)}"
d = pd.concat([pd.read_csv(f) for f in fs], ignore_index=True)
COLS = ["cand_id", "home", "variant_strain", "home_seq", "Y_allele", "mm"]   # match cross-tp file (tp is in cand_id)
d[COLS].to_csv(f"{U}/unique16/snp_variant_refinement_withintp.csv", index=False)
cross = pd.read_csv(f"{U}/unique16/snp_variant_refinement.csv", usecols=["cand_id"])
wn, cn = d.cand_id.nunique(), cross.cand_id.nunique()
print(f"within-tp: {len(d):,} records, {wn:,} distinct cand_id")
print(f"cross-tp : {len(cross):,} records, {cn:,} distinct cand_id")
print(f"SNP-variant cand_ids dropped by within-tp constraint: {cn-wn:,} ({100*(cn-wn)/cn:.1f}% of cross-tp) "
      "-> these re-enter the unique pipeline (some become genuinely-unique conserved-but-silent)")
