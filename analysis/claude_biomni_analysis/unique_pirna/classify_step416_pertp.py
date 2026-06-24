#!/usr/bin/env python3
"""WITHIN-TIMEPOINT D4 (SNP-variant refinement), per home strain X.

Mirrors classify_step416.py, but the expression test is made WITHIN the candidate's OWN developmental
timepoint instead of pooled across all three. A candidate (X, tp, seq) is a SNP-variant iff it aligns
(1-3 mm) to another strain Y's genome AND that Y-genome allele is EXPRESSED IN Y AT THE SAME tp
(pools_pertp/{Y}_{tp}.pool.txt.gz, >=2 reads, 25-32 nt at that tp). Within-tp is STRICTER than cross-tp
(same-stage expression only) -> fewer SNP-variants -> more genuinely-unique conserved-but-silent loci.

Reuses the EXISTING step4_16 cand_to_Y BAMs (mm<=3, --outSAMattributes All so MD/get_reference_sequence
works) -> NO STAR re-run. candidates16.tsv.gz `timepoints` may be multi-tp (comma-separated); each tp is
tested separately, so a candidate can be a SNP-variant at one stage but unique at another.

Memory-light: per Y, collect the small set of reference alleles seen in the BAM, then STREAM Y's per-tp
pool files (constant memory) instead of loading the (up to ~24M-sequence) pools into a set.

Output: step4_16/{X}.snp_withintp.csv.gz  cols = cand_id,home,variant_strain,home_seq,Y_allele,mm,tp
        (cand_id = X|tp|seq, matching unique16/snp_variant_refinement.csv plus a `tp` column).
Usage: classify_step416_pertp.py <strain>"""
import sys, gzip, pysam, pandas as pd
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; S4 = f"{U}/step4_16"; PP = f"{U}/pools_pertp"
ALL = ["C57BL_6NJ", "BALB_cJ", "A_J", "FVB_NJ", "C3H_HeJ", "LP_J", "129S1_SvImJ", "DBA_2J", "AKR_J",
       "CBA_J", "NZO_HlLtJ", "NOD_ShiLtJ", "WSB_EiJ", "CAST_EiJ", "PWK_PhJ", "SPRET_EiJ"]
TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]
X = sys.argv[1]; OTH = [s for s in ALL if s != X]
cand = pd.read_csv(f"{S4}/{X}.candidates16.tsv.gz")
seqs = cand.sequence.tolist(); cand_tps = [set(str(t).split(",")) for t in cand.timepoints]
id2idx = {f"c{i}": i for i in range(len(cand))}
rows = []
for Y in OTH:
    # pass 1: from the BAM, collect each 1-3 mm candidate alignment to Y (candidate idx, Y allele, mm)
    recs = []; need = set()
    bam = pysam.AlignmentFile(f"{S4}/{X}.cand_to_{Y}.Aligned.sortedByCoord.out.bam", "rb")
    for a in bam.fetch(until_eof=True):
        if a.is_unmapped: continue
        i = id2idx.get(a.query_name)
        if i is None: continue
        nm = int(a.get_tag("NM"))
        if not (1 <= nm <= 3): continue          # only 1-3 mm = SNP-variant (mm0 exact -> handled by D1)
        try: refseq = a.get_reference_sequence().upper()
        except Exception: continue
        recs.append((i, refseq, nm)); need.add(refseq)
    bam.close()
    if not recs:
        print(f"[{X}] {Y}: 0 1-3mm aln", flush=True); continue
    # pass 2: stream Y's per-tp pools, retain only the (small set of) needed alleles present at each tp
    present_at = {tp: set() for tp in TPS}
    for tp in TPS:
        with gzip.open(f"{PP}/{Y}_{tp}.pool.txt.gz", "rt") as fh:
            for line in fh:
                s = line.rstrip("\n")
                if s in need: present_at[tp].add(s)
    n0 = len(rows)
    for i, refseq, nm in recs:
        for tp in cand_tps[i]:                    # candidate's own stage(s)
            if refseq in present_at.get(tp, ()):  # Y allele expressed in Y AT that stage
                rows.append((f"{X}|{tp}|{seqs[i]}", X, Y, seqs[i], refseq, nm, tp))
    print(f"[{X}] {Y}: {len(recs)} 1-3mm aln -> {len(rows)-n0} within-tp snp records", flush=True)
out = pd.DataFrame(rows, columns=["cand_id", "home", "variant_strain", "home_seq", "Y_allele", "mm", "tp"])
out.to_csv(f"{S4}/{X}.snp_withintp.csv.gz", index=False, compression="gzip")
print(f"[{X}] DONE: {len(out):,} within-tp snp records, {out.cand_id.nunique():,} distinct cand_id", flush=True)
