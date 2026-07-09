#!/usr/bin/env python3
"""Committed producer for unique16/snp_variant_refinement.csv (the 1-3mm SNP-variant refinement set that
make_klass5.py consumes). Reconstructs the refinement TABLE from the same committed determinants as
classify_step416.py, additionally recording the variant strain / allele and expanding across timepoints.

Per home strain X: align each candidate to every OTHER strain Y (committed step4_16/{X}.cand_to_{Y}.bam,
STAR mm<=3), and where the aligned REFERENCE sequence is in Y's expressed pool (unique16/pools/{Y}.pool.txt.gz)
record the minimum-mismatch pool-matching allele. variant_strain / Y_allele = the FIRST Y (in ALL order) that
reaches the running minimum mismatch (== classify_step416.py's strict-improvement rule). A candidate is a
SNP-variant iff that minimum is in 1-3 (mm==0 would be expressed-elsewhere-exact). One output row per
tp-resolved cand_id (X|tp|seq), expanded over the candidate's `timepoints`.

CAVEAT (inherited from classify_step416.py:27-31): get_reference_sequence() is the FORWARD-strand genomic
sequence, and Y's pool is strand-aware, so a minus-strand-only expressed homolog can be missed -> mildly
under-calls SNP-variant. Kept identical to the determinant so this table reproduces that logic.

VERIFIED OVERLAP (2026-07): this faithful re-run of the committed determinants reproduces only ~50% of the
DELIVERED unique16/snp_variant_refinement.csv (pooled 49.6% across 5 strains: 129S1 66.7%, CAST 49.6%, FVB
51.2%, AKR 35.0%, C3H 38.9%; it DOES match the committed step4_classified16.csv.gz exactly, so the logic here
is faithful). The committed inputs (candidates/BAMs/pools, Jun 11) PREDATE the delivered file (Jun 13), so the
~50% gap is a LOGIC difference in the original (uncommitted, likely early tp-specific-pool best-Y) producer,
NOT input drift -- the exact original is unrecoverable, and the committed within-tp chain
(classify_step416_pertp.py) reproduces the delivered file no better (63.5% on 129S1). THEREFORE
running this does NOT regenerate the delivered klass5 SNP-variant set; it writes a SEPARATE *.RECON file and
never overwrites the delivered one. Use only if you accept re-deriving every downstream klass5 number.

Usage:  build_snp_variant_refinement.py <strain> [OUTDIR]   # writes {OUTDIR}/{X}.snp_refinement.csv.gz
        build_snp_variant_refinement.py --all [OUTDIR]      # runs all 16 then concatenates -> snp_variant_refinement.RECON.csv
"""
import sys, gzip, os, glob
import numpy as np, pandas as pd, pysam
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
S4 = f"{U}/step4_16"; POOLS = f"{U}/unique16/pools"
ALL = ["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J",
       "NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
COLS = ["cand_id","home","variant_strain","home_seq","Y_allele","mm"]


def load_pool(s):
    st = set()
    with gzip.open(f"{POOLS}/{s}.pool.txt.gz","rt") as fh:
        for line in fh: st.add(line.rstrip("\n"))
    return st


def one_strain(X):
    OTH = [s for s in ALL if s != X]
    cand = pd.read_csv(f"{S4}/{X}.candidates16.tsv.gz")            # cols: id, sequence, timepoints
    N = len(cand); id2i = {cid: i for i, cid in enumerate(cand.id)}
    seqs = cand.sequence.tolist(); tps = cand.timepoints.astype(str).tolist()
    minnm = np.full(N, 99, np.int16); bestY = [None]*N; bestAllele = [None]*N
    for Y in OTH:
        pool = load_pool(Y)
        bam = pysam.AlignmentFile(f"{S4}/{X}.cand_to_{Y}.Aligned.sortedByCoord.out.bam","rb")
        for a in bam.fetch(until_eof=True):
            if a.is_unmapped: continue
            i = id2i.get(a.query_name)
            if i is None: continue
            try: refseq = a.get_reference_sequence().upper()
            except Exception: continue
            if refseq in pool:
                nm = int(a.get_tag("NM"))
                if nm < minnm[i]:                                  # strict improvement == classify_step416.py:31
                    minnm[i] = nm; bestY[i] = Y; bestAllele[i] = refseq
        bam.close(); del pool
    rows = []
    for i in range(N):
        if bestY[i] is None or not (1 <= minnm[i] <= 3): continue  # SNP-variant = pool-match, mm in 1-3
        for tp in str(tps[i]).replace(";", ",").split(","):
            tp = tp.strip()
            if not tp or tp == "nan": continue
            rows.append((f"{X}|{tp}|{seqs[i]}", X, bestY[i], seqs[i], bestAllele[i], int(minnm[i])))
    return pd.DataFrame(rows, columns=COLS)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:]]
    if args and args[0] == "--all":
        outdir = args[1] if len(args) > 1 else U + "/unique16"
        parts = [one_strain(X) for X in ALL]
        out = pd.concat(parts, ignore_index=True)
        out.to_csv(f"{outdir}/snp_variant_refinement.RECON.csv", index=False)
        print(f"wrote {outdir}/snp_variant_refinement.RECON.csv ({len(out):,} rows, {out.cand_id.nunique():,} unique cand_id)")
    else:
        X = args[0]; outdir = args[1] if len(args) > 1 else "/tmp"
        os.makedirs(outdir, exist_ok=True)
        df = one_strain(X)
        df.to_csv(f"{outdir}/{X}.snp_refinement.csv.gz", index=False, compression="gzip")
        print(f"[{X}] {len(df):,} rows ({df.cand_id.nunique():,} unique cand_id)")
