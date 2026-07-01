#!/usr/bin/env python3
"""Classify the DESeq2 stage-peak (27/30 nt, Order-A filter-before) candidates by STRICT WITHIN-TIMEPOINT
uniqueness (user-directed 2026-06-23: uniqueness judged ONLY within a developmental stage — pre-pachytene
E16.5 vs pachytene P20.5 are different machinery, so a piRNA is unique to strain X at stage T iff NO OTHER
strain expresses it AT T). Reuses existing determinants (DESeq2 candidates are a 100% subset of the edgeR set
at strain|tp|seq -> no STAR/halLiftover recompute).

WITHIN-TP scheme (a candidate is penalised ONLY for same-stage evidence elsewhere):
  EE-same-stage : exact seq expressed in another strain AT THE SAME tp           -> NOT unique (shared at stage)
  EE-other-stage: exact seq expressed in another strain ONLY at a DIFFERENT tp   -> UNIQUE: stage-shifted
                  (heterochronic — within its own stage it is strain-specific)
  novel (exact seq in NO other strain's pool, any tp):
     SNP-variant within-tp (1-3mm of an allele expressed AT T elsewhere)         -> NOT unique (strain allele)
     ortholog locus present elsewhere (halLiftover)                              -> UNIQUE: conserved-but-silent
     no ortholog locus elsewhere                                                 -> UNIQUE: strain-private locus
                                                                                    (or low-quality if no mm0 own locus)
Genuinely-unique (within-tp) = strain-private + conserved-but-silent + stage-shifted.
Output: deseq16_lenfilt/deseq_stagepeak_classified.csv.gz"""
import pandas as pd, glob, pysam
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
D=f"{U}/deseq16_lenfilt"
STRAINS=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J","NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
EE_SAME="expressed elsewhere (same stage)"; STAGESHIFT="unique: stage-shifted (heterochronic)"
SNP="SNP-variant (1-3mm, same stage)"; CBS="unique: conserved-but-silent"; SP="unique: strain-private locus"; LQ="low-quality: no mm0 own-genome locus"
GU=[CBS,SP,STAGESHIFT]   # within-tp genuinely-unique
da=pd.read_csv(f"{D}/all_orderA_stagepeak_candidates.csv.gz"); da["cand_id"]=da.strain+"|"+da.timepoint+"|"+da.sequence
ee=pd.read_csv(f"{U}/unique16/ee_withintp_diag.csv.gz"); ee["cand_id"]=ee.strain+"|"+ee.timepoint+"|"+ee.sequence
eemap=dict(zip(ee.cand_id, ee.ee))
present={}
for Y in STRAINS:
    with open(f"{U}/unique16/loci/present_in_{Y}.bed") as fh:
        for line in fh:
            f=line.rstrip("\n").split("\t")
            if len(f)>=4: present.setdefault(f[3],set()).add(Y)
snp=set(pd.read_csv(f"{U}/unique16/snp_variant_refinement_withintp.csv",usecols=["cand_id"]).cand_id)
mapped=set()
for bam in glob.glob(f"{U}/cand_self16/*.cand_self16.bam"):
    b=pysam.AlignmentFile(bam,"rb")
    for a in b.fetch(until_eof=True):
        if not a.is_unmapped: mapped.add(a.query_name)
    b.close()
def classify(r):
    cid=r.cand_id; e=eemap.get(cid,"not-expressed-elsewhere")
    if e=="EE-same-stage": return EE_SAME              # shared at THIS stage -> NOT unique
    if e=="EE-other-stage-only": return STAGESHIFT     # within-tp: not expressed elsewhere at this stage -> UNIQUE
    gen=present.get(cid,set())-{r.strain}               # ortholog locus present elsewhere (halLiftover, primary D2)
    if gen: return SNP if cid in snp else CBS           # SNP (1-3mm same-stage allele) refines locus-present set (as make_klass5)
    return SP if cid in mapped else LQ                  # no ortholog -> strain-private (or low-quality if no mm0 own locus)
da["klass"]=[classify(r) for r in da.itertuples()]
da.to_csv(f"{D}/deseq_stagepeak_classified.csv.gz",index=False)
print("=== DESeq2 stage-peak (27/30) STRICT WITHIN-TP classification ===")
print(da.klass.value_counts().to_string())
print(f"\nGENUINELY-UNIQUE within-tp (strain-private + CBS + stage-shifted): {da.klass.isin(GU).sum():,} ({100*da.klass.isin(GU).mean():.0f}%)")
print(f"  of which stage-shifted (heterochronic): {da.klass.eq(STAGESHIFT).sum():,}")
print(f"  strict-sequence subset (CBS + strain-private only): {da.klass.isin([CBS,SP]).sum():,}")
print("\nper timepoint x klass:"); print(da.groupby(["timepoint","klass"]).size().to_string())
