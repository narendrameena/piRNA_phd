#!/usr/bin/env python3
"""WITHIN-TIMEPOINT re-derivation of klass5 (user-directed: uniqueness judged WITHIN each developmental stage,
since pre-pachytene vs pachytene piRNAs are different biology). Reclassifies the SAME adopted >=2-read, 25-32 nt
candidate set; reuses all tp-INDEPENDENT determinants + the new per-tp pools.

  D1 expressed-elsewhere (exact): seq in another strain's SAME-tp pool (pools_pertp/) -> NOT unique.
  D2 (novel within-tp): present elsewhere = seq expressed in another strain at ANY tp (any per-tp pool)  OR  the
     locus is genome-present elsewhere (loci/present_in_{Y}.bed, halLiftover) -> conserved-but-silent; else strain-private.
  D3 strain-private with NO mm0 own-genome locus (cand_self16) -> low-quality.
  D4 conserved-but-silent that is a 1-3mm SNP-variant of an allele EXPRESSED IN ANOTHER STRAIN AT THE SAME tp
     -> SNP-variant. Uses the WITHIN-tp snp set (snp_variant_refinement_withintp.csv, from
     classify_step416_pertp.py): stricter than cross-tp, so a few CBS loci move back to genuinely-unique.
Output: unique16/final_classified_clean_2read.csv.gz (within-tp); backs up the cross-tp version."""
import gzip, glob, os, pysam, pandas as pd
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
PP=f"{U}/pools_pertp"; LOCI=f"{U}/unique16/loci"; IN=f"{U}/unique16"
STRAINS=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J","NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
IDX={s:i for i,s in enumerate(STRAINS)}; TPS=["16.5dpc","12.5dpp","20.5dpp"]
tpmap={"E16.5":"16.5dpc","P12.5":"12.5dpp","P20.5":"20.5dpp"}
CBS="unique: conserved-but-silent"; SP="unique: strain-private locus"
EE="expressed elsewhere (exact)"; SNP="SNP-variant (1-3mm)"; LQ="low-quality: no mm0 own-genome locus"
print("loading per-tp masks ...",flush=True)
mask={tp:{} for tp in TPS}
for tp in TPS:
    for s in STRAINS:
        bit=1<<IDX[s]
        with gzip.open(f"{PP}/{s}_{tp}.pool.txt.gz","rt") as fh:
            for line in fh:
                q=line.strip()
                if q: mask[tp][q]=mask[tp].get(q,0)|bit
    print(f"  {tp}: {len(mask[tp]):,} distinct",flush=True)
present={}
for Y in STRAINS:
    with open(f"{LOCI}/present_in_{Y}.bed") as fh:
        for line in fh:
            f=line.rstrip("\n").split("\t")
            if len(f)>=4: present.setdefault(f[3],set()).add(Y)
mapped=set()
for bam in glob.glob(f"{U}/cand_self16/*.cand_self16.bam"):
    b=pysam.AlignmentFile(bam,"rb")
    for a in b.fetch(until_eof=True):
        if not a.is_unmapped: mapped.add(a.query_name)
    b.close()
snp=set(pd.read_csv(f"{IN}/snp_variant_refinement_withintp.csv",usecols=["cand_id"]).cand_id)
d=pd.read_csv(f"{IN}/final_classified_clean_2read.csv.gz"); d["cand_id"]=d.strain+"|"+d.timepoint+"|"+d.sequence
def classify(r):
    tp=tpmap.get(r.timepoint,r.timepoint); sb=1<<IDX[r.strain]   # final_classified tp is already dpc/dpp form
    if mask[tp].get(r.sequence,0) & ~sb: return EE
    anytp=(mask["16.5dpc"].get(r.sequence,0)|mask["12.5dpp"].get(r.sequence,0)|mask["20.5dpp"].get(r.sequence,0)) & ~sb
    gen=present.get(r.cand_id,set())-{r.strain}
    base=CBS if (anytp or gen) else SP
    if base==CBS and r.cand_id in snp: return SNP
    if base==SP and r.cand_id not in mapped: return LQ
    return base
print("reclassifying within-tp ...",flush=True)
old=d.klass5.value_counts().to_dict()
d["klass5"]=d.apply(classify,axis=1)
bk=f"{IN}/final_classified_clean_2read.crosstp_backup.csv.gz"
if os.path.exists(bk): raise SystemExit(f"ABORT: {bk} already exists -> current file may already be within-tp; refusing to overwrite the cross-tp backup")
os.rename(f"{IN}/final_classified_clean_2read.csv.gz",bk)
d.drop(columns=["cand_id"]).to_csv(f"{IN}/final_classified_clean_2read.csv.gz",index=False)
gu=[CBS,SP]
print("\n=== WITHIN-TP klass5 (vs cross-tp) ===")
for k in [EE,SNP,LQ,CBS,SP]: print(f"  {k}: {d.klass5.eq(k).sum():,}  (cross-tp {old.get(k,0):,})")
print(f"GENUINELY-UNIQUE within-tp: {d.klass5.isin(gu).sum():,} ({100*d.klass5.isin(gu).mean():.0f}%)  | cross-tp was {sum(old.get(k,0) for k in gu):,}")
print("wrote within-tp final_classified_clean_2read.csv.gz (cross-tp backed up)")
