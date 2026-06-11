#!/usr/bin/env python3
"""Uniquely-mapping refinement of the TE-driven test: restrict to strain-private piRNAs that map to
EXACTLY ONE locus in their own genome (cand_self, 1 alignment) — removing the TE-multimapping ambiguity.
A uniquely-mapping strain-private piRNA whose single locus is inside a strain-private insertion AND is
TE-annotated is an UNAMBIGUOUS 'a private TE insertion created this piRNA locus' case. Reports, per
strain: uniquely-mapping strain-private piRNAs, # whose locus is in a private insertion, and the clean
(conservative) TE-driven count + families. Compares to the all-alignment number from coord_classify."""
import os,subprocess,tempfile
from collections import Counter
import pandas as pd, pysam
S4="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/step4"
PG="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
BT="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"
rows=[]
for X in ["C57BL_6NJ","CAST_EiJ","SPRET_EiJ"]:
    d=pd.read_csv(f"{S4}/{X}.step4_classified.csv.gz")
    te=pd.read_csv(f"{S4}/{X}.private_TE_percandidate.csv.gz"); id2te=dict(zip(te.id,te.classfam))
    sp=set(d.loc[d.klass=="unique: strain-private locus","id"])
    bam=pysam.AlignmentFile(f"{S4}/{X}.cand_self.Aligned.sortedByCoord.out.bam","rb")
    cnt={}; loc={}
    for a in bam.fetch(until_eof=True):
        if a.is_unmapped: continue
        q=a.query_name; cnt[q]=cnt.get(q,0)+1; loc[q]=(a.reference_name,a.reference_start,a.reference_end)
    bam.close()
    uniq={q for q,c in cnt.items() if c==1}                # NH==1
    sp_uniq=sp & uniq
    # BED of their single loci -> intersect with private-insertion loci
    cb=tempfile.NamedTemporaryFile("w",suffix=".bed",delete=False,dir=PG)
    for q in sp_uniq:
        c,s,e=loc[q]; cb.write(f"{c}\t{s}\t{e}\t{q}\n")
    cb.close()
    out=subprocess.run([BT,"intersect","-a",cb.name,"-b",f"{PG}/{X}.ins_loci.bed","-wa","-u"],
                       capture_output=True,text=True).stdout
    os.unlink(cb.name)
    in_ins=set(l.split("\t")[3] for l in out.splitlines() if l)
    tedrv=[q for q in in_ins if q in id2te]; fam=Counter(id2te[q] for q in tedrv)
    rows.append(dict(strain=X,strain_private=len(sp),uniquely_mapping=len(sp_uniq),
                     uniq_in_priv_ins=len(in_ins),TE_driven_unambiguous=len(tedrv)))
    print(f"[{X}] strain-private={len(sp):,} | uniquely-mapping={len(sp_uniq):,} | "
          f"uniq locus in private insertion={len(in_ins):,} | UNAMBIGUOUS TE-driven={len(tedrv):,}")
    print(f"        families: {dict(fam.most_common(8))}")
pd.DataFrame(rows).to_csv(f"{PG}/TE_driven_uniquelymapping.csv",index=False)
print("\nwrote TE_driven_uniquelymapping.csv")
