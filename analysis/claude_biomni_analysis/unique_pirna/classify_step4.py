#!/usr/bin/env python3
"""Step 4 classification (per strain X): split X's strain-specific DA candidates by EXPRESSION in the
other strains, using candidate->other-genome STAR alignments (mm<=3) + each other strain's expressed-
sequence set (pool_uniq). For a candidate mapping to strain Y at a locus, the reference sequence there
(from the MD tag) is exactly the piRNA Y would make from that locus; Y EXPRESSES it iff that sequence
is in Y.pool_uniq. NM = SNP distance candidate<->Y-genome.
Classes:
  expressed elsewhere (exact, 0mm)      -> NOT unique (Y expresses the identical sequence)
  SNP-variant of expressed (1-3mm)      -> NOT unique (Y expresses a conserved-locus variant)
  unique: conserved-but-silent          -> GENUINELY unique by EXPRESSION (homolog in Y but Y silent)
  unique: strain-private locus          -> GENUINELY unique (no <=3mm homolog in any other genome)
Usage: classify_step4.py <strain>
"""
import os,sys
import numpy as np, pandas as pd, pysam
S4="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/step4"
X=sys.argv[1]
OTH={"C57BL_6NJ":["CAST_EiJ","SPRET_EiJ"],"CAST_EiJ":["C57BL_6NJ","SPRET_EiJ"],
     "SPRET_EiJ":["C57BL_6NJ","CAST_EiJ"]}[X]

def load_pool(s):
    st=set()
    with open(f"{S4}/{s}.pool_uniq.fasta") as fh:
        for line in fh:
            if line[:1]!=">": st.add(line.rstrip("\n"))
    return st

cand=pd.read_csv(f"{S4}/{X}.candidates.tsv.gz"); N=len(cand)
id2idx={f"c{i}":i for i in range(N)}
mapped   ={Y:np.zeros(N,dtype=bool)        for Y in OTH}   # homolog exists in Y (<=3mm)
expr_minnm={Y:np.full(N,99,dtype=np.int16) for Y in OTH}   # min NM over loci Y actually expresses
for Y in OTH:
    pool=load_pool(Y); print(f"[{X}] loaded {Y}.pool_uniq ({len(pool):,})",flush=True)
    bam=pysam.AlignmentFile(f"{S4}/{X}.cand_to_{Y}.Aligned.sortedByCoord.out.bam","rb")
    n=0
    for a in bam.fetch(until_eof=True):
        if a.is_unmapped: continue
        i=id2idx.get(a.query_name)
        if i is None: continue
        mapped[Y][i]=True
        try: refseq=a.get_reference_sequence().upper()
        except Exception: continue
        if refseq in pool:
            nm=int(a.get_tag("NM"))
            if nm<expr_minnm[Y][i]: expr_minnm[Y][i]=nm
        n+=1
    bam.close(); del pool
    print(f"[{X}] processed {n:,} alignments vs {Y}",flush=True)

mapped_any=np.zeros(N,dtype=bool); expr_any=np.zeros(N,dtype=bool); minnm_any=np.full(N,99,dtype=np.int16)
for Y in OTH:
    mapped_any|=mapped[Y]; expr_any|=(expr_minnm[Y]<99); minnm_any=np.minimum(minnm_any,expr_minnm[Y])
def cls(i):
    if expr_any[i]: return "expressed elsewhere (exact)" if minnm_any[i]==0 else "SNP-variant of expressed (1-3mm)"
    return "unique: conserved-but-silent" if mapped_any[i] else "unique: strain-private locus"
cand["klass"]=[cls(i) for i in range(N)]
cand["min_nm_expressed"]=[int(minnm_any[i]) if expr_any[i] else -1 for i in range(N)]
cand["homolog_in_other"]=mapped_any
cand.to_csv(f"{S4}/{X}.step4_classified.csv.gz",index=False,compression="gzip")
print(f"\n[{X}] N={N:,}")
print(cand.klass.value_counts().to_string())
g=cand.klass.str.startswith("unique").sum()
print(f"[{X}] GENUINELY UNIQUE (by expression) = {g:,}/{N:,} ({100*g/N:.1f}%)")
