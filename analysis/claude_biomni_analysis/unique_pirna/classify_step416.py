#!/usr/bin/env python3
"""Step 4 classification for the 16-strain run (per strain X): split X's candidates by EXPRESSION in the
other 15 strains via candidate->other-genome STAR alignments (mm<=3) + each other strain's expressed pool
(unique16/pools/{Y}.pool.txt.gz). Same logic as classify_step4.py, generalised to all 16 strains. Output =
step4_16/{X}.step4_classified16.csv.gz (4 classes incl. SNP-variant). Usage: classify_step416.py <strain>"""
import sys, gzip
import numpy as np, pandas as pd, pysam
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
S4=f"{U}/step4_16"; POOLS=f"{U}/unique16/pools"
ALL=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J","NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
X=sys.argv[1]; OTH=[s for s in ALL if s!=X]
def load_pool(s):
    st=set()
    with gzip.open(f"{POOLS}/{s}.pool.txt.gz","rt") as fh:
        for line in fh: st.add(line.rstrip("\n"))
    return st
cand=pd.read_csv(f"{S4}/{X}.candidates16.tsv.gz"); N=len(cand); id2idx={f"c{i}":i for i in range(N)}
mapped=np.zeros(N,bool); minnm=np.full(N,99,np.int16); expr=np.zeros(N,bool)
for Y in OTH:
    pool=load_pool(Y); print(f"[{X}] {Y}.pool ({len(pool):,})",flush=True)
    bam=pysam.AlignmentFile(f"{S4}/{X}.cand_to_{Y}.Aligned.sortedByCoord.out.bam","rb")
    for a in bam.fetch(until_eof=True):
        if a.is_unmapped: continue
        i=id2idx.get(a.query_name)
        if i is None: continue
        mapped[i]=True
        try: refseq=a.get_reference_sequence().upper()
        except Exception: continue
        if refseq in pool:
            nm=int(a.get_tag("NM"))
            if nm<minnm[i]: minnm[i]=nm; expr[i]=True
    bam.close(); del pool
def cls(i):
    if expr[i]: return "expressed elsewhere (exact)" if minnm[i]==0 else "SNP-variant of expressed (1-3mm)"
    return "unique: conserved-but-silent" if mapped[i] else "unique: strain-private locus"
cand["klass"]=[cls(i) for i in range(N)]
cand["min_nm_expressed"]=[int(minnm[i]) if expr[i] else -1 for i in range(N)]
cand["homolog_in_other"]=mapped
cand.to_csv(f"{S4}/{X}.step4_classified16.csv.gz",index=False,compression="gzip")
print(f"[{X}] N={N:,}\n"+cand.klass.value_counts().to_string())
g=cand.klass.str.startswith("unique").sum(); print(f"[{X}] genuinely unique = {g:,}/{N:,} ({100*g/N:.1f}%)")
