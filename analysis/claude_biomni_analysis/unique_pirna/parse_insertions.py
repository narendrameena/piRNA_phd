#!/usr/bin/env python3
"""Read a bcftools query stream (CHROM POS REF ALT[multi] GT_C57 GT_CAST GT_SPRET) for the 3 pilot
strains and extract STRAIN-PRIVATE INSERTIONS. ALT is the comma-list of alleles carried by the 3
strains (after `view -a --min-ac 1`); each strain's GT is the allele index (0=REF). An allele k is an
INSERTION if len(ALT[k])-len(REF) >= MINLEN. It is PRIVATE to strain X (among the 3) if X is the only
one of the three whose GT == k. Writes per-strain private-insertion sequences to FASTA (the literal
inserted allele), id = CHROM_POS_inslen. MINLEN=40 bp: large enough to embed a 24-32 nt piRNA with
flank; TE identity is NOT assumed here (it comes from Step-6a RepeatMasker of the piRNA locus).
"""
import os,sys
OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
os.makedirs(OUT,exist_ok=True)
STR=["C57BL_6NJ","CAST_EiJ","SPRET_EiJ"]; MINLEN=40
fh={s:open(f"{OUT}/{s}.private_insertions.fasta","w") for s in STR}
priv={s:0 for s in STR}; shared=0; total_ins=0; bp={s:0 for s in STR}
for line in sys.stdin:
    f=line.rstrip("\n").split("\t")
    if len(f)<7: continue
    chrom,pos,ref,alt=f[0],f[1],f[2],f[3]; gts=f[4:7]
    if "<" in alt: continue
    alts=alt.split(","); rl=len(ref)
    # allele indices actually carried by the 3 strains
    idx=set(g for g in gts if g not in ("0","."))
    for k in idx:
        ki=int(k)
        if ki>len(alts): continue
        allele=alts[ki-1]
        if len(allele)-rl < MINLEN: continue
        total_ins+=1
        carriers=[STR[i] for i in range(3) if gts[i]==k]
        if len(carriers)==1:
            x=carriers[0]; fh[x].write(f">{chrom}_{pos}_{len(allele)-rl}\n{allele}\n")
            priv[x]+=1; bp[x]+=len(allele)-rl
        else:
            shared+=1
for s in STR: fh[s].close()
print(f"insertions>={MINLEN}bp seen: {total_ins}")
print(f"private (carried by exactly 1 of 3): {priv}")
print(f"private inserted bp: {bp}")
print(f"shared (>=2 of 3 carry same allele): {shared}")
