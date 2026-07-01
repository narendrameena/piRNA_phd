#!/usr/bin/env python3
"""16-strain version of parse_insertions.py. Reads a bcftools query stream (CHROM POS REF ALT[multi] GT*16)
for ALL 16 strains (VCF column order) and extracts insertions (>=40 bp) PRIVATE to exactly one of the 16
(stricter than the 3-strain pilot). Writes per-strain private-insertion sequences to FASTA in ins16/. id =
CHROM_POS_inslen. Does NOT overwrite the pilot's 3-strain fastas (separate dir)."""
import os,sys
OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pangenome_te/ins16"
os.makedirs(OUT,exist_ok=True)
# VCF column order (bcftools query -l): MUST match the -s order in the run script
STR=["129S1_SvImJ","AKR_J","A_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
N=len(STR); MINLEN=40
fh={s:open(f"{OUT}/{s}.private_insertions.fasta","w") for s in STR}
priv={s:0 for s in STR}; bp={s:0 for s in STR}; total_ins=shared=0
for line in sys.stdin:
    f=line.rstrip("\n").split("\t")
    if len(f)<4+N: continue
    chrom,pos,ref,alt=f[0],f[1],f[2],f[3]; gts=f[4:4+N]
    if "<" in alt: continue
    alts=alt.split(","); rl=len(ref)
    for k in set(g for g in gts if g not in ("0",".")):
        ki=int(k)
        if ki>len(alts): continue
        allele=alts[ki-1]
        if len(allele)-rl < MINLEN: continue
        total_ins+=1
        carriers=[STR[i] for i in range(N) if gts[i]==k]
        if len(carriers)==1:
            x=carriers[0]; fh[x].write(f">{chrom}_{pos}_{len(allele)-rl}\n{allele}\n"); priv[x]+=1; bp[x]+=len(allele)-rl
        else: shared+=1
for s in STR: fh[s].close()
print(f"insertions>={MINLEN}bp seen: {total_ins} | shared(>=2 of 16): {shared}")
print("private (carried by exactly 1 of 16):"); [print(f"   {s}: {priv[s]:,} ({bp[s]/1e6:.1f} Mb)") for s in STR]
