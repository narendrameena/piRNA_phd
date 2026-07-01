#!/usr/bin/env python3
"""WITHIN-TIMEPOINT pools: per (strain, timepoint) expressed pool = 25-32 nt sequences with >=2 reads AT THAT tp
(pooled over the 3 reps only). Replaces the cross-tp build_pools16. Output: pools_pertp/{strain}_{tp}.pool.txt.gz.
One (strain,tp) at a time -> bounded memory."""
import gzip, collections, os
RD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/collapse"
OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/pools_pertp"
os.makedirs(OUT,exist_ok=True)
STRAINS=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J","NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
TPS=["16.5dpc","12.5dpp","20.5dpp"]; LO,HI=25,32; MIN=2
for s in STRAINS:
    for tp in TPS:
        o=f"{OUT}/{s}_{tp}.pool.txt.gz"
        if os.path.exists(o): continue
        tot=collections.Counter(); hdr=None
        for r in (1,2,3):
            f=f"{RD}/{s}-{tp}.{r}.raw.fasta.gz"
            if not os.path.exists(f): continue
            with gzip.open(f,"rt") as fh:
                for line in fh:
                    if line[0]==">": hdr=line
                    else:
                        seq=line.strip()
                        if LO<=len(seq)<=HI: tot[seq]+=int(hdr[1:].split("-")[-1])
        pool=[seq for seq,c in tot.items() if c>=MIN]
        with gzip.open(o,"wt") as w: w.write("\n".join(pool)+"\n")
        print(f"[{s} {tp}] -> {len(pool):,}",flush=True); del tot,pool
print("DONE per-tp pools")
