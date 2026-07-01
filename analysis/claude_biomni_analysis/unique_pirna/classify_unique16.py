#!/usr/bin/env python3
"""Scale-16 pangenome phase — STEP 1 of the cross-strain expression test (the 16-strain Step-4 equivalent).
SEQUENCE-LEVEL exact-expression split (certain; pangenome-free): a strain-specific candidate piRNA from
strain X (from edger16/{tp}.strain_specific_DA.csv.gz) is 'expressed elsewhere (exact)' iff its EXACT
sequence is in another strain's expressed pool (unique16/pools/{Y}.pool.txt.gz). Those are NOT genuinely
unique. The remainder = NOVEL-SEQUENCE candidates, written out per strain for the pangenome LOCUS step
(cand_loci16.sh -> lift_cand16.sh -> classify_unique16_locus.py), which splits them into SNP-variant /
unique:conserved-but-silent / unique:strain-private-locus.

Loads the 16 pools once as a seq->16-bit mask (memory-lean), then streams the per-tp candidate tables.
No timepoint arg: processes all three. Run after pools (3302627) AND edgeR DA (3302624) finish."""
import gzip, os, glob
import pandas as pd
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
POOLS=f"{U}/pools" if os.path.isdir(f"{U}/pools") else f"{U}/unique16/pools"
OUT=f"{U}/unique16"; os.makedirs(OUT,exist_ok=True)
STRAINS=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J",
         "NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
IDX={s:i for i,s in enumerate(STRAINS)}; TPS=["16.5dpc","12.5dpp","20.5dpp"]

# ---- guard: inputs ready? ----
missing=[s for s in STRAINS if not os.path.exists(f"{POOLS}/{s}.pool.txt.gz")]
if missing: raise SystemExit(f"[wait] {len(missing)} pools not built yet: {missing[:3]}... (pools job 3302627)")
da=[tp for tp in TPS if not os.path.exists(f"{U}/edger16/{tp}.strain_specific_DA.csv.gz")]
if da: raise SystemExit(f"[wait] edgeR DA not done for: {da} (job 3302624)")

# ---- load pools as seq -> 16-bit strain mask ----
print("loading 16 expressed pools -> bitmask ...",flush=True)
poolmask={}
for s in STRAINS:
    bit=1<<IDX[s]
    with gzip.open(f"{POOLS}/{s}.pool.txt.gz","rt") as fh:
        for line in fh:
            seq=line.strip()
            if seq: poolmask[seq]=poolmask.get(seq,0)|bit
    print(f"  + {s}: pool loaded ({len(poolmask):,} distinct so far)",flush=True)

# ---- classify each tp's candidates ----
allrows=[]
for tp in TPS:
    d=pd.read_csv(f"{U}/edger16/{tp}.strain_specific_DA.csv.gz")
    masks=d.sequence.map(lambda q: poolmask.get(q,0)).to_numpy()
    selfbit=d.strain.map(lambda s: 1<<IDX[s]).to_numpy()
    other=masks & ~selfbit                                    # strains OTHER than X expressing the exact seq
    d["n_other_exact"]=[bin(int(m)).count("1") for m in other]
    d["other_exact_strains"]=[",".join(STRAINS[i] for i in range(16) if (int(m)>>i)&1) for m in other]
    d["expr_class"]=["expressed elsewhere (exact)" if m else "novel-sequence (-> locus step)" for m in other]
    allrows.append(d)
    nov=int((d.expr_class.str.startswith("novel")).sum())
    print(f"[{tp}] candidates {len(d):,} | expressed-elsewhere-exact {len(d)-nov:,} | novel-sequence {nov:,}",flush=True)
    # emit novel-sequence candidates per strain (FASTA) for the locus-mapping step
    novd=d[d.expr_class.str.startswith("novel")]
    for X in STRAINS:
        sub=novd[novd.strain==X]
        if len(sub):
            with open(f"{OUT}/{X}.{tp}.novel.fasta","w") as o:
                for _,r in sub.iterrows(): o.write(f">{X}|{tp}|{r.sequence}\n{r.sequence}\n")

res=pd.concat(allrows,ignore_index=True)
res.to_csv(f"{OUT}/expr_exact_classified.csv.gz",index=False)
n=len(res); nov=int(res.expr_class.str.startswith("novel").sum())
print(f"\nTOTAL {n:,} candidates -> expressed-elsewhere-exact {n-nov:,} ({100*(n-nov)/n:.1f}%) | "
      f"novel-sequence -> pangenome locus step {nov:,} ({100*nov/n:.1f}%)")
print(f"wrote {OUT}/expr_exact_classified.csv.gz + per-strain/tp *.novel.fasta")
