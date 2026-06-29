#!/usr/bin/env python3
"""PROBE — impact of within-timepoint vs cross-timepoint 'expressed elsewhere' on the unique-piRNA classification.
Current: a candidate (strain X, tp T, seq) is 'expressed elsewhere (exact)' iff seq is in another strain's pool
POOLED OVER ALL 3 tp. Within-tp: iff another strain expresses seq AT THE SAME tp T. Within-tp is STRICTER on
'expressed elsewhere' -> fewer non-unique, MORE genuinely-unique. Quantify the flip (no re-derivation yet)."""
import gzip, glob, collections, pandas as pd, os
RD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/collapse"
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
STRAINS=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J","NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
TPS=["16.5dpc","12.5dpp","20.5dpp"]; LO,HI=25,32; MIN=2   # 25-32 (adopted) within-tp pools
print("building per-(strain,tp) pools (25-32, >=2 reads at that tp) ...",flush=True)
pool={}   # (strain,tp) -> set(seq)
for s in STRAINS:
    for tp in TPS:
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
        pool[(s,tp)]={seq for seq,c in tot.items() if c>=MIN}
    print(f"  {s} pools: "+" ".join(f"{tp}={len(pool[(s,tp)]):,}" for tp in TPS),flush=True)
tpmap={"E16.5":"16.5dpc","P12.5":"12.5dpp","P20.5":"20.5dpp"}
d=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["sequence","strain","timepoint","klass5"])
def within_tp_expressed_elsewhere(row):
    tp=tpmap.get(row.timepoint,row.timepoint)
    return any(row.sequence in pool.get((Y,tp),()) for Y in STRAINS if Y!=row.strain)
d["ee_within_tp"]=d.apply(within_tp_expressed_elsewhere,axis=1)
gu=["unique: conserved-but-silent","unique: strain-private locus"]
cur_gu=d.klass5.isin(gu).sum()
ee_cross=d.klass5=="expressed elsewhere (exact)"
flip=ee_cross & (~d.ee_within_tp)   # currently expressed-elsewhere (cross-tp) but NOT within-tp
print(f"\n=== IMPACT (25-32 set, {len(d):,} candidates) ===")
print(f"current 'expressed-elsewhere (exact)': {ee_cross.sum():,}")
print(f"  of those, NOT expressed within their own tp elsewhere (would leave 'expressed-elsewhere'): {flip.sum():,} ({100*flip.sum()/max(ee_cross.sum(),1):.0f}%)")
print(f"current genuinely-unique: {cur_gu:,} ({100*cur_gu/len(d):.0f}%)")
print(f"NOTE: within-tp is stricter on 'expressed-elsewhere' -> these {flip.sum():,} would re-enter the locus/SNP/unique pipeline (some become unique).")
print(f"upper-bound new genuinely-unique (if all flips were unique): {cur_gu+flip.sum():,} ({100*(cur_gu+flip.sum())/len(d):.0f}%)")
