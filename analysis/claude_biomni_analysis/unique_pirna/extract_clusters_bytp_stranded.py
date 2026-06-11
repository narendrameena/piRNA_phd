#!/usr/bin/env python3
"""Stranded per-timepoint PICB cluster extraction (for the sense/antisense circos). Like extract_clusters_bytp
but KEEPS the PICB cluster strand: + (sense) / - (antisense) / * (bidirectional/unstranded -> '.'), strand-aware
merge, write a 6-col BED (strand in col6) so halLiftover preserves+flips it into the GRCm39 frame ->
cluster_pav/bytp/{X}.{tp}.stranded.clusters.bed."""
import os,subprocess
import pandas as pd
RES="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/picb_result_combined"
OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/cluster_pav/bytp"
BT="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"; os.makedirs(OUT,exist_ok=True)
STRAINS=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J","NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
TPS=["16.5dpc","12.5dpp","20.5dpp"]
for X in STRAINS:
    for tp in TPS:
        f=f"{RES}/{X}/{X}-{tp}.combined.xlsx"
        if not os.path.exists(f): print(f"  [missing] {X} {tp}"); continue
        d=pd.read_excel(f,sheet_name="clusters"); rows=[]
        for i,r in d.iterrows():
            c=str(r["seqnames"]); c=c[3:] if c.lower().startswith("chr") else c
            st=str(r["strand"]); st=st if st in ("+","-") else "."
            rows.append((c,int(r["start"])-1,int(r["end"]),f"c{i}",0,st))
        raw=f"{OUT}/{X}.{tp}.s.raw.bed"; bed=f"{OUT}/{X}.{tp}.stranded.clusters.bed"
        pd.DataFrame(rows,columns=["c","s","e","n","sc","st"]).sort_values(["c","s"]).to_csv(raw,sep="\t",header=False,index=False)
        # strand-aware merge -> chrom start end strand ; reformat to 6-col (strand in col6) for halLiftover
        subprocess.run(f"{BT} merge -s -c 6 -o distinct -i {raw} | awk 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2,$3,\"c\"NR,0,$4}}' > {bed}",shell=True,check=True)
        os.remove(raw)
        from collections import Counter
        cc=Counter(l.split("\t")[5].strip() for l in open(bed))
        print(f"[{X} {tp}] {sum(cc.values())} clusters  +{cc.get('+',0)} -{cc.get('-',0)} .{cc.get('.',0)}",flush=True)
print("done stranded per-timepoint cluster BEDs")
