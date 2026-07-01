#!/usr/bin/env python3
"""Per-timepoint PICB clusters WITH expression + strand (for the sense/antisense expression circos). Each cluster
kept individually (no merge) as a 6-col BED: chrom, start, end, name=uniq_reads_FPM (expression), score=0,
strand (+ sense / - antisense / * -> '.' bidirectional). halLiftover preserves col4 (FPM) and interprets/flips
col6 (strand) into the GRCm39 frame. Output cluster_pav/bytp/{X}.{tp}.expr.clusters.bed."""
import os
import pandas as pd
RES="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/picb_result_combined"
OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/cluster_pav/bytp"; os.makedirs(OUT,exist_ok=True)
STRAINS=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J","NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
TPS=["16.5dpc","12.5dpp","20.5dpp"]
for X in STRAINS:
    for tp in TPS:
        f=f"{RES}/{X}/{X}-{tp}.combined.xlsx"
        if not os.path.exists(f): print(f"  [missing] {X} {tp}"); continue
        d=pd.read_excel(f,sheet_name="clusters"); rows=[]
        fpmcol="uniq_reads_FPM" if "uniq_reads_FPM" in d.columns else "all_reads_primary_alignments_FPM"
        for _,r in d.iterrows():
            c=str(r["seqnames"]); c=c[3:] if c.lower().startswith("chr") else c
            st=str(r["strand"]); st=st if st in ("+","-") else "."
            rows.append((c,int(r["start"])-1,int(r["end"]),round(float(r[fpmcol]),3),0,st))
        bed=f"{OUT}/{X}.{tp}.expr.clusters.bed"
        pd.DataFrame(rows,columns=["c","s","e","fpm","sc","st"]).sort_values(["c","s"]).to_csv(bed,sep="\t",header=False,index=False)
        from collections import Counter
        cc=Counter(r[5] for r in rows)
        print(f"[{X} {tp}] {len(rows)} clusters +{cc.get('+',0)} -{cc.get('-',0)} .{cc.get('.',0)} | medianFPM {pd.Series([r[3] for r in rows]).median():.1f}",flush=True)
print("done per-timepoint expr+strand cluster BEDs")
