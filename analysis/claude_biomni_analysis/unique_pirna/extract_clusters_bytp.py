#!/usr/bin/env python3
"""Per-TIMEPOINT version of extract_clusters.py: extract each strain's combined-replicate PICB clusters
SEPARATELY for each timepoint (16.5dpc/12.5dpp/20.5dpp) from the combined xlsx 'clusters' sheet, strip 'chr'
(HAL = Ensembl names), merge -> cluster_pav/bytp/{X}.{tp}.clusters.bed (HAL coords), ready for per-timepoint
halLiftover to GRCm39 (for the 3-timepoint 16-strain circos)."""
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
        d=pd.read_excel(f,sheet_name="clusters"); ivs=[]
        for _,r in d.iterrows():
            c=str(r["seqnames"]); c=c[3:] if c.lower().startswith("chr") else c
            ivs.append((c,int(r["start"])-1,int(r["end"])))
        raw=f"{OUT}/{X}.{tp}.raw.bed"; bed=f"{OUT}/{X}.{tp}.clusters.bed"
        pd.DataFrame(ivs,columns=["c","s","e"]).sort_values(["c","s"]).to_csv(raw,sep="\t",header=False,index=False)
        subprocess.run(f"{BT} merge -i {raw} > {bed}",shell=True,check=True); os.remove(raw)
        print(f"[{X} {tp}] {len(ivs)} -> {sum(1 for _ in open(bed))} merged",flush=True)
print("done per-timepoint cluster BEDs")
