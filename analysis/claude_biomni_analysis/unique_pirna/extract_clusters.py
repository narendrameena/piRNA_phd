#!/usr/bin/env python3
"""#4 cluster PAV — step 1: extract each strain's combined PICB clusters (union across the 3 timepoints)
from the combined xlsx 'clusters' sheet, strip the 'chr' prefix (PICB uses chrN; the cactus HAL uses
Ensembl names 1,2,..,X,MT), merge overlapping intervals -> per-strain cluster BED in HAL coordinates,
ready for halLiftover to GRCm39."""
import os,subprocess
import pandas as pd
RES="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/picb_result_combined"
OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"
BT="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"
os.makedirs(OUT,exist_ok=True)
STRAINS=["WSB_EiJ","CAST_EiJ","BALB_cJ","C57BL_6NJ","A_J","DBA_2J","NOD_ShiLtJ","SPRET_EiJ","AKR_J",
         "C3H_HeJ","CBA_J","PWK_PhJ","NZO_HlLtJ","FVB_NJ","129S1_SvImJ","LP_J"]
TPS=["16.5dpc","12.5dpp","20.5dpp"]
summ=[]
for X in STRAINS:
    ivs=[]
    for tp in TPS:
        f=f"{RES}/{X}/{X}-{tp}.combined.xlsx"
        if not os.path.exists(f): print(f"  [missing] {f}"); continue
        d=pd.read_excel(f,sheet_name="clusters")
        for _,r in d.iterrows():
            chrom=str(r["seqnames"])
            chrom=chrom[3:] if chrom.lower().startswith("chr") else chrom   # chr1 -> 1 ; chrM->M
            ivs.append((chrom,int(r["start"])-1,int(r["end"])))             # 1-based PICB -> 0-based BED
    raw=f"{OUT}/{X}.clusters.raw.bed"; bed=f"{OUT}/{X}.clusters.bed"
    pd.DataFrame(ivs,columns=["c","s","e"]).sort_values(["c","s"]).to_csv(raw,sep="\t",header=False,index=False)
    subprocess.run(f"{BT} merge -i {raw} > {bed}",shell=True,check=True); os.remove(raw)
    n=sum(1 for _ in open(bed))
    summ.append(dict(strain=X,clusters_merged=n)); print(f"[{X}] {len(ivs)} raw -> {n} merged clusters")
pd.DataFrame(summ).to_csv(f"{OUT}/cluster_counts_perstrain.csv",index=False)
print("wrote per-strain cluster BEDs (HAL coords) + cluster_counts_perstrain.csv")
