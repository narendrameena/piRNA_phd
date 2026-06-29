#!/usr/bin/env python3
"""FULL-SPAN (intron-included) Trinity precursor beds for the read-capture 'intron-included' comparison bar.
Same source + 100/100 (rpm&rpkm) id filter as prep_capture_all16.py, but uses each precursor's bed12 TRANSCRIPT
SPAN (cols 0,1,2 = chrom,start,end INCLUDING introns) instead of bed12tobed6 EXON blocks. This is the
intron-inflated version (the artifact we corrected by using exon blocks); the bar shows its magnitude.
Output: beds16/{st}-{tp}.trin_span.bed (per sample, union of 3 reps)."""
import pandas as pd, os
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
ALL=ROOT+"/results/trinity/filter_precursors/all_trinity_filtred_precursors.csv.gz"
BEDDIR=ROOT+"/results/filter_precursors_bed"
P=ROOT+"/analysis/claude_biomni_analysis/picb_vs_trinity"; BEDS=f"{P}/beds16"
allt=pd.read_csv(ALL,header=None,names=["ref","count","length","rpm","rpkm","sample"])
allt["s"]=allt["sample"].str.split("/").str[-1].str.replace(".500.csv","",regex=False)
ids_by=allt[(allt.rpm>100)&(allt.rpkm>100)].groupby("s").ref.apply(set).to_dict()   # thesis 100/100
sheet=pd.read_csv(f"{P}/samplesheet16.tsv",sep="\t",header=None)   # strain,tpn,tp,bam (rep1; one row per st-tp)
n=0
for _,r in sheet.iterrows():
    st,tp=r[0],r[2]; tl=[]
    for rep in (1,2,3):
        s=f"{st}-{tp}.{rep}"; bf=f"{BEDDIR}/{st}/{s}.bed"; ids=ids_by.get(s,set())
        if os.path.exists(bf) and ids:
            b=pd.read_csv(bf,sep="\t",header=None); b=b[b[3].isin(ids)]
            if len(b): tl.append(b[[0,1,2]])                       # full transcript span (incl. introns)
    if tl:
        pd.concat(tl).drop_duplicates().sort_values([0,1]).to_csv(f"{BEDS}/{st}-{tp}.trin_span.bed",sep="\t",header=False,index=False); n+=1
print(f"wrote {n} full-span (intron-included) trin_span.bed files to {BEDS}/")
