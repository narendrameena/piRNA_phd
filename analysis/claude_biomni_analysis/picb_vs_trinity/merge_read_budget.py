#!/usr/bin/env python3
"""Merge per-rep capture data for the read-budget figure:
  cap16_reps/{st}-{tp}.{rep}.csv  -> strain,tpn,rep,total_pirna(25-32),in_picb,in_trin_exon
  intron_cap/{st}-{tp}.{rep}.csv  -> strain,tpn,rep,total_all,in_trin_span(intron-included)
Output: read_budget_capture.csv (per strain,tp,rep + derived %). All counts = multimapper-weighted alignments."""
import pandas as pd, glob
P="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/picb_vs_trinity"
def load(pat, cols):
    rows=[]
    for f in glob.glob(pat):
        a=open(f).read().strip().split(",")
        if len(a)>=len(cols): rows.append(dict(zip(cols, a)))
    return pd.DataFrame(rows)
cap=load(f"{P}/cap16_reps/*.csv", ["strain","tp","rep","total_pirna","in_picb","in_trin_exon"])
intr=load(f"{P}/intron_cap/*.csv", ["strain","tp","rep","total_all","in_trin_span"])
for c in ["total_pirna","in_picb","in_trin_exon"]: cap[c]=cap[c].astype(int)
for c in ["total_all","in_trin_span"]: intr[c]=intr[c].astype(int)
d=cap.merge(intr,on=["strain","tp","rep"])
d["pirna_pct_of_sRNA"]=100*d.total_pirna/d.total_all
d["pct_picb"]=100*d.in_picb/d.total_pirna
d["pct_trin_exon"]=100*d.in_trin_exon/d.total_pirna
d["pct_trin_span"]=100*d.in_trin_span/d.total_pirna
d.to_csv(f"{P}/read_budget_capture.csv",index=False)
print(f"merged {len(d)} rows -> read_budget_capture.csv")
g=d.groupby("tp").agg(total_all_M=("total_all",lambda x:x.mean()/1e6), pirna_M=("total_pirna",lambda x:x.mean()/1e6),
    pirna_pct=("pirna_pct_of_sRNA","mean"), picb=("pct_picb","mean"), exon=("pct_trin_exon","mean"), span=("pct_trin_span","mean"))
print(g.round(2).to_string())
