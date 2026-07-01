#!/usr/bin/env python3
"""#4 cluster PAV — final: 16 strains' PICB clusters lifted to GRCm39 -> multi-strain presence/absence
-> CORE (all 16) / ACCESSORY (2-15) / PRIVATE (1 strain). Each strain's lifted intervals are merged
first; bedtools multiinter builds the per-region presence across strains. Reports region counts + bp by
class, the private-by-strain breakdown, and writes the catalogue."""
import os,subprocess,tempfile
import pandas as pd
PAV="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"
BT="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"
STR=["WSB_EiJ","CAST_EiJ","BALB_cJ","C57BL_6NJ","A_J","DBA_2J","NOD_ShiLtJ","SPRET_EiJ","AKR_J",
     "C3H_HeJ","CBA_J","PWK_PhJ","NZO_HlLtJ","FVB_NJ","129S1_SvImJ","LP_J"]
merged=[]
for X in STR:
    f=f"{PAV}/{X}.in_GRCm39.bed"
    m=f"{PAV}/{X}.GRCm39.merged.bed"
    # keep cols 1-3 (chrom,start,end), sort, merge
    subprocess.run(f"cut -f1-3 {f} | sort -k1,1 -k2,2n | {BT} merge -i - > {m}",shell=True,check=True)
    merged.append(m)
names=",".join(STR)
mi=subprocess.run([BT,"multiinter","-header","-names",*STR,"-i",*merged],capture_output=True,text=True).stdout
rows=[]
for ln in mi.splitlines():
    if ln.startswith("chrom") or not ln.strip(): continue
    f=ln.split("\t")
    chrom,s,e,num=f[0],int(f[1]),int(f[2]),int(f[3])
    rows.append((chrom,s,e,num,e-s))
df=pd.DataFrame(rows,columns=["chrom","start","end","n_strains","bp"])
def cls(n): return "core" if n==16 else ("private" if n==1 else "accessory")
df["class"]=df.n_strains.map(cls)
df.to_csv(f"{PAV}/cluster_PAV_catalogue.csv.gz",index=False,compression="gzip")
# summary
print("=== cluster PAV (GRCm39 common frame) ===")
g=df.groupby("class").agg(regions=("bp","size"),total_bp=("bp","sum"))
print(g.to_string())
print(f"\ntotal PAV regions: {len(df):,} | total bp: {df.bp.sum()/1e6:.1f} Mb")
print("\nregions by #strains present:")
print(df.groupby("n_strains").size().to_string())
# private by strain: regions present in exactly one strain -> which
priv=df[df.n_strains==1].copy()
# re-run multiinter membership for private rows: which single strain
mi2=subprocess.run([BT,"multiinter","-names",*STR,"-i",*merged],capture_output=True,text=True).stdout
from collections import Counter
pc=Counter()
for ln in mi2.splitlines():
    f=ln.split("\t")
    if len(f)<5 or int(f[3])!=1: continue
    pc[f[4]]+= (int(f[2])-int(f[1]))
print("\nprivate cluster bp by strain (top):")
for s,bp in sorted(pc.items(),key=lambda x:-x[1])[:16]: print(f"  {s}: {bp/1e3:.0f} kb")
pd.DataFrame(g.reset_index()).to_csv(f"{PAV}/cluster_PAV_summary.csv",index=False)
print("\nwrote cluster_PAV_catalogue.csv.gz + cluster_PAV_summary.csv")
