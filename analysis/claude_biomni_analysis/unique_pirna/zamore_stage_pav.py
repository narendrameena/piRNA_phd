#!/usr/bin/env python3
"""Use the Zamore stage annotation (pre-pachytene / hybrid / pachytene) for the piRNA clusters.
Stage is column 13 of the mm10 bed12 (Ozata et al. 2020); mapped onto the 214 GRCm39-lifted Zamore loci
by LOCUS NAME (build-independent, no coordinate mixing). Then test conservation-by-stage in our 16-strain
cluster PAV: for each Zamore locus, the coverage-weighted mean #strains-present (from cluster_PAV_catalogue)
-> mean by stage. Prediction (Ozata 2020): pachytene clusters deeply conserved (high #strains), pre-
pachytene/hybrid more variable."""
import os,re,subprocess,tempfile
import pandas as pd, numpy as np
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U=f"{B}/analysis/claude_biomni_analysis/unique_pirna"
PAVc=f"{U}/cluster_pav/cluster_PAV_catalogue.csv.gz"
ZB=f"{B}/analysis/claude_biomni_analysis/C57BL_6NJ_pangenome/_loci_mm39_noprefix.bed"
XLS=f"{B}/resources/zamore_piRNAs/piRNA_gene_annotation-modified-in-orange-with-Wasik-et-al-checks.xlsx"
BT="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"
def base(n): return re.sub(r"\.\d+$","",str(n))
# stage_by_locus from mm10 bed12 (col4=name idx3, col13=stage idx12)
d=pd.read_excel(XLS,sheet_name="Mus musculus",header=None,skiprows=1)
stage={}
for _,r in d.iterrows():
    nm=r[3];
    if pd.isna(nm): continue
    st=str(r[12]).strip().lower()
    stage[base(nm)]=st
print("stage label values seen:",sorted(set(stage.values()))[:10])
# GRCm39 Zamore loci -> stage by name
z=pd.read_csv(ZB,sep="\t",header=None,names=["c","s","e","name"])
z["stage"]=z["name"].map(lambda n: stage.get(base(n),"unknown"))
print("\n214 GRCm39 Zamore loci by stage:")
print(z.stage.value_counts().to_string())
# PAV conservation per locus: intersect with catalogue (n_strains in col4), coverage-weighted mean
zb=tempfile.NamedTemporaryFile("w",suffix=".bed",delete=False)
for _,r in z.sort_values(["c","s"]).iterrows(): zb.write(f"{r.c}\t{r.s}\t{r.e}\t{r['name']}\n")
zb.close()
cat=pd.read_csv(PAVc,dtype={"chrom":str}); catbed=tempfile.NamedTemporaryFile("w",suffix=".bed",delete=False)
for _,r in cat.iterrows(): catbed.write(f"{r.chrom}\t{r.start}\t{r.end}\t{r.n_strains}\n")
catbed.close()
out=subprocess.run(f"sort -k1,1 -k2,2n {catbed.name} | {BT} intersect -a {zb.name} -b - -wo",shell=True,capture_output=True,text=True).stdout
os.unlink(zb.name); os.unlink(catbed.name)
from collections import defaultdict
acc=defaultdict(lambda:[0.0,0])   # locus -> [sum(n*ov), sum(ov)]
for ln in out.splitlines():
    f=ln.split("\t"); loc=f[3]; n=int(f[7]); ov=int(f[8])
    acc[loc][0]+=n*ov; acc[loc][1]+=ov
z["mean_nstrains"]=z["name"].map(lambda l: acc[l][0]/acc[l][1] if acc.get(l) and acc[l][1]>0 else np.nan)
print("\nconservation by Zamore stage (coverage-weighted mean #strains with a cluster, /16):")
g=z.dropna(subset=["mean_nstrains"]).groupby("stage").mean_nstrains.agg(["mean","median","size"]).round(2)
print(g.to_string())
z.to_csv(f"{U}/cluster_pav/zamore_loci_stage_conservation.csv",index=False)
print("\nwrote zamore_loci_stage_conservation.csv")
