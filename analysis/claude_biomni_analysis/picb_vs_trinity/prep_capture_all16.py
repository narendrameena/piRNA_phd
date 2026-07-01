#!/usr/bin/env python3
"""Prep for the all-16-strain piRNA read-capture SLURM array: build, ONCE (no BAM access), per (strain,tp):
PICB cluster bed in PanSN coords + Trinity precursor EXON-block bed (200/200, union of 3 reps). Write to beds16/
and emit a samplesheet of samples that have a rep1 BAM. The array tasks then only run samtools counts."""
import pandas as pd, subprocess, os, io
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
TSV=ROOT+"/figures/analysis_figures/_shared_data/picb_pangenome_clusters.tsv"
ALL=ROOT+"/results/trinity/filter_precursors/all_trinity_filtred_precursors.csv.gz"
BEDDIR=ROOT+"/results/filter_precursors_bed"; BAMDIR=ROOT+"/results/STAR_srna_strain_wise"
OUT=os.path.dirname(__file__); BEDS=f"{OUT}/beds16"; os.makedirs(BEDS,exist_ok=True); BT="/usr/bin/bedtools"
TPS={"E16.5":"16.5dpc","P12.5":"12.5dpp","P20.5":"20.5dpp"}
print("loading PICB + Trinity ids ...")
picb=pd.read_csv(TSV,sep="\t",usecols=["strain","tp","own_chrom","own_start","own_end"],dtype={"own_chrom":str}).drop_duplicates()
allt=pd.read_csv(ALL,header=None,names=["ref","count","length","rpm","rpkm","sample"])
allt["s"]=allt["sample"].str.split("/").str[-1].str.replace(".500.csv","",regex=False)
ids_by=allt[(allt.rpm>100)&(allt.rpkm>100)].groupby("s").ref.apply(set).to_dict()   # thesis 100/100 (== saved all_)
sheet=[]
for st in sorted(picb.strain.unique()):
    for tpn,tp in TPS.items():
        bam=f"{BAMDIR}/{st}/{st}-{tp}.1/Aligned.sortedByCoord.out.bam"
        if not os.path.exists(bam): print("no BAM, skip",st,tpn); continue
        pb=picb[(picb.strain==st)&(picb.tp==tp)].copy()
        if not len(pb): print("no PICB, skip",st,tpn); continue
        pb["c"]=st+"#1#chr"+pb.own_chrom
        pb[["c","own_start","own_end"]].sort_values(["c","own_start"]).to_csv(f"{BEDS}/{st}-{tp}.picb.bed",sep="\t",header=False,index=False)
        tl=[]
        for rep in (1,2,3):
            s=f"{st}-{tp}.{rep}"; bf=f"{BEDDIR}/{st}/{s}.bed"; ids=ids_by.get(s,set())
            if os.path.exists(bf) and ids:
                b=pd.read_csv(bf,sep="\t",header=None); b=b[b[3].isin(ids)]
                if len(b):
                    b.to_csv(f"{BEDS}/_t12.bed",sep="\t",header=False,index=False)
                    six=subprocess.run(f"{BT} bed12tobed6 -i {BEDS}/_t12.bed",shell=True,capture_output=True,text=True).stdout
                    if six.strip(): tl.append(pd.read_csv(io.StringIO(six),sep="\t",header=None,usecols=[0,1,2]))
        if not tl: print("no Trinity, skip",st,tpn); continue
        pd.concat(tl).drop_duplicates().sort_values([0,1]).to_csv(f"{BEDS}/{st}-{tp}.trin.bed",sep="\t",header=False,index=False)
        sheet.append(f"{st}\t{tpn}\t{tp}\t{bam}")
os.path.exists(f"{BEDS}/_t12.bed") and os.remove(f"{BEDS}/_t12.bed")
open(f"{OUT}/samplesheet16.tsv","w").write("\n".join(sheet)+"\n")
print(f"wrote {len(sheet)} samples to samplesheet16.tsv + beds to beds16/")
