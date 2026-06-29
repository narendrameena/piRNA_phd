#!/usr/bin/env python3
"""piRNA-SPECIFIC coverage (user request): of TOTAL piRNA expression (piRNA-length 24-32 nt alignments, multimappers
weighted), what fraction is captured by PICB clusters vs Trinity precursor loci? Restricts the all-sRNA capture to the
piRNA read population (excludes the 21 nt miRNA peak + degradation). 24-32 nt = the project's piRNA window (phasing).
Trinity = EXON BLOCKS. BAM = results/STAR_srna_strain_wise (PanSN). total = full-BAM stream + awk length filter.
Subset: pachytene P20.5 x3 (SPRET/CAST/C57BL_6NJ) + P12.5 x2. Output: read_capture_pirna.csv."""
import pandas as pd, subprocess, os, tempfile, io
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
TSV=ROOT+"/figures/analysis_figures/_shared_data/picb_pangenome_clusters.tsv"
ALL=ROOT+"/results/trinity/filter_precursors/all_trinity_filtred_precursors.csv.gz"
BEDDIR=ROOT+"/results/filter_precursors_bed"; BAMDIR=ROOT+"/results/STAR_srna_strain_wise"
OUT=os.path.dirname(__file__); ST="/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/samtools"; BT="/usr/bin/bedtools"
SAMPLES=[("SPRET_EiJ","P20.5","20.5dpp"),("CAST_EiJ","P20.5","20.5dpp"),("C57BL_6NJ","P20.5","20.5dpp"),
         ("SPRET_EiJ","P12.5","12.5dpp"),("C57BL_6NJ","P12.5","12.5dpp")]
AWK="awk 'length($10)>=24 && length($10)<=32'"   # piRNA-length window
def cl(cmd): return int(subprocess.run(f"{cmd} | {AWK} | wc -l",shell=True,capture_output=True,text=True).stdout.strip() or 0)
def sh(c): return subprocess.run(c,shell=True,capture_output=True,text=True).stdout
picb=pd.read_csv(TSV,sep="\t",usecols=["strain","tp","own_chrom","own_start","own_end"],dtype={"own_chrom":str}).drop_duplicates()
allt=pd.read_csv(ALL,header=None,names=["ref","count","length","rpm","rpkm","sample"])
allt["s"]=allt["sample"].str.split("/").str[-1].str.replace(".500.csv","",regex=False)
ids_by=allt[(allt.rpm>100)&(allt.rpkm>100)].groupby("s").ref.apply(set).to_dict()   # thesis 100/100 (== saved all_); was 200/200
rows=[]
with tempfile.TemporaryDirectory() as td:
    for st,tpn,tp in SAMPLES:
        bam=f"{BAMDIR}/{st}/{st}-{tp}.1/Aligned.sortedByCoord.out.bam"
        if not os.path.exists(bam): print("skip",bam); continue
        total=cl(f"{ST} view {bam}")                       # total piRNA-length alignment mass (genome-wide)
        pb=picb[(picb.strain==st)&(picb.tp==tp)].copy(); pb["c"]=st+"#1#chr"+pb.own_chrom
        pbed=f"{td}/picb.bed"; pb[["c","own_start","own_end"]].sort_values(["c","own_start"]).to_csv(pbed,sep="\t",header=False,index=False)
        tl=[]
        for rep in (1,2,3):
            s=f"{st}-{tp}.{rep}"; bf=f"{BEDDIR}/{st}/{s}.bed"; ids=ids_by.get(s,set())
            if os.path.exists(bf) and ids:
                b=pd.read_csv(bf,sep="\t",header=None); b=b[b[3].isin(ids)]
                if len(b):
                    b12=f"{td}/in12.bed"; b.to_csv(b12,sep="\t",header=False,index=False)
                    six=sh(f"{BT} bed12tobed6 -i {b12}")
                    if six.strip(): tl.append(pd.read_csv(io.StringIO(six),sep="\t",header=None,usecols=[0,1,2]))
        tbed=f"{td}/trin.bed";
        if tl: pd.concat(tl).drop_duplicates().sort_values([0,1]).to_csv(tbed,sep="\t",header=False,index=False)
        in_picb=cl(f"{ST} view -L {pbed} {bam}")
        in_trin=cl(f"{ST} view -L {tbed} {bam}") if tl else 0
        rows.append(dict(strain=st,tp=tpn,total_pirna=total,in_picb=in_picb,in_trin=in_trin,
            pct_picb=round(100*in_picb/max(total,1),2),pct_trin=round(100*in_trin/max(total,1),2)))
        print(rows[-1],flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/read_capture_pirna.csv",index=False); print("wrote read_capture_pirna.csv")
