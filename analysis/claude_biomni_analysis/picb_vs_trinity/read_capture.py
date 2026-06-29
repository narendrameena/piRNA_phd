#!/usr/bin/env python3
"""piRNA read-capture: of all sRNA alignment mass (multimappers weighted, as in FPM), what fraction overlaps
PICB clusters vs Trinity precursor loci. Representative subset: 2 wild + 2 classical x P12.5/P20.5, rep1.
BAM = results/STAR_srna_strain_wise (PanSN chroms). PICB own '1' -> '{strain}#1#chr1'; Trinity BED already PanSN.
total alignment mass = sum(idxstats mapped). in_set = samtools view -c -L bed bam. Output: read_capture.csv."""
import pandas as pd, subprocess, os, tempfile
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
TSV=ROOT+"/figures/analysis_figures/_shared_data/picb_pangenome_clusters.tsv"
ALL=ROOT+"/results/trinity/filter_precursors/all_trinity_filtred_precursors.csv.gz"
BEDDIR=ROOT+"/results/filter_precursors_bed"; BAMDIR=ROOT+"/results/STAR_srna_strain_wise"
OUT=os.path.dirname(__file__); ST="/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/samtools"; BT="/usr/bin/bedtools"
import io
SAMPLES=[("SPRET_EiJ","P20.5","20.5dpp"),("CAST_EiJ","P20.5","20.5dpp"),("C57BL_6NJ","P20.5","20.5dpp"),
         ("SPRET_EiJ","P12.5","12.5dpp"),("C57BL_6NJ","P12.5","12.5dpp")]  # pachytene x3 (fast) + P12.5 x2 (developmental)
def sh(c): return subprocess.run(c,shell=True,capture_output=True,text=True).stdout
picb=pd.read_csv(TSV,sep="\t",usecols=["strain","tp","own_chrom","own_start","own_end"],dtype={"own_chrom":str}).drop_duplicates()
allt=pd.read_csv(ALL,header=None,names=["ref","count","length","rpm","rpkm","sample"])
allt["s"]=allt["sample"].str.split("/").str[-1].str.replace(".500.csv","",regex=False)
ids_by=allt[(allt.rpm>=200)&(allt.rpkm>=200)].groupby("s").ref.apply(set).to_dict()
rows=[]
with tempfile.TemporaryDirectory() as td:
    if True:
        for st,tpn,tp in SAMPLES:
            bam=f"{BAMDIR}/{st}/{st}-{tp}.1/Aligned.sortedByCoord.out.bam"
            if not os.path.exists(bam): print("skip",bam); continue
            # total mapped alignment mass (idxstats col3)
            total=sum(int(l.split("\t")[2]) for l in sh(f"{ST} idxstats {bam}").strip().split("\n") if l)
            # PICB bed in PanSN
            pb=picb[(picb.strain==st)&(picb.tp==tp)].copy()
            pb["c"]=st+"#1#chr"+pb.own_chrom
            pbed=f"{td}/picb.bed"; pb[["c","own_start","own_end"]].sort_values(["c","own_start"]).to_csv(pbed,sep="\t",header=False,index=False)
            # Trinity bed PanSN (union reps, 200/200) — EXON BLOCKS (bed12tobed6), not intron-inflated spans
            tl=[]
            for rep in (1,2,3):
                s=f"{st}-{tp}.{rep}"; bf=f"{BEDDIR}/{st}/{s}.bed"; ids=ids_by.get(s,set())
                if os.path.exists(bf) and ids:
                    b=pd.read_csv(bf,sep="\t",header=None); b=b[b[3].isin(ids)]
                    if len(b):
                        b12=f"{td}/in12.bed"; b.to_csv(b12,sep="\t",header=False,index=False)
                        six=subprocess.run(f"{BT} bed12tobed6 -i {b12}",shell=True,capture_output=True,text=True).stdout
                        if six.strip(): tl.append(pd.read_csv(io.StringIO(six),sep="\t",header=None,usecols=[0,1,2]))
            tbed=f"{td}/trin.bed"
            if tl:
                tb=pd.concat(tl).drop_duplicates(); tb.sort_values([0,1]).to_csv(tbed,sep="\t",header=False,index=False)
            in_picb=int(sh(f"{ST} view -c -L {pbed} {bam}").strip() or 0)
            in_trin=int(sh(f"{ST} view -c -L {tbed} {bam}").strip() or 0) if tl else 0
            rows.append(dict(strain=st,tp=tpn,total=total,in_picb=in_picb,in_trin=in_trin,
                pct_picb=round(100*in_picb/max(total,1),2),pct_trin=round(100*in_trin/max(total,1),2)))
            print(rows[-1])
pd.DataFrame(rows).to_csv(f"{OUT}/read_capture.csv",index=False); print("wrote read_capture.csv")
