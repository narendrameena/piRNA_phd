#!/usr/bin/env python3
"""90%-cumulative-coverage agreement test (user request): the loci producing 90% of piRNA in each method, and do
they cover the SAME genomic regions? PICB: clusters ranked by all_primary_FPM (combined run, own coords). Trinity:
100/100 precursors ranked by mean rpm (union of 3 reps); genomic loci = exon blocks. Reports per strain x tp:
n90 (how few loci give 90%), %-of-all-loci, and reciprocal genomic overlap of the two 90%-sets (+ bp Jaccard).
No BAM (FPM/rpm are the within-method abundance). Output: cov90_per_strain_tp.csv."""
import pandas as pd, numpy as np, subprocess, os, tempfile, io
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
TSV=ROOT+"/figures/analysis_figures/_shared_data/picb_pangenome_clusters.tsv"
ALL=ROOT+"/results/trinity/filter_precursors/all_trinity_filtred_precursors.csv.gz"
BEDDIR=ROOT+"/results/filter_precursors_bed"; OUT=os.path.dirname(__file__); BT="/usr/bin/bedtools"
TPS={"E16.5":"16.5dpc","P12.5":"12.5dpp","P20.5":"20.5dpp"}
def sh(c): return subprocess.run(c,shell=True,capture_output=True,text=True).stdout
def w(df,p): df.to_csv(p,sep="\t",header=False,index=False)
def sm(i,o): sh(f"sort -k1,1 -k2,2n {i} | {BT} merge -i - > {o}")
def nl(p): return sum(1 for _ in open(p)) if os.path.exists(p) and os.path.getsize(p) else 0
def bp(p): return sum(int(l.split('\t')[2])-int(l.split('\t')[1]) for l in open(p)) if os.path.exists(p) and os.path.getsize(p) else 0
def n90(v): v=np.sort(np.asarray(v,float))[::-1]; return int(np.searchsorted(np.cumsum(v)/v.sum(),0.9)+1)
picb=pd.read_csv(TSV,sep="\t",usecols=["strain","tp","own_chrom","own_start","own_end","all_primary_FPM"],dtype={"own_chrom":str})
allt=pd.read_csv(ALL,header=None,names=["ref","count","length","rpm","rpkm","sample"])
s=allt["sample"].str.split("/").str[-1].str.replace(".500.csv","",regex=False)
allt["strain"]=s.map(lambda x:[x.split(f"-{k}.")[0] for k in ("16.5dpc","12.5dpp","20.5dpp") if f"-{k}." in x][0])
allt["tpv"]=s.map(lambda x:"16.5dpc" if "16.5dpc" in x else "12.5dpp" if "12.5dpp" in x else "20.5dpp")
strip=lambda c: c.split("#")[-1].replace("chr","")
rows=[]
with tempfile.TemporaryDirectory() as td:
    for st in sorted(picb.strain.unique()):
        for tpn,tp in TPS.items():
            pb=picb[(picb.strain==st)&(picb.tp==tp)].groupby(["own_chrom","own_start","own_end"],as_index=False).all_primary_FPM.max()
            tr=allt[(allt.strain==st)&(allt.tpv==tp)].groupby("ref",as_index=False).rpm.mean()
            if not len(pb) or not len(tr): continue
            np_=n90(pb.all_primary_FPM.values); nt_=n90(tr.rpm.values)
            picb90=pb.sort_values("all_primary_FPM",ascending=False).head(np_)[["own_chrom","own_start","own_end"]]
            t90ids=set(tr.sort_values("rpm",ascending=False).head(nt_).ref)
            tl=[]
            for rep in (1,2,3):
                bf=f"{BEDDIR}/{st}/{st}-{tp}.{rep}.bed"
                if not os.path.exists(bf): continue
                b=pd.read_csv(bf,sep="\t",header=None); b=b[b[3].isin(t90ids)]
                if not len(b): continue
                b.to_csv(f"{td}/b12.bed",sep="\t",header=False,index=False)
                six=sh(f"{BT} bed12tobed6 -i {td}/b12.bed")
                if six.strip(): x=pd.read_csv(io.StringIO(six),sep="\t",header=None,usecols=[0,1,2]); x[0]=x[0].map(strip); tl.append(x)
            if not tl: continue
            pf=f"{td}/p.bed"; w(picb90.sort_values(["own_chrom","own_start"]),pf); pm=f"{td}/pm.bed"; sm(pf,pm)
            tf=f"{td}/t.bed"; pd.concat(tl).drop_duplicates().sort_values([0,1]).to_csv(tf,sep="\t",header=False,index=False); tm=f"{td}/tm.bed"; sm(tf,tm)
            npm=nl(pm); ntm=nl(tm)
            p_in_t=int(sh(f"{BT} intersect -u -a {pm} -b {tm} | wc -l").strip() or 0)
            t_in_p=int(sh(f"{BT} intersect -u -a {tm} -b {pm} | wc -l").strip() or 0)
            ov=sh(f"{BT} intersect -a {pm} -b {tm}"); ovbp=sum(int(l.split('\t')[2])-int(l.split('\t')[1]) for l in ov.strip().split("\n") if l)
            ubp=bp(pm)+bp(tm)-ovbp
            rows.append(dict(strain=st,tp=tpn,n_picb=len(pb),n_trin=len(tr),n90_picb=np_,n90_trin=nt_,
                pct_picb_for90=round(100*np_/len(pb),1),pct_trin_for90=round(100*nt_/len(tr),1),
                picb90_merged=npm,trin90_merged=ntm,pct_picb90_on_trin90=round(100*p_in_t/max(npm,1),1),
                pct_trin90_on_picb90=round(100*t_in_p/max(ntm,1),1),bp_jaccard=round(100*ovbp/max(ubp,1),1)))
        print("done",st)
r=pd.DataFrame(rows); r.to_csv(f"{OUT}/cov90_per_strain_tp.csv",index=False)
print("\n=== means per tp ===")
print(r.groupby("tp")[["n90_picb","n90_trin","pct_picb_for90","pct_trin_for90","pct_picb90_on_trin90","pct_trin90_on_picb90","bp_jaccard"]].mean().round(1).reindex(["E16.5","P12.5","P20.5"]).to_string())
print("\nwrote",f"{OUT}/cov90_per_strain_tp.csv")
