#!/usr/bin/env python3
"""PICB clusters vs Trinity precursors — genomic concordance per strain x timepoint. CORRECTED to use Trinity
EXON BLOCKS (bed12tobed6), not full genomic spans: 40% of Trinity contigs are multi-exon and their span includes
introns (5-7x inflation) — piRNA reads come from exons. Counting unit for Trinity = distinct CONTIGS (precursors),
not intron-merged intervals. Both sets in the SAME strain REL-2205 assembly (verified chr1 194,686,469). PICB chrom
'1', Trinity 'SID#1#chr1' -> strip '1'. PICB = own_chrom/own_start/own_end (merged). Trinity = 100/100 contigs
(rpm>100 & rpkm>100, thesis Ch.6), mapped, exon blocks, union of 3 reps. Output: overlap_per_strain_tp.csv."""
import pandas as pd, numpy as np, subprocess, os, tempfile, io
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
TSV=ROOT+"/figures/analysis_figures/_shared_data/picb_pangenome_clusters.tsv"
ALL=ROOT+"/results/trinity/filter_precursors/all_trinity_filtred_precursors.csv.gz"
BEDDIR=ROOT+"/results/filter_precursors_bed"; OUT=os.path.dirname(__file__); BT="/usr/bin/bedtools"
TPS=["16.5dpc","12.5dpp","20.5dpp"]
def sh(c): return subprocess.run(c,shell=True,capture_output=True,text=True).stdout
def w(df,p): df.to_csv(p,sep="\t",header=False,index=False)
def sortmerge(inp,outp): sh(f"sort -k1,1 -k2,2n {inp} | {BT} merge -i - > {outp}")
def nlines(p): return sum(1 for _ in open(p)) if os.path.exists(p) and os.path.getsize(p) else 0
def bp(p): return sum(int(l.split('\t')[2])-int(l.split('\t')[1]) for l in open(p)) if os.path.exists(p) and os.path.getsize(p) else 0
print("loading PICB tsv + Trinity 100/100 ids ...")
picb=pd.read_csv(TSV,sep="\t",usecols=["strain","tp","own_chrom","own_start","own_end"],dtype={"own_chrom":str}).drop_duplicates()
allt=pd.read_csv(ALL,header=None,names=["ref","count","length","rpm","rpkm","sample"])
allt["s"]=allt["sample"].str.split("/").str[-1].str.replace(".500.csv","",regex=False)
a2=allt[(allt.rpm>100)&(allt.rpkm>100)]   # thesis Ch.6 rule: 100 RPM & 100 RPKM (== the saved all_ file)
ids_by=a2.groupby("s").ref.apply(set).to_dict()
strip=lambda c: c.split("#")[-1].replace("chr","")
rows=[]
with tempfile.TemporaryDirectory() as td:
    for st in sorted(picb.strain.unique()):
        for tp in TPS:
            pb=picb[(picb.strain==st)&(picb.tp==tp)][["own_chrom","own_start","own_end"]]
            if not len(pb): continue
            # Trinity exon blocks (100/100) over reps, with contig id
            tl=[]
            for rep in (1,2,3):
                s=f"{st}-{tp}.{rep}"; bf=f"{BEDDIR}/{st}/{s}.bed"; ids=ids_by.get(s,set())
                if not os.path.exists(bf) or not ids: continue
                b=pd.read_csv(bf,sep="\t",header=None); b=b[b[3].isin(ids)]
                if not len(b): continue
                b12=f"{td}/in12.bed"; b.to_csv(b12,sep="\t",header=False,index=False)
                six=sh(f"{BT} bed12tobed6 -i {b12}")
                if not six.strip(): continue
                b6=pd.read_csv(io.StringIO(six),sep="\t",header=None,usecols=[0,1,2,3])
                b6[0]=b6[0].map(strip); tl.append(b6)
            if not tl: continue
            tb=pd.concat(tl).drop_duplicates()           # exon intervals + contig name (col3)
            n_trin_contigs=tb[3].nunique()
            # files: picb merged ; trin exon (4col, named) ; trin exon merged (3col)
            pf=f"{td}/p.bed"; w(pb.sort_values(["own_chrom","own_start"]),pf); pm=f"{td}/pm.bed"; sortmerge(pf,pm)
            te=f"{td}/te.bed"; w(tb[[0,1,2,3]].sort_values([0,1]),te)
            tem=f"{td}/tem.bed"; sortmerge(te,tem)
            n_picb=nlines(pm)
            # contigs (precursors) with >=1 exon overlapping a PICB cluster
            ov=sh(f"{BT} intersect -a {te} -b {pm} -u")
            contigs_in_picb=set(l.split('\t')[3] for l in ov.strip().split('\n') if l)
            n_trin_in_picb=len(contigs_in_picb)
            # PICB clusters recovered by >=1 Trinity exon
            n_picb_hit=int(sh(f"{BT} intersect -u -a {pm} -b {te} | wc -l").strip() or 0)
            # fragmentation: distinct Trinity contigs per recovered PICB cluster
            fr=sh(f"{BT} intersect -a {pm} -b {te} -wa -wb"); frag={}
            for line in fr.strip().split("\n"):
                if not line: continue
                f=line.split("\t"); frag.setdefault((f[0],f[1],f[2]),set()).add(f[6])
            fv=[len(v) for v in frag.values()]
            picb_bp=bp(pm); trin_bp=bp(tem)
            ovf=sh(f"{BT} intersect -a {pm} -b {tem}"); ovbp=sum(int(l.split(chr(9))[2])-int(l.split(chr(9))[1]) for l in ovf.strip().split("\n") if l)
            plen=int((pb.own_end-pb.own_start).median())
            tlen=int(a2[(a2.s.str.startswith(f"{st}-{tp}."))].length.median())   # precursor = Trinity CONTIG length
            rows.append(dict(strain=st,tp=tp,n_picb=n_picb,n_trin_contigs=n_trin_contigs,
                n_trin_in_picb=n_trin_in_picb,pct_trin_in_picb=round(100*n_trin_in_picb/max(n_trin_contigs,1),1),
                n_picb_hit=n_picb_hit,pct_picb_recovered=round(100*n_picb_hit/max(n_picb,1),1),
                frag_median=int(np.median(fv)) if fv else 0,frag_mean=round(float(np.mean(fv)),2) if fv else 0,
                picb_bp=picb_bp,trin_exon_bp=trin_bp,overlap_bp=ovbp,
                pct_picb_bp_by_trin=round(100*ovbp/max(picb_bp,1),1),pct_trin_bp_in_picb=round(100*ovbp/max(trin_bp,1),1),
                picb_med_len=plen,trin_med_len=tlen))
        print("  done",st)
res=pd.DataFrame(rows); res.to_csv(f"{OUT}/overlap_per_strain_tp.csv",index=False)
TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; res["TP"]=res.tp.map(TPN)
print("\n=== SPRET_EiJ (sanity, EXON-based) ===")
print(res[res.strain=="SPRET_EiJ"][["TP","n_picb","n_trin_contigs","pct_trin_in_picb","pct_picb_recovered","frag_mean","pct_picb_bp_by_trin","picb_med_len","trin_med_len"]].to_string(index=False))
print("\n=== means over 16 strains, per tp ===")
print(res.groupby("TP")[["n_picb","n_trin_contigs","pct_trin_in_picb","pct_picb_recovered","frag_mean","pct_picb_bp_by_trin"]].mean().round(1).reindex(["E16.5","P12.5","P20.5"]).to_string())
print("\nwrote",f"{OUT}/overlap_per_strain_tp.csv")
