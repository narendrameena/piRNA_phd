"""Do ragged-3' pachytene (P20.5 30nt) piRNAs come from pachytene CLUSTERS (and the high-coverage/expressed ones)?
Map piRNA 5' positions (cand_self16, self-coords) onto the strain's 20.5dpp clusters (clusters_fpm.bed), then TEST
raggedness cluster-vs-not, and within clusters high-FPM (dominant/conserved) vs low-FPM."""
import subprocess, pandas as pd, numpy as np
from scipy.stats import fisher_exact
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
ST="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/samtools"
T="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/20_pirna_3prime_length_heterogeneity"
clean=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["sequence","strain","timepoint"])
clean["L"]=clean.sequence.str.len(); clean["cand_id"]=clean.strain+"|"+clean.timepoint+"|"+clean.sequence
P=clean[clean.timepoint=="20.5dpp"]; need=set(P.cand_id); SETS={L:set(P[P.L==L].sequence) for L in range(23,35)}
rows=[]
for X in sorted(P.strain.unique()):
    out=subprocess.run([ST,"view",f"{U}/cand_self16/{X}.cand_self16.bam"],capture_output=True,text=True,timeout=900).stdout
    for ln in out.splitlines():
        f=ln.split("\t")
        if len(f)<6 or int(f[1]) not in (0,16) or f[0] not in need: continue
        seq=f[0].split("|")[2]; ch=f[2].split("#")[-1]; ch=ch[3:] if ch.startswith("chr") else ch
        rows.append((f[0],X,seq,len(seq),ch,int(f[3])))
c=pd.DataFrame(rows,columns=["cand_id","strain","seq","L","chrom","pos"]).drop_duplicates("cand_id")
print(f"P20.5 coords: {len(c):,}")
c["incluster"]=False; c["fpm"]=0.0
for X,g in c.groupby("strain"):
    bed=pd.read_csv(f"{U}/cluster_pav/{X}.clusters_fpm.bed",sep="\t",header=None,names=["chrom","start","end","fp","fm","strand","tp"],dtype={"chrom":str})
    bed=bed[bed.tp=="20.5dpp"]
    by={ch:b.sort_values("start")[["start","end","fp","fm"]].values for ch,b in bed.groupby("chrom")}
    for i in g.index.values:
        arr=by.get(c.at[i,"chrom"]);
        if arr is None: continue
        p=c.at[i,"pos"]; j=int(np.searchsorted(arr[:,0],p))-1
        for k in (j,j+1):
            if 0<=k<len(arr) and arr[k,0]<=p<=arr[k,1]: c.at[i,"incluster"]=True; c.at[i,"fpm"]=arr[k,2]+arr[k,3]; break
print(f"piRNAs in a pachytene cluster: {100*c.incluster.mean():.0f}%")
L=30; sset=SETS[L]; pre={k:{y[:L] for y in SETS[L+k]} for k in (1,2,3)}
israg=lambda x: any(x[:L-k] in SETS[L-k] for k in (1,2,3)) or any(x in pre[k] for k in (1,2,3))
g30=c[c.L==L].groupby("seq").agg(incluster=("incluster","any"),fpm=("fpm","max"))
df=pd.DataFrame([(s,israg(s),bool(g30.incluster.get(s,False)),float(g30.fpm.get(s,0.0))) for s in sset if s in g30.index],columns=["seq","rag","incluster","fpm"])
rc=df.groupby("incluster").rag.mean()*100; p=fisher_exact(pd.crosstab(df.incluster,df.rag).values)[1]
print(f"P20.5 30nt raggedness: pachytene-cluster {rc.get(True,float('nan')):.0f}% vs non-cluster {rc.get(False,float('nan')):.0f}% (Fisher p={p:.1e}); cluster fraction {100*df.incluster.mean():.0f}%")
cl=df[df.incluster].copy()
if len(cl)>20:
    q=cl.fpm.quantile([.5]).iloc[0]; hi=cl[cl.fpm>=q].rag.mean()*100; lo=cl[cl.fpm<q].rag.mean()*100
    print(f"within clusters: high-expression (FPM>={q:.0f}) {hi:.0f}% vs low-expression {lo:.0f}% ragged")
df.to_csv(f"{T}/data/SourceData_cluster_origin.csv",index=False)
print("saved SourceData_cluster_origin.csv")
