"""E16.5 ping-pong: MILI(PIWIL2)/MIWI2(PIWIL4). Compute the 10-nt 5'-5' overlap (ping-pong) signature from cand_self16
coords, then TEST whether ping-pong piRNAs (secondary, slicer-set 5' end) are more/less 3'-ragged than non-ping-pong."""
import subprocess, pandas as pd, numpy as np
from collections import Counter
from scipy.stats import fisher_exact
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
ST="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/samtools"
clean=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["sequence","strain","timepoint"])
clean["L"]=clean.sequence.str.len(); clean["cand_id"]=clean.strain+"|"+clean.timepoint+"|"+clean.sequence
E=clean[clean.timepoint=="16.5dpc"]; need=set(E.cand_id)
SETS={L:set(E[E.L==L].sequence) for L in range(23,35)}
rows=[]
for X in sorted(E.strain.unique()):
    out=subprocess.run([ST,"view",f"{U}/cand_self16/{X}.cand_self16.bam"],capture_output=True,text=True,timeout=900).stdout
    for ln in out.splitlines():
        f=ln.split("\t")
        if len(f)<6 or int(f[1]) not in (0,16) or f[0] not in need: continue
        seq=f[0].split("|")[2]; L=len(seq); pos=int(f[3]); strand="+" if int(f[1])==0 else "-"
        rows.append((f[0],f[0].split("|")[0],seq,L,f[2],pos if strand=="+" else pos+L-1,strand))
c=pd.DataFrame(rows,columns=["cand_id","strain","seq","L","chrom","five","strand"]).drop_duplicates("cand_id")
print(f"E16.5 coords: {len(c):,} candidates")
ov=Counter(); ppset=set()                              # 5'-5' overlap = q(minus) - p(plus) + 1
for (st,chrom),g in c.groupby(["strain","chrom"]):
    plus=set(g[g.strand=="+"].five); minus=set(g[g.strand=="-"].five)
    if not plus or not minus: continue
    for o in range(1,31): ov[o]+=len({p+o-1 for p in plus}&minus)
    for p in plus:
        if p+9 in minus: ppset.add((st,chrom,p,"+"))
    for q in minus:
        if q-9 in plus: ppset.add((st,chrom,q,"-"))
cnts=np.array([ov[o] for o in range(1,31)],float); oth=np.concatenate([cnts[:9],cnts[10:]])
print(f"ping-pong 10-nt-overlap z = {(cnts[9]-oth.mean())/oth.std():.1f}  (9nt={ov[9]:,} 10nt={ov[10]:,} 11nt={ov[11]:,})")
L=27; sset=SETS[L]; pre={k:{y[:L] for y in SETS[L+k]} for k in (1,2,3)}
israg=lambda x: any(x[:L-k] in SETS[L-k] for k in (1,2,3)) or any(x in pre[k] for k in (1,2,3))
c["pp"]=[(r.strain,r.chrom,r.five,r.strand) in ppset for r in c.itertuples()]
seqpp=c[c.L==L].groupby("seq").pp.any()
df=pd.DataFrame([(s,israg(s),bool(seqpp.get(s,False))) for s in sset if s in seqpp.index],columns=["seq","rag","pp"])
rt=df.groupby("pp").rag.mean()*100; p=fisher_exact(pd.crosstab(df.pp,df.rag).values)[1]
print(f"E16.5 27nt raggedness: ping-pong {rt.get(True,float('nan')):.0f}% vs non-ping-pong {rt.get(False,float('nan')):.0f}% (Fisher p={p:.1e}); ping-pong fraction {100*df.pp.mean():.0f}%")
df.to_csv("/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/20_pirna_3prime_length_heterogeneity/data/SourceData_pingpong_e16.csv",index=False)
