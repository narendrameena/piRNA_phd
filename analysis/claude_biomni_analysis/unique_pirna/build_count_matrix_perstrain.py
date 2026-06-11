#!/usr/bin/env python3
"""Per-STRAIN piRNA count matrix (24-32 nt unique seqs x 9 samples = 3 timepoints x 3 reps) for
TIMEPOINT-specific edgeR DA. Within-strain, same genome -> no cross-genome/SNP step (unlike strain-
specific): a piRNA reproducibly present at one stage and absent at the others IS timepoint-specific.
Same filterByExpr-min.total.count(=15) early pre-filter for tractability. Usage: <strain>
"""
import gzip, os, sys
import numpy as np, pandas as pd
RD ="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/collapse"
OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger_tp"
os.makedirs(OUT, exist_ok=True)
TPS=["16.5dpc","12.5dpp","20.5dpp"]; REPS=[1,2,3]; LO,HI=24,32; MIN_TOTAL=15
X=sys.argv[1]
samples=[(tp,r) for tp in TPS for r in REPS]; ncol=len(samples)
counts={}; lib_win=np.zeros(ncol,dtype=np.int64); lib_full=np.zeros(ncol,dtype=np.int64)
for j,(tp,r) in enumerate(samples):
    p=f"{RD}/{X}-{tp}.{r}.raw.fasta.gz"
    with gzip.open(p,"rt") as fh:
        cur=0
        for line in fh:
            if line[:1]==">": cur=int(line[1:].strip().rsplit("-",1)[1])
            else:
                s=line.strip()
                if not s: continue
                lib_full[j]+=cur
                if LO<=len(s)<=HI:
                    lib_win[j]+=cur
                    a=counts.get(s)
                    if a is None: a=np.zeros(ncol,dtype=np.int64); counts[s]=a
                    a[j]+=cur
    print(f"  read {X}-{tp}.{r}",flush=True)
seqs=[]; rows=[]
for s,a in counts.items():
    if a.sum()>=MIN_TOTAL: seqs.append(s); rows.append(a)
M=np.vstack(rows) if rows else np.zeros((0,ncol),dtype=np.int64)
print(f"[{X}] {len(counts):,} unique 24-32nt -> {len(seqs):,} with total>={MIN_TOTAL}")
pre=f"{OUT}/{X}"
pd.DataFrame(M,columns=[f"{tp}.{r}" for (tp,r) in samples]).to_csv(f"{pre}.counts.tsv.gz",sep="\t",index=False)
with gzip.open(f"{pre}.seqs.txt.gz","wt") as o: o.write("\n".join(seqs)+"\n")
pd.DataFrame({"sample":[f"{tp}.{r}" for (tp,r) in samples],"timepoint":[tp for (tp,r) in samples],
             "rep":[r for (tp,r) in samples],"libsize_window":lib_win,"libsize_full":lib_full}
            ).to_csv(f"{pre}.samples.tsv",sep="\t",index=False)
print(f"[{X}] wrote {pre}.counts.tsv.gz ({M.shape[0]} x {M.shape[1]}) + .seqs + .samples")
