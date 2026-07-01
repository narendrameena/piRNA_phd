#!/usr/bin/env python3
"""Build a per-timepoint piRNA count matrix (24-32 nt unique sequences x 9 samples = 3 strains x 3
reps) for edgeR differential-abundance. The reproducibility-only presence rule produced ~30M noise-
dominated candidates; the data-driven floor is edgeR filterByExpr + FDR-controlled DA (user-approved).
Here we only PARSE + assemble the matrix and apply filterByExpr's own documented min.total.count(=15)
early, purely for tractability (it removes nothing filterByExpr would keep). Full filterByExpr (CPM +
min group size) and the QL F-test run in edgeR downstream.
Usage: build_count_matrix.py <timepoint>   e.g. 16.5dpc | 12.5dpp | 20.5dpp
"""
import gzip, os, sys
import numpy as np, pandas as pd
RD ="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/collapse"
OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger"
os.makedirs(OUT, exist_ok=True)
PILOT=["C57BL_6NJ","CAST_EiJ","SPRET_EiJ"]; REPS=[1,2,3]; LO,HI=24,32
MIN_TOTAL=15   # = edgeR filterByExpr default min.total.count; applied early for tractability ONLY
tp=sys.argv[1]
samples=[(s,r) for s in PILOT for r in REPS]; ncol=len(samples)
counts={}                                   # seq -> np.array(ncol) int64
lib_win =np.zeros(ncol,dtype=np.int64)      # 24-32 nt total per sample (piRNA library size)
lib_full=np.zeros(ncol,dtype=np.int64)      # all-length total per sample
for j,(s,r) in enumerate(samples):
    p=f"{RD}/{s}-{tp}.{r}.raw.fasta.gz"
    with gzip.open(p,"rt") as fh:
        cur=0
        for line in fh:
            if line[:1]==">": cur=int(line[1:].strip().rsplit("-",1)[1])
            else:
                seq=line.strip()
                if not seq: continue
                lib_full[j]+=cur
                if LO<=len(seq)<=HI:
                    lib_win[j]+=cur
                    a=counts.get(seq)
                    if a is None: a=np.zeros(ncol,dtype=np.int64); counts[seq]=a
                    a[j]+=cur
    print(f"  read {s}-{tp}.{r}: lib_full={lib_full[j]:,} lib_window={lib_win[j]:,}",flush=True)

seqs=[]; rows=[]
for seq,a in counts.items():
    if a.sum()>=MIN_TOTAL: seqs.append(seq); rows.append(a)
M=np.vstack(rows) if rows else np.zeros((0,ncol),dtype=np.int64)
print(f"[{tp}] {len(counts):,} unique 24-32nt seqs -> {len(seqs):,} with total>={MIN_TOTAL}")
pre=f"{OUT}/{tp}"
pd.DataFrame(M,columns=[f"{s}.{r}" for (s,r) in samples]).to_csv(f"{pre}.counts.tsv.gz",sep="\t",index=False)
with gzip.open(f"{pre}.seqs.txt.gz","wt") as o: o.write("\n".join(seqs)+"\n")
pd.DataFrame({"sample":[f"{s}.{r}" for (s,r) in samples],"strain":[s for (s,r) in samples],
             "rep":[r for (s,r) in samples],"libsize_window":lib_win,"libsize_full":lib_full}
            ).to_csv(f"{pre}.samples.tsv",sep="\t",index=False)
print(f"[{tp}] wrote {pre}.counts.tsv.gz ({M.shape[0]} x {M.shape[1]}) + .seqs.txt.gz + .samples.tsv")
