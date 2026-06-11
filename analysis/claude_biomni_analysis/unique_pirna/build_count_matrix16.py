#!/usr/bin/env python3
"""Scale-16 — per-timepoint piRNA count matrix across ALL 16 sRNA strains (16 strains x 3 reps = up to
48 samples). 24-32 nt unique sequences; filterByExpr min.total.count(=15) early pre-filter for
tractability. Robust to missing samples (uses whatever collapse files exist). Usage: <timepoint>"""
import gzip, os, sys
import numpy as np, pandas as pd
RD ="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/collapse"
OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger16"
os.makedirs(OUT, exist_ok=True)
STRAINS=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J",
         "NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]   # 16 sRNA strains
REPS=[1,2,3]; LO,HI=24,32; MIN_TOTAL=15
tp=sys.argv[1]
samples=[(s,r) for s in STRAINS for r in REPS]
present=[(s,r) for (s,r) in samples if os.path.exists(f"{RD}/{s}-{tp}.{r}.raw.fasta.gz")]
print(f"[{tp}] {len(present)}/{len(samples)} samples present",flush=True)
ncol=len(present); counts={}; lib_win=np.zeros(ncol,dtype=np.int64); lib_full=np.zeros(ncol,dtype=np.int64)
for j,(s,r) in enumerate(present):
    with gzip.open(f"{RD}/{s}-{tp}.{r}.raw.fasta.gz","rt") as fh:
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
    print(f"  read {s}-{tp}.{r}",flush=True)
seqs=[]; rows=[]
for seq,a in counts.items():
    if a.sum()>=MIN_TOTAL: seqs.append(seq); rows.append(a)
M=np.vstack(rows) if rows else np.zeros((0,ncol),dtype=np.int64)
print(f"[{tp}] {len(counts):,} unique 24-32nt -> {len(seqs):,} total>={MIN_TOTAL}")
pre=f"{OUT}/{tp}"
pd.DataFrame(M,columns=[f"{s}.{r}" for (s,r) in present]).to_csv(f"{pre}.counts.tsv.gz",sep="\t",index=False)
with gzip.open(f"{pre}.seqs.txt.gz","wt") as o: o.write("\n".join(seqs)+"\n")
pd.DataFrame({"sample":[f"{s}.{r}" for (s,r) in present],"strain":[s for (s,r) in present],
             "rep":[r for (s,r) in present],"libsize_window":lib_win,"libsize_full":lib_full}
            ).to_csv(f"{pre}.samples.tsv",sep="\t",index=False)
print(f"[{tp}] wrote {pre}.counts.tsv.gz ({M.shape[0]} x {M.shape[1]}) + .seqs + .samples")
