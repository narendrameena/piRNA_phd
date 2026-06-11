#!/usr/bin/env python3
"""Step 4 prep (per strain X). (1) Union of X's strain-specific DA candidate sequences across the 3
timepoints -> X.candidates.fasta (+ .tsv recording which timepoints each is specific at). (2) X's own
expressed-sequence pool = union of X's 24-32 nt collapse unique sequences (all 9 samples) ->
X.pool_uniq.fasta. The 'others pool' for the expression test is later built by concatenating the OTHER
strains' pool_uniq.fasta (each collapse file read once here). Usage: prep_step4.py <strain>
"""
import gzip, os, sys
import pandas as pd
ED ="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/edger"
RD ="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/collapse"
OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/step4"
os.makedirs(OUT, exist_ok=True)
TPS=[("16.5dpc","E16.5"),("12.5dpp","P12.5"),("20.5dpp","P20.5")]; REPS=[1,2,3]; LO,HI=24,32
X=sys.argv[1]

# (1) candidate union across timepoints
cand={}                                   # seq -> set(timepoint labels)
for f,lab in TPS:
    d=pd.read_csv(f"{ED}/{f}.strain_specific_DA.csv.gz")
    for s in d.loc[d.strain==X,"sequence"]: cand.setdefault(s,set()).add(lab)
cseqs=sorted(cand)
with open(f"{OUT}/{X}.candidates.fasta","w") as o:
    for i,s in enumerate(cseqs): o.write(f">c{i}\n{s}\n")
pd.DataFrame({"id":[f"c{i}" for i in range(len(cseqs))],"sequence":cseqs,
             "timepoints":[",".join(sorted(cand[s])) for s in cseqs],
             "n_tp":[len(cand[s]) for s in cseqs]}).to_csv(f"{OUT}/{X}.candidates.tsv.gz",index=False,compression="gzip")
print(f"[{X}] {len(cseqs):,} unique strain-specific candidate seqs -> {X}.candidates.fasta",flush=True)

# (2) X's own expressed-sequence pool (24-32 nt unique), all 9 samples
pool=set()
for f,_ in TPS:
    for r in REPS:
        with gzip.open(f"{RD}/{X}-{f}.{r}.raw.fasta.gz","rt") as fh:
            for line in fh:
                if line[:1]!=">":
                    s=line.strip()
                    if s and LO<=len(s)<=HI: pool.add(s)
        print(f"  [{X}] pooled {X}-{f}.{r} (running uniq={len(pool):,})",flush=True)
with open(f"{OUT}/{X}.pool_uniq.fasta","w") as o:
    for i,s in enumerate(sorted(pool)): o.write(f">p{i}\n{s}\n")
print(f"[{X}] {len(pool):,} unique 24-32nt expressed seqs -> {X}.pool_uniq.fasta")
