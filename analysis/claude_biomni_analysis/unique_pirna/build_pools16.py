#!/usr/bin/env python3
"""Scale-16 pangenome phase — PREP: per-strain EXPRESSED piRNA POOL.
The pool = the set of 24-32 nt unique sequences a strain DETECTABLY expresses (>=2 reads total = non-
singleton), pooled over its 3 timepoints x 3 reps collapse files. This is the reference for "is this
sequence expressed in strain Y" in the cross-strain expression test (the 16-strain Step-4 equivalent).
Built directly from collapse (NOT from the min.total>=15 union matrix) so a sequence expressed only in
one strain is not lost. Usage: <strain>"""
import gzip, os, sys
RD="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/collapse"
OUT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna/unique16/pools"
os.makedirs(OUT,exist_ok=True)
TPS=["16.5dpc","12.5dpp","20.5dpp"]; REPS=[1,2,3]; LO,HI=24,32; MIN=2
s=sys.argv[1]; tot={}
for tp in TPS:
    for r in REPS:
        f=f"{RD}/{s}-{tp}.{r}.raw.fasta.gz"
        if not os.path.exists(f): continue
        with gzip.open(f,"rt") as fh:
            cur=0
            for line in fh:
                if line[:1]==">": cur=int(line[1:].strip().rsplit("-",1)[1])
                else:
                    seq=line.strip()
                    if LO<=len(seq)<=HI: tot[seq]=tot.get(seq,0)+cur
pool=[seq for seq,c in tot.items() if c>=MIN]
with gzip.open(f"{OUT}/{s}.pool.txt.gz","wt") as o: o.write("\n".join(pool)+"\n")
print(f"[{s}] {len(tot):,} distinct 24-32nt -> pool {len(pool):,} (>={MIN} reads non-singleton)")
