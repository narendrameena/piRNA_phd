#!/usr/bin/env python3
"""Step 1: data-driven mismatch cutoff (no magic number) + extract genuinely-novel SPRET piRNAs.
Logic: among SPRET candidates that DO match an expressed piRNA in C57/CAST, true SNP-variants of
the SAME piRNA follow a Poisson(lambda) over mismatches (lambda = SNP-divergence x piRNA length).
Estimate lambda from the data as n(mm=1)/n(mm=0) (= Poisson P(1)/P(0) = lambda). Cutoff k* = smallest
k with Poisson CDF(k; lambda) >= 0.999 (captures >=99.9% of true SNP-variants). genuinely-novel =
candidates with NO expressed match within k* mismatches (min-mm > k* or no match at all).
"""
import numpy as np, pandas as pd
from scipy.stats import poisson
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
sN2seq={}
with open(f"{U}/SPRET_candidates.fasta") as fh:
    sid=None
    for ln in fh:
        if ln[0]==">": sid=ln[1:].strip()
        else: sN2seq[sid]=ln.strip()
mm=pd.read_csv(f"{U}/SPRET_candidate_minmm_to_others_expressed.tsv",sep="\t",header=None,names=["sN","mm"])
seq2mm={sN2seq[s]:int(m) for s,m in zip(mm.sN,mm.mm)}
allseq=set(sN2seq.values()); n_total=len(allseq)
dist={k:sum(1 for s in allseq if seq2mm.get(s,99)==k) for k in (0,1,2,3)}
n_nomatch=n_total-sum(dist.values())
lam = dist[1]/dist[0]                      # Poisson lambda (SNP mean at piRNA loci)
div = lam/27.0                              # implied per-base divergence (piRNA length ~27)
kstar = next(k for k in range(0,6) if poisson.cdf(k,lam)>=0.999)
print(f"unique SPRET candidate seqs: {n_total}")
print(f"min-mm distribution: mm0={dist[0]} mm1={dist[1]} mm2={dist[2]} mm3={dist[3]} no_match<=3mm={n_nomatch}")
print(f"lambda (data-derived SNP mean at piRNA loci) = n1/n0 = {lam:.3f}  -> implied divergence ~{div*100:.2f}%/bp")
print(f"Poisson CDF: P(<=1)={poisson.cdf(1,lam):.4f} P(<=2)={poisson.cdf(2,lam):.4f} P(<=3)={poisson.cdf(3,lam):.4f}")
print(f"==> data-driven cutoff k* = {kstar} mismatches (captures >=99.9% of true SNP-variants)")
# genuinely novel = no expressed match within k* mismatches
def novel(seq): return seq2mm.get(seq,99) > kstar
novel_seqs=[s for s in allseq if novel(s)]
print(f"\nGENUINELY-NOVEL SPRET sequences (no expressed match <= {kstar}mm): {len(novel_seqs)} / {n_total}")
# per timepoint
cand=pd.read_csv(f"{U}/strain_specific_presenceAbsence_candidates.csv.gz")
sp=cand[cand.strain=='SPRET_EiJ'].copy(); sp['novel']=sp.sequence.map(novel)
TPMAP={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}; sp['tp']=sp.timepoint.map(TPMAP)
print("\nper timepoint (candidates -> genuinely novel):")
for tp in ["E16.5","P12.5","P20.5"]:
    g=sp[sp.tp==tp]; print(f"  {tp}: {g.novel.sum()} novel / {len(g)} candidates ({100*g.novel.mean():.1f}%)")
# write novel FASTA for genome mapping (Step 2)
with open(f"{U}/SPRET_novel.fasta","w") as out:
    for i,s in enumerate(sorted(novel_seqs)): out.write(f">n{i}\n{s}\n")
print(f"\nwrote {len(novel_seqs)} novel seqs -> SPRET_novel.fasta (for genome-origin mapping)")
pd.DataFrame({"sensitivity_cutoff_mm":[1,2,3],
   "n_novel":[sum(1 for s in allseq if seq2mm.get(s,99)>k) for k in (1,2,3)]}).to_csv(f"{U}/novel_cutoff_sensitivity.csv",index=False)
