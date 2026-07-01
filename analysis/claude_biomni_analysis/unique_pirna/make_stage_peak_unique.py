#!/usr/bin/env python3
"""STAGE-PEAK-LENGTH unique-piRNA SUBSET (user-requested 2026-06-23, ADDITIONAL / non-destructive).

From the genuinely-unique piRNAs (conserved-but-silent + strain-private), keep ONLY the canonical
stage-characteristic piRNA length(s) per developmental timepoint -- the pre-pachytene (27 nt) vs pachytene
(30 nt) peaks, verified against our data (E16.5 & P12.5 peak 27 nt, P20.5 peak 30 nt):

   E16.5 (16.5dpc) -> 27 nt          (pre-pachytene)
   P12.5 (12.5dpp) -> 27 nt AND 30 nt (transition: both populations)
   P20.5 (20.5dpp) -> 30 nt          (pachytene)

Length = exact sequence length (user-specified single lengths, not windows). The FULL unique set is left
untouched; this writes a separate subset CSV only.

Input : unique16/final_classified_clean_2read.csv.gz  (the adopted >=2-read, 25-32 nt, WITHIN-TP set)
Output : unique16/unique_stage_peak_27_30.csv.gz
Usage : make_stage_peak_unique.py [INPUT_CSV] [OUTPUT_CSV]"""
import sys, pandas as pd
U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
INP = sys.argv[1] if len(sys.argv) > 1 else f"{U}/unique16/final_classified_clean_2read.csv.gz"
OUT = sys.argv[2] if len(sys.argv) > 2 else f"{U}/unique16/unique_stage_peak_27_30.csv.gz"
KEEP = {"16.5dpc": {27}, "12.5dpp": {27, 30}, "20.5dpp": {30}}
GU = ["unique: conserved-but-silent", "unique: strain-private locus"]
d = pd.read_csv(INP)
g = d[d.klass5.isin(GU)].copy()
g["L"] = g.sequence.str.len()
g = g[[r.L in KEEP.get(r.timepoint, set()) for r in g.itertuples()]].drop(columns=["L"])
g.to_csv(OUT, index=False)
tot = d.klass5.isin(GU).sum()
print(f"input: {INP}")
print(f"genuinely-unique in input: {tot:,}  ->  stage-peak subset: {len(g):,} ({100*len(g)/tot:.1f}%)")
print(g.groupby(["timepoint", "klass5"]).size().to_string())
print(f"wrote {OUT}")
