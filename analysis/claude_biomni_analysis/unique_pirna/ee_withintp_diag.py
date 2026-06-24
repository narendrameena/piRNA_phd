#!/usr/bin/env python3
"""Diagnose the within-tp 'expressed elsewhere' situation that produced EE=0 + CBS explosion.

For each adopted candidate (strain X, tp T, seq), determine whether the EXACT seq is expressed (>=2 reads)
in ANOTHER strain at the SAME stage T (EE-same-stage) and/or ONLY at a DIFFERENT stage (EE-other-stage-only).
One streaming pass over the 48 per-tp pools, keeping only candidate-seq hits (memory-light).

Reads the CROSS-TP backup (the well-defined baseline) so the cross-tab is against the established klass5."""
import gzip, collections, pandas as pd
U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
PP = f"{U}/pools_pertp"
STRAINS = ["C57BL_6NJ", "BALB_cJ", "A_J", "FVB_NJ", "C3H_HeJ", "LP_J", "129S1_SvImJ", "DBA_2J", "AKR_J",
           "CBA_J", "NZO_HlLtJ", "NOD_ShiLtJ", "WSB_EiJ", "CAST_EiJ", "PWK_PhJ", "SPRET_EiJ"]
TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]
d = pd.read_csv(f"{U}/unique16/final_classified_clean_2read.crosstp_backup.csv.gz",
                usecols=["sequence", "strain", "timepoint", "klass5"])
cand_seqs = set(d.sequence)
print(f"{len(d):,} candidates, {len(cand_seqs):,} distinct sequences", flush=True)
seq_in = collections.defaultdict(set)   # seq -> {(strain, tp)}
for s in STRAINS:
    for tp in TPS:
        with gzip.open(f"{PP}/{s}_{tp}.pool.txt.gz", "rt") as fh:
            for line in fh:
                q = line.rstrip("\n")
                if q in cand_seqs:
                    seq_in[q].add((s, tp))
        print(f"  scanned {s} {tp}", flush=True)
print(f"distinct candidate seqs seen in some pool: {len(seq_in):,}", flush=True)
def ee_status(r):
    locs = seq_in.get(r.sequence, ())
    same = any(Y != r.strain and (Y, r.timepoint) in locs for Y in STRAINS)
    if same: return "EE-same-stage"
    other = any(Y != r.strain for (Y, _tp) in locs)
    return "EE-other-stage-only" if other else "not-expressed-elsewhere"
d["ee"] = d.apply(ee_status, axis=1)
print("\n=== within-tp EE status (all candidates) ===")
print(d.ee.value_counts().to_string())
print("\n=== cross-tab: cross-tp klass5 (rows) x within-tp EE status (cols) ===")
print(pd.crosstab(d.klass5, d.ee).to_string())
d[["sequence", "strain", "timepoint", "klass5", "ee"]].to_csv(f"{U}/unique16/ee_withintp_diag.csv.gz", index=False)
print(f"\nwrote {U}/unique16/ee_withintp_diag.csv.gz")
