#!/usr/bin/env python3
"""
PILOT: strain-specific & timepoint-specific UNIQUE piRNA identification (sequence level).
User-approved strategy:
  - Unit = collapsed unique sequences per sample (tstk collapse: >seqid-nreads), RPM-normalized,
    3 replicates kept as replicates.
  - presence(strain) = >= MIN_RPM in >= MIN_REPS of that strain's reps (occupancy).
  - strain-specific (within a timepoint) = present in strain X, ABSENT in all other strains.
  - timepoint-specific (within a strain) = present at one stage, absent at the others.
  - (edgeR/DESeq2 differential-abundance criterion + SNP-variant-vs-private-locus split are
    SEPARATE downstream steps; this builds the matrix, QC, and the presence/absence layer.)
QC FIRST (so thresholds are calibrated from data, not hard-coded): per-sample read totals,
unique-sequence counts, count-weighted length distribution, piRNA-window fraction.
Robust to partial input (smoke-test on whatever collapse files exist).
Usage: python build_unique_pirna_pilot.py [MIN_RPM] [MIN_REPS]
"""
import gzip, os, sys
from collections import defaultdict
import numpy as np, pandas as pd

RD   = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/collapse"
OUT  = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
os.makedirs(OUT, exist_ok=True)
PILOT = ["C57BL_6NJ", "CAST_EiJ", "SPRET_EiJ"]
TPS   = ["16.5dpc", "12.5dpp", "20.5dpp"]
REPS  = [1, 2, 3]
LO, HI   = 24, 32   # piRNA length window (data-driven, user-locked 2026-06-10)
MIN_REPS = int(sys.argv[1]) if len(sys.argv) > 1 else 2   # majority of 3 replicates
# Presence = detected (>=1 read) in >= MIN_REPS reps (reproducibility); NO magic RPM threshold.
# >= MIN_REPS(=2) reps implies >= 2 reads => non-singleton automatically. RPM kept for quantification.

def read_collapse(path):
    """seq -> nreads ; tstk format '>seqid-nreads' then sequence line."""
    d = {}
    with gzip.open(path, "rt") as fh:
        cur = 0
        for line in fh:
            if line and line[0] == ">":
                cur = int(line[1:].strip().rsplit("-", 1)[1])
            else:
                s = line.strip()
                if s:
                    d[s] = d.get(s, 0) + cur
    return d

def sample_path(strain, tp, rep):
    return f"{RD}/{strain}-{tp}.{rep}.raw.fasta.gz"

# ---------- load what exists + QC ----------
data = {}          # (strain,tp,rep) -> {seq:count}
qc_rows, len_rows = [], []
for s in PILOT:
    for tp in TPS:
        for r in REPS:
            p = sample_path(s, tp, r)
            if not os.path.exists(p):
                continue
            d = read_collapse(p)
            data[(s, tp, r)] = d
            tot = sum(d.values()); nuniq = len(d)
            qc_rows.append(dict(strain=s, timepoint=tp, rep=r, total_reads=tot,
                                n_unique_seq=nuniq))
            lc = defaultdict(int)
            for seq, c in d.items():
                lc[len(seq)] += c
            for L, c in lc.items():
                len_rows.append(dict(strain=s, timepoint=tp, rep=r, length=L,
                                     reads=c, frac=c / tot if tot else 0))
            print(f"  loaded {s}-{tp}.{r}: {tot:,} reads, {nuniq:,} unique", flush=True)

qc = pd.DataFrame(qc_rows)
lend = pd.DataFrame(len_rows)
if len(qc):
    qc.to_csv(f"{OUT}/QC_per_sample.csv", index=False)
    lend.to_csv(f"{OUT}/QC_length_distribution.csv", index=False)
    # data-driven piRNA window: count-weighted length profile (pooled)
    prof = lend.groupby("length")["reads"].sum()
    prof = (prof / prof.sum() * 100).round(2)
    print("\n=== QC: count-weighted length distribution (% of reads), pooled ===")
    print(prof[(prof.index >= 18) & (prof.index <= 40)].to_string())
    print("\n=== QC: per sample ===")
    print(qc.to_string(index=False))
else:
    print("No collapse files found yet."); sys.exit(0)

# ---------- presence (reproducibility-based, 24-32 nt; no magic RPM) ----------
# present(strain,tp) = sequence detected (>=1 read) in >= MIN_REPS of that strain's reps.
# mean RPM over the 3 reps retained for quantification/reporting only (absent rep counts as 0).
presence = defaultdict(lambda: defaultdict(int))     # (strain,tp) -> seq -> #reps detected
meanrpm  = defaultdict(lambda: defaultdict(float))   # (strain,tp) -> seq -> mean RPM over reps
have     = defaultdict(set)                           # (strain,tp) -> set of reps loaded
for (s, tp, r), d in data.items():
    tot = sum(d.values()); have[(s, tp)].add(r)
    for seq, c in d.items():
        if LO <= len(seq) <= HI:
            presence[(s, tp)][seq] += 1
            meanrpm[(s, tp)][seq]  += (c / tot * 1e6) / len(REPS)

def n_reps_loaded(s, tp): return len(have[(s, tp)])

# ---------- strain-specific (within a timepoint): reproducible in exactly ONE strain ----------
strain_calls = []
for tp in TPS:
    if not all(n_reps_loaded(s, tp) == len(REPS) for s in PILOT):
        print(f"\n[{tp}] not all 3 reps x 3 strains loaded -> skip strain-specific calls.")
        continue
    seqs = set().union(*[set(presence[(s, tp)]) for s in PILOT])
    for seq in seqs:
        pres = {s: presence[(s, tp)].get(seq, 0) for s in PILOT}
        xs = [s for s in PILOT if pres[s] >= MIN_REPS]   # reproducible strains
        if len(xs) == 1:                                  # present in one strain, <MIN_REPS in others
            x = xs[0]
            strain_calls.append(dict(timepoint=tp, strain=x, length=len(seq), sequence=seq,
                                     reps_present=pres[x], mean_rpm=round(meanrpm[(x, tp)][seq], 3)))
    nby = pd.Series([c["strain"] for c in strain_calls if c["timepoint"] == tp]).value_counts().to_dict()
    print(f"[{tp}] strain-specific candidates (24-32 nt, >={MIN_REPS}/3 reps) by strain: {nby}")

# ---------- timepoint-specific (within a strain): reproducible at exactly ONE stage ----------
tp_calls = []
for s in PILOT:
    if not all(n_reps_loaded(s, tp) == len(REPS) for tp in TPS):
        print(f"\n[{s}] not all 3 timepoints loaded -> skip timepoint-specific calls.")
        continue
    seqs = set().union(*[set(presence[(s, tp)]) for tp in TPS])
    for seq in seqs:
        pres = {tp: presence[(s, tp)].get(seq, 0) for tp in TPS}
        ts = [tp for tp in TPS if pres[tp] >= MIN_REPS]
        if len(ts) == 1:
            t = ts[0]
            tp_calls.append(dict(strain=s, timepoint=t, length=len(seq), sequence=seq,
                                 reps_present=pres[t], mean_rpm=round(meanrpm[(s, t)][seq], 3)))
    nby = pd.Series([c["timepoint"] for c in tp_calls if c["strain"] == s]).value_counts().to_dict()
    print(f"[{s}] timepoint-specific candidates (24-32 nt, >={MIN_REPS}/3 reps) by stage: {nby}")

if strain_calls:
    sdf = pd.DataFrame(strain_calls)
    sdf.to_csv(f"{OUT}/strain_specific_presenceAbsence_candidates.csv.gz", index=False, compression="gzip")
    print(f"\nWROTE {len(sdf)} strain-specific candidates -> strain_specific_presenceAbsence_candidates.csv.gz")
if tp_calls:
    tdf = pd.DataFrame(tp_calls)
    tdf.to_csv(f"{OUT}/timepoint_specific_presenceAbsence_candidates.csv.gz", index=False, compression="gzip")
    print(f"WROTE {len(tdf)} timepoint-specific candidates -> timepoint_specific_presenceAbsence_candidates.csv.gz")
if not strain_calls and not tp_calls:
    print("\nNo full-set calls yet (waiting for all 27 collapse files).")
else:
    print("\n(next: STAR expression-match split [strain-specific] + genomic origin for both candidate sets)")
