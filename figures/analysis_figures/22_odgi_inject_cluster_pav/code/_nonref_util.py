#!/usr/bin/env python3
"""Shared helper for theme-22 non-reference cluster analyses.

BUG it fixes: `nonref.bed` and `clusters.bed` are bedtools-MERGED (one interval per cluster, min-start/max-end),
but `clusters_fpm.bed` is UNMERGED (one row per timepoint per strand, original PICB coords). The previous code
flagged non-reference clusters with an EXACT (chrom,start,end) match of an fpm row against nonref.bed — which fails
whenever the merged extent equals no single fpm row. That silently dropped ~12% (171/1393) of non-ref clusters and
lost the non-matching strand/timepoint rows of even the matched clusters, understating non-ref expression and
producing an internal n_nonref inconsistency (1222 vs the canonical 1393).

FIX: aggregate the unmerged fpm rows up to the MERGED clusters.bed intervals (per-chrom searchsorted; merged
clusters are non-overlapping and each fpm row is a merge input, so it falls inside exactly one), then flag
non-reference by exact match to nonref.bed — valid at the merged level because nonref.bed is a subset of
clusters.bed. Returns one row per MERGED cluster with summed all-primary + unique FPM and the nonref flag.
"""
import pandas as pd, numpy as np


def merged_cluster_expr(CP, T22, X):
    """Per-MERGED-cluster expression for strain X.
    Returns a DataFrame: chrom, start, end, nonref (bool), allF (Σ all-primary FPM), uniqF (Σ unique FPM)."""
    fpm = pd.read_csv(f"{CP}/{X}.clusters_fpm.bed", sep="\t", header=None,
                      names=["chrom", "start", "end", "allF", "uniqF", "strand", "tp"], dtype={"chrom": str})
    mc = pd.read_csv(f"{CP}/{X}.clusters.bed", sep="\t", header=None, usecols=[0, 1, 2],
                     names=["chrom", "start", "end"], dtype={"chrom": str}).drop_duplicates().reset_index(drop=True)
    nr = pd.read_csv(f"{T22}/nonref/{X}.nonref.bed", sep="\t", header=None,
                     names=["chrom", "start", "end", "id"], dtype={"chrom": str})
    nrset = set(zip(nr.chrom.astype(str), nr.start.astype(int), nr.end.astype(int)))   # nonref ⊂ clusters (both MERGED) → exact match valid here
    mc["nonref"] = [(c, int(s), int(e)) in nrset for c, s, e in zip(mc.chrom.astype(str), mc.start, mc.end)]
    allF = np.zeros(len(mc)); uniqF = np.zeros(len(mc))
    for c, g in fpm.groupby(fpm.chrom.astype(str)):                      # assign each unmerged fpm row to its merged cluster
        sub = mc[mc.chrom.astype(str) == c]
        if not len(sub): continue
        st = sub.start.values; en = sub.end.values; idx = sub.index.values
        o = np.argsort(st); st = st[o]; en = en[o]; idx = idx[o]
        pos = np.searchsorted(st, g.start.values, side="right") - 1      # largest cluster-start <= fpm.start
        gs = g.start.values; ga = g.allF.values; gu = g.uniqF.values
        for j, k in enumerate(pos):
            if k >= 0 and gs[j] < en[k]:                                 # verify overlap (fpm row lies inside its merged cluster)
                allF[idx[k]] += ga[j]; uniqF[idx[k]] += gu[j]
    mc["allF"] = allF; mc["uniqF"] = uniqF
    return mc
