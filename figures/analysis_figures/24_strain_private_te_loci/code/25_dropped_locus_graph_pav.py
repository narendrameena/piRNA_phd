#!/usr/bin/env python3
"""THEME 24 — recover the graph PAV for a DROPPED (fragment-boundary-spanning) private cluster that odgi inject skipped
(e.g. the 82-kb / 31,188-FPM WSB chr7 private cluster). No re-inject: map the native cluster onto the strain's graph
FRAGMENTS (STRAIN#0#chr#offset, offsets from graph_path_lengths.txt) + local coords, run odgi pav on those fragment
regions IN THE FRAGMENT FRAME (mode bed -> dropped_locus_frag.bed; run step-20 odgi pav -> dropped_locus_pav.tsv), then
length-weighted-merge the per-fragment coverage into the cluster's per-strain coverage and append it to
graph_check_pav3.tsv so make_nonref_locus.py can draw panel A.
Usage: 25_dropped_locus_graph_pav.py <bed|merge> <strain> <chrom> <start> <end> <cid>"""
import sys, pandas as pd
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
D = f"{ROOT}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"
mode, strain, chrom, start, end, cid = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
if mode == "bed":
    frags = []
    for l in open(f"{D}/graph_path_lengths.txt"):
        p, ln = l.rstrip("\n").split("\t"); pr = p.split("#")
        if len(pr) >= 4 and pr[0] == strain and pr[2] == chrom: frags.append((int(pr[3]), int(ln), p))
    frags.sort()
    with open(f"{D}/dropped_locus_frag.bed", "w") as o:
        n = 0
        for off, ln, p in frags:
            if off + ln <= start or off >= end: continue
            ls, le = max(0, start - off), min(ln, end - off)
            if le > ls: o.write(f"{p}\t{ls}\t{le}\tSEG{n}\n"); n += 1
    print(f"{cid} ({strain} chr{chrom}:{start:,}-{end:,}) -> {n} fragment segment(s); wrote dropped_locus_frag.bed")
elif mode == "merge":
    df = pd.read_csv(f"{D}/dropped_locus_pav.tsv", sep="\t"); df["len"] = df.end - df.start
    hdr = open(f"{D}/graph_check_pav3.tsv").readline().rstrip("\n").split("\t"); samples = hdr[4:]
    tot = df.len.sum(); merged = {s: float((df[s] * df.len).sum() / tot) for s in samples}
    vals = [str(chrom), "0", str(end - start), cid] + [f"{merged.get(s, 0.0):.5f}" for s in samples]
    with open(f"{D}/graph_check_pav3.tsv", "a") as o: o.write("\t".join(vals) + "\n")
    print(f"appended {cid}: carrier({strain})={merged.get(strain,0):.2f}  max_other={max(v for s,v in merged.items() if s!=strain and s!='GRCm39'):.2f}  GRCm39={merged.get('GRCm39',0):.2f}  (merged {len(df)} segments)")
