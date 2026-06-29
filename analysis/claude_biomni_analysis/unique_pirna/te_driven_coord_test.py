#!/usr/bin/env python3
"""COORDINATE-BASED TE-driven test (NOT sequence containment). A piRNA cluster is TE-INSERTION-DRIVEN iff its
mapping COORDINATES fall inside a strain-PRIVATE insertion's COORDINATES in the strain genome. Private insertions
(g39 breakpoints, from the pangenome VCF, parse_insertions16.py) are lifted g39->strain (one-time pangenome step)
to get strain-coordinate intervals; these are intersected with the strain's PICB clusters (clusters_fpm.bed, own
coords). Output: clusters that COORDINATE-overlap a private TE insertion. Usage: te_driven_coord_test.py <STRAIN>"""
import os, sys, glob, subprocess
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"; CP = f"{U}/cluster_pav"
SIF = f"{ROOT}/cactus_v2.9.3.sif"; HAL = f"{ROOT}/results/pangenome/output/mouse_17strain_pangenome.full.hal"
BT = "/mnt/home3/miska/nm667/miniconda3/envs/biomni_e1/bin/bedtools"
X = sys.argv[1] if len(sys.argv) > 1 else "CAST_EiJ"; MININS = 150
MODE = sys.argv[2] if len(sys.argv) > 2 else "pilot"   # "ins16" => 16-strain-private (reads ins16/, writes TE_driven_COORDINATE16_{X}.csv); preserves the 3-strain pilot outputs
INS = f"{PG}/ins16/{X}.private_insertions.fasta" if MODE == "ins16" else f"{PG}/{X}.private_insertions.fasta"
# 1) parse private insertions (g39 breakpoints) >= MININS bp
ins = []
for ln in open(INS):
    if ln[0] != ">": continue
    c, p, L = ln[1:].strip().rsplit("_", 2)
    if int(L) >= MININS: ins.append((c, int(p), int(L)))
print(f"[{X}] private insertions >= {MININS}bp: {len(ins)}", flush=True)
# 2) lift g39 breakpoint windows -> strain (one halLiftover call); name carries id+inslen
inb = f"{CP}/.tdc_{X}_g39.bed"
with open(inb, "w") as o:
    for c, p, L in ins: o.write(f"{c}\t{max(0,p-25)}\t{p+25}\t{c}:{p}:{L}\n")
outb = f"{CP}/.tdc_{X}_strain.bed"
subprocess.run(f"singularity exec --bind /mnt {SIF} halLiftover {HAL} GRCm39 {inb} {X} {outb}", shell=True, capture_output=True)
# 3) strain-coordinate insertion intervals (expand by inslen to cover the inserted body)
reg = {}
for ln in open(outb):
    f = ln.rstrip("\n").split("\t")
    if len(f) < 4: continue
    ch = f[0][3:] if f[0].startswith("chr") else f[0]; nm = f[3]; L = int(nm.split(":")[2])
    s = max(0, int(f[1]) - L - 50); e = int(f[2]) + L + 50
    if nm in reg: reg[nm] = (ch, min(reg[nm][1], s), max(reg[nm][2], e), L)
    else: reg[nm] = (ch, s, e, L)
insbed = f"{CP}/.tdc_{X}_insreg.bed"
with open(insbed, "w") as o:
    for nm, (ch, s, e, L) in reg.items(): o.write(f"{ch}\t{s}\t{e}\t{nm}|{L}\n")
print(f"[{X}] insertions placed in strain coords: {len(reg)}", flush=True)
# 4) intersect strain PICB clusters with strain-coordinate private insertions
clb = f"{CP}/.tdc_{X}_clusters.bed"
with open(clb, "w") as o:
    for ln in open(f"{CP}/{X}.clusters_fpm.bed"):
        f = ln.split("\t"); o.write(f"{f[0]}\t{f[1]}\t{f[2]}\t{f[3]}:{f[5]}:{f[6].strip()}\n")  # chrom start end allFPM:strand:tp
hits = subprocess.run(f"sort -k1,1 -k2,2n {clb} > {clb}.s; sort -k1,1 -k2,2n {insbed} > {insbed}.s; {BT} intersect -a {clb}.s -b {insbed}.s -wo", shell=True, capture_output=True, text=True).stdout
# 5) report clusters overlapping a private insertion (coordinate-based TE-driven), ranked by FPM
rows = []
for ln in hits.splitlines():
    f = ln.split("\t"); fpm, strand, tp = f[3].split(":"); inslen = f[7].split("|")[1]; ov = int(f[8])
    rows.append((float(fpm), f[0], int(f[1]), int(f[2]), strand, tp, int(inslen), ov))
rows.sort(reverse=True)
print(f"[{X}] PICB clusters coordinate-overlapping a PRIVATE insertion: {len(set((r[1],r[2]) for r in rows))}", flush=True)
print(f"\n{'cluster (own)':<26}{'FPM':>8}{'strand':>7}{'tp':>9}{'ins_bp':>8}{'overlap':>8}")
seen = set()
for fpm, c, s, e, st, tp, L, ov in rows:
    if (c, s) in seen: continue
    seen.add((c, s))
    print(f"chr{c}:{s:,}-{e:,}".ljust(26) + f"{fpm:>8.1f}{st:>7}{tp:>9}{L:>8}{ov:>8}")
    if len(seen) >= 15: break
import pandas as pd
_outcsv = f"{PG}/TE_driven_COORDINATE16_{X}.csv" if MODE == "ins16" else f"{PG}/TE_driven_COORDINATE_{X}.csv"
pd.DataFrame(rows, columns=["FPM","chrom","start","end","strand","tp","ins_bp","overlap_bp"]).to_csv(_outcsv, index=False)
print(f"\nwrote {_outcsv}")
for f in [inb, outb, insbed, clb, clb+".s", insbed+".s"]:
    try: os.remove(f)
    except OSError: pass
