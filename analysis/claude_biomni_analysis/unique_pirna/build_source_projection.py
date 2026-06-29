#!/usr/bin/env python3
"""Batched per-carrier HAL projection of selected source loci into all 16 strains (genome presence + own coords),
so make_source_pav_multi reads it instead of doing per-locus liftovers. 16 carriers x 15 targets = 240 calls.
Usage: build_source_projection.py <master.tsv> <out_projection.tsv>"""
import sys, os, subprocess, tempfile, pandas as pd
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; PG = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
SIF = f"{ROOT}/cactus_v2.9.3.sif"; HAL = f"{ROOT}/results/pangenome/output/mouse_17strain_pangenome.full.hal"
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]
master = sys.argv[1] if len(sys.argv) > 1 else f"{PG}/source_loci_master_insertion.tsv"
outf = sys.argv[2] if len(sys.argv) > 2 else f"{PG}/source_projection_insertion.tsv"
M = pd.read_csv(master, sep="\t")
recs = []
for carrier, grp in M.groupby("carrier"):
    bed = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False); locinfo = {}
    for _, r in grp.iterrows():
        ch = r.chrom_own.split("chr")[-1]; lid = f"{carrier}:{ch}:{r.start}"
        bed.write(f"{ch}\t{r.start}\t{r.end}\t{lid}\n"); locinfo[lid] = (ch, int(r.start), int(r.end))
    bed.close()
    for lid, (ch, s, e) in locinfo.items(): recs.append((lid, carrier, carrier, True, ch, s, e))   # carrier present at its own locus
    for X in ORDER:
        if X == carrier: continue
        out = subprocess.run(f"singularity exec --bind /mnt {SIF} halLiftover {HAL} {carrier} {bed.name} {X} /dev/stdout", shell=True, capture_output=True, text=True).stdout
        agg = {}
        for ln in out.splitlines():
            f = ln.split("\t")
            if len(f) < 4: continue
            lid = f[3]; ch = f[0]; a = int(f[1]); b = int(f[2])
            agg.setdefault(lid, {}).setdefault(ch, [10**18, 0, 0])
            g = agg[lid][ch]; g[0] = min(g[0], a); g[1] = max(g[1], b); g[2] += (b - a)
        for lid, (ch0, s0, e0) in locinfo.items():
            win = max(1, e0 - s0); done = False
            if lid in agg:
                bestch = max(agg[lid], key=lambda c: agg[lid][c][2]); g = agg[lid][bestch]
                if g[2] >= 0.30 * win: recs.append((lid, carrier, X, True, bestch, g[0], g[1])); done = True
            if not done: recs.append((lid, carrier, X, False, "", 0, 0))
    os.unlink(bed.name)
    print(f"[{carrier}] projected {len(locinfo)} loci -> 15 strains", flush=True)
pd.DataFrame(recs, columns=["locus_id", "carrier", "target", "present", "tgt_chrom", "tgt_start", "tgt_end"]).to_csv(outf, sep="\t", index=False)
print(f"wrote {outf}: {len(recs)} records ({M.carrier.nunique()} carriers)")
