#!/usr/bin/env python3
"""
Pangenome RETENTION metric for piRNA loci  [cactus halLiftover, NOT UCSC liftOver].

For each reference locus, bin its span into 100-bp windows (in the source genome
frame) and halLiftover every window into each strain through the 17-strain cactus
pangenome graph. The fraction of a locus's windows that project =
`fraction_lifted` = a graded, source-span-based RETENTION score.

This REPLACES the binary UCSC-liftOver 'not_lifted' call, which conflated rapid
sequence divergence with structural absence (UCSC liftOver needs ~95% high-identity
mapping; pachytene piRNA loci sit at conserved positions but diverge in sequence,
Yu et al. 2021 Nat Commun PMID 33397987 — so they project through the pangenome
but fail chain-based liftOver).

Usage:
  pangenome_fraction_lifted.py  <ref_bed>  <src_genome>  <out_csv>  <workdir>
    ref_bed     : chr<TAB>start<TAB>end<TAB>locus_id  (source-frame, bare-numeric chr)
    src_genome  : HAL genome name of the source frame ("GRCm39" or "C57BL_6NJ")
    out_csv     : output table  (locus, strain, n_bins, n_bins_lifted, fraction_lifted)
    workdir     : scratch dir for the binned + lifted BEDs
"""
import os, sys, subprocess
import pandas as pd

BASE = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
SIF  = f"{BASE}/cactus_v2.9.3.sif"
HAL  = f"{BASE}/results/pangenome/output/mouse_17strain_pangenome.full.hal"
STRAINS = ["129S1_SvImJ","AKR_J","A_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J",
           "DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
BIN = 100

ref_bed, src_genome, out_csv, workdir = sys.argv[1:5]
os.makedirs(workdir, exist_ok=True)

ref = pd.read_csv(ref_bed, sep="\t", header=None, names=["chr","start","end","locus"])
ref["chr"] = ref["chr"].astype(str)

# ── bin every locus into BIN-bp windows (source frame) ────────────────────────
rows, totbins = [], {}
for _, r in ref.iterrows():
    s, e = int(r.start), int(r.end)
    for i, bs in enumerate(range(s, e, BIN)):
        be = min(bs + BIN, e)
        rows.append((r.chr, bs, be, f"{r.locus}||b{i}"))
        totbins[r.locus] = totbins.get(r.locus, 0) + 1
binbed = f"{workdir}/_bins_src.bed"
pd.DataFrame(rows, columns=["chr","start","end","name"]).to_csv(binbed, sep="\t", header=False, index=False)
print(f"binned {len(totbins)} loci into {len(rows)} windows ({BIN} bp)")

# ── halLiftover bins src -> each strain; count windows that project ───────────
res = []
for strain in STRAINS:
    if strain == src_genome:
        lifted = dict(totbins)                       # identity: all windows present
    else:
        out = f"{workdir}/_bins_{strain}.bed"
        subprocess.run(["singularity","exec","--bind","/mnt",SIF,"halLiftover",
                        HAL, src_genome, binbed, strain, out], capture_output=True)
        present = set()
        if os.path.exists(out) and os.path.getsize(out) > 0:
            d = pd.read_csv(out, sep="\t", header=None, usecols=[3], names=["name"])
            present = set(d["name"].unique())
        lifted = {}
        for nm in present:
            loc = nm.rsplit("||b", 1)[0]
            lifted[loc] = lifted.get(loc, 0) + 1
    for loc, nb in totbins.items():
        lb = lifted.get(loc, 0)
        res.append({"locus": loc, "strain": strain, "n_bins": nb,
                    "n_bins_lifted": lb, "fraction_lifted": lb / nb})
    print(f"  {strain}: mean fraction_lifted = "
          f"{sum(lifted.get(l,0)/totbins[l] for l in totbins)/len(totbins):.3f}")

df = pd.DataFrame(res)
df.to_csv(out_csv, index=False)
print(f"\nwrote {out_csv}: {len(df)} rows")
print(f"overall RETAINED (fraction>=0.5): {(df.fraction_lifted>=0.5).mean()*100:.1f}%  | "
      f"DISRUPTED (<0.5): {(df.fraction_lifted<0.5).mean()*100:.1f}%  | "
      f"fully absent (==0): {(df.fraction_lifted==0).mean()*100:.1f}%")
