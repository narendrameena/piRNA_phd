#!/usr/bin/env python3
"""cis vs trans test for DIVERGENCE loci: is the strain-restricted expression due to TE-SEQUENCE divergence (cis)
or REGULATORY divergence (trans)? For each divergence locus, halSnps gives carrier-vs-silent SNP rate AT the locus,
compared to a matched random-region baseline (same carrier→silent pairs, same chrom, same length). locus_rate ~
baseline => the TE copy is as conserved as the genome average => REGULATORY/trans; locus_rate >> baseline => cis
TE-sequence divergence at the locus."""
import sys, subprocess, random
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
import pav_clusters as pc, pandas as pd, numpy as np, pysam
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; PG = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"; CP = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"
SIF = f"{ROOT}/cactus_v2.9.3.sif"; HAL = f"{ROOT}/results/pangenome/output/mouse_17strain_pangenome.full.hal"; ASM = lambda X: f"{ROOT}/resources/REL-2205-Assembly/{X}_chromosomes_MT.fasta"
random.seed(3)
cand = pd.read_csv(f"{PG}/divergence_candidates.tsv", sep="\t", dtype={"g39_chrom": str})
sup = pd.read_csv(f"{CP}/locus_genome_pav_divergence.tsv", sep="\t", dtype={"g39_chrom": str})
gn = sup[sup.present].groupby(["g39_chrom", "g39_start"]).size().rename("genome_n").reset_index()
div = cand.merge(gn, on=["g39_chrom", "g39_start"], how="left"); div = div[div.genome_n >= 10]
def snps(carrier, targets, oc, s, L):
    out = subprocess.run(f"singularity exec --bind /mnt {SIF} halSnps {HAL} {carrier} {','.join(targets)} --refSequence {oc} --start {s} --length {L}", shell=True, capture_output=True, text=True).stdout
    ts = tp = 0
    for ln in out.splitlines():
        f = ln.split()
        if len(f) >= 3 and f[0] in targets:
            try: ts += int(f[1]); tp += int(f[2])
            except ValueError: pass
    return ts, tp
rows = []; samp = div.sample(min(40, len(div)), random_state=3)
for _, r in samp.iterrows():
    sub = pc.clusters_at(r.g39_chrom, int(r.g39_start), int(r.g39_end)); cc = sub[sub.strain == r.carrier]
    if len(cc) == 0: continue
    oc = str(cc.iloc[0].own_chrom); s = int(cc.own_start.min()); e = int(cc.own_end.max()); L = max(500, e - s)
    expr = set(str(r.expressing).split(",")); pres = sup[(sup.g39_chrom == r.g39_chrom) & (sup.g39_start == int(r.g39_start)) & (sup.present)].strain.tolist()
    silent = [X for X in pres if X not in expr and X != r.carrier]
    if not silent: continue
    ls, lp = snps(r.carrier, silent, oc, s, L)
    try:
        fa = pysam.FastaFile(ASM(r.carrier)); clen = fa.get_reference_length(f"{r.carrier}#1#chr{oc}"); fa.close()
    except Exception: clen = None
    bs = bp = 0
    if clen and clen > L + 2_000_000:
        rs = random.randint(1_000_000, clen - L - 1_000_000); bs, bp = snps(r.carrier, silent, oc, rs, L)
    rows.append((r.g39_chrom, int(r.g39_start), r.carrier, len(silent), ls, lp, (ls / lp if lp else np.nan), bs, bp, (bs / bp if bp else np.nan)))
    print(f"  {r.carrier} chr{oc}:{s} locus_rate={ls/lp if lp else float('nan'):.4f} base_rate={bs/bp if bp else float('nan'):.4f}", flush=True)
res = pd.DataFrame(rows, columns=["g39_chrom", "g39_start", "carrier", "n_silent", "locus_snps", "locus_pairs", "locus_rate", "base_snps", "base_pairs", "base_rate"]).dropna()
res.to_csv(f"{PG}/divergence_cis_vs_trans.csv", index=False)
ratio = (res.locus_rate / res.base_rate.replace(0, np.nan)).dropna()
print(f"\nloci tested: {len(res)}")
print(f"median locus SNP rate: {res.locus_rate.median():.4f} | median baseline: {res.base_rate.median():.4f}")
print(f"locus/baseline ratio median: {ratio.median():.2f}")
print(f"REGULATORY/trans (locus_rate <= 1.5x baseline): {(res.locus_rate <= 1.5*res.base_rate).mean()*100:.0f}%")
print(f"cis-divergence (locus_rate > 2x baseline): {(res.locus_rate > 2*res.base_rate).mean()*100:.0f}%")
