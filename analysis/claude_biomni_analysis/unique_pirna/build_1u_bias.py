#!/usr/bin/env python3
"""Per-2Mb-bin 1U bias (5'-Uridine fraction) of piRNA-sized small RNAs, per strain x timepoint — input for the
1U-signature circos. PARALLEL across BAMs (16 workers): the sRNA BAMs are huge (22x multimapping), so reading is
the bottleneck. Per BAM: samtools view -F 0x104 (primary mapped only) -s 0.1 (10% subsample; 1U FRACTION is
unaffected, only depth scales) -> 5' base (reverse = complement of last base) -> count 1U vs total for 24-32 nt
reads per bin. Output active_1u_bias_tp.tsv (strain, tp, chrom, bin, 1Ufrac, subsampled_depth)."""
import glob, os, subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"
SAMTOOLS = sorted(glob.glob("/mnt/home3/miska/nm667/miniconda3/envs/*/bin/samtools"))[0]
CHROMS = set([str(i) for i in range(1, 20)] + ["X"]); TPMAP = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}
comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
def process(args):
    strain, tp, bam = args; cnt = defaultdict(lambda: [0, 0, 0])   # 1U, total, 10A
    p = subprocess.Popen([SAMTOOLS, "view", "-@", "2", "-F", "0x104", "-s", "0.1", bam], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        f = line.split("\t", 11)
        if len(f) < 10: continue
        seq = f[9]; L = len(seq)
        if not (24 <= L <= 32): continue
        c = f[2].split("#")[-1]; c = c[3:] if c.startswith("chr") else c
        if c not in CHROMS: continue
        rev = int(f[1]) & 16
        b5 = comp.get(seq[-1].upper(), "N") if rev else seq[0].upper()                       # 5' base (1U)
        b10 = (comp.get(seq[-10].upper(), "N") if rev else seq[9].upper()) if L >= 10 else "N" # 10th base (10A = ping-pong)
        key = (strain, tp, c, (int(f[3]) - 1) // 2_000_000); cnt[key][1] += 1
        if b5 == "T": cnt[key][0] += 1
        if b10 == "A": cnt[key][2] += 1
    p.wait()
    return {k: v for k, v in cnt.items()}
if __name__ == "__main__":
    bams = []
    for strain in os.listdir(f"{ROOT}/results/STAR_srna_strain_wise"):
        for tpk, tpv in TPMAP.items():
            b = sorted(glob.glob(f"{ROOT}/results/STAR_srna_strain_wise/{strain}/{strain}-{tpk}.1/Aligned.sortedByCoord.out.bam"))
            if b: bams.append((strain, tpv, b[0]))
    print(f"{len(bams)} BAMs, 16 parallel workers, primary-only + 10% subsample", flush=True)
    merged = defaultdict(lambda: [0, 0, 0]); done = 0
    with ProcessPoolExecutor(max_workers=16) as ex:
        for d in ex.map(process, bams):
            for k, (u, t, a) in d.items(): merged[k][0] += u; merged[k][1] += t; merged[k][2] += a
            done += 1; print(f"  {done}/{len(bams)} BAMs merged", flush=True)
    with open(f"{U}/active_1u_bias_tp.tsv", "w") as o:   # cols: strain tp chrom bin 1Ufrac depth 10Afrac
        for (s, tp, c, b), (u1, tot, a10) in sorted(merged.items()):
            if tot >= 5: o.write(f"{s}\t{tp}\t{c}\t{b}\t{u1/tot:.4f}\t{tot}\t{a10/tot:.4f}\n")
    print(f"wrote active_1u_bias_tp.tsv (1U+10A): {sum(1 for v in merged.values() if v[1]>=5)} bins", flush=True)
