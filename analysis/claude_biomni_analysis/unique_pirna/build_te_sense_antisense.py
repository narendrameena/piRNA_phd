#!/usr/bin/env python3
"""Per-strain, per-timepoint, per-family, per-2Mb-bin SENSE-to-TE vs ANTISENSE-to-TE small-RNA on TEs, using the
VERIFIED read-vs-TE-strand logic (antisense = ((TE_strand=='-') != read.is_reverse); sense = complement) — NOT
genomic +/- strand. PRIMARY sRNA reads only (24-32 nt; no multimapping spillover). TE annotation = RepeatMasker
own genome. This is the small-RNA basis for the circos: SENSE-to-TE = TE expression; ANTISENSE-to-TE = piRNA-on-TE.
Output cluster_pav/te_sense_antisense/{strain}.tsv (strain, chrom, bin, family, tp, sense, antisense).
Usage: build_te_sense_antisense.py <strain>   (or SLURM array index)."""
import sys, os, glob, bisect
from collections import defaultdict
import pysam
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; CP = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/cluster_pav"
STR = ["C57BL_6NJ", "BALB_cJ", "A_J", "FVB_NJ", "C3H_HeJ", "LP_J", "129S1_SvImJ", "DBA_2J", "AKR_J", "CBA_J", "NZO_HlLtJ", "NOD_ShiLtJ", "WSB_EiJ", "CAST_EiJ", "PWK_PhJ", "SPRET_EiJ"]
X = sys.argv[1] if len(sys.argv) > 1 else STR[int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))]
TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]; CHROMS = set([str(i) for i in range(1, 20)] + ["X"])
def fam_of(cf):                                   # RepeatMasker class/family -> coarse family for the circos
    if cf.startswith("LINE/L1"): return "L1"
    if "ERVK" in cf or "IAP" in cf: return "ERVK"
    if cf.startswith("LTR/ERVL-MaLR") or cf.startswith("LTR/ERVL"): return "ERVL"
    if cf.startswith("SINE"): return "SINE"
    if cf.startswith("LTR"): return "LTR_other"
    if cf.startswith("LINE"): return "LINE_other"
    return None
# load RepeatMasker -> per chrom sorted (start, end, te_strand, family)
outf = glob.glob(f"{ROOT}/resources/repeatMasker/{X}_*.out")
assert outf, f"no RepeatMasker for {X}"
te = defaultdict(list)
with open(outf[0]) as fh:
    for i, ln in enumerate(fh):
        if i < 3: continue
        p = ln.split()
        if len(p) < 11: continue
        chrom = p[4].split("#")[-1]; chrom = chrom[3:] if chrom.startswith("chr") else chrom
        if chrom not in CHROMS: continue
        fam = fam_of(p[10])
        if fam is None: continue
        st = "-" if p[8] == "C" else "+"
        te[chrom].append((int(p[5]), int(p[6]), st, fam))
for c in te: te[c].sort()
starts = {c: [t[0] for t in te[c]] for c in te}
def te_at(chrom, pos):                            # TE overlapping pos (nearest by start, check end)
    arr = te.get(chrom);
    if not arr: return None
    i = bisect.bisect_right(starts[chrom], pos) - 1
    for j in range(i, max(-1, i - 6), -1):        # a few back (TEs can nest/overlap)
        if j < 0: break
        s, e, stt, fam = arr[j]
        if s <= pos < e: return stt, fam
    return None
out = defaultdict(lambda: [0, 0])                 # (chrom,bin,fam,tp) -> [sense, antisense]
for tp in TPS:
    for r in (1, 2, 3):
        b = f"{ROOT}/results/STAR_srna_strain_wise/{X}/{X}-{tp}.{r}/Aligned.sortedByCoord.out.bam"
        if not os.path.exists(b): continue
        try: bam = pysam.AlignmentFile(b, "rb")
        except Exception: continue
        for a in bam.fetch():
            if a.is_unmapped or a.is_secondary or not a.query_sequence: continue   # PRIMARY only
            if not (25 <= a.reference_end - a.reference_start <= 32): continue
            chrom = a.reference_name.split("#")[-1]; chrom = chrom[3:] if chrom.startswith("chr") else chrom
            if chrom not in CHROMS: continue
            mid = (a.reference_start + a.reference_end) // 2
            hit = te_at(chrom, mid)
            if hit is None: continue
            st, fam = hit
            anti = ((st == "-") != a.is_reverse)   # VERIFIED antisense-to-TE
            out[(chrom, mid // 2_000_000, fam, tp)][1 if anti else 0] += 1
        bam.close()
        print(f"[{X}] {tp}.{r} done", flush=True)
os.makedirs(f"{CP}/te_sense_antisense", exist_ok=True)
with open(f"{CP}/te_sense_antisense/{X}.tsv", "w") as o:
    for (c, bn, fam, tp), (se, an) in sorted(out.items()):
        o.write(f"{X}\t{c}\t{bn}\t{fam}\t{tp}\t{se}\t{an}\n")
print(f"[{X}] wrote {len(out)} (chrom,bin,fam,tp) rows", flush=True)
