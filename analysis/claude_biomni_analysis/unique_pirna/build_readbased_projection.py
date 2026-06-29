#!/usr/bin/env python3
"""Read-based genome-PAV projection for source-locus figures (accurate ●/○ for CREATION loci, where the whole-window
HAL projection wrongly counts conserved FLANKS as 'present'). For each locus: carrier present; each other strain is
'present' (the piRNA SOURCE exists there) only if its orthologous position has primary-piRNA reads >= max(ABS, REL x
carrier primary) — i.e. it actually PRODUCES piRNA, not merely shares flanking sequence. Output matches
build_source_projection.py columns so make_source_pav_multi reads it directly.
Usage: build_readbased_projection.py <master.tsv> <out_projection.tsv>"""
import sys, os, subprocess, tempfile
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis"); sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna")
import pav_clusters as pc, pandas as pd
from strain_order import STRAIN_ORDER
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; SIF = f"{ROOT}/cactus_v2.9.3.sif"; HAL = f"{ROOT}/results/pangenome/output/mouse_17strain_pangenome.full.hal"; PG = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna/pangenome_te"
ORDER = [s for s in STRAIN_ORDER if s != "C57BL_6"]; TPS = ["16.5dpc", "12.5dpp", "20.5dpp"]
REL = 0.10; ABS = 50
master = sys.argv[1]; outf = sys.argv[2]
M = pd.read_csv(master, sep="\t")
def reads_at(X, ch, s, e): return sum((pc.fetch_primary(X, ch, s, e, tp, 50) or {}).get("ntot", 0) for tp in TPS)
recs = []
for r in M.itertuples():
    carrier = r.carrier; ch = r.chrom_own.split("chr")[-1]; S, E = int(r.start), int(r.end); lid = f"{carrier}:{ch}:{S}"
    cr = reads_at(carrier, ch, S, E); thr = max(ABS, REL * cr)
    recs.append((lid, carrier, carrier, True, ch, S, E))   # carrier always present
    bed = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False); bed.write(f"{ch}\t{S}\t{E}\tl\n"); bed.close()
    npres = 1
    for Y in ORDER:
        if Y == carrier: recs[-1]; continue
        out = subprocess.run(f"singularity exec --bind /mnt {SIF} halLiftover {HAL} {carrier} {bed.name} {Y} /dev/stdout", shell=True, capture_output=True, text=True).stdout
        frs = [(f[0], int(f[1]), int(f[2])) for f in (l.split("\t") for l in out.splitlines()) if len(f) >= 3]
        if not frs: recs.append((lid, carrier, Y, False, "", 0, 0)); continue
        ch2 = max(set(f[0] for f in frs), key=lambda c: sum(f[2] - f[1] for f in frs if f[0] == c))
        fr = [f for f in frs if f[0] == ch2]; s2 = min(f[1] for f in fr); e2 = max(f[2] for f in fr)
        if reads_at(Y, ch2, s2, e2) >= thr: recs.append((lid, carrier, Y, True, ch2, s2, e2)); npres += 1
        else: recs.append((lid, carrier, Y, False, ch2, s2, e2))
    os.unlink(bed.name)
    print(f"[{lid}] carrier_primary={cr} thr={thr:.0f} -> read-based present in {npres}/16", flush=True)
pd.DataFrame(recs, columns=["locus_id", "carrier", "target", "present", "tgt_chrom", "tgt_start", "tgt_end"]).to_csv(outf, sep="\t", index=False)
print(f"wrote {outf}: {len(recs)} records ({M.shape[0]} loci)")
