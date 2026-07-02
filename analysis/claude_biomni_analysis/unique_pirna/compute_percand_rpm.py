#!/usr/bin/env python3
"""Per-piRNA-SEQUENCE expression (RPM) for the genuinely-unique candidates, the biologically-standard
piRNA abundance metric: RPM = read count / library size (24-32 nt window) x 1e6, per (sequence, strain,
timepoint), averaged over the strain's replicates. Reads the edger16 per-sequence count matrices
(6.7M seqs x 48 samples per timepoint) in chunks, keeping only candidate rows."""
import pandas as pd, numpy as np, gzip
U = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
GU = {"unique: strain-private locus", "unique: conserved-but-silent"}
k = pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz", usecols=["sequence","strain","timepoint","klass5"])
gu = k[k.klass5.isin(GU)][["sequence","strain","timepoint"]].drop_duplicates()
allrows = []
for tp in ["16.5dpc", "12.5dpp", "20.5dpp"]:
    cand = gu[gu.timepoint == tp]
    cand_seqs = set(cand.sequence)
    if not cand_seqs:
        continue
    samp = pd.read_csv(f"{U}/edger16/{tp}.samples.tsv", sep="\t")
    lib = dict(zip(samp["sample"], samp["libsize_window"]))
    seqs = np.array([l.rstrip("\n") for l in gzip.open(f"{U}/edger16/{tp}.seqs.txt.gz", "rt")], dtype=object)
    kept_seq, kept, off, cols = [], [], 0, None
    for chunk in pd.read_csv(f"{U}/edger16/{tp}.counts.tsv.gz", sep="\t", chunksize=1_000_000):
        if cols is None:
            cols = list(chunk.columns)
        sc = seqs[off:off+len(chunk)]
        mask = pd.Series(sc).isin(cand_seqs).values
        if mask.any():
            kept_seq.append(sc[mask]); kept.append(chunk.values[mask])
        off += len(chunk)
    ks = np.concatenate(kept_seq); kc = np.vstack(kept).astype(float)
    libv = np.array([lib[c] for c in cols], dtype=float); libv[libv == 0] = 1.0
    rpm = pd.DataFrame(kc / libv * 1e6, columns=cols, index=ks)
    colstrain = [c.rsplit(".", 1)[0] for c in cols]                 # strain.rep -> strain
    rpm = rpm.T.groupby(colstrain).mean().T                          # mean RPM over the strain's replicates
    dfl = rpm.stack().rename("rpm").reset_index(); dfl.columns = ["sequence","strain","rpm"]
    cc = cand.merge(dfl, on=["sequence","strain"], how="left"); cc["rpm"] = cc.rpm.fillna(0.0)
    allrows.append(cc)
    print(f"[{tp}] {len(cand):,} candidates; median RPM {cc.rpm.median():.3f}, max {cc.rpm.max():.1f}")
out = pd.concat(allrows, ignore_index=True)
out.to_csv(f"{U}/unique16/percand_rpm_expr.csv.gz", index=False)
print(f"wrote percand_rpm_expr.csv.gz: {len(out):,} candidate rows")
