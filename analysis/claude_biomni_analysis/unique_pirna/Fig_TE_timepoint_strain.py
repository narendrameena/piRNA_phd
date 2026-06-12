#!/usr/bin/env python3
"""At EACH timepoint, which TE families drive the strain-private piRNAs? Per strain, intersect strain-private-locus
piRNA loci (cand_self16 BAM, own genome) with the per-strain RepeatMasker BED, assign each locus its largest-overlap
TE family AND keep its timepoint (locus id = strain|timepoint|sequence). Output: per-locus table (strain, tp, family)
+ three faceted heatmaps (E16.5 / P12.5 / P20.5), TE family x 16 strains, cell = # strain-private piRNAs (log colour,
white = none). Reveals the developmental TE-silencing dynamic: pre-pachytene (E16.5) IAP/L1/ERVK vs later windows.
Descriptive; biological framing queued for BioMNI. (Same intersect method as Fig_TE_private_families16, timepoint-split.)"""
import warnings; warnings.filterwarnings("ignore")
import sys, subprocess, tempfile, os
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
import pandas as pd, numpy as np, pysam
from collections import Counter, defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"
BT = "/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"; CANON = [s for s in STRAIN_ORDER if s != "C57BL_6"]
TPMAP = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}; TPS = ["E16.5", "P12.5", "P20.5"]
d = pd.read_csv(f"{U}/unique16/final_classified.csv.gz")
rows = []   # (strain, tp, family)
for X in CANON:
    g = d[(d.strain == X) & (d.klass == "unique: strain-private locus")].copy(); g["id"] = X + "|" + g.timepoint + "|" + g.sequence
    sp = set(g.id); rm = f"{ROOT}/resources/repeatMasker/{X}_repeatmasker.bed"
    if not sp or not os.path.exists(rm): continue
    bam = pysam.AlignmentFile(f"{U}/cand_self16/{X}.cand_self16.bam", "rb")
    tb = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False, dir=PG); seen = set()
    for a in bam.fetch(until_eof=True):
        if a.is_unmapped or a.query_name not in sp or a.query_name in seen: continue
        seen.add(a.query_name); tb.write(f"{a.reference_name.split('#')[-1]}\t{a.reference_start}\t{a.reference_end}\t{a.query_name}\n")
    tb.close(); bam.close()
    out = subprocess.run(f"sort -k1,1 -k2,2n {tb.name} | {BT} intersect -a - -b {rm} -wa -wb", shell=True, capture_output=True, text=True).stdout
    os.unlink(tb.name)
    best = {}
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 8: continue
        ov = min(int(f[2]), int(f[6])) - max(int(f[1]), int(f[5])); fam = f[7].split("|")[-1] if "|" in f[7] else f[7]
        if f[3] not in best or ov > best[f[3]][0]: best[f[3]] = (ov, fam)
    for locid, (ov, fam) in best.items():
        rows.append((X, TPMAP.get(locid.split("|")[1], locid.split("|")[1]), fam))
    print(f"{X}: TE-annotated loci={len(best):,}")
T = pd.DataFrame(rows, columns=["strain", "tp", "family"])
T.to_csv(f"{PG}/SourceData_TE_timepoint_strain.csv", index=False)
topfam = T.family.value_counts().head(14).index.tolist()   # shared family rows across panels
mats = {tp: T[T.tp == tp].groupby(["family", "strain"]).size().unstack(fill_value=0).reindex(index=topfam, columns=CANON).fillna(0) for tp in TPS}
vmax = np.log10(max(m.values.max() for m in mats.values()) + 1)
plt.rcParams.update({"font.family": "Liberation Sans"})
fig, axes = plt.subplots(1, 3, figsize=(18, 6.2), dpi=300)
cmap = plt.get_cmap("magma").copy(); cmap.set_bad("white")
for ax, tp in zip(axes, TPS):
    Mt = mats[tp]; im = ax.imshow(np.ma.masked_where(Mt.values == 0, np.log10(Mt.values + 1)), aspect="auto", cmap=cmap, vmin=0, vmax=vmax)
    ax.set_xticks(range(len(CANON))); ax.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=90, fontsize=6.5)
    for t, X in zip(ax.get_xticklabels(), CANON):
        if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
    ax.set_yticks(range(len(topfam))); ax.set_yticklabels(topfam if tp == "E16.5" else [], fontsize=7)
    ax.set_xticks(np.arange(-.5, len(CANON), 1), minor=True); ax.set_yticks(np.arange(-.5, len(topfam), 1), minor=True)
    ax.grid(which="minor", color="#e6e6e6", linewidth=0.5); ax.tick_params(which="minor", length=0)
    for i in range(len(topfam)):
        for j in range(len(CANON)):
            v = int(Mt.values[i, j])
            if v > 0: ax.text(j, i, f"{v:,}" if v < 1000 else f"{v // 1000}k", ha="center", va="center", fontsize=3.8, color="white" if np.log10(v + 1) < vmax * 0.6 else "black")
    ax.set_title(tp, fontsize=11, fontweight="bold")
fig.colorbar(im, ax=axes, fraction=0.012, pad=0.01).set_label("log10(# strain-private piRNAs +1)", fontsize=8)
fig.suptitle("At each timepoint, which TE families drive strain-private piRNAs — TE family × 16 strains, faceted by developmental window\n"
             "wild strains (red) carry far more; watch the TE-family mix shift across E16.5 → P12.5 → P20.5", fontsize=11.5, fontweight="bold", y=1.02)
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_TE_timepoint_strain.{e}", bbox_inches="tight")
print("wrote Fig_TE_timepoint_strain + per-locus table")
