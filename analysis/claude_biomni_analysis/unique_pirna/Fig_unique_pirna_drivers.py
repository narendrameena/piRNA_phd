#!/usr/bin/env python3
"""What DRIVES each strain's UNIQUE (strain-private-locus) piRNAs — TE vs ncRNA(lncRNA) vs protein-coding gene vs
intergenic? Per strain, each strain-private piRNA locus (cand_self16 BAM, own genome) is intersected with the
per-strain RepeatMasker BED (TE) and the per-strain Ensembl gene models (protein_coding gene / lncRNA gene). Each
locus gets ONE mutually-exclusive driver by priority TE > lncRNA > protein-coding > intergenic (TE wins because a
TE inside an intron is still TE-driven). Keeps timepoint (locus id = strain|tp|seq). Two panels (canonical strain
order, wild red): (A) driver composition per strain; (B) driver composition per developmental window (does the
TE-driven prepachytene programme give way to lncRNA/other later?). Biology framing queued for BioMNI.
Linked: [[project_te_driven_finding]] [[project_ncrna_driven_finding]]."""
import warnings; warnings.filterwarnings("ignore")
import sys, subprocess, tempfile, os
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
from strain_order import STRAIN_ORDER, WILD
import pandas as pd, numpy as np, pysam
from collections import defaultdict
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"; U = f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; PG = f"{U}/pangenome_te"
BT = "/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/bedtools"; GB = f"{PG}/gene_beds"
CANON = [s for s in STRAIN_ORDER if s != "C57BL_6"]; TPMAP = {"16.5dpc": "E16.5", "12.5dpp": "P12.5", "20.5dpp": "P20.5"}; TPS = ["E16.5", "P12.5", "P20.5"]
DRIVERS = [("TE", "#1B7837"), ("lncRNA", "#762A83"), ("protein-coding", "#4575B4"), ("intergenic", "#BBBBBB")]
d = pd.read_csv(f"{U}/unique16/final_classified.csv.gz")
def ids_overlapping(locibed, bed):
    if not os.path.exists(bed) or os.path.getsize(bed) == 0: return set()
    out = subprocess.run(f"sort -k1,1 -k2,2n {locibed} | {BT} intersect -a - -b {bed} -u", shell=True, capture_output=True, text=True).stdout
    return {ln.split('\t')[3] for ln in out.splitlines() if ln}
rows = []
for X in CANON:
    g = d[(d.strain == X) & (d.klass == "unique: strain-private locus")].copy(); g["id"] = X + "|" + g.timepoint + "|" + g.sequence
    sp = set(g.id); rm = f"{ROOT}/resources/repeatMasker/{X}_repeatmasker.bed"
    if not sp: continue
    bam = pysam.AlignmentFile(f"{U}/cand_self16/{X}.cand_self16.bam", "rb")
    tb = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False, dir=PG); seen = {}
    for a in bam.fetch(until_eof=True):
        if a.is_unmapped or a.query_name not in sp or a.query_name in seen: continue
        seen[a.query_name] = 1; tb.write(f"{a.reference_name.split('#')[-1]}\t{a.reference_start}\t{a.reference_end}\t{a.query_name}\n")
    tb.close(); bam.close()
    TEset = ids_overlapping(tb.name, rm); LNCset = ids_overlapping(tb.name, f"{GB}/{X}.lnc.bed"); PCset = ids_overlapping(tb.name, f"{GB}/{X}.pc.bed")
    os.unlink(tb.name)
    for locid in seen:
        tp = TPMAP.get(locid.split("|")[1], locid.split("|")[1])
        drv = "TE" if locid in TEset else "lncRNA" if locid in LNCset else "protein-coding" if locid in PCset else "intergenic"
        rows.append((X, tp, drv))
    print(f"{X}: n={len(seen):,} TE={len(TEset):,} lnc={len(LNCset):,} pc={len(PCset):,}")
T = pd.DataFrame(rows, columns=["strain", "tp", "driver"])
T.to_csv(f"{PG}/SourceData_unique_pirna_drivers.csv", index=False)
A = T.groupby(["strain", "driver"]).size().unstack(fill_value=0).reindex(index=CANON, columns=[k for k, _ in DRIVERS]).fillna(0)
Afrac = A.div(A.sum(1).replace(0, np.nan), axis=0) * 100
B = T.groupby(["tp", "driver"]).size().unstack(fill_value=0).reindex(index=TPS, columns=[k for k, _ in DRIVERS]).fillna(0)
Bfrac = B.div(B.sum(1), axis=0) * 100
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 8, "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(13, 5.6), dpi=300); gs = fig.add_gridspec(1, 3, width_ratios=[2.6, 1, 0.02], wspace=0.28)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
x = np.arange(len(CANON)); bottom = np.zeros(len(CANON))
for k, col in DRIVERS:
    axA.bar(x, Afrac[k].values, bottom=bottom, width=0.82, color=col, edgecolor="white", linewidth=0.3, label=k)
    bottom = bottom + np.nan_to_num(Afrac[k].values)
axA.set_xticks(x); axA.set_xticklabels([s.replace("_", "/") for s in CANON], rotation=90, fontsize=7)
for t, X in zip(axA.get_xticklabels(), CANON):
    if X in WILD: t.set_color("#C0392B"); t.set_fontweight("bold")
axA.set_ylim(0, 100); axA.set_ylabel("driver composition (%)", fontsize=9)
axA.set_title("A  Driver of strain-private piRNAs per strain (TE > lncRNA > protein-coding > intergenic)", fontsize=9.6, fontweight="bold", loc="left")
axA.legend(handles=[Patch(facecolor=c, label=k) for k, c in DRIVERS], ncol=4, fontsize=8, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.42))
xb = np.arange(len(TPS)); bottom = np.zeros(len(TPS))
for k, col in DRIVERS:
    axB.bar(xb, Bfrac[k].values, bottom=bottom, width=0.7, color=col, edgecolor="white", linewidth=0.3)
    bottom = bottom + Bfrac[k].values
axB.set_xticks(xb); axB.set_xticklabels(TPS, fontsize=8.5); axB.set_ylim(0, 100); axB.set_ylabel("driver composition (%)", fontsize=9)
axB.set_title("B  by timepoint\n(pooled)", fontsize=9.6, fontweight="bold", loc="left")
fig.suptitle("What drives UNIQUE (strain-private) piRNAs — TE vs ncRNA(lncRNA) vs protein-coding vs intergenic, 16 strains × 3 windows",
             fontsize=11, fontweight="bold", y=1.0)
fig.tight_layout(rect=[0, 0.02, 1, 0.97])
for e in ("pdf", "svg", "png"): fig.savefig(f"{PG}/Fig_unique_pirna_drivers.{e}", bbox_inches="tight")
print("wrote Fig_unique_pirna_drivers + source data"); print("by-timepoint %:\n", Bfrac.round(1).to_string())
