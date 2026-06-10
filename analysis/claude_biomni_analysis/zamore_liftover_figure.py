#!/usr/bin/env python3
"""
Liftover Zamore piRNA gene annotation (mm10/GRCm38) to C57BL_6NJ REL-2205 assembly,
then intersect with PICB P12.5/P20.5 clusters and produce a comparison figure.

Liftover chain:
  Step 1: mm10 → GRCm39   (mm10ToMm39.over.chain.gz)
  Step 2: GRCm39 → C57BL_6NJ  (GRCm39_C57BL_6NJ_chromosomes_MT_unplaced.chain)
  Note: GRCm39 chain uses numeric chr names (no "chr" prefix)
  Note: C57BL_6NJ output uses "C57BL_6NJ#1#chrN" → strip to "chrN"
"""
import os, subprocess, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

BASE    = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
OUT_DIR = f"{BASE}/analysis/claude_biomni_analysis/C57BL_6NJ_pangenome"

ZAMORE_XLSX = f"{BASE}/resources/zamore_piRNAs/piRNA_gene_annotation-modified-in-orange-with-Wasik-et-al-checks.xlsx"
CHAIN_MM10_MM39  = "/mnt/beegfs/scratch/miska/nm667/inProgress/mice_PiRNA/results/liftOverChainFiles/mm10ToMm39.over.chain.gz"
CHAIN_MM39_6NJ   = f"{BASE}/resources/REL-2205-Assembly/GRCm39_chains/GRCm39_C57BL_6NJ_chromosomes_MT_unplaced.chain"
LIFTOVER_BIN     = f"{BASE}/workflow/scripts/liftOver"

P12_BED  = os.path.join(OUT_DIR, "C57BL_6NJ_P12.5_merged.bed")
P20_BED  = os.path.join(OUT_DIR, "C57BL_6NJ_P20.5_merged.bed")

TMP = OUT_DIR  # write intermediate files here

# ── Wong 2011 colourblind palette ─────────────────────────────────────────────
C_PREP  = "#0072B2"   # blue — Prepachytene
C_PACH  = "#E69F00"   # amber — Pachytene
C_HYB   = "#009E73"   # green — Hybrid
C_P12   = "#56B4E9"   # sky blue — P12.5 PICB
C_P20   = "#CC79A7"   # pink — P20.5 PICB

plt.rcParams.update({
    "font.family": "Liberation Sans",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Parse Zamore xlsx → BED4 in mm10
# ═══════════════════════════════════════════════════════════════════════════════
print("=== STEP 1: Parse Zamore xlsx ===")

# Row 0 is a description/header row (col name in pandas); data starts at row 1
# Use header=None and skiprows=1 to get all data cleanly
df = pd.read_excel(ZAMORE_XLSX, sheet_name="Mus musculus", header=None, skiprows=1)
# BED12 cols: 0=chr 1=start 2=end 3=name 4=score 5=strand 6=thickStart 7=thickEnd
#             8=itemRgb 9=blockCount 10=blockSizes 11=blockStarts 12=stage
df.columns = ["chr","start","end","name","score","strand",
              "thickStart","thickEnd","itemRgb","blockCount",
              "blockSizes","blockStarts","stage"]

# Normalise stage labels
df["stage"] = df["stage"].astype(str).str.strip()
stage_map = {"Prepachytene": "Prepachytene", "Pachytene": "Pachytene",
             "Hybrid": "Hybrid"}
df["stage"] = df["stage"].map(stage_map).fillna("Other")

print(f"  Total Zamore genes: {len(df)}")
print(f"  Stage breakdown: {df['stage'].value_counts().to_dict()}")

# Write mm10 BED4
zamore_mm10 = os.path.join(TMP, "zamore_mm10.bed")
df[["chr","start","end","name","score","strand"]].to_csv(
    zamore_mm10, sep="\t", index=False, header=False)
# Also save stage info for later lookup
stage_lookup = df.set_index("name")["stage"].to_dict()


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Liftover mm10 → GRCm39
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== STEP 2: Liftover mm10 → GRCm39 ===")

zamore_mm39_raw   = os.path.join(TMP, "zamore_mm39_raw.bed")
zamore_mm39_unmap = os.path.join(TMP, "zamore_mm39_unmap1.bed")

cmd = (f"{LIFTOVER_BIN} {zamore_mm10} {CHAIN_MM10_MM39} "
       f"{zamore_mm39_raw} {zamore_mm39_unmap}")
r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
if r.returncode not in (0, 1):
    raise RuntimeError(f"liftOver step1 failed:\n{r.stderr}")

mapped1   = sum(1 for _ in open(zamore_mm39_raw))
unmapped1 = sum(1 for l in open(zamore_mm39_unmap) if not l.startswith("#"))
print(f"  Mapped to GRCm39:   {mapped1}")
print(f"  Unmapped in step 1: {unmapped1}")

# Strip "chr" prefix so names match numeric names used in GRCm39→6NJ chain
zamore_mm39 = os.path.join(TMP, "zamore_mm39.bed")
with open(zamore_mm39_raw) as fin, open(zamore_mm39, "w") as fout:
    for line in fin:
        cols = line.split("\t")
        cols[0] = re.sub(r"^chr", "", cols[0])
        fout.write("\t".join(cols))


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Liftover GRCm39 → C57BL_6NJ
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== STEP 3: Liftover GRCm39 → C57BL_6NJ ===")

zamore_6nj_raw   = os.path.join(TMP, "zamore_6nj_raw.bed")
zamore_6nj_unmap = os.path.join(TMP, "zamore_6nj_unmap2.bed")

cmd = (f"{LIFTOVER_BIN} {zamore_mm39} {CHAIN_MM39_6NJ} "
       f"{zamore_6nj_raw} {zamore_6nj_unmap}")
r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
if r.returncode not in (0, 1):
    raise RuntimeError(f"liftOver step2 failed:\n{r.stderr}")

mapped2   = sum(1 for _ in open(zamore_6nj_raw))
unmapped2 = sum(1 for l in open(zamore_6nj_unmap) if not l.startswith("#"))
print(f"  Mapped to C57BL_6NJ: {mapped2}")
print(f"  Unmapped in step 2:  {unmapped2}")

# Convert "C57BL_6NJ#1#chr1" → "chr1"
zamore_6nj = os.path.join(TMP, "zamore_6nj.bed")
with open(zamore_6nj_raw) as fin, open(zamore_6nj, "w") as fout:
    for line in fin:
        cols = line.split("\t")
        m = re.search(r"chr\S+", cols[0])
        if m:
            cols[0] = m.group(0)
        fout.write("\t".join(cols))

# Load final lifted-over annotation (all isoforms)
zamore_iso = pd.read_csv(zamore_6nj, sep="\t", header=None,
                         names=["chr","start","end","name","score","strand"],
                         usecols=[0,1,2,3,4,5])
zamore_iso["stage"] = zamore_iso["name"].map(stage_lookup).fillna("Other")
zamore_iso = zamore_iso[zamore_iso["chr"].str.startswith("chr")]

print(f"\n  Lifted isoforms in C57BL_6NJ coords: {len(zamore_iso)}")
print(f"  Stage breakdown (isoforms): {zamore_iso['stage'].value_counts().to_dict()}")

# Collapse isoforms → one row per piRNA LOCUS
# Base gene name = strip trailing .N suffix (e.g. pi-Zbtb37.1 → pi-Zbtb37)
import re as _re
zamore_iso["base_gene"] = zamore_iso["name"].str.replace(r'\.\d+$', '', regex=True)

# Stage is consistent across isoforms; take first occurrence
stage_by_locus = zamore_iso.groupby("base_gene")["stage"].first()

# Union span per locus: same chr (assert), min start, max end
zamore_final = (
    zamore_iso.groupby("base_gene")
    .agg(chr=("chr","first"), start=("start","min"),
         end=("end","max"), strand=("strand","first"))
    .reset_index()
    .rename(columns={"base_gene": "name"})
)
zamore_final["stage"] = zamore_final["name"].map(stage_by_locus)

print(f"\n  Unique piRNA loci after collapsing isoforms: {len(zamore_final)}")
print(f"  Stage breakdown (loci): {zamore_final['stage'].value_counts().to_dict()}")

# Save locus-level BED (used for intersection)
zamore_loci_bed = os.path.join(TMP, "zamore_C57BL_6NJ_loci.bed")
zamore_final[["chr","start","end","name","strand"]].to_csv(
    zamore_loci_bed, sep="\t", index=False, header=False)
print("  Saved zamore_C57BL_6NJ_loci.bed")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Intersect with PICB P12.5 and P20.5 clusters
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== STEP 4: Intersect with PICB clusters ===")

zamore_sorted = os.path.join(TMP, "_zamore_sorted.bed")
subprocess.run(f"sort -k1,1 -k2,2n {zamore_loci_bed} > {zamore_sorted}",
               shell=True, check=True)

def intersect_zamore_picb(zamore_bed, picb_bed, label):
    """Zamore genes overlapping (any bp) PICB clusters at a timepoint."""
    cmd = (f"bedtools intersect -a {zamore_bed} -b {picb_bed} -wa -u | "
           f"awk '{{print $4}}'")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
    hits = set(x for x in r.stdout.strip().split("\n") if x)
    print(f"  Zamore genes overlapping {label}: {len(hits)}")
    return hits

hits_p12 = intersect_zamore_picb(zamore_sorted, P12_BED, "P12.5 PICB")
hits_p20 = intersect_zamore_picb(zamore_sorted, P20_BED, "P20.5 PICB")

# Classify each lifted Zamore gene by PICB overlap
def picb_class(row):
    n  = row["name"]
    p12 = n in hits_p12
    p20 = n in hits_p20
    if p12 and p20:   return "P12.5 & P20.5"
    elif p12:         return "P12.5 only"
    elif p20:         return "P20.5 only"
    else:             return "No overlap"

zamore_final["picb_class"] = zamore_final.apply(picb_class, axis=1)

print("\n  Cross-tabulation (Zamore stage × PICB overlap):")
ct = pd.crosstab(zamore_final["stage"], zamore_final["picb_class"])
print(ct.to_string())

ct.to_csv(os.path.join(TMP, "zamore_PICB_crosstab.csv"))
print("  Saved zamore_PICB_crosstab.csv")


# ═══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Figure
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== STEP 5: Figure ===")

PICB_COLS  = {"P12.5 & P20.5": "#8856a7",  # purple
              "P12.5 only":    C_P12,       # sky blue
              "P20.5 only":    C_P20,       # pink
              "No overlap":    "#cccccc"}   # grey

STAGE_COLS = {"Prepachytene": C_PREP,
              "Pachytene":    C_PACH,
              "Hybrid":       C_HYB,
              "Other":        "#aaaaaa"}

stages_order = ["Prepachytene", "Pachytene", "Hybrid"]
picb_order   = ["P12.5 only", "P12.5 & P20.5", "P20.5 only", "No overlap"]

fig = plt.figure(figsize=(12, 10))
gs  = GridSpec(2, 2, figure=fig, hspace=0.65, wspace=0.42)

# ── Panel A: stacked bar — per Zamore stage, % PICB overlap class ─────────────
ax_a = fig.add_subplot(gs[0, 0])

totals = zamore_final["stage"].value_counts()
pct_data = {}
for pc in picb_order:
    vals = []
    for st in stages_order:
        sub = zamore_final[zamore_final["stage"] == st]
        n_overlap = (sub["picb_class"] == pc).sum()
        pct = 100.0 * n_overlap / len(sub) if len(sub) else 0
        vals.append(pct)
    pct_data[pc] = vals

x = np.arange(len(stages_order))
bottoms = np.zeros(len(stages_order))
for pc in picb_order:
    bars = ax_a.bar(x, pct_data[pc], bottom=bottoms,
                    color=PICB_COLS[pc], label=pc, width=0.55, edgecolor="white", lw=0.3)
    bottoms += np.array(pct_data[pc])

ax_a.set_xticks(x)
ax_a.set_xticklabels(stages_order, fontsize=8)
ax_a.set_ylabel("% of Zamore genes", fontsize=8)
ax_a.set_title("A  Zamore piRNA genes:\nPICB cluster overlap by stage", fontsize=9, loc="left", pad=4)
ax_a.set_ylim(0, 108)
# n= labels just above bars, well within ylim
n_labels = [f"n={totals.get(s,0)}" for s in stages_order]
for i, nl in enumerate(n_labels):
    ax_a.text(i, 102, nl, ha="center", va="bottom", fontsize=7)
# Legend placed below the panel, outside axes
ax_a.legend(title="PICB overlap", fontsize=7, title_fontsize=7,
            loc="upper center", bbox_to_anchor=(0.5, -0.18),
            ncol=2, frameon=False,
            handles=[mpatches.Patch(color=PICB_COLS[pc], label=pc) for pc in picb_order])

# ── Panel B: stacked bar — per PICB class, composition of Zamore stages ──────
ax_b = fig.add_subplot(gs[0, 1])

stage_in_class = {}
for pc in picb_order:
    sub = zamore_final[zamore_final["picb_class"] == pc]
    stage_in_class[pc] = {st: (sub["stage"] == st).sum() for st in stages_order}

x2 = np.arange(len(picb_order))
bottoms2 = np.zeros(len(picb_order))
for st in stages_order:
    vals2 = [stage_in_class[pc][st] for pc in picb_order]
    ax_b.bar(x2, vals2, bottom=bottoms2,
             color=STAGE_COLS[st], label=st, width=0.55,
             edgecolor="white", lw=0.3)
    bottoms2 += np.array(vals2, dtype=float)

ax_b.set_xticks(x2)
ax_b.set_xticklabels(picb_order, fontsize=7, rotation=35, ha="right")
ax_b.set_ylabel("Number of Zamore genes", fontsize=8)
ax_b.set_title("B  PICB overlap class:\ncomposition of Zamore stages", fontsize=9, loc="left", pad=4)
y_max_b = max(bottoms2) * 1.15
ax_b.set_ylim(0, y_max_b)
for i, pc in enumerate(picb_order):
    n = int(bottoms2[i])
    ax_b.text(i, n + y_max_b * 0.01, str(n), ha="center", va="bottom", fontsize=7)
# Legend placed below the panel, outside axes
ax_b.legend(title="Zamore stage", fontsize=7, title_fontsize=7,
            loc="upper center", bbox_to_anchor=(0.5, -0.25),
            ncol=3, frameon=False,
            handles=[mpatches.Patch(color=STAGE_COLS[st], label=st) for st in stages_order])

# ── Panel C: genomic span (bp) of lifted Zamore genes by stage ───────────────
ax_c = fig.add_subplot(gs[1, 0])

span_data = []
for st in stages_order:
    sub = zamore_final[zamore_final["stage"] == st]
    spans = (sub["end"] - sub["start"]).values
    span_data.append(spans)

bp = ax_c.boxplot(span_data, positions=range(len(stages_order)),
                  patch_artist=True, notch=False,
                  widths=0.45, showfliers=True,
                  flierprops=dict(marker=".", markersize=2, alpha=0.5),
                  medianprops=dict(color="black", lw=1.5))
for patch, st in zip(bp["boxes"], stages_order):
    patch.set_facecolor(STAGE_COLS[st])
    patch.set_alpha(0.8)
ax_c.set_xticks(range(len(stages_order)))
ax_c.set_xticklabels(stages_order, fontsize=8)
ax_c.set_ylabel("Cluster span (bp)", fontsize=8)
ax_c.set_yscale("log")
ax_c.set_title("C  Genomic span of Zamore genes\nin C57BL_6NJ coords", fontsize=9, loc="left", pad=4)

# ── Panel D: chromosomal distribution of P20.5-only Zamore genes ─────────────
ax_d = fig.add_subplot(gs[1, 1])

# Show chr distribution for each stage across all autosomes + X
chr_order = [f"chr{i}" for i in list(range(1,20))] + ["chrX"]

def chr_counts(df, stage):
    sub = df[df["stage"] == stage]
    counts = sub["chr"].value_counts().reindex(chr_order, fill_value=0)
    return counts.values

x3  = np.arange(len(chr_order))
bot3 = np.zeros(len(chr_order))
for st in stages_order:
    vals3 = chr_counts(zamore_final, st)
    ax_d.bar(x3, vals3, bottom=bot3, color=STAGE_COLS[st],
             label=st, width=0.8, edgecolor="none")
    bot3 += vals3

ax_d.set_xticks(x3)
ax_d.set_xticklabels([c.replace("chr","") for c in chr_order],
                     fontsize=6, rotation=0)
ax_d.set_ylabel("Number of Zamore genes", fontsize=8)
ax_d.set_title("D  Chromosomal distribution\nof lifted Zamore genes", fontsize=9, loc="left", pad=4)
ax_d.legend(title="Zamore stage", fontsize=7, title_fontsize=7,
            loc="upper right", frameon=False,
            handles=[mpatches.Patch(color=STAGE_COLS[st], label=st) for st in stages_order])

# ── Save ──────────────────────────────────────────────────────────────────────
fig.suptitle(
    "Zamore piRNA gene annotation (mm10→C57BL_6NJ REL-2205)\nvs. PICB P12.5/P20.5 clusters",
    fontsize=10, y=1.01)

for ext in ("pdf", "svg", "png"):
    path = os.path.join(OUT_DIR, f"Fig6_Zamore_liftover.{ext}")
    fig.savefig(path, dpi=300, bbox_inches="tight",
                format=ext if ext != "png" else None)
plt.close(fig)
print("  Saved Fig6_Zamore_liftover.{pdf,svg,png}")
print(f"\nOutput: {OUT_DIR}")
