#!/usr/bin/env python3
"""
C57BL_6NJ_TE_annotation.py

Intersect C57BL_6NJ piRNA cluster BEDs with RepeatMasker annotation.
Report: per-developmental-class TE content (% coverage, dominant TE families).

Outputs in: analysis/claude_biomni_analysis/C57BL_6NJ_pangenome/
"""

import os, subprocess
import numpy as np
import pandas as pd

BASE   = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
RM_BED = f"{BASE}/resources/repeatMasker/C57BL_6NJ_repeatmasker.bed"
OUT    = f"{BASE}/analysis/claude_biomni_analysis/C57BL_6NJ_pangenome"

# ── Step 1: per-class TE coverage ──────────────────────────────────────────
print("=== TE annotation of C57BL_6NJ piRNA clusters ===")
print(f"RepeatMasker BED: {RM_BED}")
print()

CLASSES = ["shared_postnatal", "pachytene", "P12.5_only"]

class_bed = {c: os.path.join(OUT, f"C57BL_6NJ_{c}.bed") for c in CLASSES}

results = []

for dc in CLASSES:
    cbed = class_bed[dc]

    # Count total bp in class
    tot_cmd = f"awk '{{sum += $3-$2}} END {{print sum}}' {cbed}"
    r = subprocess.run(tot_cmd, shell=True, capture_output=True, text=True)
    total_bp = int(r.stdout.strip() or 0)

    n_loci = int(subprocess.run(f"wc -l < {cbed}", shell=True,
                                 capture_output=True, text=True).stdout.strip())

    # bedtools intersect -wo to get overlapping TE base pairs per TE type
    # Column 4 of RM BED: TE_name|TE_class (e.g. "L1MdTf_II|LINE/L1")
    int_cmd = (
        f"bedtools intersect -a {cbed} -b {RM_BED} -wo 2>/dev/null | "
        f"awk '{{te_class=$8; overlap=$NF; split(te_class,a,\"|\"); "
        f"print a[2], overlap}}'"
    )
    r = subprocess.run(int_cmd, shell=True, capture_output=True, text=True)

    # Accumulate bp per TE class
    te_bp = {}
    total_te_bp = 0
    for line in r.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        te_cls = parts[0]
        try:
            bp = int(parts[1])
        except ValueError:
            continue
        te_bp[te_cls] = te_bp.get(te_cls, 0) + bp
        total_te_bp += bp

    # TE coverage fraction
    te_frac = total_te_bp / total_bp if total_bp > 0 else 0

    # Top TE classes (by bp)
    top_te = sorted(te_bp.items(), key=lambda x: x[1], reverse=True)[:6]

    print(f"  {dc}:")
    print(f"    n_loci   = {n_loci:,}")
    print(f"    total_bp = {total_bp:,}")
    print(f"    TE_bp    = {total_te_bp:,}  ({100*te_frac:.1f}% of cluster bp)")
    print(f"    Top TE classes:")
    for te_cls, bp in top_te:
        pct = 100 * bp / total_bp if total_bp > 0 else 0
        print(f"      {te_cls:35s}: {bp:8,} bp  ({pct:.1f}%)")
    print()

    results.append({
        "dev_class":        dc,
        "n_loci":           n_loci,
        "total_bp":         total_bp,
        "TE_bp":            total_te_bp,
        "TE_fraction_pct":  round(100*te_frac, 2),
        **{f"bp_{te_cls}": bp for te_cls, bp in top_te},
    })

# ── Step 2: Save results ───────────────────────────────────────────────────
te_df = pd.DataFrame(results)
te_csv = os.path.join(OUT, "C57BL_6NJ_TE_annotation.csv")
te_df.to_csv(te_csv, index=False)
print(f"Saved: {te_csv}")

# ── Step 3: TE class breakdown per developmental class ─────────────────────
print("\n=== Detailed TE class breakdown per cluster class ===")
print("(Percentage of total cluster bp covered by each TE superfamily)\n")

TE_SUPERFAMILIES = {
    "LINE/L1":          "LINE-1 (L1MdTf, L1MdA, etc.)",
    "LTR/ERVK":         "IAP / ETn / MusD (ERVK LTR)",
    "LTR/ERV1":         "MLV, VL30 (ERV1 LTR)",
    "LTR/ERVL-MaLR":    "MaLR (ancient LTR)",
    "SINE/Alu":         "B1/B2 SINEs (Alu-like)",
    "SINE/B2":          "B2/B3 SINEs",
    "SINE/B4":          "ID/B4 SINEs",
    "DNA":              "DNA transposons",
    "Simple_repeat":    "Simple repeats / satellites",
    "Low_complexity":   "Low-complexity regions",
}

for dc in CLASSES:
    cbed = class_bed[dc]
    int_cmd = (
        f"bedtools intersect -a {cbed} -b {RM_BED} -wo 2>/dev/null | "
        f"awk '{{split($8,a,\"|\"); print a[2], $NF}}'"
    )
    r = subprocess.run(int_cmd, shell=True, capture_output=True, text=True)

    te_bp = {}
    total_te_bp = 0
    for line in r.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.strip().split()
        if len(parts) < 2:
            continue
        te_cls = parts[0]
        try:
            bp = int(parts[1])
        except ValueError:
            continue
        te_bp[te_cls] = te_bp.get(te_cls, 0) + bp
        total_te_bp += bp

    tot_bp_cmd = f"awk '{{sum += $3-$2}} END {{print sum}}' {cbed}"
    total_bp = int(subprocess.run(tot_bp_cmd, shell=True, capture_output=True,
                                   text=True).stdout.strip() or 0)

    print(f"  {dc} ({total_bp:,} bp total):")
    for fam, desc in TE_SUPERFAMILIES.items():
        bp = te_bp.get(fam, 0)
        pct = 100 * bp / total_bp if total_bp > 0 else 0
        if pct > 0.1:
            print(f"    {fam:25s}: {pct:5.1f}%  [{desc}]")
    # non-TE fraction
    non_te_pct = 100 - (100 * total_te_bp / total_bp if total_bp > 0 else 0)
    print(f"    {'Non-TE (genic)':25s}: {non_te_pct:5.1f}%")
    print()

    # Save per-class TE breakdown CSV
    breakdown = sorted(te_bp.items(), key=lambda x: x[1], reverse=True)
    pd.DataFrame({"te_class": [x[0] for x in breakdown],
                  "bp": [x[1] for x in breakdown],
                  "pct_of_total_cluster_bp": [100*x[1]/total_bp for x in breakdown]
                  }).to_csv(os.path.join(OUT, f"C57BL_6NJ_{dc}_TE_breakdown.csv"), index=False)

print("=== DONE ===")
