#!/usr/bin/env python3
"""
Re-run Steps 3-4 of the pangenome starter analysis with the fixed overlap logic
(-F 0.5 instead of -f 0.5 -r for timepoint presence classification), then
re-run TE annotation and all figures. Steps 1-2 (xlsx loading, 2-of-3 filter,
merge) are skipped because they are verified correct.
"""
import os, subprocess
import numpy as np
import pandas as pd

BASE    = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
OUT_DIR = f"{BASE}/analysis/claude_biomni_analysis/C57BL_6NJ_pangenome"
RM_BED  = f"{BASE}/resources/repeatMasker/C57BL_6NJ_repeatmasker.bed"

p12_path = os.path.join(OUT_DIR, "C57BL_6NJ_P12.5_merged.bed")
p20_path = os.path.join(OUT_DIR, "C57BL_6NJ_P20.5_merged.bed")
master_bed = os.path.join(OUT_DIR, "C57BL_6NJ_master.bed")

# ── Verify inputs exist ───────────────────────────────────────────────────────
for f in [p12_path, p20_path]:
    assert os.path.exists(f), f"MISSING: {f}"
    print(f"  Input verified: {os.path.basename(f)} ({sum(1 for _ in open(f))} loci)")

# ── Re-build master union (unchanged) ────────────────────────────────────────
print("\n=== STEP 3a: Re-build master union ===")
cmd = (f"cat {p12_path} {p20_path} | sort -k1,1 -k2,2n | "
       f"bedtools merge -d 1000 > {master_bed}")
subprocess.run(cmd, shell=True, check=True)
master = pd.read_csv(master_bed, sep="\t", header=None, names=["chr","start","end"])
print(f"  Master: {len(master)} loci")

# ── Classification with fixed -F 0.5 ────────────────────────────────────────
print("\n=== STEP 3b: Classification (-F 0.5) ===")

def timepoint_present(master_path, tp_merged_path, label):
    """Master loci that contain a timepoint cluster (-F 0.5: ≥50% of tp cluster in master)."""
    cmd = (f"bedtools intersect -a {master_path} -b {tp_merged_path} -F 0.5 -wa | sort -u")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
    hits = set()
    for line in r.stdout.strip().split("\n"):
        if not line: continue
        c = line.split("\t")
        hits.add((c[0], c[1], c[2]))
    print(f"  Present at {label}: {len(hits)}")
    return hits

flag_p12 = os.path.join(OUT_DIR, "C57BL_6NJ_master_present_P12.5.bed")
flag_p20 = os.path.join(OUT_DIR, "C57BL_6NJ_master_present_P20.5.bed")
hits_p12 = timepoint_present(master_bed, p12_path, "P12.5")
hits_p20 = timepoint_present(master_bed, p20_path, "P20.5")

def classify(row):
    key = (row["chr"], str(row["start"]), str(row["end"]))
    p12 = key in hits_p12
    p20 = key in hits_p20
    if p12 and p20:   return "shared_postnatal"
    elif p20:         return "pachytene"
    elif p12:         return "P12.5_only"
    else:             return "unclassified"

master["dev_class"] = master.apply(classify, axis=1)

print("\n  Classification summary:")
counts = master["dev_class"].value_counts()
for dc, n in counts.items():
    print(f"    {dc:25s}: {n:6d}  ({100*n/len(master):.1f}%)")

# Save
classified_bed = os.path.join(OUT_DIR, "C57BL_6NJ_classified.bed")
master.to_csv(classified_bed, sep="\t", index=False, header=False)
for dc in ("shared_postnatal", "pachytene", "P12.5_only", "unclassified"):
    sub = master[master["dev_class"] == dc]
    sub[["chr","start","end","dev_class"]].to_csv(
        os.path.join(OUT_DIR, f"C57BL_6NJ_{dc}.bed"), sep="\t", index=False, header=False)
print(f"  Saved classified BED and per-class BEDs")

# ── STEP 4: Per-class statistics ─────────────────────────────────────────────
print("\n=== STEP 4: Per-class statistics ===")

p12_merged = pd.read_csv(p12_path, sep="\t", header=None,
                         names=["chr","start","end","mean_FPM","type"])
p20_merged = pd.read_csv(p20_path, sep="\t", header=None,
                         names=["chr","start","end","mean_FPM","type"])

def get_fpm_for_class(class_bed_path, source_df):
    """FPM values for source clusters with ≥50% of their length in a class locus."""
    tmp_src = os.path.join(OUT_DIR, "_tmp_src.bed")
    source_df[["chr","start","end","mean_FPM"]].to_csv(
        tmp_src, sep="\t", index=False, header=False)
    cmd = (f"bedtools intersect -a {tmp_src} -b {class_bed_path} -f 0.5 -wa "
           f"| awk '{{print $4}}'")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    os.remove(tmp_src)
    return [float(x) for x in r.stdout.strip().split("\n") if x]

stats_rows = []
for dc in ("shared_postnatal", "pachytene", "P12.5_only"):
    sub  = master[master["dev_class"] == dc]
    n    = len(sub)
    sizes = sub["end"] - sub["start"]
    class_bed = os.path.join(OUT_DIR, f"C57BL_6NJ_{dc}.bed")

    fpm12 = get_fpm_for_class(class_bed, p12_merged)
    fpm20 = get_fpm_for_class(class_bed, p20_merged)

    row = {
        "dev_class":         dc,
        "n_loci":            n,
        "median_size_bp":    int(sizes.median()),
        "mean_size_bp":      int(sizes.mean()),
        "max_size_bp":       int(sizes.max()),
        "n_P12.5_expressed": len(fpm12),
        "median_FPM_P12.5":  round(float(np.median(fpm12)), 2) if fpm12 else None,
        "n_P20.5_expressed": len(fpm20),
        "median_FPM_P20.5":  round(float(np.median(fpm20)), 2) if fpm20 else None,
    }
    stats_rows.append(row)
    print(f"  {dc}: n={n}, med_size={row['median_size_bp']:,} bp, "
          f"P12.5 FPM_med={row['median_FPM_P12.5']}, P20.5 FPM_med={row['median_FPM_P20.5']}")

stats_df = pd.DataFrame(stats_rows)
stats_df.to_csv(os.path.join(OUT_DIR, "C57BL_6NJ_class_stats.csv"), index=False)
print("  Saved class_stats.csv")

# ── STEP 7: TE annotation (re-run on new per-class BEDs) ────────────────────
print("\n=== STEP 7: TE annotation ===")

te_rows = []
for dc in ("shared_postnatal", "pachytene", "P12.5_only"):
    class_bed = os.path.join(OUT_DIR, f"C57BL_6NJ_{dc}.bed")
    n_loci    = int(stats_df[stats_df.dev_class==dc]["n_loci"].iloc[0])
    total_bp  = int((master[master.dev_class==dc]["end"] -
                     master[master.dev_class==dc]["start"]).sum())

    cmd = (f"bedtools intersect -a {class_bed} -b {RM_BED} -wo 2>/dev/null | "
           f"awk '{{split($8,a,\"|\"); split(a[2],b,\"/\"); "
           f"print b[1]\"/\"b[2], $NF}}' | "
           f"awk '{{s[$1]+=$2}} END{{for(k in s) print k, s[k]}}'")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    te_bp_total = 0
    te_by_class = {}
    for line in r.stdout.strip().split("\n"):
        if not line: continue
        parts = line.split()
        if len(parts) != 2: continue
        tc, bp = parts[0], int(parts[1])
        te_by_class[tc] = bp
        if not tc.startswith("Simple") and not tc.startswith("Low") and tc != "Non-TE":
            te_bp_total += bp

    row = {
        "dev_class":       dc,
        "n_loci":          n_loci,
        "total_bp":        total_bp,
        "TE_bp":           te_bp_total,
        "TE_fraction_pct": round(100.0 * te_bp_total / total_bp, 2) if total_bp else 0,
    }
    # Top TE classes
    for tc, bp in sorted(te_by_class.items(), key=lambda x: -x[1])[:6]:
        row[f"bp_{tc}"] = bp
    te_rows.append(row)
    print(f"  {dc}: {n_loci} loci, {total_bp:,} bp, {row['TE_fraction_pct']}% TE")

    # Save per-class breakdown
    breakdown = sorted(te_by_class.items(), key=lambda x: -x[1])
    pd.DataFrame(breakdown, columns=["te_class","bp"]).assign(
        pct_of_total_cluster_bp=lambda d: 100.0 * d["bp"] / total_bp
    ).to_csv(os.path.join(OUT_DIR, f"C57BL_6NJ_{dc}_TE_breakdown.csv"), index=False)

te_df = pd.DataFrame(te_rows)
te_df.to_csv(os.path.join(OUT_DIR, "C57BL_6NJ_TE_annotation.csv"), index=False)
print("  Saved C57BL_6NJ_TE_annotation.csv")

print("\n=== ALL STEPS COMPLETE ===")
print(f"Output: {OUT_DIR}")
