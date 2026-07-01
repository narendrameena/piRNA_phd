#!/usr/bin/env python3
"""
Verify combined-PICB xlsx are COMPLETE (ground truth, not file existence).
Complete = clusters span all autosomes chr1-19 + chrX. Reports incomplete/corrupt ones.
Usage: python verify_picb_completeness.py [glob]   (default: all combined xlsx)
"""
import glob, os, sys, pandas as pd, warnings
warnings.simplefilter("ignore")
EXPECT = {str(i) for i in range(1, 20)} | {"X"}
pat = sys.argv[1] if len(sys.argv) > 1 else \
    "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/results/picb_result_combined/*/*.combined.xlsx"
files = sorted(glob.glob(pat))
rows = []
for f in files:
    name = os.path.basename(f).replace(".combined.xlsx", "")
    try:
        xl = pd.ExcelFile(f); best = None
        for sh in xl.sheet_names:
            d = xl.parse(sh)
            c = next((c for c in d.columns if str(c).lower() in
                      ("seqnames", "seqname", "chr", "chrom", "chromosome")), None)
            if c is not None and len(d):
                ch = {str(x).replace("chr", "") for x in d[c].dropna().unique()}
                if best is None or len(ch) > len(best): best = ch
        miss = sorted(EXPECT - (best or set()), key=lambda s: (s == "X", int(s) if s.isdigit() else 99))
        rows.append((name, len(best or []), "COMPLETE" if not miss else "MISSING:" + ",".join(miss)))
    except Exception as e:
        rows.append((name, -1, f"ERR:{e}"[:50]))
df = pd.DataFrame(rows, columns=["sample", "n_chrom", "status"]).sort_values(["status", "sample"])
print(df.to_string(index=False))
nC = (df.status == "COMPLETE").sum()
print(f"\nCOMPLETE: {nC} / {len(df)}")
bad = list(df[df.status != "COMPLETE"]["sample"])
if bad: print("INCOMPLETE ->", bad)
