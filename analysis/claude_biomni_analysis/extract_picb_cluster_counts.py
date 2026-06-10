#!/usr/bin/env python3
"""
Extract PICB cluster counts from every per-replicate and combined xlsx, for:
  (1) combined clusters across all 16 strains
  (2) the effect of single-replicate vs combined-replicate PICB.
Counts the final `clusters` sheet (+ seeds/cores), cluster total width, and a
chromosome-completeness QC (clusters should span chr1-19 + X). Writes a tidy CSV.
Per-rep:  results/picb_result/{strain}/{strain}-{tp}.{rep}.xlsx
Combined: results/picb_result_combined/{strain}/{strain}-{tp}.combined.xlsx
"""
import glob, os, sys, warnings
import pandas as pd
warnings.simplefilter("ignore")

ROOT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
ISV = {"129S1_SvImJ","AKR_J","A_J","BALB_cJ","C3H_HeJ","CAST_EiJ","CBA_J","DBA_2J",
       "FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ","C57BL_6NJ"}
TPS = ["16.5dpc","12.5dpp","20.5dpp"]
EXPECT = {str(i) for i in range(1,20)} | {"X"}

def parse(base):
    for tp in TPS:
        key = f"-{tp}."
        if key in base:
            strain = base.split(key)[0]
            rep = base.split(key)[1].rsplit(".xlsx",1)[0]   # '1'/'2'/'3' or 'combined'
            return strain, tp, rep
    return None, None, None

files = sorted(glob.glob(f"{ROOT}/results/picb_result/*/*.xlsx") +
               glob.glob(f"{ROOT}/results/picb_result_combined/*/*.combined.xlsx"))
rows = []
for i, f in enumerate(files):
    strain, tp, rep = parse(os.path.basename(f))
    if strain not in ISV or tp is None:
        continue
    try:
        xl = pd.ExcelFile(f)
        n = {sh: 0 for sh in ("seeds","cores","clusters")}
        width = chrom = 0
        for sh in xl.sheet_names:
            d = xl.parse(sh)
            if sh in n: n[sh] = len(d)
            if sh == "clusters" and len(d):
                width = int(d["width"].sum())
                chrom = d["seqnames"].astype(str).str.replace("chr","").nunique()
        rows.append(dict(strain=strain, timepoint=tp, replicate=rep,
                         n_seeds=n["seeds"], n_cores=n["cores"], n_clusters=n["clusters"],
                         cluster_total_width_bp=width, cluster_n_chrom=chrom,
                         complete=(chrom>=20)))
    except Exception as e:
        rows.append(dict(strain=strain, timepoint=tp, replicate=rep, n_seeds=-1,
                         n_cores=-1, n_clusters=-1, cluster_total_width_bp=-1,
                         cluster_n_chrom=-1, complete=False, error=str(e)[:60]))
    if (i+1) % 20 == 0:
        print(f"  ...{i+1}/{len(files)} read", flush=True)

df = pd.DataFrame(rows).sort_values(["strain","timepoint","replicate"])
out = f"{ROOT}/analysis/claude_biomni_analysis/source_data/SourceData_PICB_cluster_counts.csv"
df.to_csv(out, index=False)
print(f"\nWROTE {out}  ({len(df)} rows)")
print("  per-rep rows:", (df.replicate!='combined').sum(), "| combined rows:", (df.replicate=='combined').sum())
print("  any incomplete (n_chrom<20):", list(df[~df.complete][['strain','timepoint','replicate','cluster_n_chrom']].itertuples(index=False)) or "none")
print("\n  combined clusters per strain x tp:")
piv = df[df.replicate=='combined'].pivot_table(index='strain', columns='timepoint', values='n_clusters', aggfunc='first')
print(piv[TPS].to_string())
