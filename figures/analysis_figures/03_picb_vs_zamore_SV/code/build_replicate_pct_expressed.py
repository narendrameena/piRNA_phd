#!/usr/bin/env python3
"""
Theme 03 — per-REPLICATE % of Zamore loci expressed (for Panel A error bars).

The combined expression matrix (all_strains_expression_matrix.csv) is built from the
COMBINED-replicate PICB run (one value per strain x timepoint). To draw replicate
error bars we recompute % expressed SEPARATELY for each replicate's PICB clusters,
reusing the SAME projected Zamore loci (_zamore_{strain}.bed, cactus halLiftover) and
the SAME overlap rule (bedtools intersect -wa -u) as rebuild_zamore_COMBINED.py Step 4.

Per-rep PICB xlsx 'clusters' sheet uses chr-prefixed seqnames ('chr1'); projected loci
+ combined BEDs are bare-numeric ('1') -> strip 'chr' before intersect (same as the
combined run). Completeness of each replicate xlsx is verified (chr1-19 + X present)
to guard against the known PICB concurrency-bug corruption.

Output: combined_rebuild/replicate_pct_expressed.csv
        (strain,timepoint,replicate,n_clusters,n_chroms,n_expressed,pct_expressed)
"""
import os, glob, subprocess, tempfile
import pandas as pd, numpy as np

# bedtools from the snakemake env (the base miniconda3/bin has none)
os.environ["PATH"] = "/mnt/home3/miska/nm667/miniconda3/envs/snakemake/bin:" + os.environ["PATH"]

BASE = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
PICB = f"{BASE}/results/picb_result"
CR   = f"{BASE}/analysis/claude_biomni_analysis/all_strains_pangenome/combined_rebuild"
CB   = f"{BASE}/analysis/claude_biomni_analysis/all_strains_pangenome/combined_beds"

STRAINS = ["129S1_SvImJ","AKR_J","A_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ",
           "CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ",
           "SPRET_EiJ","WSB_EiJ"]
TPMAP = {"E16.5":"16.5dpc","P12.5":"12.5dpp","P20.5":"20.5dpp"}
MAJOR = {f"chr{c}" for c in list(range(1,20))+["X"]}

# denominator = total Zamore loci (same as combined matrix)
em = pd.read_csv(f"{CR}/all_strains_expression_matrix.csv")
N_TOT = em.locus.nunique()
print(f"denominator N_TOT (Zamore loci) = {N_TOT}")

# combined % expressed per strain x tp (for the sanity comparison combined vs rep-mean)
comb_pct = (em.assign(e=(em.status=="expressed"))
              .pivot_table(index="strain", columns="timepoint", values="e", aggfunc="mean")*100)

rows = []
for strain in STRAINS:
    loci_bed = f"{CR}/_zamore_{strain}.bed"
    if not os.path.exists(loci_bed):
        print(f"  !! missing projected loci {loci_bed}"); continue
    # sort projected loci once
    dl = pd.read_csv(loci_bed, sep="\t", header=None, names=["chr","start","end","name"])
    dl["chr"] = dl["chr"].astype(str)
    sorted_loci = f"{CR}/_tmp_repbuild_loci_{strain}.bed"
    dl.sort_values(["chr","start"]).to_csv(sorted_loci, sep="\t", index=False, header=False)
    for tp, dpc in TPMAP.items():
        for rep in (1,2,3):
            xlsx = f"{PICB}/{strain}/{strain}-{dpc}.{rep}.xlsx"
            if not os.path.exists(xlsx):
                continue
            cl = pd.read_excel(xlsx, sheet_name="clusters")
            seq = cl["seqnames"].astype(str)
            n_major = len(set(seq) & MAJOR)
            if n_major < 20:
                # PICB concurrency-bug corruption: a whole chromosome is missing -> this replicate's cluster
                # set undercounts loci on that chromosome. Exclude it from the replicate SD (the combined-run
                # bar height is unaffected; verified chr1-19+X complete there).
                print(f"  !! {strain}-{dpc}.{rep}: only {n_major}/20 major chroms (missing {sorted(MAJOR-set(seq))}) -> EXCLUDED from replicate SD")
                continue
            bed = pd.DataFrame({"chr": seq.str.replace(r"^chr","",regex=True),
                                "start": cl["start"].astype(int)-1,   # PICB/GRanges 1-based -> BED 0-based
                                "end":   cl["end"].astype(int)})
            cbed = f"{CR}/_tmp_repbuild_cl_{strain}_{dpc}_{rep}.bed"
            bed.sort_values(["chr","start"]).to_csv(cbed, sep="\t", index=False, header=False)
            r = subprocess.run(["bedtools","intersect","-a",sorted_loci,"-b",cbed,"-wa","-u"],
                               capture_output=True, text=True)
            expr = {l.split("\t")[3] for l in r.stdout.strip().split("\n") if l}
            n_expr = len(expr)
            rows.append({"strain":strain,"timepoint":tp,"replicate":rep,
                         "n_clusters":len(cl),"n_chroms":n_major,
                         "n_expressed":n_expr,"pct_expressed":round(100*n_expr/N_TOT,3)})
            os.remove(cbed)
    os.remove(sorted_loci)
    print(f"  {strain}: done ({sum(1 for x in rows if x['strain']==strain)} reps)")

out = pd.DataFrame(rows)
out.to_csv(f"{CR}/replicate_pct_expressed.csv", index=False)
print(f"\nwrote {CR}/replicate_pct_expressed.csv  ({len(out)} rows)")

# sanity: combined vs replicate-mean (should be close; combined >= rep-mean as pooled reads detect more)
print("\nstrain x tp : combined%  rep-mean%  rep-SD%  (n reps)")
for strain in STRAINS:
    for tp in TPMAP:
        sub = out[(out.strain==strain)&(out.timepoint==tp)]
        if not len(sub): continue
        cb = comb_pct.loc[strain, tp] if strain in comb_pct.index else float("nan")
        print(f"  {strain:13s} {tp}: {cb:5.1f}   {sub.pct_expressed.mean():5.1f}   "
              f"{sub.pct_expressed.std(ddof=1):4.1f}   (n={len(sub)})")
d = out.pct_expressed
print(f"\noverall rep pct range: {d.min():.1f}-{d.max():.1f}%")
