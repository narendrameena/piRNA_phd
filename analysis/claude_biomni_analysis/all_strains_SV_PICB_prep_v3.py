#!/usr/bin/env python3
"""
Multi-strain SV × piRNA expression — v3: PICB C57BL_6NJ clusters as reference loci.

Replaces Zamore 214-locus annotation with all PICB-detected clusters from C57BL_6NJ
across all 3 timepoints (E16.5, P12.5, P20.5), using 2-of-3 consensus per timepoint
merged to a non-redundant union.

Stage assignment:
  - Clusters overlapping a Zamore locus: inherit Zamore stage
  - Novel clusters (no Zamore overlap):
    E16.5-only        → Prepachytene
    P12.5-only/P20.5  → Pachytene
    Multi-timepoint   → Pachytene

SVs from pangenome VCF (same as v2).
Outputs separate CSVs: all_strains_expression_matrix_picb.csv, all_strains_SV_matrix_picb.csv
"""
import os, re, subprocess, tempfile
import pandas as pd
import numpy as np

BASE      = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
VCF       = ("/mnt/beegfs/scratch/miska/nm667/inProgress/mice_PiRNA"
             "/results/pangenome/output/mouse_17strain_pangenome.raw.vcf.gz")
PICB_DIR  = f"{BASE}/results/picb_result"
CHAIN_DIR = f"{BASE}/resources/REL-2205-Assembly/GRCm39_chains"
LIFTOVER  = f"{BASE}/workflow/scripts/liftOver"
OUT       = f"{BASE}/analysis/claude_biomni_analysis/all_strains_pangenome"
ZAMORE_XLS = (f"{BASE}/resources/zamore_piRNAs/"
              "piRNA_gene_annotation-modified-in-orange-with-Wasik-et-al-checks.xlsx")
os.makedirs(OUT, exist_ok=True)

SV_WIN  = 50000
SV_MIN  = 300

VCF_SAMPLES = ["129S1_SvImJ","AKR_J","A_J","BALB_cJ","C3H_HeJ","C57BL_6NJ",
               "CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ",
               "NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
STRAINS    = VCF_SAMPLES[:]
TIMEPOINTS = {"E16.5":"16.5dpc","P12.5":"12.5dpp","P20.5":"20.5dpp"}

# ── Zamore stage lookup ───────────────────────────────────────────────────────
df_xlsx = pd.read_excel(ZAMORE_XLS, sheet_name='Mus musculus', header=None, skiprows=1)
df_xlsx.columns = ['chr','start','end','name','score','strand',
                   'thickStart','thickEnd','itemRgb','blockCount',
                   'blockSizes','blockStarts','stage']
df_xlsx['base_gene'] = df_xlsx['name'].str.replace(r'\.\d+$','',regex=True)
df_xlsx['stage'] = df_xlsx['stage'].astype(str).str.strip()
stage_by_locus = df_xlsx.groupby('base_gene')['stage'].first().to_dict()

# Zamore loci BED (GRCm39, numeric chr) for stage intersection
zamore_noprefix = f"{OUT}/_zamore_loci_noprefix.bed"
if not os.path.exists(zamore_noprefix):
    raise FileNotFoundError(f"Missing {zamore_noprefix} — run v2 prep first")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 0 — Build C57BL_6NJ PICB union reference (all 3 timepoints)
# ══════════════════════════════════════════════════════════════════════════════
print("Step 0: Building C57BL_6NJ PICB reference loci...")

# Existing 2-of-3 BEDs for C57BL_6NJ (chr-prefixed)
tp_beds = {
    "E16.5": f"{OUT}/C57BL_6NJ_E16.5_2of3.bed",
    "P12.5": f"{OUT}/C57BL_6NJ_P12.5_2of3.bed",
    "P20.5": f"{OUT}/C57BL_6NJ_P20.5_2of3.bed",
}

# Strip chr prefix and tag each cluster with its detection timepoint(s)
tp_dfs = {}
for tp, bed in tp_beds.items():
    if os.path.exists(bed):
        df = pd.read_csv(bed, sep='\t', header=None, names=['chr','start','end'])
        df['chr'] = df['chr'].str.replace('^chr','', regex=True)
        # Keep only standard numeric chromosomes + X/Y/M
        df = df[df['chr'].str.match(r'^\d+$|^X$|^Y$|^M$')].copy()
        tp_dfs[tp] = df
        print(f"  {tp}: {len(df)} clusters")
    else:
        print(f"  WARNING: {bed} not found")

# Merge across timepoints: union with bedtools merge
all_tp_bed = f"{OUT}/_picb_all_tp.bed"
merged_frames = []
for tp, df in tp_dfs.items():
    df2 = df.copy(); df2['tp'] = tp
    merged_frames.append(df2)
all_tp = pd.concat(merged_frames, ignore_index=True)
# Write and sort
all_tp_sorted = all_tp.sort_values(['chr', 'start'], key=lambda c: c if c.name!='chr' else c.str.zfill(3))
all_tp_sorted[['chr','start','end']].to_csv(all_tp_bed, sep='\t', index=False, header=False)

# Union merge
union_bed = f"{OUT}/_picb_union.bed"
subprocess.run(
    f"sort -k1,1 -k2,2n {all_tp_bed} | bedtools merge > {union_bed}",
    shell=True, check=True
)
os.remove(all_tp_bed)

picb_union = pd.read_csv(union_bed, sep='\t', header=None, names=['chr','start','end'])
print(f"  Union: {len(picb_union)} clusters after merge")

# ── Stage assignment ──────────────────────────────────────────────────────────
# 1. Intersect with Zamore noprefix BED to inherit stage
zamore_df = pd.read_csv(zamore_noprefix, sep='\t', header=None,
                         names=['chr','start','end','name'])
zamore_df['stage'] = zamore_df['name'].str.replace(r'\.\d+$','',regex=True).map(stage_by_locus)
zamore_bed_for_stage = f"{OUT}/_zamore_stage.bed"
zamore_df[['chr','start','end','stage']].to_csv(
    zamore_bed_for_stage, sep='\t', index=False, header=False)

# Sort picb_union for intersection
picb_sorted_bed = f"{OUT}/_picb_union_sorted.bed"
picb_union.sort_values(['chr','start'], key=lambda c: c if c.name!='chr' else c.str.zfill(3)).to_csv(
    picb_sorted_bed, sep='\t', index=False, header=False)

# Intersect picb_union with zamore to get stage for overlapping clusters
r_stage = subprocess.run(
    ['bedtools','intersect','-a',picb_sorted_bed,'-b',zamore_bed_for_stage,
     '-wa','-wb','-f','0.10'],   # 10% overlap threshold
    capture_output=True, text=True
)

stage_from_zamore = {}  # (chr,start,end) → stage
for line in r_stage.stdout.strip().split('\n'):
    if not line: continue
    f = line.split('\t')
    if len(f) < 7: continue
    key = (f[0], int(f[1]), int(f[2]))
    stage_val = f[6]
    if stage_val in ('Pachytene','Prepachytene','Hybrid'):
        stage_from_zamore[key] = stage_val

# For non-overlapping: assign by timepoint of detection
# Which timepoints detect each cluster?
tp_detection = {}
for tp, df in tp_dfs.items():
    for _, row in df.iterrows():
        # Find which union cluster this belongs to (rough match)
        pass  # we'll handle this via bedtools

# Simpler approach: bedtools intersect picb_union with each timepoint BED
for tp, df in tp_dfs.items():
    tp_bed = f"{OUT}/_picb_tp_{tp}.bed"
    df.sort_values(['chr','start'], key=lambda c: c if c.name!='chr' else c.str.zfill(3)).to_csv(
        tp_bed, sep='\t', index=False, header=False)
    r_tp = subprocess.run(
        ['bedtools','intersect','-a',picb_sorted_bed,'-b',tp_bed,'-wa','-u'],
        capture_output=True, text=True
    )
    for line in r_tp.stdout.strip().split('\n'):
        if not line: continue
        f = line.split('\t')
        if len(f) < 3: continue
        key = (f[0], int(f[1]), int(f[2]))
        if key not in tp_detection:
            tp_detection[key] = set()
        tp_detection[key].add(tp)
    os.remove(tp_bed)

def assign_stage(key, zamore_stage, tp_set):
    if key in zamore_stage:
        return zamore_stage[key]
    if tp_set == {'E16.5'}:
        return 'Prepachytene'
    if tp_set == {'P12.5'} or tp_set == {'P20.5'}:
        return 'Pachytene'
    return 'Pachytene'  # multi-timepoint → pachytene

picb_union['stage'] = [
    assign_stage((row.chr, row.start, row.end), stage_from_zamore,
                 tp_detection.get((row.chr, row.start, row.end), {'E16.5'}))
    for _, row in picb_union.iterrows()
]

# Assign unique cluster IDs
picb_union['cluster_id'] = (picb_union['chr'].astype(str) + '_' +
                             picb_union['start'].astype(str) + '_' +
                             picb_union['end'].astype(str))
picb_union = picb_union[picb_union['stage'].isin(['Pachytene','Prepachytene','Hybrid'])].copy()
n_loci = len(picb_union)
print(f"  Final reference: {n_loci} clusters (stage-annotated)")
print(picb_union['stage'].value_counts().to_string())

# Write reference BED (4 col: chr, start, end, cluster_id)
ref_bed = f"{OUT}/_picb_ref_loci.bed"
picb_union[['chr','start','end','cluster_id']].to_csv(
    ref_bed, sep='\t', index=False, header=False)

# Clean up temp files
for f in [zamore_bed_for_stage, picb_sorted_bed, union_bed]:
    if os.path.exists(f):
        os.remove(f)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Extract SVs from pangenome VCF (same as v2)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Step 1: Extracting SVs from pangenome VCF ──")

regions_bed = f"{OUT}/_picb_regions_50kb.bed"
reg = picb_union.copy()
reg['rs'] = (reg['start'] - SV_WIN).clip(lower=0).astype(int)
reg['re'] = (reg['end']   + SV_WIN).astype(int)
reg[['chr','rs','re']].sort_values(
    ['chr','rs'], key=lambda c: c if c.name!='chr' else c.str.zfill(3)
).drop_duplicates().to_csv(regions_bed, sep='\t', index=False, header=False)

cmd = ["bcftools","view","-H","-R",regions_bed,VCF]
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1<<20)

sv_records = {s: [] for s in STRAINS}
n_sv_total = 0
for line in proc.stdout:
    f = line.rstrip('\n').split('\t')
    if len(f) < 10: continue
    chrom=f[0]; pos=int(f[1]); ref=f[3]; alts=f[4].split(','); gts=f[9:]
    best_diff=0; best_type=None; best_ai=None
    for ai, alt in enumerate(alts):
        diff=len(alt)-len(ref)
        if abs(diff)>abs(best_diff):
            best_diff=diff; best_type='INS' if diff>0 else 'DEL'; best_ai=str(ai+1)
    if abs(best_diff)<SV_MIN: continue
    sv_size=abs(best_diff); bstart=pos-1; bend=pos-1+len(ref)
    n_sv_total+=1
    for i, strain in enumerate(VCF_SAMPLES):
        gt=gts[i].split(':')[0] if i<len(gts) else '.'
        if gt=='.' or gt=='0': continue
        sv_records[strain].append((chrom, bstart, bend, best_type, sv_size))
proc.wait()
print(f"  Total SV records in regions: {n_sv_total}")

# Write per-strain SV BEDs (tagged _picb to avoid overwriting v2)
sv_beds = {}
for strain in STRAINS:
    sv_path = f"{OUT}/_sv_picb_{strain}.bed"
    df_sv = pd.DataFrame(sv_records[strain],
                         columns=['chr','start','end','sv_type','sv_size'])
    df_sv = df_sv.drop_duplicates().sort_values(['chr','start'])
    df_sv.to_csv(sv_path, sep='\t', index=False, header=False)
    sv_beds[strain] = sv_path
    print(f"  {strain}: {len(df_sv)} SVs")

os.remove(regions_bed)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Merge PICB 2-of-3 per strain × timepoint (already done in v2, reuse)
# ══════════════════════════════════════════════════════════════════════════════
def xlsx_to_bed(xlsx_path, tmp_path):
    df = pd.read_excel(xlsx_path, usecols=['seqnames','start','end'])
    df = df.rename(columns={'seqnames':'chr'})
    df['start'] = df['start'] - 1
    df = df[['chr','start','end']].dropna().astype({'start':int,'end':int})
    df.to_csv(tmp_path, sep='\t', index=False, header=False)

def merge_2of3(rep_paths, out_path):
    sorted_beds = []
    for p in rep_paths:
        sp = p + '.sorted'
        subprocess.run(f"sort -k1,1 -k2,2n {p} > {sp}", shell=True, check=True)
        sorted_beds.append(sp)
    cmd = (f"bedtools multiinter -i {' '.join(sorted_beds)} "
           f"| awk '$4>=2' | cut -f1-3 "
           f"| sort -k1,1 -k2,2n | bedtools merge > {out_path}")
    subprocess.run(cmd, shell=True, check=True)
    for sp in sorted_beds:
        os.remove(sp)

print("\n── Step 2: Checking/building PICB 2-of-3 BEDs ──")
picb_merged = {}
with tempfile.TemporaryDirectory() as tmp:
    for strain in STRAINS:
        for tp_label, tp_code in TIMEPOINTS.items():
            out_bed = f"{OUT}/{strain}_{tp_label}_2of3.bed"
            if os.path.exists(out_bed):
                picb_merged[(strain, tp_label)] = out_bed
                continue
            rep_paths = []
            for rep in [1, 2, 3]:
                xls = f"{PICB_DIR}/{strain}/{strain}-{tp_code}.{rep}.xlsx"
                if not os.path.exists(xls): continue
                bed = f"{tmp}/{strain}_{tp_label}_rep{rep}.bed"
                xlsx_to_bed(xls, bed)
                rep_paths.append(bed)
            if len(rep_paths) < 2: continue
            if len(rep_paths) == 2:
                sp1, sp2 = rep_paths[0]+'.s', rep_paths[1]+'.s'
                subprocess.run(f"sort -k1,1 -k2,2n {rep_paths[0]} > {sp1}", shell=True)
                subprocess.run(f"sort -k1,1 -k2,2n {rep_paths[1]} > {sp2}", shell=True)
                subprocess.run(f"bedtools intersect -a {sp1} -b {sp2} -u | bedtools merge > {out_bed}",
                               shell=True)
            else:
                merge_2of3(rep_paths, out_bed)
            picb_merged[(strain, tp_label)] = out_bed
print(f"  Total PICB merged BEDs: {len(picb_merged)}")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — LiftOver PICB reference loci → each strain
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Step 3: LiftOver PICB reference loci ──")
lifted = {}
lift_stats = {}
for strain in STRAINS:
    chain = f"{CHAIN_DIR}/GRCm39_{strain}_chromosomes_MT_unplaced.chain"
    out_lifted   = f"{OUT}/_picb_{strain}.bed"
    out_unmapped = f"{OUT}/_picb_{strain}_unmapped.bed"
    subprocess.run([LIFTOVER, ref_bed, chain, out_lifted, out_unmapped],
                   capture_output=True)
    df_lift = pd.DataFrame(columns=['chr','start','end','name'])
    if os.path.exists(out_lifted) and os.path.getsize(out_lifted) > 0:
        df_lift = pd.read_csv(out_lifted, sep='\t', header=None,
                              names=['chr','start','end','name'])
        df_lift['chr'] = df_lift['chr'].str.replace(r'^[^#]+#\d+#','',regex=True)
    lift_stats[strain] = {'n_lifted':len(df_lift),'n_total':n_loci}
    lifted[strain] = df_lift.set_index('name')[['chr','start','end']].to_dict('index')
    pct = 100*len(df_lift)/n_loci
    print(f"  {strain}: {len(df_lift)}/{n_loci} ({pct:.0f}%)")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Expression matrix
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Step 4: Expression matrix ──")
expr_rows = []
for strain in STRAINS:
    lf_path = f"{OUT}/_picb_{strain}.bed"
    df_lift = pd.DataFrame(columns=['chr','start','end','name'])
    if os.path.exists(lf_path) and os.path.getsize(lf_path) > 0:
        df_lift = pd.read_csv(lf_path, sep='\t', header=None,
                              names=['chr','start','end','name'])
        df_lift['chr'] = df_lift['chr'].str.replace(r'^[^#]+#\d+#','',regex=True)

    for tp_label in TIMEPOINTS:
        key = (strain, tp_label)
        if key not in picb_merged:
            continue
        sorted_lift = f"{OUT}/_tmp_picb_{strain}_{tp_label}.bed"
        df_lift.sort_values(['chr','start']).to_csv(
            sorted_lift, sep='\t', index=False, header=False)
        r = subprocess.run(
            ['bedtools','intersect','-a',sorted_lift,'-b',picb_merged[key],'-wa','-u'],
            capture_output=True, text=True)
        expressed = {l.split('\t')[3] for l in r.stdout.strip().split('\n') if l}
        os.remove(sorted_lift)

        for _, row in picb_union.iterrows():
            cid = row['cluster_id']
            if cid not in lifted[strain]:
                status = 'not_lifted'
            elif cid in expressed:
                status = 'expressed'
            else:
                status = 'not_expressed'
            expr_rows.append({'locus': cid, 'stage': row['stage'],
                               'strain': strain, 'timepoint': tp_label,
                               'status': status})

expr_df = pd.DataFrame(expr_rows)
expr_df.to_csv(f"{OUT}/all_strains_expression_matrix_picb.csv", index=False)
print(f"  Saved: {len(expr_df)} rows")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — SV matrix
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Step 5: SV matrix ──")
sv_rows = []
for strain in STRAINS:
    sv_bed = sv_beds[strain]
    for label, win in [('direct',0),('10kb',10000),('50kb',50000)]:
        if win == 0:
            r = subprocess.run(['bedtools','intersect','-a',ref_bed,'-b',sv_bed,'-wo'],
                               capture_output=True, text=True)
        else:
            r = subprocess.run(['bedtools','window','-a',ref_bed,'-b',sv_bed,'-w',str(win)],
                               capture_output=True, text=True)
        hits = {}
        for line in r.stdout.strip().split('\n'):
            if not line: continue
            c = line.split('\t')
            gene=c[3]; sv_type=c[7]; sv_size=int(c[8])
            if gene not in hits:
                hits[gene] = {'INS_bp':0,'DEL_bp':0,'INS_n':0,'DEL_n':0}
            hits[gene][f'{sv_type}_bp'] += sv_size
            hits[gene][f'{sv_type}_n']  += 1
        for _, row in picb_union.iterrows():
            cid = row['cluster_id']
            h = hits.get(cid, {})
            sv_rows.append({'locus':cid,'stage':row['stage'],'strain':strain,
                            'window':label,
                            'INS_bp':h.get('INS_bp',0),'INS_n':h.get('INS_n',0),
                            'DEL_bp':h.get('DEL_bp',0),'DEL_n':h.get('DEL_n',0),
                            'has_SV':(h.get('INS_bp',0)+h.get('DEL_bp',0))>0})

sv_df = pd.DataFrame(sv_rows)
sv_df.to_csv(f"{OUT}/all_strains_SV_matrix_picb.csv", index=False)
print(f"  Saved: {len(sv_df)} rows")

# Quick summary
print("\n── Expression vs SV (direct, P20.5, Pachytene, M.m. domesticus) ──")
MMD = ["C57BL_6NJ","DBA_2J","BALB_cJ","A_J","CBA_J","129S1_SvImJ",
       "NOD_ShiLtJ","AKR_J","C3H_HeJ","NZO_HlLtJ","LP_J","FVB_NJ"]
e = expr_df[(expr_df['timepoint']=='P20.5') & (expr_df['stage']=='Pachytene') &
            expr_df['strain'].isin(MMD)]
s = sv_df[sv_df['window']=='direct'][['locus','strain','has_SV']]
m = e.merge(s, on=['locus','strain'], how='left')
m['has_SV'] = m['has_SV'].fillna(False)
print(m.groupby(['has_SV','status']).size().unstack(fill_value=0))
print(f"\nDone — outputs in {OUT}/")
print(f"  Reference: {n_loci} PICB C57BL_6NJ clusters")
