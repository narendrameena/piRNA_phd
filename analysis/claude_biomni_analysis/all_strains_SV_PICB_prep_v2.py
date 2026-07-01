#!/usr/bin/env python3
"""
Multi-strain SV × piRNA expression — v2: SVs from pangenome VCF.

Replaces chain-file gap parsing with per-sample genotype calls from the
17-strain pangenome VCF. Uses tabix index to query only Zamore locus
regions (± 50 kb), so scans <1% of the 2.8 GB file.

Steps:
  1. Extract per-strain SVs (≥300 bp INS/DEL) from pangenome VCF
  2. Merge PICB clusters 2-of-3 per strain × timepoint
  3. LiftOver Zamore loci (GRCm39) → each strain
  4. Expression matrix
  5. SV matrix
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
os.makedirs(OUT, exist_ok=True)

SV_WIN    = 50000     # query window around loci for VCF lookup
SV_MIN    = 300       # minimum SV size (bp)

# Sample order in VCF (from header)
VCF_SAMPLES = ["129S1_SvImJ","AKR_J","A_J","BALB_cJ","C3H_HeJ","C57BL_6NJ",
               "CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ",
               "NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]

STRAINS    = VCF_SAMPLES[:]
TIMEPOINTS = {"E16.5": "16.5dpc", "P12.5": "12.5dpp", "P20.5": "20.5dpp"}

ZAMORE_BED = f"{BASE}/analysis/claude_biomni_analysis/C57BL_6NJ_pangenome/_loci_mm39_noprefix.bed"
ZAMORE_XLS = (f"{BASE}/resources/zamore_piRNAs/"
              "piRNA_gene_annotation-modified-in-orange-with-Wasik-et-al-checks.xlsx")

# ── Stage lookup ──────────────────────────────────────────────────────────────
df_xlsx = pd.read_excel(ZAMORE_XLS, sheet_name='Mus musculus',
                        header=None, skiprows=1)
df_xlsx.columns = ['chr','start','end','name','score','strand',
                   'thickStart','thickEnd','itemRgb','blockCount',
                   'blockSizes','blockStarts','stage']
df_xlsx['base_gene'] = df_xlsx['name'].str.replace(r'\.\d+$','',regex=True)
df_xlsx['stage'] = df_xlsx['stage'].astype(str).str.strip()
stage_by_locus = df_xlsx.groupby('base_gene')['stage'].first().to_dict()

zamore_raw = pd.read_csv(ZAMORE_BED, sep='\t', header=None,
                         names=['chr','start','end','name'])
zamore_raw['base_gene'] = zamore_raw['name'].str.replace(r'\.\d+$','',regex=True)
zamore_loci = (zamore_raw.groupby('base_gene')
               .agg(chr=('chr','first'), start=('start','min'), end=('end','max'))
               .reset_index())
zamore_loci['stage'] = zamore_loci['base_gene'].map(stage_by_locus)
zamore_loci = zamore_loci[zamore_loci['stage'].isin(
    ['Pachytene','Prepachytene','Hybrid'])].copy().reset_index(drop=True)
n_loci = len(zamore_loci)
print(f"Zamore loci: {n_loci}")

zamore_bed_path = f"{OUT}/_zamore_loci_noprefix.bed"
zamore_loci[['chr','start','end','base_gene']].to_csv(
    zamore_bed_path, sep='\t', index=False, header=False)


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Extract SVs from pangenome VCF (tabix-indexed, per-strain genotypes)
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Step 1: Extracting SVs from pangenome VCF ──")

# Build padded regions BED for bcftools -R
# Genome chromsizes not needed: bcftools clips to contig boundaries
regions_bed = f"{OUT}/_zamore_regions_50kb.bed"
reg = zamore_loci.copy()
reg['rs'] = (reg['start'] - SV_WIN).clip(lower=0).astype(int)
reg['re'] = (reg['end']   + SV_WIN).astype(int)
reg[['chr','rs','re']].to_csv(regions_bed, sep='\t', index=False, header=False)

# Run bcftools view restricted to those regions, pipe to Python
cmd = ["bcftools", "view", "-H", "-R", regions_bed, VCF]
proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True,
                        bufsize=1 << 20)

sv_records = {s: [] for s in STRAINS}
n_sv_total = 0

for line in proc.stdout:
    f = line.rstrip('\n').split('\t')
    if len(f) < 10:
        continue
    chrom = f[0]
    pos   = int(f[1])
    ref   = f[3]
    alts  = f[4].split(',')
    gts   = f[9:]          # per-sample GT (same order as VCF_SAMPLES)

    # Find the largest SV among all ALT alleles
    best_diff = 0
    best_type = None
    best_alt_idx = None      # 1-based ALT index for GT matching
    for ai, alt in enumerate(alts):
        diff = len(alt) - len(ref)
        if abs(diff) > abs(best_diff):
            best_diff = diff
            best_type = 'INS' if diff > 0 else 'DEL'
            best_alt_idx = str(ai + 1)

    if abs(best_diff) < SV_MIN:
        continue

    sv_size = abs(best_diff)
    # BED coords: (POS-1, POS-1+len(REF))  — 0-based half-open
    bstart = pos - 1
    bend   = pos - 1 + len(ref)   # for INS: bend=pos (1-bp insertion point)

    n_sv_total += 1
    for i, strain in enumerate(VCF_SAMPLES):
        gt = gts[i].split(':')[0] if i < len(gts) else '.'
        # GT = "1" or "2" etc. means alt allele; "0" = ref; "." = missing
        if gt == '.' or gt == '0':
            continue
        sv_records[strain].append((chrom, bstart, bend, best_type, sv_size))

proc.wait()
print(f"  Total SV records in regions: {n_sv_total}")

# Write per-strain SV BEDs
sv_beds = {}
for strain in STRAINS:
    sv_path = f"{OUT}/_sv_{strain}.bed"
    df_sv = pd.DataFrame(sv_records[strain],
                         columns=['chr','start','end','sv_type','sv_size'])
    df_sv = df_sv.drop_duplicates().sort_values(['chr','start'])
    df_sv.to_csv(sv_path, sep='\t', index=False, header=False)
    sv_beds[strain] = sv_path
    print(f"  {strain}: {len(df_sv)} SVs "
          f"(INS={(df_sv.sv_type=='INS').sum()}, "
          f"DEL={(df_sv.sv_type=='DEL').sum()})")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Merge PICB XLSXs 2-of-3 per strain × timepoint
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
        subprocess.run(f"sort -k1,1 -k2,2n {p} > {sp}",
                       shell=True, check=True)
        sorted_beds.append(sp)
    cmd = (f"bedtools multiinter -i {' '.join(sorted_beds)} "
           f"| awk '$4>=2' | cut -f1-3 "
           f"| sort -k1,1 -k2,2n | bedtools merge > {out_path}")
    subprocess.run(cmd, shell=True, check=True)
    for sp in sorted_beds:
        os.remove(sp)


print("\n── Step 2: Merging PICB XLSXs (2-of-3) ──")
picb_merged = {}
with tempfile.TemporaryDirectory() as tmp:
    for strain in STRAINS:
        for tp_label, tp_code in TIMEPOINTS.items():
            rep_paths = []
            for rep in [1, 2, 3]:
                xls = f"{PICB_DIR}/{strain}/{strain}-{tp_code}.{rep}.xlsx"
                if not os.path.exists(xls):
                    continue
                bed = f"{tmp}/{strain}_{tp_label}_rep{rep}.bed"
                xlsx_to_bed(xls, bed)
                rep_paths.append(bed)
            if len(rep_paths) < 2:
                continue
            out_bed = f"{OUT}/{strain}_{tp_label}_2of3.bed"
            if len(rep_paths) == 2:
                sp1, sp2 = rep_paths[0]+'.s', rep_paths[1]+'.s'
                subprocess.run(f"sort -k1,1 -k2,2n {rep_paths[0]} > {sp1}", shell=True, check=True)
                subprocess.run(f"sort -k1,1 -k2,2n {rep_paths[1]} > {sp2}", shell=True, check=True)
                subprocess.run(f"bedtools intersect -a {sp1} -b {sp2} -u | bedtools merge > {out_bed}",
                               shell=True, check=True)
            else:
                merge_2of3(rep_paths, out_bed)
            picb_merged[(strain, tp_label)] = out_bed
        print(f"  {strain}: merged")
print(f"  Total merged BEDs: {len(picb_merged)}")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — LiftOver Zamore loci → each strain
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Step 3: LiftOver Zamore loci ──")
lifted = {}
lift_stats = {}

for strain in STRAINS:
    chain = f"{CHAIN_DIR}/GRCm39_{strain}_chromosomes_MT_unplaced.chain"
    out_lifted   = f"{OUT}/_zamore_{strain}.bed"
    out_unmapped = f"{OUT}/_zamore_{strain}_unmapped.bed"
    subprocess.run([LIFTOVER, zamore_bed_path, chain, out_lifted, out_unmapped],
                   capture_output=True)
    df_lift = pd.DataFrame(columns=['chr','start','end','name'])
    if os.path.exists(out_lifted) and os.path.getsize(out_lifted) > 0:
        df_lift = pd.read_csv(out_lifted, sep='\t', header=None,
                              names=['chr','start','end','name'])
        df_lift['chr'] = df_lift['chr'].str.replace(r'^[^#]+#\d+#','',regex=True)
    lift_stats[strain] = {'n_lifted': len(df_lift), 'n_total': n_loci}
    lifted[strain] = df_lift.set_index('name')[['chr','start','end']].to_dict('index')
    pct = 100 * len(df_lift) / n_loci
    print(f"  {strain}: {len(df_lift)}/{n_loci} ({pct:.0f}%)")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Expression matrix
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Step 4: Expression matrix ──")
expr_rows = []
for strain in STRAINS:
    lf_path = f"{OUT}/_zamore_{strain}.bed"
    df_lift = pd.DataFrame(columns=['chr','start','end','name'])
    if os.path.exists(lf_path) and os.path.getsize(lf_path) > 0:
        df_lift = pd.read_csv(lf_path, sep='\t', header=None,
                              names=['chr','start','end','name'])
        df_lift['chr'] = df_lift['chr'].str.replace(r'^[^#]+#\d+#','',regex=True)

    for tp_label in TIMEPOINTS:
        key = (strain, tp_label)
        if key not in picb_merged:
            continue
        sorted_lift = f"{OUT}/_tmp_{strain}_{tp_label}.bed"
        df_lift.sort_values(['chr','start']).to_csv(
            sorted_lift, sep='\t', index=False, header=False)
        r = subprocess.run(
            ['bedtools','intersect','-a',sorted_lift,'-b',picb_merged[key],'-wa','-u'],
            capture_output=True, text=True)
        expressed = {l.split('\t')[3] for l in r.stdout.strip().split('\n') if l}
        os.remove(sorted_lift)

        for _, row in zamore_loci.iterrows():
            gene = row['base_gene']
            if gene not in lifted[strain]:
                status = 'not_lifted'
            elif gene in expressed:
                status = 'expressed'
            else:
                status = 'not_expressed'
            expr_rows.append({'locus': gene, 'stage': row['stage'],
                               'strain': strain, 'timepoint': tp_label,
                               'status': status})

expr_df = pd.DataFrame(expr_rows)
expr_df.to_csv(f"{OUT}/all_strains_expression_matrix.csv", index=False)
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
            r = subprocess.run(
                ['bedtools','intersect','-a',zamore_bed_path,'-b',sv_bed,'-wo'],
                capture_output=True, text=True)
        else:
            r = subprocess.run(
                ['bedtools','window','-a',zamore_bed_path,'-b',sv_bed,'-w',str(win)],
                capture_output=True, text=True)
        hits = {}
        for line in r.stdout.strip().split('\n'):
            if not line:
                continue
            c = line.split('\t')
            gene = c[3]; sv_type = c[7]; sv_size = int(c[8])
            if gene not in hits:
                hits[gene] = {'INS_bp':0,'DEL_bp':0,'INS_n':0,'DEL_n':0}
            hits[gene][f'{sv_type}_bp'] += sv_size
            hits[gene][f'{sv_type}_n']  += 1
        for _, row in zamore_loci.iterrows():
            gene = row['base_gene']
            h = hits.get(gene, {})
            sv_rows.append({'locus':gene,'stage':row['stage'],'strain':strain,
                            'window':label,
                            'INS_bp':h.get('INS_bp',0),'INS_n':h.get('INS_n',0),
                            'DEL_bp':h.get('DEL_bp',0),'DEL_n':h.get('DEL_n',0),
                            'has_SV':(h.get('INS_bp',0)+h.get('DEL_bp',0))>0})

sv_df = pd.DataFrame(sv_rows)
sv_df.to_csv(f"{OUT}/all_strains_SV_matrix.csv", index=False)
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
