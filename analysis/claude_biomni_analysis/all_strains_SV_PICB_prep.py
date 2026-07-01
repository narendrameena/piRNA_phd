#!/usr/bin/env python3
"""
Multi-strain SV × piRNA expression analysis.

For each of 16 strains:
  1. Extract SVs (≥300 bp INS/DEL) from GRCm39→strain chain file gaps
  2. Merge PICB clusters 2-of-3 per strain × timepoint
  3. Lift Zamore piRNA loci (GRCm39) → each strain via liftOver
  4. Intersect: expression matrix (loci × strains × timepoints)
  5. Intersect: SV matrix (loci × strains, at direct/10kb/50kb windows)

Outputs:
  all_strains_expression_matrix.csv  — per-locus × strain × timepoint expression
  all_strains_SV_matrix.csv          — per-locus × strain SV summary
"""
import os, re, subprocess, tempfile
import pandas as pd
import numpy as np

BASE       = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
PICB_DIR   = f"{BASE}/results/picb_result"
CHAIN_DIR  = f"{BASE}/resources/REL-2205-Assembly/GRCm39_chains"
LIFTOVER   = f"{BASE}/workflow/scripts/liftOver"
OUT        = f"{BASE}/analysis/claude_biomni_analysis/all_strains_pangenome"
os.makedirs(OUT, exist_ok=True)

SV_MIN     = 300   # bp threshold for SV

STRAINS = [
    "129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ",
    "CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ",
    "NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ",
]
TIMEPOINTS = {"E16.5": "16.5dpc", "P12.5": "12.5dpp", "P20.5": "20.5dpp"}

# Zamore loci in GRCm39 (numeric chr names, no prefix)
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

# ── Zamore loci: unique base genes in GRCm39 ──────────────────────────────────
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
print(f"Zamore loci (Pachytene/Prepachytene/Hybrid): {n_loci}")

# Write loci BED (no-chr-prefix for chain compatibility)
zamore_bed_path = f"{OUT}/_zamore_loci_noprefix.bed"
zamore_loci[['chr','start','end','base_gene']].to_csv(
    zamore_bed_path, sep='\t', index=False, header=False)


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Extract SVs per strain from chain file gaps
# ══════════════════════════════════════════════════════════════════════════════
def parse_chain_svs(chain_path, sv_min=300):
    """Parse chain file; return DataFrame of SVs in GRCm39 coords.
    dt > 0 → gap in GRCm39 target = DEL in strain
    dq > 0 → gap in strain query  = INS in strain
    """
    records = []
    chrom = None
    t_pos = 0
    with open(chain_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('chain'):
                parts = line.split()
                # tName is parts[2] (e.g. "1", "2", ..., "X")
                chrom = parts[2]
                t_pos = int(parts[5])   # tStart
            else:
                fields = line.split('\t')
                size = int(fields[0])
                t_pos += size
                if len(fields) == 3:
                    dt, dq = int(fields[1]), int(fields[2])
                    if dt >= sv_min:
                        records.append((chrom, t_pos, t_pos + dt, 'DEL', dt))
                    if dq >= sv_min:
                        records.append((chrom, t_pos, t_pos + 1, 'INS', dq))
                    t_pos += dt
    return pd.DataFrame(records,
                        columns=['chr','start','end','sv_type','sv_size'])


print("\n── Step 1: Extracting SVs from chain files ──")
sv_beds = {}
for strain in STRAINS:
    chain = f"{CHAIN_DIR}/GRCm39_{strain}_chromosomes_MT_unplaced.chain"
    df_sv = parse_chain_svs(chain)
    # Write BED
    sv_path = f"{OUT}/_sv_{strain}.bed"
    df_sv[['chr','start','end','sv_type','sv_size']].to_csv(
        sv_path, sep='\t', index=False, header=False)
    sv_beds[strain] = sv_path
    print(f"  {strain}: {len(df_sv)} SVs "
          f"(INS={( df_sv.sv_type=='INS').sum()}, "
          f"DEL={(df_sv.sv_type=='DEL').sum()})")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Merge PICB XLSXs 2-of-3 per strain × timepoint
# ══════════════════════════════════════════════════════════════════════════════
def xlsx_to_bed(xlsx_path, tmp_path):
    """Read PICB xlsx, write 0-based BED (chr, start-1, end)."""
    df = pd.read_excel(xlsx_path, usecols=['seqnames','start','end'])
    df = df.rename(columns={'seqnames':'chr'})
    df['start'] = df['start'] - 1          # 1-based → 0-based
    df = df[['chr','start','end']].dropna().astype({'start':int,'end':int})
    df.to_csv(tmp_path, sep='\t', index=False, header=False)


def merge_2of3(rep_paths, out_path):
    """Use bedtools multiinter to find regions in ≥2 of 3 replicates."""
    # Sort each rep first
    sorted_beds = []
    for i, p in enumerate(rep_paths):
        sp = p + '.sorted'
        subprocess.run(f"sort -k1,1 -k2,2n {p} > {sp}",
                       shell=True, check=True)
        sorted_beds.append(sp)
    # multiinter then filter ≥2 then merge
    cmd = (f"bedtools multiinter -i {' '.join(sorted_beds)} "
           f"| awk '$4>=2' | cut -f1-3 "
           f"| sort -k1,1 -k2,2n | bedtools merge > {out_path}")
    subprocess.run(cmd, shell=True, check=True)
    for sp in sorted_beds:
        os.remove(sp)


print("\n── Step 2: Merging PICB XLSXs (2-of-3) ──")
picb_merged = {}   # (strain, tp) → bed path
with tempfile.TemporaryDirectory() as tmp:
    for strain in STRAINS:
        for tp_label, tp_code in TIMEPOINTS.items():
            xls_pattern = f"{PICB_DIR}/{strain}/{strain}-{tp_code}."
            rep_paths = []
            for rep in [1, 2, 3]:
                xls = f"{xls_pattern}{rep}.xlsx"
                if not os.path.exists(xls):
                    continue
                bed = f"{tmp}/{strain}_{tp_label}_rep{rep}.bed"
                xlsx_to_bed(xls, bed)
                rep_paths.append(bed)
            if len(rep_paths) < 2:
                print(f"  SKIP {strain} {tp_label}: only {len(rep_paths)} reps")
                continue
            out_bed = f"{OUT}/{strain}_{tp_label}_2of3.bed"
            if len(rep_paths) == 2:
                # Only 2 reps: just intersect both
                sp1 = rep_paths[0] + '.s'
                sp2 = rep_paths[1] + '.s'
                subprocess.run(f"sort -k1,1 -k2,2n {rep_paths[0]} > {sp1}",
                               shell=True, check=True)
                subprocess.run(f"sort -k1,1 -k2,2n {rep_paths[1]} > {sp2}",
                               shell=True, check=True)
                subprocess.run(
                    f"bedtools intersect -a {sp1} -b {sp2} -u "
                    f"| bedtools merge > {out_bed}",
                    shell=True, check=True)
            else:
                merge_2of3(rep_paths, out_bed)
            n = sum(1 for _ in open(out_bed))
            picb_merged[(strain, tp_label)] = out_bed
        print(f"  {strain}: merged")

print(f"  Total merged BEDs: {len(picb_merged)}")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — LiftOver Zamore loci GRCm39 → each strain
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Step 3: LiftOver Zamore loci to each strain ──")
lifted = {}    # strain → {base_gene: (chr, start, end)} in strain coords
lift_stats = {}

for strain in STRAINS:
    chain = f"{CHAIN_DIR}/GRCm39_{strain}_chromosomes_MT_unplaced.chain"
    out_lifted   = f"{OUT}/_zamore_{strain}.bed"
    out_unmapped = f"{OUT}/_zamore_{strain}_unmapped.bed"
    r = subprocess.run(
        [LIFTOVER, zamore_bed_path, chain, out_lifted, out_unmapped],
        capture_output=True, text=True)
    # Parse lifted BED
    df_lift = pd.DataFrame(columns=['chr','start','end','name'])
    if os.path.exists(out_lifted) and os.path.getsize(out_lifted) > 0:
        df_lift = pd.read_csv(out_lifted, sep='\t', header=None,
                              names=['chr','start','end','name'])
        # Strip STRAIN#1# prefix from chr names
        df_lift['chr'] = df_lift['chr'].str.replace(
            r'^[^#]+#\d+#', '', regex=True)
    n_lifted   = len(df_lift)
    n_total    = n_loci
    n_unmapped = n_total - n_lifted
    lift_stats[strain] = {'n_lifted': n_lifted, 'n_total': n_total}
    lifted[strain] = df_lift.set_index('name')[['chr','start','end']].to_dict('index')
    print(f"  {strain}: {n_lifted}/{n_total} loci lifted "
          f"({100*n_lifted/n_total:.0f}%)")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Expression matrix: intersect lifted loci × PICB 2-of-3
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Step 4: Expression matrix ──")
expr_rows = []
for strain in STRAINS:
    df_lift = pd.DataFrame(columns=['chr','start','end','name'])
    lf_path = f"{OUT}/_zamore_{strain}.bed"
    if os.path.exists(lf_path) and os.path.getsize(lf_path) > 0:
        df_lift = pd.read_csv(lf_path, sep='\t', header=None,
                              names=['chr','start','end','name'])
        df_lift['chr'] = df_lift['chr'].str.replace(
            r'^[^#]+#\d+#', '', regex=True)

    for tp_label in TIMEPOINTS:
        key = (strain, tp_label)
        if key not in picb_merged or not os.path.exists(picb_merged[key]):
            continue
        # Sort lifted BED and intersect
        sorted_lift = f"{OUT}/_tmp_{strain}_{tp_label}.bed"
        df_lift_sorted = df_lift.sort_values(['chr','start'])
        df_lift_sorted.to_csv(sorted_lift, sep='\t', index=False, header=False)

        r = subprocess.run(
            ['bedtools', 'intersect',
             '-a', sorted_lift, '-b', picb_merged[key], '-wa', '-u'],
            capture_output=True, text=True)
        expressed_genes = set()
        for line in r.stdout.strip().split('\n'):
            if line:
                expressed_genes.add(line.split('\t')[3])
        os.remove(sorted_lift)

        for _, row in zamore_loci.iterrows():
            gene = row['base_gene']
            if gene not in lifted[strain]:
                status = 'not_lifted'
            elif gene in expressed_genes:
                status = 'expressed'
            else:
                status = 'not_expressed'
            expr_rows.append({
                'locus': gene, 'stage': row['stage'],
                'strain': strain, 'timepoint': tp_label, 'status': status
            })

expr_df = pd.DataFrame(expr_rows)
expr_df.to_csv(f"{OUT}/all_strains_expression_matrix.csv", index=False)
print(f"  Expression matrix: {len(expr_df)} rows saved")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — SV matrix: intersect Zamore loci (GRCm39) × per-strain SVs
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Step 5: SV matrix ──")
sv_rows = []
for strain in STRAINS:
    sv_bed = sv_beds[strain]
    for label, win in [('direct', 0), ('10kb', 10000), ('50kb', 50000)]:
        if win == 0:
            r = subprocess.run(
                ['bedtools', 'intersect',
                 '-a', zamore_bed_path, '-b', sv_bed, '-wo'],
                capture_output=True, text=True)
        else:
            r = subprocess.run(
                ['bedtools', 'window',
                 '-a', zamore_bed_path, '-b', sv_bed, '-w', str(win)],
                capture_output=True, text=True)
        hits = {}
        for line in r.stdout.strip().split('\n'):
            if not line:
                continue
            c = line.split('\t')
            gene    = c[3]
            sv_type = c[7]    # INS or DEL
            sv_size = int(c[8])
            if gene not in hits:
                hits[gene] = {'INS_bp': 0, 'DEL_bp': 0,
                              'INS_n': 0, 'DEL_n': 0}
            hits[gene][f'{sv_type}_bp'] += sv_size
            hits[gene][f'{sv_type}_n']  += 1
        for _, row in zamore_loci.iterrows():
            gene = row['base_gene']
            h    = hits.get(gene, {})
            sv_rows.append({
                'locus': gene, 'stage': row['stage'],
                'strain': strain, 'window': label,
                'INS_bp': h.get('INS_bp', 0), 'INS_n': h.get('INS_n', 0),
                'DEL_bp': h.get('DEL_bp', 0), 'DEL_n': h.get('DEL_n', 0),
                'has_SV': (h.get('INS_bp', 0) + h.get('DEL_bp', 0)) > 0,
            })

sv_df = pd.DataFrame(sv_rows)
sv_df.to_csv(f"{OUT}/all_strains_SV_matrix.csv", index=False)
print(f"  SV matrix: {len(sv_df)} rows saved")

# ── Liftover summary ──────────────────────────────────────────────────────────
lift_df = pd.DataFrame(lift_stats).T
lift_df.index.name = 'strain'
lift_df['pct_lifted'] = 100 * lift_df['n_lifted'] / lift_df['n_total']
lift_df.to_csv(f"{OUT}/liftover_stats.csv")
print("\n── Liftover summary ──")
print(lift_df[['n_lifted','pct_lifted']].to_string())

# ── Quick cross-tab ───────────────────────────────────────────────────────────
print("\n── Expression vs SV (direct overlap, P20.5, Pachytene) ──")
e = expr_df[(expr_df['timepoint'] == 'P20.5') & (expr_df['stage'] == 'Pachytene')].copy()
s = sv_df[(sv_df['window'] == 'direct')][['locus','strain','has_SV']]
merged = e.merge(s, on=['locus','strain'], how='left')
merged['has_SV'] = merged['has_SV'].fillna(False)
tab = merged.groupby(['has_SV','status']).size().unstack(fill_value=0)
print(tab)

print(f"\nDone. Outputs in {OUT}/")
