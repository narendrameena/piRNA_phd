#!/usr/bin/env python3
"""Compute per-locus SV + PICB status for all Zamore piRNA loci and save to CSV."""
import pandas as pd, numpy as np, subprocess, os

OUT = '/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/C57BL_6NJ_pangenome'

# Stage lookup
df_xlsx = pd.read_excel(
    '/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/resources/zamore_piRNAs/'
    'piRNA_gene_annotation-modified-in-orange-with-Wasik-et-al-checks.xlsx',
    sheet_name='Mus musculus', header=None, skiprows=1)
df_xlsx.columns = ['chr','start','end','name','score','strand',
                   'thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts','stage']
df_xlsx['base_gene'] = df_xlsx['name'].str.replace(r'\.\d+$','',regex=True)
df_xlsx['stage'] = df_xlsx['stage'].astype(str).str.strip()
stage_by_locus = df_xlsx.groupby('base_gene')['stage'].first().to_dict()

# Loci in GRCm39 (after mm10→GRCm39 liftover, collapse isoforms)
iso = pd.read_csv(f'{OUT}/zamore_mm39_raw.bed', sep='\t', header=None,
                  names=['chr','start','end','name','score','strand'], usecols=[0,1,2,3])
iso['base_gene'] = iso['name'].str.replace(r'\.\d+$','',regex=True)
loci_mm39 = iso.groupby('base_gene').agg(
    chr=('chr','first'), start=('start','min'), end=('end','max')).reset_index()
loci_mm39['stage'] = loci_mm39['base_gene'].map(stage_by_locus)
loci_mm39 = loci_mm39[loci_mm39['stage'].isin(['Pachytene','Prepachytene','Hybrid'])].copy().reset_index(drop=True)
print(f'Zamore loci in GRCm39: {len(loci_mm39)}')

# Write BED without chr prefix (matches SV chr naming)
loci_noprefix = f'{OUT}/_loci_mm39_noprefix.bed'
loci_mm39.assign(chr=loci_mm39['chr'].str.replace('^chr','',regex=True))[
    ['chr','start','end','base_gene']].to_csv(loci_noprefix, sep='\t', index=False, header=False)

# SVs in GRCm39
sv = pd.read_csv(f'{OUT}/C57BL_6NJ_TE_sized_SVs_GRCm39.csv')
sv_bed = f'{OUT}/_sv_mm39.bed'
sv.assign(chr=sv['chr'].astype(str))[['chr','start','end','id','sv_size_bp','sv_type']].to_csv(
    sv_bed, sep='\t', index=False, header=False)

# Loci in C57BL_6NJ — PICB overlap
loci_6nj = pd.read_csv(f'{OUT}/zamore_C57BL_6NJ_loci.bed', sep='\t', header=None,
                       names=['chr','start','end','name','strand'])
loci_6nj['stage'] = loci_6nj['name'].map(stage_by_locus)
loci_6nj_sorted = f'{OUT}/_loci_6nj_sorted2.bed'
loci_6nj.sort_values(['chr','start'])[['chr','start','end','name']].to_csv(
    loci_6nj_sorted, sep='\t', index=False, header=False)

def picb_hits(bed_b):
    r = subprocess.run(['bedtools','intersect','-a',loci_6nj_sorted,'-b',bed_b,'-wa','-u'],
                       capture_output=True, text=True, check=True)
    return set(l.split('\t')[3] for l in r.stdout.strip().split('\n') if l)

hits_p12 = picb_hits(f'{OUT}/C57BL_6NJ_P12.5_merged.bed')
hits_p20 = picb_hits(f'{OUT}/C57BL_6NJ_P20.5_merged.bed')
print(f'P12.5 hits: {len(hits_p12)},  P20.5 hits: {len(hits_p20)}')

def picb_cls(name):
    p12 = name in hits_p12; p20 = name in hits_p20
    if p12 and p20:  return 'both'
    elif p12:        return 'P12.5_only'
    elif p20:        return 'P20.5_only'
    else:            return 'none'

loci_6nj['picb_class'] = loci_6nj['name'].apply(picb_cls)
lifted_names = set(loci_6nj['name'])
loci_mm39['lifted']     = loci_mm39['base_gene'].isin(lifted_names)
loci_mm39['picb_class'] = loci_mm39['base_gene'].map(
    loci_6nj.set_index('name')['picb_class']).fillna('not_lifted')

# SV intersection at three windows
def sv_window(window):
    if window == 0:
        r = subprocess.run(['bedtools','intersect','-a',loci_noprefix,'-b',sv_bed,'-wo'],
                           capture_output=True, text=True)
    else:
        r = subprocess.run(['bedtools','window','-a',loci_noprefix,'-b',sv_bed,'-w',str(window)],
                           capture_output=True, text=True)
    hits = {}
    for line in r.stdout.strip().split('\n'):
        if not line: continue
        c = line.split('\t')
        gene=c[3]; sv_type=c[9]; sv_size=int(c[8])
        if gene not in hits: hits[gene]={'INS_bp':0,'DEL_bp':0}
        hits[gene][f'{sv_type}_bp'] += sv_size
    return hits

for label, win in [('direct',0),('10kb',10000),('50kb',50000)]:
    h = sv_window(win)
    loci_mm39[f'INS_{label}'] = loci_mm39['base_gene'].map(lambda g,h=h: h.get(g,{}).get('INS_bp',0))
    loci_mm39[f'DEL_{label}'] = loci_mm39['base_gene'].map(lambda g,h=h: h.get(g,{}).get('DEL_bp',0))
    loci_mm39[f'SV_{label}']  = (loci_mm39[f'INS_{label}']+loci_mm39[f'DEL_{label}'])>0

loci_mm39.to_csv(f'{OUT}/zamore_SV_PICB_analysis.csv', index=False)
print(f'\nSaved: {OUT}/zamore_SV_PICB_analysis.csv')

print('\n=== Direct SV × PICB class ===')
print(loci_mm39.groupby(['SV_direct','picb_class']).size().unstack(fill_value=0))
print('\n=== 10kb window × PICB class ===')
print(loci_mm39.groupby(['SV_10kb','picb_class']).size().unstack(fill_value=0))
print('\n=== 50kb window × PICB class ===')
print(loci_mm39.groupby(['SV_50kb','picb_class']).size().unstack(fill_value=0))
print('\n=== picb_class totals ===')
print(loci_mm39['picb_class'].value_counts())
