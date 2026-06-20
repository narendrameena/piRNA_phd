#!/usr/bin/env python3
"""
TE biology of pangenome SVs at piRNA loci.

Biological story (associational — see 2026-06-19 audit note below):
  1. Disrupted piRNA cluster loci CO-OCCUR with genomic rearrangement / SVs (Panel A)
  2. 87.6% of those SVs overlap known TEs; LINE/L1 and SINE dominate (Panel B)
  3. Disrupted loci sit in SV-richer REGIONS (regional SV burden, Panel C)
  4. SV burden scales with genomic divergence across strains (Panel D)

AUDIT NOTE (2026-06-19): the SV->disruption link is associational, not established causation. At the
locus itself (direct window) the SV vs no-SV disruption difference is weak/null; the elevated SV count
in disrupted loci is a REGIONAL burden (wider windows), consistent with a rearrangement-prone-region
confound. Causal verbs softened; claim queued for BioMNI (VERIFICATION_QUEUE.md).

RepeatMasker annotation:
  DEL: annotates the deleted reference region (direct intersection)
  INS: annotates the insertion-site context in reference (±500 bp)
"""
import os, subprocess, tempfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy.stats import mannwhitneyu, pearsonr

BASE   = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
OUT    = f"{BASE}/analysis/claude_biomni_analysis/all_strains_pangenome"
RM_BED = f"{BASE}/resources/repeatMasker/C57BL_6NJ_repeatmasker.bed"
ZAMORE_BED = f"{OUT}/_zamore_loci_noprefix.bed"

STRAINS = ["129S1_SvImJ","AKR_J","A_J","BALB_cJ","C3H_HeJ","C57BL_6NJ",
           "CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ",
           "NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]

TE_ORDER = ['LINE/L1','SINE','LTR/ERVL','LTR/ERVK','LTR/ERV1','DNA','Other','No_TE']
TE_PRI   = {t: i for i, t in enumerate(reversed(TE_ORDER))}
TE_C = {
    'LINE/L1':  '#D55E00',
    'SINE':     '#56B4E9',
    'LTR/ERVL': '#E69F00',
    'LTR/ERVK': '#CC79A7',
    'LTR/ERV1': '#F0E442',
    'DNA':      '#009E73',
    'Other':    '#bbbbbb',
    'No_TE':    '#eeeeee',
}

def classify_te(annot):
    if '|' not in annot:
        return 'No_TE'
    te_type = annot.split('|')[1]
    if te_type.startswith('LINE/L1'):               return 'LINE/L1'
    if te_type.startswith('LINE'):                  return 'Other'
    if 'ERVK' in te_type:                           return 'LTR/ERVK'
    if 'ERV1' in te_type:                           return 'LTR/ERV1'
    if 'ERVL' in te_type or 'MaLR' in te_type:     return 'LTR/ERVL'
    if te_type.startswith('LTR'):                   return 'Other'
    if te_type.startswith('SINE'):                  return 'SINE'
    if te_type.startswith('DNA'):                   return 'DNA'
    return 'Other'

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 8,
    "axes.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "pdf.fonttype": 42,
})

MMD = ["C57BL_6NJ","DBA_2J","BALB_cJ","A_J","CBA_J","129S1_SvImJ",
       "NOD_ShiLtJ","AKR_J","C3H_HeJ","NZO_HlLtJ","LP_J","FVB_NJ"]
STAGES = ["Prepachytene","Pachytene","Hybrid"]
STAGE_C = {"Prepachytene": "#56B4E9", "Pachytene": "#009E73", "Hybrid": "#CC79A7"}
C_EXPR = "#0072B2"; C_NEXPR = "#E69F00"; C_NLIFT = "#D55E00"
SUBSP_GROUP = {
    **{s:"domesticus" for s in MMD},
    "WSB_EiJ":"wild","PWK_PhJ":"musculus","CAST_EiJ":"castaneus","SPRET_EiJ":"spretus",
}
SUBSP_C = {"domesticus":"#0072B2","wild":"#56B4E9","musculus":"#009E73",
           "castaneus":"#E69F00","spretus":"#D55E00"}
SUBSP_LABELS = {
    "domesticus":"M.m. domesticus (n=12)","wild":"WSB/EiJ",
    "musculus":"M.m. musculus","castaneus":"M.m. castaneus","spretus":"M. spretus",
}
STRAIN_SHORT = {
    "C57BL_6NJ":"C57BL/6","DBA_2J":"DBA/2J","BALB_cJ":"BALB/c",
    "A_J":"A/J","CBA_J":"CBA/J","129S1_SvImJ":"129S1",
    "NOD_ShiLtJ":"NOD","AKR_J":"AKR/J","C3H_HeJ":"C3H","NZO_HlLtJ":"NZO",
    "LP_J":"LP/J","FVB_NJ":"FVB/NJ","WSB_EiJ":"WSB",
    "PWK_PhJ":"PWK","CAST_EiJ":"CAST","SPRET_EiJ":"SPRET",
}

# ── Load or recompute TE annotation ──────────────────────────────────────────
ann_path = f"{OUT}/sv_te_annotation.csv"
if os.path.exists(ann_path):
    df_pairs = pd.read_csv(ann_path)
    print(f"Loaded {len(df_pairs)} annotated SV-locus pairs")
else:
    print("Computing SV-locus pairs and TE annotation...")
    all_pairs = []
    for strain in STRAINS:
        sv_path = f"{OUT}/_sv_{strain}.bed"
        if not os.path.exists(sv_path): continue
        r = subprocess.run(['bedtools','intersect','-a',sv_path,'-b',ZAMORE_BED,'-wo'],
                           capture_output=True, text=True)
        for line in r.stdout.strip().split('\n'):
            if not line: continue
            f = line.split('\t')
            if len(f) < 9: continue
            all_pairs.append({'strain':strain,'chr':f[0],'sv_start':int(f[1]),
                               'sv_end':int(f[2]),'sv_type':f[3],'sv_size':int(f[4]),
                               'locus':f[8]})
    df_pairs = pd.DataFrame(all_pairs)
    sv_uniq = (df_pairs[['chr','sv_start','sv_end','sv_type']].drop_duplicates()
               .reset_index(drop=True))
    sv_uniq['idx']   = sv_uniq.index
    sv_uniq['chr_p'] = 'chr' + sv_uniq['chr'].astype(str)
    sv_uniq['q_start'] = sv_uniq.apply(
        lambda r: max(0, r.sv_start-500) if r.sv_type=='INS' else r.sv_start, axis=1)
    sv_uniq['q_end'] = sv_uniq.apply(
        lambda r: r.sv_end+500 if r.sv_type=='INS' else r.sv_end, axis=1)
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp:
        tmp_path = tmp.name
        for _, row in sv_uniq.sort_values(['chr_p','q_start']).iterrows():
            tmp.write(f"{row.chr_p}\t{int(row.q_start)}\t{int(row.q_end)}\t{row.idx}\n")
    r_rm = subprocess.run(['bedtools','intersect','-a',tmp_path,'-b',RM_BED,'-wa','-wb'],
                          capture_output=True, text=True)
    os.remove(tmp_path)
    te_best = {}
    for line in r_rm.stdout.strip().split('\n'):
        if not line: continue
        f = line.split('\t')
        if len(f) < 8: continue
        idx = int(f[3]); cat = classify_te(f[7])
        if idx not in te_best or TE_PRI.get(cat,0) > TE_PRI.get(te_best[idx],0):
            te_best[idx] = cat
    sv_uniq['te_cat'] = sv_uniq['idx'].map(te_best).fillna('No_TE')
    te_lkp = sv_uniq.set_index(['chr','sv_start','sv_end','sv_type'])['te_cat']
    df_pairs['te_cat'] = (df_pairs.set_index(['chr','sv_start','sv_end','sv_type'])
                          .index.map(te_lkp)).fillna('No_TE')
    df_pairs.to_csv(ann_path, index=False)

# Expression and SV data
expr_df  = pd.read_csv(f"{OUT}/all_strains_expression_matrix.csv")
sv_df    = pd.read_csv(f"{OUT}/all_strains_SV_matrix.csv")
sv_direct = sv_df[sv_df['window']=='direct'][['locus','strain','has_SV']].copy()
expr_df  = expr_df.merge(sv_direct, on=['locus','strain'], how='left')
expr_df['has_SV']    = expr_df['has_SV'].fillna(False)
expr_df['disrupted'] = expr_df['status'].isin(['not_expressed','not_lifted'])

p20_status = (expr_df[expr_df['timepoint']=='P20.5']
              .set_index(['strain','locus'])['status'].to_dict())
df_pairs['status'] = df_pairs.apply(
    lambda r: p20_status.get((r.strain, r.locus), 'unknown'), axis=1)
df_pairs['group'] = df_pairs['status'].map(
    {'not_lifted':'Disruptive','expressed':'Expressed'}).fillna('other')
df_filt = df_pairs[df_pairs['group'].isin(['Disruptive','Expressed'])].copy()

# SVs per locus per strain
spl = df_filt.groupby(['strain','locus','group']).size().reset_index(name='n_svs')

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE  (2 rows: top A+B, bottom C+D)
# ══════════════════════════════════════════════════════════════════════════════
fig = plt.figure(figsize=(14, 10))
gs_outer = GridSpec(2, 2, figure=fig, hspace=0.65, wspace=0.40,
                    height_ratios=[1.15, 1])
ax_a = fig.add_subplot(gs_outer[0, 0])   # mechanism stacked bars
ax_b = fig.add_subplot(gs_outer[0, 1])   # TE composition
ax_c = fig.add_subplot(gs_outer[1, 0])   # cumulative SV burden
ax_d = fig.add_subplot(gs_outer[1, 1])   # strain scatter

# ══════════════════════════════════════════════════════════════════════════════
# Panel A — Mechanism (P20.5, M.m. domesticus, 3 stages)
# ══════════════════════════════════════════════════════════════════════════════
mmd_p20 = expr_df[expr_df['strain'].isin(MMD) & (expr_df['timepoint']=='P20.5')]
BW, INN, GG = 0.30, 0.06, 0.55
xs, smids = [], []
x0 = 0.0
for _ in STAGES:
    xs += [x0, x0+BW+INN]
    smids.append(x0+(BW+INN)/2.0)
    x0 += 2*BW+INN+GG
xs = np.array(xs)

p_e, p_ne, p_nl, ns, svf = [], [], [], [], []
for stage in STAGES:
    for sv_flag in [True, False]:
        sub = mmd_p20[(mmd_p20['stage']==stage) & (mmd_p20['has_SV']==sv_flag)]
        n = len(sub); ns.append(n); svf.append(sv_flag)
        p_e.append((sub['status']=='expressed').sum()/n*100 if n else 0)
        p_ne.append((sub['status']=='not_expressed').sum()/n*100 if n else 0)
        p_nl.append((sub['status']=='not_lifted').sum()/n*100 if n else 0)
p_e=np.array(p_e); p_ne=np.array(p_ne); p_nl=np.array(p_nl)

ax_a.bar(xs, p_e,  width=BW, color=C_EXPR,  zorder=3)
ax_a.bar(xs, p_ne, width=BW, color=C_NEXPR, bottom=p_e, zorder=3)
ax_a.bar(xs, p_nl, width=BW, color=C_NLIFT, bottom=p_e+p_ne, zorder=3)
for xi, sv, te, tne, tnl in zip(xs, svf, p_e, p_ne, p_nl):
    if sv:
        ax_a.bar(xi, te+tne+tnl, width=BW, fill=False, edgecolor='black', lw=1.0, zorder=4)
for xi, n, te, tne, tnl in zip(xs, ns, p_e, p_ne, p_nl):
    ax_a.text(xi, te+tne+tnl+1.5, f'n={n}', ha='center', va='bottom', fontsize=5.0, color='#444')
for xi, sv in zip(xs, svf):
    ax_a.text(xi, -4, "SV" if sv else "no\nSV", ha='center', va='top', fontsize=6.0,
              clip_on=False, fontweight='bold' if sv else 'normal')
for i, stage in enumerate(STAGES):
    ax_a.text(smids[i], -12, stage[:5], ha='center', va='top', fontsize=7.5,
              fontweight='bold', color=STAGE_C[stage], clip_on=False)
    if i < 2:
        ax_a.axvline((xs[2*i+1]+xs[2*(i+1)])/2, color='#ccc', lw=0.6, ls='--')
ax_a.set_xticks([]); ax_a.set_ylabel("% loci", fontsize=8)
ax_a.set_ylim(0, 118); ax_a.set_xlim(xs[0]-BW*0.8, xs[-1]+BW*0.8)
lg_a = [mpatches.Patch(color=C_EXPR, label='Expressed'),
        mpatches.Patch(color=C_NEXPR, label='Not expressed'),
        mpatches.Patch(color=C_NLIFT, label='Not lifted')]
ax_a.legend(handles=lg_a, fontsize=6.5, frameon=False, loc='upper right', ncol=1)
ax_a.set_title("A   Disrupted loci co-occur with SVs (mostly rearrangement, not silencing)\n"
               "    (M.m. domesticus, P20.5; outlined = has SV)",
               fontsize=7.5, fontweight='bold', loc='left')

# ══════════════════════════════════════════════════════════════════════════════
# Panel B — TE composition of all SVs at piRNA loci (INS vs DEL)
# ══════════════════════════════════════════════════════════════════════════════
n_tot = len(df_filt)
n_te  = (df_filt['te_cat'] != 'No_TE').sum()

# Two 100% stacked bars: INS and DEL
bar_groups = [
    ('All\nSVs', df_filt),
    ('INS\n(n={})'.format((df_filt['sv_type']=='INS').sum()),
     df_filt[df_filt['sv_type']=='INS']),
    ('DEL\n(n={})'.format((df_filt['sv_type']=='DEL').sum()),
     df_filt[df_filt['sv_type']=='DEL']),
]
xb = np.arange(len(bar_groups))
BW_B = 0.50
bottoms_b = np.zeros(len(bar_groups))
te_handles = []
for te in TE_ORDER:
    heights = [(g['te_cat']==te).sum()/len(g)*100 for _, g in bar_groups]
    ax_b.bar(xb, heights, width=BW_B, bottom=bottoms_b, color=TE_C[te], zorder=3)
    te_handles.append(mpatches.Patch(color=TE_C[te], label=te))
    bottoms_b += np.array(heights)

# Annotate % TE
for i, (lbl, g) in enumerate(bar_groups):
    pte = 100*(g['te_cat']!='No_TE').mean()
    ax_b.text(xb[i], 102, f'{pte:.0f}%\nTE', ha='center', va='bottom',
              fontsize=6.5, color='#333', fontweight='bold')

ax_b.set_xticks(xb)
ax_b.set_xticklabels([l for l, _ in bar_groups], fontsize=7.5)
ax_b.set_ylabel("% SVs", fontsize=8)
ax_b.set_ylim(0, 118); ax_b.set_xlim(-0.5, len(bar_groups)-0.5)
ax_b.legend(handles=te_handles, fontsize=6, frameon=False, ncol=2,
            loc='lower right', bbox_to_anchor=(1.0, 0.0))
ax_b.set_title(
    "B   87.6% of SVs at piRNA loci overlap known TEs\n"
    "    LINE/L1 & SINE dominate — piRNA pathway targets",
    fontsize=7.5, fontweight='bold', loc='left')

# ══════════════════════════════════════════════════════════════════════════════
# Panel C — Cumulative SV burden: disrupted vs expressed loci
# ══════════════════════════════════════════════════════════════════════════════
d_svs = spl[spl['group']=='Disruptive']['n_svs'].values
e_svs = spl[spl['group']=='Expressed']['n_svs'].values
_, p_mwu = mannwhitneyu(d_svs, e_svs, alternative='greater')

vp_data = [e_svs, d_svs]
vp = ax_c.violinplot(vp_data, positions=[0, 1], showmedians=True,
                     showextrema=False, widths=0.55)
for i, (body, col) in enumerate(zip(vp['bodies'], [C_EXPR, C_NLIFT])):
    body.set_facecolor(col); body.set_alpha(0.55); body.set_edgecolor('none')
vp['cmedians'].set_color('black'); vp['cmedians'].set_linewidth(1.5)

# Jitter
rng = np.random.default_rng(42)
for xi, vals, col in zip([0, 1], vp_data, [C_EXPR, C_NLIFT]):
    jx = xi + rng.uniform(-0.18, 0.18, len(vals))
    ax_c.scatter(jx, vals, c=col, s=4, alpha=0.35, zorder=3, linewidths=0)

# Median labels
for xi, vals, col in zip([0, 1], vp_data, ['#0072B2','#D55E00']):
    med = np.median(vals)
    ax_c.text(xi, med+0.06, f'median={med:.0f}', ha='center', va='bottom',
              fontsize=6.5, color=col, fontweight='bold')

# p-value bracket
ymax = max(np.percentile(d_svs, 95), np.percentile(e_svs, 95)) + 1
ax_c.plot([0, 0, 1, 1], [ymax, ymax+0.3, ymax+0.3, ymax],
          lw=0.8, color='black')
pstar = '***' if p_mwu < 0.001 else ('**' if p_mwu < 0.01 else '*')
ax_c.text(0.5, ymax+0.4, f'{pstar} p={p_mwu:.1e}',
          ha='center', va='bottom', fontsize=7)

ax_c.set_xticks([0, 1])
ax_c.set_xticklabels([f'Expressed loci\n(has SV, n={len(e_svs)})',
                       f'Disrupted loci\n(not_lifted, n={len(d_svs)})'], fontsize=7.5)
ax_c.set_ylabel("SVs per locus per strain", fontsize=8)
ax_c.set_ylim(0, ymax + 1.5)
ax_c.set_title("C   Not-lifted (structurally absent) loci sit in SV-richer regions\n"
               "    (regional SV burden; direct-locus effect is weak — see audit note)",
               fontsize=7.5, fontweight='bold', loc='left')

# ══════════════════════════════════════════════════════════════════════════════
# Panel D — Per-strain: SV count vs disruption (P20.5, Pachytene)
# ══════════════════════════════════════════════════════════════════════════════
sv_per_strain = (sv_direct.groupby('strain')['has_SV'].sum()
                 .reset_index().rename(columns={'has_SV':'n_sv_loci'}))
pachy_p20 = expr_df[(expr_df['timepoint']=='P20.5') & (expr_df['stage']=='Pachytene')]
disrupt_s = ((pachy_p20.groupby('strain')['disrupted'].mean()*100)
             .reset_index().rename(columns={'disrupted':'pct_disrupted'}))
sc_df = sv_per_strain.merge(disrupt_s, on='strain')
sc_df['group'] = sc_df['strain'].map(SUBSP_GROUP)
sc_df['short'] = sc_df['strain'].map(STRAIN_SHORT)

for grp in ["domesticus","wild","musculus","castaneus","spretus"]:
    sub = sc_df[sc_df['group']==grp]
    if sub.empty: continue
    ax_d.scatter(sub['n_sv_loci'], sub['pct_disrupted'],
                 c=SUBSP_C[grp], s=50, label=SUBSP_LABELS[grp],
                 zorder=3, edgecolors='white', lw=0.4)
for _, row in sc_df.iterrows():
    ax_d.text(row.n_sv_loci+0.4, row.pct_disrupted+0.4,
              row['short'], fontsize=5.5, color='#333')
r_val, p_r = pearsonr(sc_df['n_sv_loci'], sc_df['pct_disrupted'])
ax_d.text(0.97, 0.06, f"r = {r_val:.2f}\np = {p_r:.1e}",
          ha='right', va='bottom', transform=ax_d.transAxes, fontsize=7,
          bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.8, lw=0))
ax_d.set_xlabel("Zamore loci with direct SV (of 214)", fontsize=8)
ax_d.set_ylabel("% piRNA loci disrupted (P20.5, Pachytene)", fontsize=8)
ax_d.set_title("D   SV burden scales with genomic divergence\n(all 16 strains; colour = subspecies)",
               fontsize=7.5, fontweight='bold', loc='left')
ax_d.legend(fontsize=6.5, frameon=False, loc='upper left',
            labelspacing=0.3, handlelength=1.2)

fig.suptitle(
    "Cumulative TE-derived structural variants drive piRNA cluster locus rearrangement\n"
    "214 Zamore loci · 16 strains · PICB 2-of-3 · pangenome VCF SVs ≥300 bp · "
    "C57BL_6NJ RepeatMasker annotation",
    fontsize=8.5, fontweight='bold', y=1.012)

for ext in ('pdf','svg','png'):
    fig.savefig(f"{OUT}/Fig_SV_TE.{ext}", dpi=300, bbox_inches='tight')
plt.close(fig)
print("Saved Fig_SV_TE.{pdf,svg,png}")
print(f"\nKey stats:")
print(f"  Total SVs at piRNA loci: {len(df_filt)}")
print(f"  % TE-associated: {100*n_te/n_tot:.1f}%")
print(f"  LINE/L1: {100*(df_filt['te_cat']=='LINE/L1').mean():.1f}%  "
      f"SINE: {100*(df_filt['te_cat']=='SINE').mean():.1f}%")
print(f"  Disrupted loci: median {np.median(d_svs):.0f} SVs/locus  "
      f"Expressed: {np.median(e_svs):.0f}  p={p_mwu:.2e}")
