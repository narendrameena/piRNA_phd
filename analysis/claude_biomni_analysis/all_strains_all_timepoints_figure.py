#!/usr/bin/env python3
"""
Figure: all strains × all timepoints — SV effect on piRNA locus expression.

Panels A/B/C:  Heatmap for each timepoint (E16.5, P12.5, P20.5)
               Rows = strains, Cols = all 214 Zamore loci
               Loci sorted by P20.5 conservation; stage bands annotated above
               SV-frequency bar annotated above top panel only
Panel D:       SV vs no-SV disruption rate — all stage × timepoint combos
               M. musculus domesticus only (12 strains; SPRET/CAST/PWK/WSB excluded)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy.stats import fisher_exact

OUT = ("/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
       "/analysis/claude_biomni_analysis/all_strains_pangenome")

# ── colours ──────────────────────────────────────────────────────────────────
C_EXPR   = "#0072B2"   # blue    — expressed
C_NEXPR  = "#E69F00"   # amber   — lifted, no PICB
C_NLIFT  = "#D55E00"   # orange  — not lifted (genomic rearrangement)
C_NA     = "#dddddd"   # grey    — no data
STAGE_C  = {"Prepachytene": "#56B4E9", "Pachytene": "#009E73", "Hybrid": "#CC79A7"}

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 8,
    "axes.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "pdf.fonttype": 42,
})

STRAIN_ORDER = [
    "C57BL_6NJ","DBA_2J","BALB_cJ","A_J","CBA_J","129S1_SvImJ",
    "NOD_ShiLtJ","AKR_J","C3H_HeJ","NZO_HlLtJ","LP_J","FVB_NJ",
    "WSB_EiJ","PWK_PhJ","CAST_EiJ","SPRET_EiJ",
]
STRAIN_SHORT = {
    "C57BL_6NJ":"C57BL/6NJ","DBA_2J":"DBA/2J","BALB_cJ":"BALB/cJ",
    "A_J":"A/J","CBA_J":"CBA/J","129S1_SvImJ":"129S1",
    "NOD_ShiLtJ":"NOD","AKR_J":"AKR/J","C3H_HeJ":"C3H","NZO_HlLtJ":"NZO",
    "LP_J":"LP/J","FVB_NJ":"FVB/NJ","WSB_EiJ":"WSB/EiJ",
    "PWK_PhJ":"PWK/PhJ","CAST_EiJ":"CAST/EiJ","SPRET_EiJ":"SPRET/EiJ",
}
# domesticus-only: exclude wild-derived and M. spretus
MMD = ["C57BL_6NJ","DBA_2J","BALB_cJ","A_J","CBA_J","129S1_SvImJ",
       "NOD_ShiLtJ","AKR_J","C3H_HeJ","NZO_HlLtJ","LP_J","FVB_NJ"]
TIMEPOINTS = ["E16.5","P12.5","P20.5"]
TP_LABEL   = {"E16.5":"E16.5 (embryonic)", "P12.5":"P12.5 (postnatal)",
              "P20.5":"P20.5 (postnatal)"}
STAGES     = ["Prepachytene","Pachytene","Hybrid"]

# ── Load ──────────────────────────────────────────────────────────────────────
expr_df = pd.read_csv(f"{OUT}/all_strains_expression_matrix.csv")
sv_df   = pd.read_csv(f"{OUT}/all_strains_SV_matrix.csv")
sv_direct = sv_df[sv_df['window'] == 'direct'][['locus','strain','has_SV']].copy()
expr_df = expr_df.merge(sv_direct, on=['locus','strain'], how='left')
expr_df['has_SV']     = expr_df['has_SV'].fillna(False)
expr_df['disrupted']  = expr_df['status'].isin(['not_expressed','not_lifted'])

# ── Consistent locus ordering (by P20.5 conservation, all stages) ─────────────
all_loci_stages = (expr_df[['locus','stage']]
                   .drop_duplicates()
                   .set_index('locus')['stage'].to_dict())
p20_cons = (expr_df[expr_df['timepoint']=='P20.5']
            .groupby('locus')
            .apply(lambda x: (x['status']=='expressed').mean(), include_groups=False)
            .sort_values(ascending=False))
# Sort by stage then by P20.5 conservation within stage
stage_order = {"Prepachytene": 0, "Pachytene": 1, "Hybrid": 2}
locus_sort = pd.DataFrame({'locus': p20_cons.index,
                            'p20_cons': p20_cons.values})
locus_sort['stage'] = locus_sort['locus'].map(all_loci_stages)
locus_sort['stage_order'] = locus_sort['stage'].map(stage_order)
locus_sort = locus_sort.sort_values(['stage_order','p20_cons'],
                                     ascending=[True, False])
ALL_LOCI   = locus_sort['locus'].tolist()
n_loci     = len(ALL_LOCI)
locus_idx  = {l: i for i, l in enumerate(ALL_LOCI)}

# Stage boundary indices for colour bands
stage_bounds = {}
for stage in STAGES:
    idxs = [i for i, l in enumerate(ALL_LOCI)
            if all_loci_stages.get(l) == stage]
    if idxs:
        stage_bounds[stage] = (min(idxs), max(idxs))

# SV frequency per locus (across all strains, direct overlap)
sv_freq = sv_direct.groupby('locus')['has_SV'].mean().reindex(ALL_LOCI, fill_value=0)

# ── Build heatmap matrices ────────────────────────────────────────────────────
STATUS_NUM = {'expressed': 0, 'not_expressed': 1, 'not_lifted': 2}
matrices = {}
for tp in TIMEPOINTS:
    mat = np.full((len(STRAIN_ORDER), n_loci), 3, dtype=float)
    sub = expr_df[expr_df['timepoint'] == tp]
    for _, row in sub.iterrows():
        if row['strain'] not in STRAIN_ORDER:
            continue
        si = STRAIN_ORDER.index(row['strain'])
        li = locus_idx.get(row['locus'], -1)
        if li >= 0:
            mat[si, li] = STATUS_NUM.get(row['status'], 3)
    matrices[tp] = mat

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 14))
gs_main = GridSpec(4, 1, figure=fig, height_ratios=[1, 1, 1, 1.2],
                   hspace=0.55)

cmap = matplotlib.colors.ListedColormap([C_EXPR, C_NEXPR, C_NLIFT, C_NA])
norm = matplotlib.colors.BoundaryNorm([0, 1, 2, 3, 4], cmap.N)

axes_hm = []
for panel_idx, tp in enumerate(TIMEPOINTS):
    ax = fig.add_subplot(gs_main[panel_idx])
    axes_hm.append(ax)
    mat = matrices[tp]
    ax.imshow(mat, aspect='auto', cmap=cmap, norm=norm,
              interpolation='nearest')

    # Stage colour bands below x-axis
    for stage, (s, e) in stage_bounds.items():
        ax.add_patch(mpatches.FancyBboxPatch(
            (s - 0.5, len(STRAIN_ORDER) - 0.5),
            e - s + 1, 0.7,
            boxstyle="square,pad=0",
            facecolor=STAGE_C[stage], alpha=0.85,
            clip_on=False, zorder=4))
    # Stage labels
    for stage, (s, e) in stage_bounds.items():
        ax.text((s + e) / 2, len(STRAIN_ORDER) + 0.3, stage[:5],
                ha='center', va='bottom', fontsize=5.5,
                color=STAGE_C[stage], fontweight='bold', clip_on=False)

    # SV bar above top panel only
    if panel_idx == 0:
        ax_sv = ax.inset_axes([0, 1.12, 1, 0.06])
        ax_sv.bar(range(n_loci), sv_freq.values, width=1.0,
                  color="#CC79A7", linewidth=0)
        ax_sv.set_xlim(-0.5, n_loci - 0.5)
        ax_sv.set_ylim(0, 1)
        ax_sv.axis('off')
        ax_sv.text(-2, 0.5, "SV\nfreq.", fontsize=6, ha='right',
                   va='center', color="#CC79A7", clip_on=False)

    ax.set_yticks(range(len(STRAIN_ORDER)))
    ax.set_yticklabels([STRAIN_SHORT[s] for s in STRAIN_ORDER], fontsize=6.5)
    ax.set_xticks([])
    panel_letter = chr(ord('A') + panel_idx)
    ax.set_title(f"{panel_letter}   {TP_LABEL[tp]}  (all {n_loci} loci × 16 strains)",
                 fontsize=8.5, fontweight='bold', loc='left')

# Shared legend
legend_patches = [
    mpatches.Patch(color=C_EXPR,  label='Expressed'),
    mpatches.Patch(color=C_NEXPR, label='Not expressed (lifted)'),
    mpatches.Patch(color=C_NLIFT, label='Not lifted (rearranged)'),
    mpatches.Patch(color=C_NA,    label='No data'),
]
axes_hm[2].legend(handles=legend_patches, fontsize=7, frameon=False,
                  loc='lower right', bbox_to_anchor=(1.0, -0.28), ncol=4)

# ══════════════════════════════════════════════════════════════════════════════
# Panel D — SV vs no-SV disruption by stage × timepoint (M.m. domesticus only)
# ══════════════════════════════════════════════════════════════════════════════
ax_d = fig.add_subplot(gs_main[3])

mmd_e = expr_df[expr_df['strain'].isin(MMD)].copy()

# Build bar data: for each stage × timepoint, % disrupted with and without SV
groups = []
for tp in TIMEPOINTS:
    for stage in STAGES:
        sub = mmd_e[(mmd_e['timepoint'] == tp) & (mmd_e['stage'] == stage)]
        for sv_flag, sv_label in [(True, 'SV'), (False, 'no SV')]:
            s = sub[sub['has_SV'] == sv_flag]
            n = len(s)
            d = s['disrupted'].mean() * 100 if n > 0 else 0
            # Fisher p vs complement
            a = s['disrupted'].sum()
            b = n - a
            c_ = sub[sub['has_SV'] != sv_flag]['disrupted'].sum()
            d2 = len(sub[sub['has_SV'] != sv_flag]) - c_
            try:
                _, p = fisher_exact([[a, b], [c_, d2]])
            except Exception:
                p = float('nan')
            groups.append({'tp': tp, 'stage': stage, 'sv': sv_label,
                           'pct': d, 'n': n, 'p': p})

gdf = pd.DataFrame(groups)

# Grouped bars: x = stage × timepoint (9 groups), two bars each (SV / no-SV)
n_groups = len(STAGES) * len(TIMEPOINTS)
x_pos = np.arange(n_groups)
w = 0.35
labels = []
for stage in STAGES:
    for tp in TIMEPOINTS:
        labels.append(f"{tp[:3]}\n{stage[:5]}")

for sv_flag, sv_label, col, offset in [
        ('SV',    'Has SV',   '#CC79A7', -w/2),
        ('no SV', 'No SV',    '#888888', +w/2)]:
    pcts = []
    for stage in STAGES:
        for tp in TIMEPOINTS:
            row = gdf[(gdf['stage']==stage) & (gdf['tp']==tp) &
                      (gdf['sv']==sv_flag)]
            pcts.append(row.iloc[0]['pct'] if len(row) else 0)
    ax_d.bar(x_pos + offset, pcts, width=w, color=col,
             label=sv_label, alpha=0.85, edgecolor='white', lw=0.3)

# Significance markers
for i, (stage, tp) in enumerate([(s, t) for s in STAGES for t in TIMEPOINTS]):
    row = gdf[(gdf['stage']==stage) & (gdf['tp']==tp) & (gdf['sv']=='SV')]
    if len(row) and row.iloc[0]['p'] < 0.001:
        pct_sv = row.iloc[0]['pct']
        ax_d.text(i, pct_sv + 2, '***', ha='center', va='bottom', fontsize=7)
    elif len(row) and row.iloc[0]['p'] < 0.05:
        ax_d.text(i, row.iloc[0]['pct'] + 2, '*', ha='center',
                  va='bottom', fontsize=7)

# Stage background shading
stage_ranges = [(0, 2, STAGES[0]), (3, 5, STAGES[1]), (6, 8, STAGES[2])]
for s, e, stage in stage_ranges:
    ax_d.axvspan(s - 0.5, e + 0.5, alpha=0.07,
                 color=STAGE_C[stage], zorder=0)
    ax_d.text((s + e) / 2, -10, stage, ha='center', va='top',
              fontsize=7, color=STAGE_C[stage], fontweight='bold',
              clip_on=False)

ax_d.set_xticks(x_pos)
ax_d.set_xticklabels(labels, fontsize=6.5)
ax_d.set_ylabel("% loci disrupted\n(M. musculus domesticus, n=12 strains)", fontsize=8)
ax_d.set_ylim(0, 95)
ax_d.set_title("D   SV vs no-SV disruption by stage × timepoint  (*** = p<0.001)",
               fontsize=8.5, fontweight='bold', loc='left')
ax_d.legend(fontsize=7.5, frameon=False, loc='upper right')

fig.suptitle(
    "Pangenome SVs and piRNA cluster expression — all 16 strains, all timepoints\n"
    "214 Zamore loci · chain-file SVs ≥300 bp · PICB 2-of-3 consensus",
    fontsize=9.5, y=1.01)

for ext in ('pdf', 'svg', 'png'):
    fig.savefig(f"{OUT}/Fig_allstrains_all_timepoints.{ext}",
                dpi=300, bbox_inches='tight')
plt.close(fig)
print("Saved Fig_allstrains_all_timepoints.{pdf,svg,png}")

# ── Print summary table ───────────────────────────────────────────────────────
print("\n=== Disruption rate: SV vs no-SV (M.m. domesticus) ===")
pivot = gdf.pivot_table(index=['stage','tp'], columns='sv', values=['pct','n'],
                        aggfunc='first')
print(pivot.to_string())
