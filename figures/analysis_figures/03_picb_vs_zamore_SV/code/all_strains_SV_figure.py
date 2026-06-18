#!/usr/bin/env python3
"""
Figure: Multi-strain SV × piRNA locus expression across 16 mouse strains.

Panels:
  A  Heatmap — P20.5 pachytene loci × strains (expressed/not_expressed/not_lifted/no_PICB_data)
     with SV annotation bar at right
  B  Bar — per-strain: % loci disrupted (not_expressed + not_lifted) at each timepoint
  C  Scatter — per-locus: SV frequency (% strains) vs disruption rate (% strains);
     colour by stage; labelled outliers
  D  Bar — for M. musculus domesticus only: SVs vs no-SVs, expression rate by stage × timepoint
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy.stats import fisher_exact

OUT = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/all_strains_pangenome"

# ── colours ──────────────────────────────────────────────────────────────────
C_EXPR   = "#0072B2"   # blue     — expressed
C_NEXPR  = "#E69F00"   # amber    — not expressed (lifted, no PICB)
C_NLIFT  = "#D55E00"   # vermil.  — not lifted (major rearrangement)
C_NA     = "#eeeeee"   # grey     — no data
C_SV     = "#CC79A7"   # pink     — has SV
C_NOSV   = "#bbbbbb"   # grey     — no SV

PLT_STAGES = ["Prepachytene", "Pachytene", "Hybrid"]
STAGE_C    = {"Prepachytene": "#56B4E9", "Pachytene": "#009E73", "Hybrid": "#E69F00"}

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 8,
    "axes.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "pdf.fonttype": 42,
})

# Strain ordering: by liftover success (divergence from GRCm39)
STRAIN_ORDER = [
    "C57BL_6NJ","DBA_2J","BALB_cJ","A_J","CBA_J","129S1_SvImJ",
    "NOD_ShiLtJ","AKR_J","C3H_HeJ","NZO_HlLtJ","LP_J",
    "FVB_NJ","WSB_EiJ","PWK_PhJ","CAST_EiJ","SPRET_EiJ",
]
STRAIN_SHORT = {
    "C57BL_6NJ":"C57BL/6NJ","DBA_2J":"DBA/2J","BALB_cJ":"BALB/cJ",
    "A_J":"A/J","CBA_J":"CBA/J","129S1_SvImJ":"129S1","NOD_ShiLtJ":"NOD",
    "AKR_J":"AKR/J","C3H_HeJ":"C3H","NZO_HlLtJ":"NZO","LP_J":"LP/J",
    "FVB_NJ":"FVB/NJ","WSB_EiJ":"WSB/EiJ","PWK_PhJ":"PWK/PhJ",
    "CAST_EiJ":"CAST/EiJ","SPRET_EiJ":"SPRET/EiJ",
}
# M. musculus domesticus-only strains (exclude wild-derived and M. spretus)
MM_DOMESTICUS = ["C57BL_6NJ","DBA_2J","BALB_cJ","A_J","CBA_J","129S1_SvImJ",
                 "NOD_ShiLtJ","AKR_J","C3H_HeJ","NZO_HlLtJ","LP_J","FVB_NJ"]

# ── Load data ─────────────────────────────────────────────────────────────────
expr_df = pd.read_csv(f"{OUT}/all_strains_expression_matrix.csv")
sv_df   = pd.read_csv(f"{OUT}/all_strains_SV_matrix.csv")

sv_direct = sv_df[sv_df['window'] == 'direct'][
    ['locus','strain','has_SV']].copy()

expr_df = expr_df.merge(sv_direct, on=['locus','strain'], how='left')
expr_df['has_SV'] = expr_df['has_SV'].fillna(False)

# Disrupted = not_expressed OR not_lifted
expr_df['disrupted'] = expr_df['status'].isin(['not_expressed','not_lifted'])
expr_df['expressed_bin'] = (expr_df['status'] == 'expressed').astype(int)

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 12))
gs  = GridSpec(2, 2, figure=fig, hspace=0.50, wspace=0.42)

# ══════════════════════════════════════════════════════════════════════════════
# Panel A — Heatmap: P20.5 × Pachytene loci × strains
# ══════════════════════════════════════════════════════════════════════════════
ax_a = fig.add_subplot(gs[0, :])   # span both columns

pachy_p20 = expr_df[
    (expr_df['timepoint'] == 'P20.5') &
    (expr_df['stage'] == 'Pachytene')].copy()

# Locus ordering: by conservation (% strains expressing, descending)
pachy_loci = pachy_p20.groupby('locus').apply(
    lambda x: x[x['status'] == 'expressed'].shape[0] / len(STRAIN_ORDER)
).sort_values(ascending=False).index.tolist()

# Build numeric matrix
STATUS_NUM = {'expressed': 0, 'not_expressed': 1, 'not_lifted': 2}
mat = np.full((len(pachy_loci), len(STRAIN_ORDER)), 3)   # 3 = no data
for i, locus in enumerate(pachy_loci):
    sub = pachy_p20[pachy_p20['locus'] == locus]
    for j, strain in enumerate(STRAIN_ORDER):
        row = sub[sub['strain'] == strain]
        if len(row):
            mat[i, j] = STATUS_NUM.get(row.iloc[0]['status'], 3)

# Draw heatmap
cmap = matplotlib.colors.ListedColormap([C_EXPR, C_NEXPR, C_NLIFT, C_NA])
norm = matplotlib.colors.BoundaryNorm([0, 1, 2, 3, 4], cmap.N)
im = ax_a.imshow(mat.T, aspect='auto', cmap=cmap, norm=norm,
                 interpolation='nearest')
ax_a.set_xticks([])
ax_a.set_yticks(range(len(STRAIN_ORDER)))
ax_a.set_yticklabels(
    [STRAIN_SHORT[s] for s in STRAIN_ORDER], fontsize=7)
ax_a.set_xlabel(
    f"Pachytene piRNA loci (n={len(pachy_loci)}, sorted by expression breadth)",
    fontsize=8)

# SV annotation bar above heatmap
sv_pachy = sv_direct[sv_direct['locus'].isin(pachy_loci)]
sv_per_locus = sv_pachy.groupby('locus')['has_SV'].sum()
sv_frac = np.array([sv_per_locus.get(l, 0) / len(STRAIN_ORDER)
                    for l in pachy_loci])
# Plot as a thin strip above the main heatmap
ax_sv = ax_a.inset_axes([0, 1.02, 1, 0.04])
ax_sv.bar(range(len(pachy_loci)), sv_frac, width=1.0,
          color=C_SV, linewidth=0)
ax_sv.set_xlim(-0.5, len(pachy_loci) - 0.5)
ax_sv.set_ylim(0, 1)
ax_sv.axis('off')
ax_sv.text(-1.5, 0.5, "% strains\nwith SV", fontsize=6,
           ha='right', va='center', color=C_SV)

legend_patches = [
    mpatches.Patch(color=C_EXPR,  label='Expressed'),
    mpatches.Patch(color=C_NEXPR, label='Not expressed (lifted)'),
    mpatches.Patch(color=C_NLIFT, label='Not lifted (rearranged)'),
    mpatches.Patch(color=C_NA,    label='No data'),
]
ax_a.legend(handles=legend_patches, fontsize=7, frameon=False,
            loc='lower right', bbox_to_anchor=(1.0, -0.18), ncol=4)
ax_a.set_title("A   Pachytene piRNA loci — P20.5 expression across 16 strains",
               fontsize=9, fontweight='bold', loc='left')

# ══════════════════════════════════════════════════════════════════════════════
# Panel B — Per-strain disruption rate at P20.5 (Pachytene loci)
# ══════════════════════════════════════════════════════════════════════════════
ax_b = fig.add_subplot(gs[1, 0])

disruption = (pachy_p20.groupby('strain')['disrupted'].mean() * 100).reindex(STRAIN_ORDER)
liftover_fail = (pachy_p20[pachy_p20['status'] == 'not_lifted']
                 .groupby('strain').size().reindex(STRAIN_ORDER, fill_value=0)
                 / pachy_p20.groupby('strain').size().reindex(STRAIN_ORDER) * 100)
nexpr = (pachy_p20[pachy_p20['status'] == 'not_expressed']
         .groupby('strain').size().reindex(STRAIN_ORDER, fill_value=0)
         / pachy_p20.groupby('strain').size().reindex(STRAIN_ORDER) * 100)

x = np.arange(len(STRAIN_ORDER))
ax_b.bar(x, liftover_fail.fillna(0), color=C_NLIFT, label='Not lifted')
ax_b.bar(x, nexpr.fillna(0),
         bottom=liftover_fail.fillna(0), color=C_NEXPR, label='Not expressed')
ax_b.set_xticks(x)
ax_b.set_xticklabels([STRAIN_SHORT[s] for s in STRAIN_ORDER],
                     rotation=40, ha='right', fontsize=6.5)
ax_b.set_ylabel("% loci disrupted", fontsize=8)
ax_b.set_ylim(0, 100)
ax_b.set_title("B   Disruption rate per strain\n(Pachytene loci, P20.5)",
               fontsize=8.5, fontweight='bold', loc='left')
ax_b.legend(fontsize=7, frameon=False, loc='upper left')

# Divergence label on right
ax_b.annotate('← M. musculus  |  M. spretus →',
              xy=(0.5, -0.38), xycoords='axes fraction',
              ha='center', fontsize=7, color='#555555')

# ══════════════════════════════════════════════════════════════════════════════
# Panel C — Scatter: SV frequency vs disruption rate per locus
# ══════════════════════════════════════════════════════════════════════════════
ax_c = fig.add_subplot(gs[1, 1])

# Compute per-locus stats (all strains, direct SVs, P20.5 only → combine stages)
locus_stats = []
for stage in PLT_STAGES:
    sub_e = expr_df[(expr_df['timepoint'] == 'P20.5') &
                    (expr_df['stage'] == stage)]
    for locus, grp in sub_e.groupby('locus'):
        n = len(grp)
        if n < 3:
            continue
        sv_freq   = grp['has_SV'].mean() * 100
        disrupt   = grp['disrupted'].mean() * 100
        locus_stats.append({'locus': locus, 'stage': stage,
                             'sv_freq': sv_freq, 'disrupt': disrupt, 'n': n})

ls_df = pd.DataFrame(locus_stats)

for stage in PLT_STAGES:
    sub = ls_df[ls_df['stage'] == stage]
    ax_c.scatter(sub['sv_freq'], sub['disrupt'],
                 c=STAGE_C[stage], s=20, alpha=0.7,
                 label=stage, edgecolors='none')

# Label top disrupted loci with manual y-offsets to avoid overlap
top = ls_df.nlargest(4, 'disrupt').reset_index(drop=True)
y_shifts = [2, -4, 2, -4]
for i, (_, r) in enumerate(top.iterrows()):
    ax_c.text(r['sv_freq'] + 1, r['disrupt'] + y_shifts[i],
              r['locus'], fontsize=5.5, color='#333333',
              bbox=dict(boxstyle='round,pad=0.1', fc='white', alpha=0.6, lw=0))

ax_c.set_xlabel("% strains with direct SV", fontsize=8)
ax_c.set_ylabel("% strains disrupted (P20.5)", fontsize=8)
ax_c.set_title("C   SV frequency vs disruption rate per locus",
               fontsize=8.5, fontweight='bold', loc='left')
ax_c.legend(fontsize=7, frameon=False)

# Fisher exact for M. musculus domesticus
mmd_e = expr_df[expr_df['strain'].isin(MM_DOMESTICUS) &
                (expr_df['timepoint'] == 'P20.5')]
a = mmd_e[mmd_e['has_SV']]['disrupted'].sum()
b = mmd_e[mmd_e['has_SV']].shape[0] - a
c = mmd_e[~mmd_e['has_SV']]['disrupted'].sum()
d = mmd_e[~mmd_e['has_SV']].shape[0] - c
_, p = fisher_exact([[a, b], [c, d]])
ax_c.text(0.98, 0.02, f"Fisher p={p:.2e}\n(M.m. domesticus)",
          ha='right', va='bottom', transform=ax_c.transAxes,
          fontsize=7, color='#555555')

fig.suptitle(
    "Pangenome SVs and piRNA cluster expression across 16 inbred mouse strains\n"
    "GRCm39 reference  ·  Zamore piRNA annotation  ·  PICB cluster calls",
    fontsize=9.5, y=1.01)

for ext in ('pdf', 'svg', 'png'):
    fig.savefig(f"{OUT}/Fig_allstrains_SV_expression.{ext}",
                dpi=300, bbox_inches='tight')
plt.close(fig)
print("Saved Fig_allstrains_SV_expression.{pdf,svg,png}")

# ── Print key stats ───────────────────────────────────────────────────────────
print("\n=== Summary stats ===")
for tp in ['E16.5','P12.5','P20.5']:
    for stage in PLT_STAGES:
        sub = expr_df[(expr_df['timepoint']==tp) & (expr_df['stage']==stage) &
                      (expr_df['strain'].isin(MM_DOMESTICUS))]
        n_sv  = sub[sub['has_SV']].shape[0]
        n_nosv = sub[~sub['has_SV']].shape[0]
        if n_sv + n_nosv == 0:
            continue
        d_sv   = sub[sub['has_SV']]['disrupted'].mean()*100 if n_sv else 0
        d_nosv = sub[~sub['has_SV']]['disrupted'].mean()*100 if n_nosv else 0
        a = sub[sub['has_SV']]['disrupted'].sum()
        b = n_sv - a
        c_ = sub[~sub['has_SV']]['disrupted'].sum()
        d2 = n_nosv - c_
        try:
            _, p = fisher_exact([[a,b],[c_,d2]])
        except Exception:
            p = float('nan')
        print(f"  {tp} {stage:15s}: SV={n_sv} disrupted={d_sv:.1f}%  "
              f"noSV={n_nosv} disrupted={d_nosv:.1f}%  p={p:.3f}")
