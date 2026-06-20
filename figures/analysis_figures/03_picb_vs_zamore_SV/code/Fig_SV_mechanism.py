#!/usr/bin/env python3
"""
Pangenome SVs CO-OCCUR with piRNA cluster disruption (association, not established causation).
NOTE (2026-06-19 data-integrity audit): at the locus itself (direct window, used here) the SV vs
no-SV disruption difference is weak/null (Pachytene ~24.9% vs ~24.8%); a positive association appears
only at wider windows (10-50 kb), consistent with a regional rearrangement confound rather than a
direct causal effect. Causal framing softened; SV->disruption claim queued for BioMNI (VERIFICATION_QUEUE.md).

Panel A:  Stacked proportion bars — expressed/not_expressed/not_lifted,
          SV vs no-SV loci, by stage (M.m. domesticus, P20.5)
Panel B:  SV vs no-SV disruption scatter — all 9 stage × timepoint groups
Panel C:  Per-strain SV count vs disruption rate (all 16 strains, P20.5 Pachytene)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy.stats import pearsonr

OUT = ("/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
       "/analysis/claude_biomni_analysis/all_strains_pangenome")

C_EXPR  = "#0072B2"
C_NEXPR = "#E69F00"
C_NLIFT = "#D55E00"
STAGE_C = {"Prepachytene": "#56B4E9", "Pachytene": "#009E73", "Hybrid": "#CC79A7"}
MARKER  = {"E16.5": "o", "P12.5": "s", "P20.5": "^"}

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 8,
    "axes.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "pdf.fonttype": 42,
})

MMD = ["C57BL_6NJ","DBA_2J","BALB_cJ","A_J","CBA_J","129S1_SvImJ",
       "NOD_ShiLtJ","AKR_J","C3H_HeJ","NZO_HlLtJ","LP_J","FVB_NJ"]
STAGES     = ["Prepachytene","Pachytene","Hybrid"]
TIMEPOINTS = ["E16.5","P12.5","P20.5"]
STRAIN_SHORT = {
    "C57BL_6NJ":"C57BL/6","DBA_2J":"DBA/2J","BALB_cJ":"BALB/c",
    "A_J":"A/J","CBA_J":"CBA/J","129S1_SvImJ":"129S1",
    "NOD_ShiLtJ":"NOD","AKR_J":"AKR/J","C3H_HeJ":"C3H","NZO_HlLtJ":"NZO",
    "LP_J":"LP/J","FVB_NJ":"FVB/NJ","WSB_EiJ":"WSB",
    "PWK_PhJ":"PWK","CAST_EiJ":"CAST","SPRET_EiJ":"SPRET",
}
SUBSP_GROUP = {
    **{s: "domesticus" for s in MMD},
    "WSB_EiJ": "wild", "PWK_PhJ": "musculus",
    "CAST_EiJ": "castaneus", "SPRET_EiJ": "spretus",
}
SUBSP_C = {
    "domesticus": "#0072B2", "wild": "#56B4E9",
    "musculus": "#009E73", "castaneus": "#E69F00", "spretus": "#D55E00",
}
SUBSP_LABELS = {
    "domesticus": "M.m. domesticus (n=12)",
    "wild":       "WSB/EiJ (wild dom.)",
    "musculus":   "M.m. musculus (PWK)",
    "castaneus":  "M.m. castaneus (CAST)",
    "spretus":    "M. spretus (SPRET)",
}

# ── Load ──────────────────────────────────────────────────────────────────────
expr_df = pd.read_csv(f"{OUT}/all_strains_expression_matrix.csv")
sv_df   = pd.read_csv(f"{OUT}/all_strains_SV_matrix.csv")
sv_direct = sv_df[sv_df['window'] == 'direct'][['locus','strain','has_SV']].copy()
expr_df = expr_df.merge(sv_direct, on=['locus','strain'], how='left')
expr_df['has_SV']    = expr_df['has_SV'].fillna(False)
expr_df['disrupted'] = expr_df['status'].isin(['not_expressed','not_lifted'])

# ── Layout ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(13, 9))
gs  = GridSpec(2, 2, figure=fig, height_ratios=[1.2, 1],
               hspace=0.65, wspace=0.42)
ax_a = fig.add_subplot(gs[0, :])
ax_b = fig.add_subplot(gs[1, 0])
ax_c = fig.add_subplot(gs[1, 1])

# ══════════════════════════════════════════════════════════════════════════════
# Panel A — Stacked bars: mechanism (M.m. domesticus, P20.5)
# ══════════════════════════════════════════════════════════════════════════════
mmd_p20 = expr_df[expr_df['strain'].isin(MMD) & (expr_df['timepoint'] == 'P20.5')]

BW, INN, GG = 0.30, 0.06, 0.55

xs, smids = [], []
x0 = 0.0
for _ in STAGES:
    xs += [x0, x0 + BW + INN]
    smids.append(x0 + (BW + INN) / 2.0)
    x0 += 2 * BW + INN + GG
xs = np.array(xs)

p_e, p_ne, p_nl, ns, svf = [], [], [], [], []
for stage in STAGES:
    for sv_flag in [True, False]:
        sub = mmd_p20[(mmd_p20['stage'] == stage) & (mmd_p20['has_SV'] == sv_flag)]
        n = len(sub)
        ns.append(n); svf.append(sv_flag)
        if n:
            p_e.append((sub['status'] == 'expressed').sum() / n * 100)
            p_ne.append((sub['status'] == 'not_expressed').sum() / n * 100)
            p_nl.append((sub['status'] == 'not_lifted').sum() / n * 100)
        else:
            p_e.append(0); p_ne.append(0); p_nl.append(0)

p_e  = np.array(p_e)
p_ne = np.array(p_ne)
p_nl = np.array(p_nl)

ax_a.bar(xs, p_e,  width=BW, color=C_EXPR,  zorder=3)
ax_a.bar(xs, p_ne, width=BW, color=C_NEXPR, bottom=p_e,        zorder=3)
ax_a.bar(xs, p_nl, width=BW, color=C_NLIFT, bottom=p_e + p_ne, zorder=3)

# Black outline on SV bars to distinguish them
for xi, sv, te, tne, tnl in zip(xs, svf, p_e, p_ne, p_nl):
    if sv:
        ax_a.bar(xi, te + tne + tnl, width=BW, fill=False,
                 edgecolor='black', linewidth=1.0, zorder=4)

# n= labels above each bar
for xi, n, te, tne, tnl in zip(xs, ns, p_e, p_ne, p_nl):
    ax_a.text(xi, te + tne + tnl + 1.5, f'n={n}',
              ha='center', va='bottom', fontsize=5.5, color='#444444')

# SV / no-SV sub-labels
for xi, sv in zip(xs, svf):
    ax_a.text(xi, -4, "SV" if sv else "no SV",
              ha='center', va='top', fontsize=6.5, clip_on=False,
              fontweight='bold' if sv else 'normal')

# Stage labels and vertical dividers
for i, stage in enumerate(STAGES):
    ax_a.text(smids[i], -11, stage, ha='center', va='top',
              fontsize=8, fontweight='bold',
              color=STAGE_C[stage], clip_on=False)
    if i < len(STAGES) - 1:
        xd = (xs[2*i + 1] + xs[2*(i+1)]) / 2
        ax_a.axvline(xd, color='#cccccc', lw=0.6, ls='--', zorder=1)

ax_a.set_xticks([])
ax_a.set_ylabel("% loci", fontsize=8)
ax_a.set_ylim(0, 115)
ax_a.set_xlim(xs[0] - BW * 0.8, xs[-1] + BW * 0.8)

legend_patches = [
    mpatches.Patch(color=C_EXPR,  label='Expressed'),
    mpatches.Patch(color=C_NEXPR, label='Not expressed (silenced)'),
    mpatches.Patch(color=C_NLIFT, label='Not lifted (rearranged)'),
]
ax_a.legend(handles=legend_patches, fontsize=7.5, frameon=False,
            loc='upper right', ncol=3)
ax_a.set_title(
    "A   Disrupted piRNA loci CO-OCCUR with structural variants — disruption is mostly genomic rearrangement, not silencing  "
    "(M. musculus domesticus, n=12 strains, P20.5)\n"
    "    Outlined bar = loci with SV · among disrupted loci, most are liftover-failure/rearrangement (orange), "
    "not expression loss in intact loci (amber). Direct-window SV vs no-SV difference is weak (see docstring).",
    fontsize=7.5, fontweight='bold', loc='left')

# ══════════════════════════════════════════════════════════════════════════════
# Panel B — Disruption scatter: SV vs no-SV, all stage × timepoint
# ══════════════════════════════════════════════════════════════════════════════
mmd_e = expr_df[expr_df['strain'].isin(MMD)].copy()

pts = []
for tp in TIMEPOINTS:
    for stage in STAGES:
        sub = mmd_e[(mmd_e['timepoint'] == tp) & (mmd_e['stage'] == stage)]
        d_sv   = sub[sub['has_SV']]['disrupted'].mean() * 100
        d_nosv = sub[~sub['has_SV']]['disrupted'].mean() * 100
        pts.append({'tp': tp, 'stage': stage, 'd_sv': d_sv, 'd_nosv': d_nosv})
ptsdf = pd.DataFrame(pts)

maxv = float(ptsdf[['d_sv','d_nosv']].max().max()) + 8
ax_b.plot([0, maxv], [0, maxv], color='#aaaaaa', lw=0.8, ls='--', zorder=1)
ax_b.fill_between([0, maxv], [0, maxv], [maxv, maxv],
                  alpha=0.06, color='#D55E00', zorder=0)

for tp in TIMEPOINTS:
    for stage in STAGES:
        row = ptsdf[(ptsdf['tp'] == tp) & (ptsdf['stage'] == stage)].iloc[0]
        ax_b.scatter(row.d_nosv, row.d_sv, c=STAGE_C[stage],
                     marker=MARKER[tp], s=65, zorder=3,
                     edgecolors='white', linewidths=0.5)

ax_b.set_xlim(-1, maxv); ax_b.set_ylim(-1, maxv)
ax_b.set_xlabel("% disrupted — no SV", fontsize=8)
ax_b.set_ylabel("% disrupted — has SV", fontsize=8)
ax_b.set_title("B   SV vs no-SV disruption (direct window): weak,\nlargely on the diagonal across stage × timepoint groups",
               fontsize=8, fontweight='bold', loc='left')
ax_b.text(maxv * 0.96, maxv * 0.04, "SV worse\n(above line)",
          ha='right', va='bottom', fontsize=6.5, color='#D55E00')

for stage in STAGES:
    ax_b.scatter([], [], c=STAGE_C[stage], s=30, label=stage)
for tp, mk in MARKER.items():
    ax_b.scatter([], [], c='#888888', marker=mk, s=30, label=tp)
ax_b.legend(fontsize=6.5, frameon=False, ncol=2, loc='upper left',
            labelspacing=0.3, handlelength=1.2)

# ══════════════════════════════════════════════════════════════════════════════
# Panel C — Per-strain: SV locus count vs disruption (P20.5 Pachytene)
# ══════════════════════════════════════════════════════════════════════════════
sv_per_strain = (sv_direct.groupby('strain')['has_SV'].sum()
                 .reset_index().rename(columns={'has_SV': 'n_sv_loci'}))
pachy_p20 = expr_df[(expr_df['timepoint'] == 'P20.5') &
                    (expr_df['stage'] == 'Pachytene')]
disrupt_s = ((pachy_p20.groupby('strain')['disrupted'].mean() * 100)
             .reset_index().rename(columns={'disrupted': 'pct_disrupted'}))
sc_df = sv_per_strain.merge(disrupt_s, on='strain')
sc_df['group'] = sc_df['strain'].map(SUBSP_GROUP)
sc_df['short'] = sc_df['strain'].map(STRAIN_SHORT)

for grp in ["domesticus","wild","musculus","castaneus","spretus"]:
    sub = sc_df[sc_df['group'] == grp]
    if sub.empty:
        continue
    ax_c.scatter(sub['n_sv_loci'], sub['pct_disrupted'],
                 c=SUBSP_C[grp], s=55, label=SUBSP_LABELS[grp],
                 zorder=3, edgecolors='white', linewidths=0.4)

for _, row in sc_df.iterrows():
    ax_c.text(row.n_sv_loci + 0.4, row.pct_disrupted + 0.4,
              row['short'], fontsize=5.5, color='#333333')

r_val, p_val = pearsonr(sc_df['n_sv_loci'], sc_df['pct_disrupted'])
ax_c.text(0.97, 0.06, f"Pearson r = {r_val:.2f}\np = {p_val:.1e}",
          ha='right', va='bottom', transform=ax_c.transAxes, fontsize=7,
          bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.8, lw=0))

ax_c.set_xlabel("Zamore loci with direct SV (of 214)", fontsize=8)
ax_c.set_ylabel("% piRNA loci disrupted (P20.5, Pachytene)", fontsize=8)
ax_c.set_title("C   SV burden scales with genomic divergence\n(all 16 strains; colour = subspecies)",
               fontsize=8, fontweight='bold', loc='left')
ax_c.legend(fontsize=6.5, frameon=False, loc='upper left',
            labelspacing=0.3, handlelength=1.2)

fig.suptitle(
    "Pangenome structural variants predict piRNA cluster disruption through genomic rearrangement\n"
    "214 Zamore loci · 16 inbred mouse strains · PICB 2-of-3 consensus · pangenome VCF SVs ≥300 bp",
    fontsize=9, fontweight='bold', y=1.015)

for ext in ('pdf', 'svg', 'png'):
    fig.savefig(f"{OUT}/Fig_SV_mechanism.{ext}", dpi=300, bbox_inches='tight')
plt.close(fig)
print("Saved Fig_SV_mechanism.{pdf,svg,png}")

# Print Panel A data
print("\n=== Panel A data (P20.5, M.m. domesticus) ===")
for stage in STAGES:
    for sv_flag in [True, False]:
        sub = mmd_p20[(mmd_p20['stage'] == stage) & (mmd_p20['has_SV'] == sv_flag)]
        n = len(sub)
        if n:
            print(f"  {stage:16s} {'SV' if sv_flag else 'noSV':5s}: "
                  f"expressed={100*(sub['status']=='expressed').mean():.1f}%  "
                  f"not_expr={100*(sub['status']=='not_expressed').mean():.1f}%  "
                  f"not_lift={100*(sub['status']=='not_lifted').mean():.1f}%  n={n}")
