#!/usr/bin/env python3
"""
PICB all-clusters vs Zamore conserved loci: SV analysis comparison.

Panel A: Scale and expression overview — 28k PICB clusters (C57BL_6NJ reference)
         vs 214 Zamore conserved loci; expression rates across all strains
Panel B: Status breakdown with/without SV — PICB clusters (all 16 strains, all timepoints)
         vs Zamore loci; focus on not_lifted (structural disruption)
Panel C: LiftOver failure rate with/without SV — PICB vs Zamore
Panel D: Per-strain SV count vs disruption — PICB clusters (r=0.33) vs Zamore (r=0.80)
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy.stats import pearsonr

BASE = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
OUT  = f"{BASE}/analysis/claude_biomni_analysis/all_strains_pangenome"

plt.rcParams.update({
    "font.family": "Liberation Sans", "font.size": 8,
    "axes.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 300, "savefig.dpi": 300, "pdf.fonttype": 42,
})

C_EXPR  = "#0072B2"
C_NEXPR = "#E69F00"
C_NLIFT = "#D55E00"
MMD = ["C57BL_6NJ","DBA_2J","BALB_cJ","A_J","CBA_J","129S1_SvImJ",
       "NOD_ShiLtJ","AKR_J","C3H_HeJ","NZO_HlLtJ","LP_J","FVB_NJ"]
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
    "domesticus": "M.m. domesticus (n=12)", "wild": "WSB/EiJ",
    "musculus": "M.m. musculus (PWK)", "castaneus": "M.m. castaneus (CAST)",
    "spretus": "M. spretus (SPRET)",
}
STRAIN_SHORT = {
    "C57BL_6NJ": "C57BL/6", "DBA_2J": "DBA/2J", "BALB_cJ": "BALB/c",
    "A_J": "A/J", "CBA_J": "CBA/J", "129S1_SvImJ": "129S1",
    "NOD_ShiLtJ": "NOD", "AKR_J": "AKR/J", "C3H_HeJ": "C3H",
    "NZO_HlLtJ": "NZO", "LP_J": "LP/J", "FVB_NJ": "FVB/NJ",
    "WSB_EiJ": "WSB", "PWK_PhJ": "PWK", "CAST_EiJ": "CAST", "SPRET_EiJ": "SPRET",
}

# ── Load Zamore data ──────────────────────────────────────────────────────────
z_expr = pd.read_csv(f"{OUT}/all_strains_expression_matrix.csv")
z_sv   = pd.read_csv(f"{OUT}/all_strains_SV_matrix.csv")
z_sv_d = z_sv[z_sv['window'] == 'direct'][['locus','strain','has_SV']].copy()
z_expr = z_expr.merge(z_sv_d, on=['locus','strain'], how='left')
z_expr['has_SV']    = z_expr['has_SV'].fillna(False)
z_expr['disrupted'] = z_expr['status'].isin(['not_expressed','not_lifted'])

# ── Load PICB data ─────────────────────────────────────────────────────────────
p_expr = pd.read_csv(f"{OUT}/all_strains_expression_matrix_picb.csv")
p_sv   = pd.read_csv(f"{OUT}/all_strains_SV_matrix_picb.csv")
p_sv_d = p_sv[p_sv['window'] == 'direct'][['locus','stage','strain','has_SV']].copy()
p_expr = p_expr.merge(p_sv_d, on=['locus','stage','strain'], how='left')
p_expr['has_SV']    = p_expr['has_SV'].fillna(False)
p_expr['disrupted'] = p_expr['status'].isin(['not_expressed','not_lifted'])

# Summary stats for print
print("=== PICB overview ===")
print(f"  Unique loci: {p_expr['locus'].nunique()}")
print(f"  Total rows: {len(p_expr)}")
print(f"  Status counts: {p_expr['status'].value_counts().to_dict()}")
print(f"  Stage: {p_expr['stage'].value_counts().to_dict()}")
print(f"  has_SV rows: {p_expr['has_SV'].sum()} / {len(p_expr)}")

for sv_flag in [True, False]:
    sub = p_expr[p_expr['has_SV'] == sv_flag]
    vc = sub['status'].value_counts(normalize=True)
    print(f"  has_SV={sv_flag}: {dict(vc.round(3))}  n={len(sub)}")

print("\n=== Zamore overview ===")
print(f"  Unique loci: {z_expr['locus'].nunique()}")
for sv_flag in [True, False]:
    sub = z_expr[z_expr['has_SV'] == sv_flag]
    vc = sub['status'].value_counts(normalize=True)
    print(f"  has_SV={sv_flag}: {dict(vc.round(3))}  n={len(sub)}")

# ── Figure layout ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 11))
gs = GridSpec(2, 2, figure=fig, hspace=0.60, wspace=0.42,
              height_ratios=[1, 1])
ax_a = fig.add_subplot(gs[0, 0])
ax_b = fig.add_subplot(gs[0, 1])
ax_c = fig.add_subplot(gs[1, 0])
ax_d = fig.add_subplot(gs[1, 1])

BW = 0.24
INN = 0.07
GRP_GAP = 0.65

# ══════════════════════════════════════════════════════════════════════════════
# Panel A — Scale & expression overview: PICB vs Zamore
# ══════════════════════════════════════════════════════════════════════════════
grp_data = [
    (p_expr[p_expr['has_SV']],   p_expr[~p_expr['has_SV']]),
    (z_expr[z_expr['has_SV']],   z_expr[~z_expr['has_SV']]),
]
grp_labels = ["PICB\n(28,058 loci)", "Zamore\n(214 loci)"]

# 2 pairs: within each pair SV left, no-SV right
# formula: pair k bar j  →  k*(2*(BW+INN)+GRP_GAP) + j*(BW+INN)
xpos = np.array([0, BW + INN,
                 2*(BW + INN) + GRP_GAP,
                 3*(BW + INN) + GRP_GAP])

p_e, p_ne, p_nl, ns = [], [], [], []
all_subs = [grp_data[0][0], grp_data[0][1], grp_data[1][0], grp_data[1][1]]
for sub in all_subs:
    n = len(sub); ns.append(n)
    p_e.append(100 * (sub['status'] == 'expressed').mean())
    p_ne.append(100 * (sub['status'] == 'not_expressed').mean())
    p_nl.append(100 * (sub['status'] == 'not_lifted').mean())
p_e = np.array(p_e); p_ne = np.array(p_ne); p_nl = np.array(p_nl)

ax_a.bar(xpos, p_e,  width=BW, color=C_EXPR,  zorder=3)
ax_a.bar(xpos, p_ne, width=BW, color=C_NEXPR, bottom=p_e,        zorder=3)
ax_a.bar(xpos, p_nl, width=BW, color=C_NLIFT, bottom=p_e + p_ne, zorder=3)

# Outline SV (first of each pair) bars
for xi in [xpos[0], xpos[2]]:
    ax_a.bar(xi, 100, width=BW, fill=False, edgecolor='black', lw=0.9, zorder=4)

# x-tick labels: "SV\n(n=207K)" / "no SV\n(n=1140K)" — carry n value to avoid above-bar text
tick_lbls = []
for (lbl_sv, lbl_nosv), (sv_sub, nosv_sub) in zip(
        [("SV", "no SV"), ("SV", "no SV")],
        [(all_subs[0], all_subs[1]), (all_subs[2], all_subs[3])]):
    for sub, lbl in [(sv_sub, lbl_sv), (nosv_sub, lbl_nosv)]:
        n = len(sub)
        n_str = f"{n/1e3:.0f}K" if n > 9999 else str(n)
        tick_lbls.append(f"{lbl}\n(n={n_str})")

ax_a.set_xticks(xpos)
ax_a.set_xticklabels(tick_lbls, fontsize=6.0)

# Group labels centred under each pair
for i, lbl in enumerate(grp_labels):
    mid = (xpos[2*i] + xpos[2*i + 1]) / 2
    ax_a.text(mid, -18, lbl, ha='center', va='top', fontsize=7.5,
              fontweight='bold', color='#333', clip_on=False)

# Vertical divider
x_div = (xpos[1] + xpos[2]) / 2
ax_a.axvline(x_div, color='#cccccc', lw=0.8, ls='--', zorder=1)

lg = [mpatches.Patch(color=C_EXPR,  label='Expressed'),
      mpatches.Patch(color=C_NEXPR, label='Not expressed'),
      mpatches.Patch(color=C_NLIFT, label='Not lifted')]
ax_a.legend(handles=lg, fontsize=6.5, frameon=False,
            loc='upper center', bbox_to_anchor=(0.5, -0.30), ncol=3)
ax_a.set_ylabel("% loci", fontsize=8)
ax_a.set_ylim(0, 106)
ax_a.set_xlim(xpos[0] - BW * 0.7, xpos[-1] + BW * 0.7)
ax_a.set_title(
    "A   C57BL_6NJ PICB clusters are mostly strain-specific\n"
    "    only ~2% expressed in other strains (vs ~40% Zamore conserved loci)",
    fontsize=7.5, fontweight='bold', loc='left')

# ══════════════════════════════════════════════════════════════════════════════
# Panel B — LiftOver failure rate: the SV signature in both datasets
#   Bar chart: % not_lifted for PICB (has_SV / no_SV) and Zamore (has_SV / no_SV)
#   Shows SVs predict liftover failure in both datasets
# ══════════════════════════════════════════════════════════════════════════════
pairs = [
    ("PICB\n(has SV)",    p_expr[p_expr['has_SV']]),
    ("PICB\n(no SV)",     p_expr[~p_expr['has_SV']]),
    ("Zamore\n(has SV)",  z_expr[z_expr['has_SV']]),
    ("Zamore\n(no SV)",   z_expr[~z_expr['has_SV']]),
]
colors_b = [C_NLIFT, '#f5a672', C_NLIFT, '#f5a672']
nl_pcts = [100 * (sub['status'] == 'not_lifted').mean() for _, sub in pairs]
xb = np.arange(len(pairs))

bars = ax_b.bar(xb, nl_pcts, width=0.55, color=colors_b, zorder=3,
                edgecolor='white', linewidth=0.5)
for xi, pct, (lbl, sub) in zip(xb, nl_pcts, pairs):
    ax_b.text(xi, pct + 0.4, f'{pct:.1f}%', ha='center', va='bottom',
              fontsize=7, fontweight='bold')

# Fold-change labels between SV / no-SV pairs
for i_sv, i_nosv in [(0, 1), (2, 3)]:
    fc = nl_pcts[i_sv] / nl_pcts[i_nosv]
    ym = max(nl_pcts[i_sv], nl_pcts[i_nosv]) + 2.5
    ax_b.annotate('', xy=(xb[i_nosv], ym), xytext=(xb[i_sv], ym),
                  arrowprops=dict(arrowstyle='<->', color='#555', lw=0.8))
    ax_b.text((xb[i_sv] + xb[i_nosv]) / 2, ym + 0.4,
              f'{fc:.1f}x', ha='center', va='bottom', fontsize=7, color='#333')

ax_b.set_xticks(xb)
ax_b.set_xticklabels([lbl for lbl, _ in pairs], fontsize=7.5)
ax_b.set_ylabel("% loci not lifted (structural rearrangement)", fontsize=8)
ax_b.set_ylim(0, max(nl_pcts) * 1.35)
ax_b.axvline(1.5, color='#cccccc', lw=0.8, ls='--', zorder=1)
ax_b.set_title(
    "B   SVs predict liftOver failure in both locus sets\n"
    "    3x enrichment in PICB, 5x in Zamore (pangenome SV = rearrangement)",
    fontsize=7.5, fontweight='bold', loc='left')

# ══════════════════════════════════════════════════════════════════════════════
# Panel C — Per-strain SV count vs liftOver failure rate: PICB
# ══════════════════════════════════════════════════════════════════════════════
p_sv_strain = (p_sv_d.groupby('strain')['has_SV'].sum()
               .reset_index().rename(columns={'has_SV': 'n_sv_loci'}))
p_nl_strain = (p_expr.groupby('strain')['status']
               .apply(lambda g: 100 * (g == 'not_lifted').mean())
               .reset_index(name='pct_not_lifted'))
pc_df = p_sv_strain.merge(p_nl_strain, on='strain')
pc_df['group'] = pc_df['strain'].map(SUBSP_GROUP)
pc_df['short'] = pc_df['strain'].map(STRAIN_SHORT)

r_p, pv_p = pearsonr(pc_df['n_sv_loci'], pc_df['pct_not_lifted'])
for grp in ["domesticus", "wild", "musculus", "castaneus", "spretus"]:
    sub = pc_df[pc_df['group'] == grp]
    if sub.empty: continue
    ax_c.scatter(sub['n_sv_loci'], sub['pct_not_lifted'],
                 c=SUBSP_C[grp], s=55, label=SUBSP_LABELS[grp],
                 zorder=3, edgecolors='white', linewidths=0.4)

# Label only clear outliers; annotate domesticus cluster as a group
OUTLIERS_C = {'SPRET_EiJ', 'PWK_PhJ', 'CAST_EiJ', 'WSB_EiJ', 'C57BL_6NJ'}
for _, row in pc_df.iterrows():
    if row['strain'] in OUTLIERS_C:
        ax_c.text(row.n_sv_loci + 55, row.pct_not_lifted,
                  row['short'], fontsize=6, color='#333', va='center')

# Annotate the dense domesticus cluster
dom = pc_df[~pc_df['strain'].isin(OUTLIERS_C)]
ax_c.annotate("M.m. domesticus\n(11 strains)",
              xy=(dom['n_sv_loci'].mean(), dom['pct_not_lifted'].mean()),
              xytext=(dom['n_sv_loci'].mean() - 1100, dom['pct_not_lifted'].mean() + 6),
              fontsize=6, color=SUBSP_C['domesticus'],
              arrowprops=dict(arrowstyle='->', color=SUBSP_C['domesticus'], lw=0.8),
              ha='center')
ax_c.text(0.97, 0.06, f"r = {r_p:.2f}\np = {pv_p:.1e}",
          ha='right', va='bottom', transform=ax_c.transAxes, fontsize=7,
          bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.8, lw=0))
ax_c.set_xlabel("PICB loci with direct SV", fontsize=8)
ax_c.set_ylabel("% PICB loci not lifted", fontsize=8)
ax_c.set_title(
    "C   PICB: SV burden predicts locus rearrangement\n"
    f"    (all 16 strains; r={r_p:.2f}; colour = subspecies)",
    fontsize=7.5, fontweight='bold', loc='left')
ax_c.legend(fontsize=6.5, frameon=False, loc='upper left',
            labelspacing=0.3, handlelength=1.2)

# ══════════════════════════════════════════════════════════════════════════════
# Panel D — Per-strain Zamore: SV count vs disruption (P20.5, Pachytene)
#   Reference comparison from the Zamore conserved analysis
# ══════════════════════════════════════════════════════════════════════════════
z_sv_strain = (z_sv_d.groupby('strain')['has_SV'].sum()
               .reset_index().rename(columns={'has_SV': 'n_sv_loci'}))
pachy_p20 = z_expr[(z_expr['timepoint'] == 'P20.5') & (z_expr['stage'] == 'Pachytene')]
z_dis_strain = ((pachy_p20.groupby('strain')['disrupted'].mean() * 100)
                .reset_index().rename(columns={'disrupted': 'pct_disrupted'}))
zc_df = z_sv_strain.merge(z_dis_strain, on='strain')
zc_df['group'] = zc_df['strain'].map(SUBSP_GROUP)
zc_df['short'] = zc_df['strain'].map(STRAIN_SHORT)

r_z, pv_z = pearsonr(zc_df['n_sv_loci'], zc_df['pct_disrupted'])
for grp in ["domesticus", "wild", "musculus", "castaneus", "spretus"]:
    sub = zc_df[zc_df['group'] == grp]
    if sub.empty: continue
    ax_d.scatter(sub['n_sv_loci'], sub['pct_disrupted'],
                 c=SUBSP_C[grp], s=55, label=SUBSP_LABELS[grp],
                 zorder=3, edgecolors='white', linewidths=0.4)
OUTLIERS_D = {'SPRET_EiJ', 'PWK_PhJ', 'CAST_EiJ', 'WSB_EiJ', 'C57BL_6NJ'}
for _, row in zc_df.iterrows():
    if row['strain'] in OUTLIERS_D:
        ax_d.text(row.n_sv_loci + 1.0, row.pct_disrupted,
                  row['short'], fontsize=6, color='#333', va='center')

dom_d = zc_df[~zc_df['strain'].isin(OUTLIERS_D)]
ax_d.annotate("M.m. domesticus\n(11 strains)",
              xy=(dom_d['n_sv_loci'].mean(), dom_d['pct_disrupted'].mean()),
              xytext=(dom_d['n_sv_loci'].mean() - 20, dom_d['pct_disrupted'].mean() + 5),
              fontsize=6, color=SUBSP_C['domesticus'],
              arrowprops=dict(arrowstyle='->', color=SUBSP_C['domesticus'], lw=0.8),
              ha='center')
ax_d.text(0.97, 0.06, f"r = {r_z:.2f}\np = {pv_z:.1e}",
          ha='right', va='bottom', transform=ax_d.transAxes, fontsize=7,
          bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.8, lw=0))
ax_d.set_xlabel("Zamore loci with direct SV (of 214)", fontsize=8)
ax_d.set_ylabel("% Zamore loci disrupted (P20.5, Pachytene)", fontsize=8)
ax_d.set_title(
    "D   Zamore conserved: stronger SV-disruption correlation\n"
    f"    (all 16 strains; r={r_z:.2f}; colour = subspecies)",
    fontsize=7.5, fontweight='bold', loc='left')
ax_d.legend(fontsize=6.5, frameon=False, loc='upper left',
            labelspacing=0.3, handlelength=1.2)

fig.suptitle(
    "C57BL_6NJ PICB clusters (n=28,058) vs conserved Zamore loci (n=214): "
    "SV impact on piRNA locus integrity across 16 inbred mouse strains\n"
    "Pangenome VCF SVs >=300 bp · PICB 2-of-3 consensus · liftOver chain files",
    fontsize=8.5, fontweight='bold', y=1.012)

for ext in ('pdf', 'svg', 'png'):
    fig.savefig(f"{OUT}/Fig_PICB_vs_Zamore.{ext}", dpi=300, bbox_inches='tight')
plt.close(fig)
print("Saved Fig_PICB_vs_Zamore.{pdf,svg,png}")
