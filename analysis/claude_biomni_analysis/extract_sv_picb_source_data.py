#!/usr/bin/env python3
"""
Extract publication source-data for the SV/PICB figures, reproducing the EXACT
plotted values from Fig_PICB_vs_Zamore.py and Fig_SV_TE.py (logic copied verbatim).
Non-destructive: reads the same matrices, writes only SourceData_*.csv. Per-strain
tables are ordered by the canonical figure strain order (strain_order.py).
"""
import os, sys
sys.path.insert(0, "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis")
import numpy as np, pandas as pd
from scipy.stats import pearsonr, mannwhitneyu
from strain_order import strain_rank, WILD

BASE = "/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
OUT  = f"{BASE}/analysis/claude_biomni_analysis/all_strains_pangenome"
SD   = f"{BASE}/analysis/claude_biomni_analysis/source_data"

MMD = ["C57BL_6NJ","DBA_2J","BALB_cJ","A_J","CBA_J","129S1_SvImJ",
       "NOD_ShiLtJ","AKR_J","C3H_HeJ","NZO_HlLtJ","LP_J","FVB_NJ"]
SUBSP_GROUP = {**{s:"M.m.domesticus" for s in MMD},
    "WSB_EiJ":"M.m.domesticus(wild,WSB)","PWK_PhJ":"M.m.musculus(PWK)",
    "CAST_EiJ":"M.m.castaneus(CAST)","SPRET_EiJ":"M.spretus(SPRET)"}
TE_ORDER = ['LINE/L1','SINE','LTR/ERVL','LTR/ERVK','LTR/ERV1','DNA','Other','No_TE']
STAGES = ["Prepachytene","Pachytene","Hybrid"]

def by_strain(df, col="strain"):
    """sort a per-strain table into canonical figure order"""
    d = df.copy(); d["_r"] = d[col].map(strain_rank)
    d = d.sort_values("_r").drop(columns="_r").reset_index(drop=True)
    d.insert(0, "plot_order", range(1, len(d)+1))
    return d

manifest = []   # (figure, panel, file, note)

# ════════════════════════════════════════════════════════════════════════════
# Fig_PICB_vs_Zamore.py
# ════════════════════════════════════════════════════════════════════════════
z_expr = pd.read_csv(f"{OUT}/all_strains_expression_matrix.csv")
z_sv   = pd.read_csv(f"{OUT}/all_strains_SV_matrix.csv")
z_sv_d = z_sv[z_sv['window']=='direct'][['locus','strain','has_SV']].copy()
z_expr = z_expr.merge(z_sv_d, on=['locus','strain'], how='left')
z_expr['has_SV']=z_expr['has_SV'].fillna(False)
z_expr['disrupted']=z_expr['status'].isin(['not_expressed','not_lifted'])

p_expr = pd.read_csv(f"{OUT}/all_strains_expression_matrix_picb.csv")
p_sv   = pd.read_csv(f"{OUT}/all_strains_SV_matrix_picb.csv")
p_sv_d = p_sv[p_sv['window']=='direct'][['locus','stage','strain','has_SV']].copy()
p_expr = p_expr.merge(p_sv_d, on=['locus','stage','strain'], how='left')
p_expr['has_SV']=p_expr['has_SV'].fillna(False)
p_expr['disrupted']=p_expr['status'].isin(['not_expressed','not_lifted'])

# --- Panels A & B: status breakdown for 4 groups -----------------------------
rows=[]
for dataset, sub_sv, sub_nosv in [("PICB", p_expr[p_expr.has_SV], p_expr[~p_expr.has_SV]),
                                   ("Zamore", z_expr[z_expr.has_SV], z_expr[~z_expr.has_SV])]:
    for svlab, sub in [("has_SV", sub_sv), ("no_SV", sub_nosv)]:
        n=len(sub)
        rows.append(dict(dataset=dataset, sv_status=svlab, n_loci_rows=n,
            pct_expressed=100*(sub.status=='expressed').mean(),
            pct_not_expressed=100*(sub.status=='not_expressed').mean(),
            pct_not_lifted=100*(sub.status=='not_lifted').mean()))
ab=pd.DataFrame(rows)
ab.to_csv(f"{SD}/SourceData_Fig_PICB_vs_Zamore_panelAB.csv", index=False)
manifest.append(("Fig_PICB_vs_Zamore","A (stacked status) + B (%not_lifted)",
                 "SourceData_Fig_PICB_vs_Zamore_panelAB.csv",
                 "4 groups (PICB/Zamore x has_SV/no_SV); panel B = pct_not_lifted column"))

# --- Panel C: PICB per-strain SV count vs %not_lifted ------------------------
pC = (p_sv_d.groupby('strain')['has_SV'].sum().reset_index()
      .rename(columns={'has_SV':'n_sv_loci'})).merge(
      p_expr.groupby('strain')['status'].apply(lambda g:100*(g=='not_lifted').mean())
      .reset_index(name='pct_not_lifted'), on='strain')
rC,pvC = pearsonr(pC['n_sv_loci'], pC['pct_not_lifted'])
pC['subspecies']=pC['strain'].map(SUBSP_GROUP)
pC=by_strain(pC)[['plot_order','strain','subspecies','n_sv_loci','pct_not_lifted']]
pC.to_csv(f"{SD}/SourceData_Fig_PICB_vs_Zamore_panelC.csv", index=False)
manifest.append(("Fig_PICB_vs_Zamore","C (PICB per-strain scatter)",
                 "SourceData_Fig_PICB_vs_Zamore_panelC.csv",
                 f"16 strains; Pearson r={rC:.3f}, p={pvC:.2e}"))

# --- Panel D: Zamore per-strain SV count vs %disrupted (P20.5 Pachytene) ------
pD = (z_sv_d.groupby('strain')['has_SV'].sum().reset_index()
      .rename(columns={'has_SV':'n_sv_loci'}))
pachy = z_expr[(z_expr.timepoint=='P20.5') & (z_expr.stage=='Pachytene')]
pD = pD.merge((pachy.groupby('strain')['disrupted'].mean()*100).reset_index()
              .rename(columns={'disrupted':'pct_disrupted'}), on='strain')
rD,pvD = pearsonr(pD['n_sv_loci'], pD['pct_disrupted'])
pD['subspecies']=pD['strain'].map(SUBSP_GROUP)
pD=by_strain(pD)[['plot_order','strain','subspecies','n_sv_loci','pct_disrupted']]
pD.to_csv(f"{SD}/SourceData_Fig_PICB_vs_Zamore_panelD.csv", index=False)
manifest.append(("Fig_PICB_vs_Zamore","D (Zamore per-strain scatter)",
                 "SourceData_Fig_PICB_vs_Zamore_panelD.csv",
                 f"16 strains; Pearson r={rD:.3f}, p={pvD:.2e}"))

# ════════════════════════════════════════════════════════════════════════════
# Fig_SV_TE.py
# ════════════════════════════════════════════════════════════════════════════
df_pairs = pd.read_csv(f"{OUT}/sv_te_annotation.csv")
expr_df  = pd.read_csv(f"{OUT}/all_strains_expression_matrix.csv")
sv_df    = pd.read_csv(f"{OUT}/all_strains_SV_matrix.csv")
sv_direct= sv_df[sv_df.window=='direct'][['locus','strain','has_SV']].copy()
expr_df  = expr_df.merge(sv_direct, on=['locus','strain'], how='left')
expr_df['has_SV']=expr_df['has_SV'].fillna(False)
expr_df['disrupted']=expr_df['status'].isin(['not_expressed','not_lifted'])
p20_status=(expr_df[expr_df.timepoint=='P20.5'].set_index(['strain','locus'])['status'].to_dict())
df_pairs['status']=df_pairs.apply(lambda r:p20_status.get((r.strain,r.locus),'unknown'),axis=1)
df_pairs['group']=df_pairs['status'].map({'not_lifted':'Disruptive','expressed':'Expressed'}).fillna('other')
df_filt=df_pairs[df_pairs['group'].isin(['Disruptive','Expressed'])].copy()
spl=df_filt.groupby(['strain','locus','group']).size().reset_index(name='n_svs')

# --- Panel A: mechanism (M.m.domesticus, P20.5, 3 stages x SV/no_SV) ---------
rowsA=[]
mmd_p20=expr_df[expr_df.strain.isin(MMD) & (expr_df.timepoint=='P20.5')]
for stage in STAGES:
    for sv_flag in [True,False]:
        sub=mmd_p20[(mmd_p20.stage==stage)&(mmd_p20.has_SV==sv_flag)]; n=len(sub)
        rowsA.append(dict(stage=stage, sv_status="has_SV" if sv_flag else "no_SV", n_loci_rows=n,
            pct_expressed=(sub.status=='expressed').sum()/n*100 if n else 0,
            pct_not_expressed=(sub.status=='not_expressed').sum()/n*100 if n else 0,
            pct_not_lifted=(sub.status=='not_lifted').sum()/n*100 if n else 0))
pd.DataFrame(rowsA).to_csv(f"{SD}/SourceData_Fig_SV_TE_panelA.csv", index=False)
manifest.append(("Fig_SV_TE","A (mechanism stacked bars)","SourceData_Fig_SV_TE_panelA.csv",
                 "M.m.domesticus, P20.5; 3 stages x has_SV/no_SV"))

# --- Panel B: TE composition (All / INS / DEL) -------------------------------
rowsB=[]
groupsB=[("All_SVs",df_filt),("INS",df_filt[df_filt.sv_type=='INS']),("DEL",df_filt[df_filt.sv_type=='DEL'])]
for lbl,g in groupsB:
    for te in TE_ORDER:
        rowsB.append(dict(sv_group=lbl, n_svs=len(g), te_category=te,
                          pct_of_SVs=(g.te_cat==te).sum()/len(g)*100 if len(g) else 0))
pd.DataFrame(rowsB).to_csv(f"{SD}/SourceData_Fig_SV_TE_panelB.csv", index=False)
manifest.append(("Fig_SV_TE","B (TE composition)","SourceData_Fig_SV_TE_panelB.csv",
                 "100% stacked TE category fractions for All/INS/DEL"))

# --- Panel C: SV burden per locus (violin points), canonical strain order ----
splC=spl.copy(); splC['_r']=splC['strain'].map(strain_rank)
splC=splC.sort_values(['group','_r','locus']).drop(columns='_r')
splC.to_csv(f"{SD}/SourceData_Fig_SV_TE_panelC.csv", index=False)
d_svs=spl[spl.group=='Disruptive']['n_svs'].values; e_svs=spl[spl.group=='Expressed']['n_svs'].values
_,p_mwu=mannwhitneyu(d_svs,e_svs,alternative='greater')
manifest.append(("Fig_SV_TE","C (SV burden violin)","SourceData_Fig_SV_TE_panelC.csv",
    f"per (strain,locus,group) n_svs; Expressed median={np.median(e_svs):.0f} "
    f"(n={len(e_svs)}), Disruptive median={np.median(d_svs):.0f} (n={len(d_svs)}), MWU p={p_mwu:.2e}"))

# --- Panel D: per-strain SV count vs %disrupted (P20.5 Pachytene) ------------
sD=(sv_direct.groupby('strain')['has_SV'].sum().reset_index().rename(columns={'has_SV':'n_sv_loci'}))
pachy2=expr_df[(expr_df.timepoint=='P20.5')&(expr_df.stage=='Pachytene')]
sD=sD.merge((pachy2.groupby('strain')['disrupted'].mean()*100).reset_index()
            .rename(columns={'disrupted':'pct_disrupted'}), on='strain')
rD2,pvD2=pearsonr(sD['n_sv_loci'],sD['pct_disrupted'])
sD['subspecies']=sD['strain'].map(SUBSP_GROUP)
sD=by_strain(sD)[['plot_order','strain','subspecies','n_sv_loci','pct_disrupted']]
sD.to_csv(f"{SD}/SourceData_Fig_SV_TE_panelD.csv", index=False)
manifest.append(("Fig_SV_TE","D (per-strain scatter)","SourceData_Fig_SV_TE_panelD.csv",
                 f"16 strains; Pearson r={rD2:.3f}, p={pvD2:.2e}"))

print("WROTE SOURCE DATA:")
for fig,panel,f,note in manifest:
    print(f"  [{fig}] {panel}\n      -> {f}\n      {note}")
# machine-readable manifest for README assembly
pd.DataFrame(manifest, columns=["figure","panel","source_data_file","note"]).to_csv(
    f"{SD}/_sv_picb_manifest.csv", index=False)
