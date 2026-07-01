#!/usr/bin/env python3
"""THEME 22 step 5 — compare the GRAPH-inject cluster-PAV (odgi pav of injected cluster-paths, graph node coverage)
against the HAL-liftover cluster-PAV (theme 21), at the same 42,384 master loci. Both ask 'which strains have a piRNA
CLUSTER here', via two methods: graph = injected cluster-paths' node coverage (reference-free placement); liftover =
halLiftover coordinate overlap. KEEP BOTH. Expected: broad agreement, with the graph recovering clusters HAL
mis-lifts/fragments (esp. divergent strains)."""
import pandas as pd, numpy as np
B="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
T21=f"{B}/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav/data"
T22=f"{B}/figures/analysis_figures/22_odgi_inject_cluster_pav/data"
S=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
WILD={"SPRET_EiJ","CAST_EiJ","PWK_PhJ","WSB_EiJ"}
THRESH=0.1
lift=pd.read_csv(f"{T21}/liftover_pav_matrix.tsv",sep="\t",dtype={"chrom":str})
g=pd.read_csv(f"{T22}/graph_cluster_pav_matrix.tsv",sep="\t").rename(columns={"name":"locus_id"})
m=lift.merge(g[["locus_id"]+S].rename(columns={s:"g_"+s for s in S}),on="locus_id")
print(f"merged {len(m):,} loci (lift {len(lift):,}, graph {len(g):,})")
lp=m[S].values.astype(int); gr=m[["g_"+s for s in S]].values
cls=lambda n:"core" if n>=16 else("private" if n==1 else("dispensable" if n>=2 else "absent"))
for TH in [0.05,0.1,0.2]:
    gp=(gr>=TH).astype(int); ln=lp.sum(1); gn=gp.sum(1)
    agree=100*(lp==gp).mean()
    lc=pd.Series([cls(n) for n in ln]).value_counts(); gc=pd.Series([cls(n) for n in gn]).value_counts()
    print(f"\nTHRESH={TH}: per strain×locus agreement={agree:.0f}% | median lift_n={int(np.median(ln))} graph_n={int(np.median(gn))} | graph extra calls={int(gp.sum()-lp.sum()):+,}")
    print(f"  liftover classes: core={lc.get('core',0)} disp={lc.get('dispensable',0)} private={lc.get('private',0)}")
    print(f"  graph    classes: core={gc.get('core',0)} disp={gc.get('dispensable',0)} private={gc.get('private',0)} absent={gc.get('absent',0)}")
# main comparison at THRESH=0.1
gp=(gr>=THRESH).astype(int); ln=lp.sum(1); gn=gp.sum(1)
print("\n=== confusion: liftover class (rows) x graph class (cols) @0.1 ===")
print(pd.crosstab(pd.Series([cls(n) for n in ln],name="liftover"),pd.Series([cls(n) for n in gn],name="graph")))
print("\n=== per-strain agreement + cluster calls (divergent strains differ more?) ===")
for i,s in enumerate(S):
    ag=100*(lp[:,i]==gp[:,i]).mean()
    print(f"  {s:13s} {'(wild)' if s in WILD else '      '} agree={ag:.0f}%  liftover+={int(lp[:,i].sum()):6d}  graph+={int(gp[:,i].sum()):6d}  graph-only={int(((gp[:,i]==1)&(lp[:,i]==0)).sum()):5d}  lift-only={int(((gp[:,i]==0)&(lp[:,i]==1)).sum()):5d}")
out=m[["locus_id","chrom","start","end"]].copy(); out["liftover_n"]=ln; out["graph_n"]=gn
out["liftover_class"]=[cls(n) for n in ln]; out["graph_class"]=[cls(n) for n in gn]
out.to_csv(f"{T22}/cluster_pav_comparison.tsv",sep="\t",index=False)
print(f"\nwrote cluster_pav_comparison.tsv ({len(out):,} loci)")
