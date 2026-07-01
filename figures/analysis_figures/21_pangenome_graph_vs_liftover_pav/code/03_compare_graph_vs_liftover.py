#!/usr/bin/env python3
"""THEME 21 step 3 — compare GRAPH sequence-PAV (odgi pav) vs HAL-LIFTOVER cluster-PAV at the 42,384 master loci.
Liftover = which strains have a piRNA CLUSTER lifting here (regulatory presence). Graph = PAV ratio (fraction of the
locus's graph nodes a strain traverses) >= THRESH = genetic sequence presence. Disagreements:
  liftover-absent & graph-present  -> SILENCING (sequence there, no cluster; epigenetic)
  liftover-present & graph-present  -> expressed (concordant)
  both absent                       -> genetic LOSS / never-present
A liftover-PRIVATE locus the graph shows as shared = reference-bias rescue. odgi pav matrix: chrom start end name
<17 samples>; merge to the liftover matrix by locus_id (= 'name')."""
import pandas as pd, numpy as np
D="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/figures/analysis_figures/21_pangenome_graph_vs_liftover_pav/data"
THRESH=0.5
STRAINS=["129S1_SvImJ","A_J","AKR_J","BALB_cJ","C3H_HeJ","C57BL_6NJ","CAST_EiJ","CBA_J","DBA_2J","FVB_NJ","LP_J","NOD_ShiLtJ","NZO_HlLtJ","PWK_PhJ","SPRET_EiJ","WSB_EiJ"]
lift=pd.read_csv(f"{D}/liftover_pav_matrix.tsv",sep="\t",dtype={"chrom":str})
g=pd.read_csv(f"{D}/graph_pav_matrix.tsv",sep="\t").rename(columns={"name":"locus_id"})
print("graph PAV-ratio distribution (strain cells):",np.percentile(g[STRAINS].values.ravel(),[5,25,50,75,95]).round(2))
gren={s:f"g_{s}" for s in STRAINS}; g=g.rename(columns=gren)
m=lift.merge(g[["locus_id"]+[f"g_{s}" for s in STRAINS]],on="locus_id",how="inner")
print(f"merged {len(m):,} loci (liftover {len(lift):,}, graph {len(g):,})")
lp=m[STRAINS].values.astype(int)                                   # liftover cluster presence
gp=(m[[f"g_{s}" for s in STRAINS]].values>=THRESH).astype(int)     # graph sequence presence
cls=lambda n: "core" if n>=16 else ("private" if n==1 else ("dispensable" if n>=2 else "absent"))
ln=lp.sum(1); gn=gp.sum(1)
out=m[["locus_id","chrom","start","end"]].copy()
out["liftover_n"]=ln; out["graph_n"]=gn
out["liftover_class"]=[cls(n) for n in ln]; out["graph_class"]=[cls(n) for n in gn]
out["n_silenced_strains"]=((lp==0)&(gp==1)).sum(1)
out["n_loss_strains"]=((lp==0)&(gp==0)).sum(1)
out.to_csv(f"{D}/graph_vs_liftover_comparison.tsv",sep="\t",index=False)
print("\n=== confusion: liftover class (rows) x graph class (cols) ===")
print(pd.crosstab(out.liftover_class,out.graph_class))
sil=int(((lp==0)&(gp==1)).sum()); loss=int(((lp==0)&(gp==0)).sum()); both=int(((lp==1)&(gp==1)).sum())
tot=lp.size
print(f"\nstrain-locus events ({tot:,}): concordant-present={both:,} ({100*both/tot:.0f}%)  SILENCING(seq+/cluster-)={sil:,} ({100*sil/tot:.0f}%)  loss/absent={loss:,} ({100*loss/tot:.0f}%)")
resc=int(((out.liftover_class=="private")&out.graph_class.isin(["dispensable","core"])).sum()); npriv=int((out.liftover_class=="private").sum())
print(f"liftover-PRIVATE loci RESCUED to shared by graph: {resc:,}/{npriv:,} ({100*resc/max(npriv,1):.0f}%) = reference bias")
print(f"loci with >=1 silenced strain (seq present, cluster absent): {(out.n_silenced_strains>0).sum():,} ({100*(out.n_silenced_strains>0).mean():.0f}%)")
print(f"median graph_n vs liftover_n: {np.median(gn):.0f} vs {np.median(ln):.0f}  (graph sees more sequence than liftover sees clusters)")
