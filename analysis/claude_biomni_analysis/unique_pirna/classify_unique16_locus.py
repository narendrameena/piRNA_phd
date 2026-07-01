#!/usr/bin/env python3
"""Scale-16 pangenome phase — STEP 3 (final): combine the sequence-level exact-expression split
(classify_unique16.py -> expr_exact_classified.csv.gz) with the pangenome cross-strain LOCUS-PRESENCE
result (lift_presence16.sh -> loci/present_in_{Y}.bed) into the 4-way Step-4 classification across 16 strains:

  expressed elsewhere (exact)      -> NOT unique          (exact sequence in another strain's pool)
  unique: strain-private locus     -> UNIQUE (locus gain) (novel seq; orthologous locus in NO other strain)
  unique: conserved-but-silent     -> UNIQUE (expression) (novel seq; orthologous locus EXISTS elsewhere,
                                                            but not exact-expressed there)
The SNP-variant (1-3 mm) refinement of conserved-but-silent (fetch the other strain's allele at the
orthologous locus, compare <=3 mm, check it is in that strain's pool) is a documented add-on (needs the
per-strain genome FASTA + reverse-lift) — see PANGENOME_PHASE_README.md. Run after lift_presence16.sh."""
import gzip, os
import pandas as pd
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
OUT=f"{U}/unique16"; LOCI=f"{OUT}/loci"
STRAINS=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J",
         "NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
ec=f"{OUT}/expr_exact_classified.csv.gz"
if not os.path.exists(ec): raise SystemExit("[wait] run classify_unique16.py first (expr_exact_classified.csv.gz)")
miss=[Y for Y in STRAINS if not os.path.exists(f"{LOCI}/present_in_{Y}.bed")]
if miss: raise SystemExit(f"[wait] presence lifts missing for {miss[:3]}... (lift_presence16.sh)")

# candidate id -> set of strains whose genome contains the orthologous locus (any of its loci)
present={}
for Y in STRAINS:
    with open(f"{LOCI}/present_in_{Y}.bed") as fh:
        for line in fh:
            f=line.rstrip("\n").split("\t")
            if len(f)>=4: present.setdefault(f[3],set()).add(Y)

d=pd.read_csv(ec)
d["cand_id"]=d.strain+"|"+d.timepoint+"|"+d.sequence
def final_klass(row):
    if not row.expr_class.startswith("novel"): return "expressed elsewhere (exact)"
    others=present.get(row.cand_id,set())-{row.strain}
    return "unique: conserved-but-silent" if others else "unique: strain-private locus"
d["klass"]=d.apply(final_klass,axis=1)
d["homolog_strains"]=d.cand_id.map(lambda c: ",".join(sorted(present.get(c,set()))))
d.drop(columns=["cand_id"]).to_csv(f"{OUT}/final_classified.csv.gz",index=False)

print("Scale-16 cross-strain classification (per strain x timepoint candidates):")
print(d.groupby(["timepoint","klass"]).size().unstack(fill_value=0).to_string())
print("\nby strain:"); print(d.groupby(["strain","klass"]).size().unstack(fill_value=0).to_string())
print(f"\nwrote {OUT}/final_classified.csv.gz ({len(d):,} rows)")
