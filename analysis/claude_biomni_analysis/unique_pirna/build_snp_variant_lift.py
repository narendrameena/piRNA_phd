#!/usr/bin/env python3
"""DEFINITIVE lift-anchored SNP-variant producer -- rescues the STAR-align blind-spot via the cactus lift.

For each candidate, its ORTHOLOGOUS locus in each strain Y is given by unique16/loci/present_in_{Y}.bed
(halLiftover of the candidate's GRCm39 locus onto Y through the cactus HAL -- this captures divergent /
off-assembly orthologs that STAR's read-alignment cannot). Extract Y's genome at that locus (= the ortholog
sequence G), strand-corrected; if G is EXPRESSED in Y (G in Y's expressed pool) and 1-3 substitutions from the
candidate, the candidate is a 1-3mm SNP-variant of the piRNA expressed at its OWN orthologous locus in Y.
Take the min mismatch across strains (0 = expressed-exact -> excluded; 1-3 = SNP-variant).

This needs NO STAR alignment and NO bowtie: the lift + per-strain genome + expressed pool give "expressed AT the
ortholog, within 1-3 substitutions" directly, which is the biologically-correct criterion with no STAR blind-spot
(the lift's only gap is candidates whose GRCm39 projection failed -- those are strain-private, not base-CBS).
Same-length only (SNP = substitution; indel/fragmented lifts skipped). Writes snp_variant_refinement.lift.csv.
"""
import os, gzip, pandas as pd, pysam
from concurrent.futures import ProcessPoolExecutor
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; LOCI=f"{U}/unique16/loci"; POOLS=f"{U}/unique16/pools"
GEN=f"{ROOT}/results/pangenome/prepared"; WORK="/tmp/liftsnp"; os.makedirs(WORK, exist_ok=True)
ALL=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J",
     "NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
_C={"A":"T","T":"A","C":"G","G":"C","N":"N"}
def rc(s): return "".join(_C.get(b,"N") for b in reversed(s))
def ham(a,b): return sum(1 for x,y in zip(a,b) if x!=y)


def lift_Y(Y):
    """For every candidate lifted onto Y, is Y's genome at the ortholog EXPRESSED and 1-3mm from the candidate?"""
    pool=set()
    with gzip.open(f"{POOLS}/{Y}.pool.txt.gz","rt") as fh:
        for l in fh: pool.add(l.rstrip("\n"))
    fa=pysam.FastaFile(f"{GEN}/{Y}.fa")
    best={}                                        # cand_id -> (min mm, ortholog G) where the ortholog is expressed
    with open(f"{LOCI}/present_in_{Y}.bed") as fh:
        for line in fh:
            f=line.rstrip("\n").split("\t")
            if len(f)<6: continue
            chrom,start,end,cid,strand=f[0],int(f[1]),int(f[2]),f[3],f[5]
            seq=cid.split("|")[-1]
            if end-start != len(seq): continue     # same-length substitution only
            if cid.split("|")[0]==Y: continue       # skip self (candidate's own strain)
            try: g=fa.fetch(chrom,start,end).upper()
            except Exception: continue
            if len(g)!=len(seq): continue
            G = g if strand=="+" else rc(g)
            if G not in pool: continue              # ortholog must be EXPRESSED in Y
            mm=ham(seq,G)
            if mm>3: continue
            if mm<best.get(cid,(99,))[0]: best[cid]=(mm,G)
    fa.close(); del pool
    with open(f"{WORK}/{Y}.tsv","w") as out:
        for cid,(mm,G) in best.items(): out.write(f"{cid}\t{mm}\t{G}\t{Y}\n")
    return Y, sum(1 for v in best.values() if 1<=v[0]<=3)


if __name__=="__main__":
    with ProcessPoolExecutor(max_workers=16) as ex:
        for Y,n in ex.map(lift_Y, ALL): print(f"[{Y}] {n:,} lift-confirmed SNP-variants (ortholog expressed, 1-3mm)",flush=True)
    gmin={}; gY={}; gAl={}
    for Y in ALL:
        with open(f"{WORK}/{Y}.tsv") as fh:
            for line in fh:
                cid,m_s,G,yy=line.rstrip("\n").split("\t"); m=int(m_s)
                if m<gmin.get(cid,99): gmin[cid]=m; gY[cid]=yy; gAl[cid]=G    # min across strains (0=exact excludes)
    rows=[(cid,cid.split("|")[0],gY[cid],cid.split("|")[-1],gAl[cid],m) for cid,m in gmin.items() if 1<=m<=3]
    out=pd.DataFrame(rows,columns=["cand_id","home","variant_strain","home_seq","Y_allele","mm"])
    out.to_csv(f"{U}/unique16/snp_variant_refinement.lift.csv",index=False)
    print(f"lift SNP: {len(out):,} rows; mm={out.mm.value_counts().sort_index().to_dict()}",flush=True)
    D=set(pd.read_csv(f"{U}/unique16/snp_variant_refinement.delivered_orig.csv",usecols=["cand_id"]).cand_id)
    fc=pd.read_csv(f"{U}/unique16/final_classified.csv.gz",usecols=["strain","timepoint","sequence","klass"])
    fc["cand_id"]=fc.strain+"|"+fc.timepoint+"|"+fc.sequence
    CBS=set(fc[fc.klass=="unique: conserved-but-silent"].cand_id)
    RC=set(out.cand_id)&CBS
    print(f"lift∩base-CBS = {len(RC):,}  delivered = {len(D):,}  overlap {len(RC&D):,} ({100*len(RC&D)/len(D):.1f}% of delivered)  extra {len(RC-D):,}  missed {len(D-RC):,}",flush=True)
    print(f"ADOPTION: SNP-variant 217,559 -> {len(RC):,};  CBS 86,115 -> {303674-len(RC):,};  genuinely-unique 106,961 -> {303674-len(RC)+20846:,}",flush=True)
