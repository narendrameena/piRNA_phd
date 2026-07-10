#!/usr/bin/env python3
"""DEFINITIVE (locus-anchored) SNP-variant producer: bowtie EXPRESSION + genomic CO-LOCATION.

Of bowtie's exhaustive 1-3mm expressed matches (build_snp_variant_bowtie.py, per-pool matches in
/tmp/snpbt/partial_*.tsv), keep ONLY those where the matched Y_allele IS the sequence at the candidate's OWN
orthologous locus in that strain -- i.e. the candidate aligns to strain Y's genome (committed cand_to_Y BAM)
and the aligned-orientation genomic sequence there == Y_allele. This is the biologically-correct criterion for
"1-3mm SNP-variant of the SAME piRNA expressed elsewhere", and it drops the coincidental pool matches (candidate
1-3mm from an expressed piRNA at a DIFFERENT locus) that pure pool-matching cannot exclude.

RESIDUAL BLIND-SPOT: a candidate whose ortholog is off-assembly or too divergent for STAR to align (no
cand_to_Y alignment) cannot be co-location-confirmed and is dropped (~14% of the pool-only extras), so this
UNDER-counts by that blind-spot -- it is the lower/locus-strict end of the SNP/CBS boundary band, complementing
the delivered (middle) and pure-pool bowtie (upper) estimates. See SNP_VARIANT_METHOD_TEST.md.

PREREQ: run build_snp_variant_bowtie.py first (produces /tmp/snpbt/partial_*.tsv). Writes
unique16/snp_variant_refinement.coloc.csv (does NOT overwrite delivered) + reports vs delivered.
"""
import os, pandas as pd, pysam
from concurrent.futures import ProcessPoolExecutor
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; S4=f"{U}/step4_16"; WORK="/tmp/snpbt"
ALL=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J",
     "NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]
_C={"A":"T","T":"A","C":"G","G":"C","N":"N"}
def rc(s): return "".join(_C.get(b,"N") for b in reversed(s))
SEQ=[]; HOME=[]; CID=[]; TPS=[]
for X in ALL:
    c=pd.read_csv(f"{S4}/{X}.candidates16.tsv.gz")
    SEQ+=c.sequence.tolist(); HOME+=[X]*len(c); CID+=c.id.tolist(); TPS+=c.timepoints.astype(str).tolist()
N=len(SEQ)


def coloc_Y(Y):
    """Keep pool-Y bowtie matches whose Y_allele == the candidate's own aligned-orientation genomic seq in Y."""
    part={}
    with open(f"{WORK}/partial_{Y}.tsv") as f:
        for line in f:
            gi_s,m_s,al=line.rstrip("\n").split("\t"); part[int(gi_s)]=(int(m_s),al)
    byhome={}
    for gi in part: byhome.setdefault(HOME[gi],[]).append(gi)
    out={}
    for home,gis in byhome.items():
        if home==Y: continue
        want={CID[gi]:gi for gi in gis}
        genom={}
        try: bam=pysam.AlignmentFile(f"{S4}/{home}.cand_to_{Y}.Aligned.sortedByCoord.out.bam","rb")
        except Exception: continue
        for a in bam.fetch(until_eof=True):
            if a.query_name in want and not a.is_unmapped:
                try: ref=a.get_reference_sequence().upper()
                except Exception: continue
                genom.setdefault(a.query_name,set()).add(rc(ref) if a.is_reverse else ref)
        bam.close()
        for cid,gi in want.items():
            mm,al=part[gi]
            if al in genom.get(cid,()): out[gi]=(mm,al)      # co-located: expressed match IS the candidate's own-locus genome
    with open(f"{WORK}/coloc_{Y}.tsv","w") as f:
        for gi,(mm,al) in out.items(): f.write(f"{gi}\t{mm}\t{al}\t{Y}\n")
    return Y, len(out)


if __name__=="__main__":
    if not os.path.exists(f"{WORK}/partial_{ALL[0]}.tsv"):
        raise SystemExit("run build_snp_variant_bowtie.py first (need /tmp/snpbt/partial_*.tsv)")
    with ProcessPoolExecutor(max_workers=16) as ex:
        for Y,n in ex.map(coloc_Y, ALL): print(f"[{Y}] {n:,} co-located matches",flush=True)
    gmin=[99]*N; gY=[None]*N; gAl=[None]*N
    for Y in ALL:
        with open(f"{WORK}/coloc_{Y}.tsv") as f:
            for line in f:
                gi_s,m_s,al,yy=line.rstrip("\n").split("\t"); gi=int(gi_s); m=int(m_s)
                if m<gmin[gi]: gmin[gi]=m; gY[gi]=yy; gAl[gi]=al   # min across strains (0=exact excluded below)
    rows=[]
    for gi in range(N):
        if 1<=gmin[gi]<=3:
            for tp in TPS[gi].replace(";",",").split(","):
                tp=tp.strip()
                if tp and tp!="nan": rows.append((f"{HOME[gi]}|{tp}|{SEQ[gi]}",HOME[gi],gY[gi],SEQ[gi],gAl[gi],gmin[gi]))
    out=pd.DataFrame(rows,columns=["cand_id","home","variant_strain","home_seq","Y_allele","mm"])
    out.to_csv(f"{U}/unique16/snp_variant_refinement.coloc.csv",index=False)
    print(f"co-located SNP: {len(out):,} rows, {out.cand_id.nunique():,} uniq; mm={out.mm.value_counts().sort_index().to_dict()}",flush=True)
    D=set(pd.read_csv(f"{U}/unique16/snp_variant_refinement.delivered_orig.csv",usecols=["cand_id"]).cand_id)
    fc=pd.read_csv(f"{U}/unique16/final_classified.csv.gz",usecols=["strain","timepoint","sequence","klass"])
    fc["cand_id"]=fc.strain+"|"+fc.timepoint+"|"+fc.sequence
    CBS=set(fc[fc.klass=="unique: conserved-but-silent"].cand_id)
    RC=set(out.cand_id)&CBS
    print(f"coloc∩base-CBS = {len(RC):,}  delivered = {len(D):,}  overlap {len(RC&D):,} ({100*len(RC&D)/len(D):.1f}% of delivered)  extra {len(RC-D):,}  missed {len(D-RC):,}",flush=True)
    print(f"ADOPTION: SNP-variant 217,559 -> {len(RC):,};  CBS 86,115 -> {303674-len(RC):,};  genuinely-unique 106,961 -> {303674-len(RC)+20846:,}",flush=True)
