#!/usr/bin/env python3
"""Committed producer for unique16/snp_variant_refinement.csv via DIRECT expressed-pool matching (bowtie1 -v3).

THE CORRECT METHOD (see 07_unique_piRNA_identification/SNP_VARIANT_METHOD_TEST.md). A candidate is a SNP-variant
iff it is within 1-3 SUBSTITUTIONS (Hamming, same length) of a piRNA EXPRESSED in ANOTHER strain (that strain's
expressed pool). This replaces the classify_step416.py GENOMIC-alignment proxy, which asks "does it align to a
genomic locus that is expressed" (only ~50-86% reproducible; misses candidates that don't align to the genome
or whose genomic sequence != the expressed allele). Head-to-head: direct pool-search = 100% recall of the
delivered set; genome proxy = 84.5%.

METHOD: build ONE bowtie index over all candidates; stream each of the 16 expressed pools through it once
(`bowtie --norc -v3 -a`), keeping same-length forward alignments; per candidate track the MIN mismatch across
OTHER strains' pools (exclude self). min==0 = expressed-exact (excluded); min in 1-3 = SNP-variant. Substitution-
only same-length Hamming (delivered rows are 100% same-length with hamming==mm; SNP != indel/length-isoform).

Writes unique16/snp_variant_refinement.bowtie.csv (does NOT overwrite the delivered file) and reports
reproduction vs delivered. To adopt: cp over snp_variant_refinement.csv, then re-run make_klass5.py.
Usage: build_snp_variant_bowtie.py
"""
import os, subprocess
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
ROOT="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA"
U=f"{ROOT}/analysis/claude_biomni_analysis/unique_pirna"; S4=f"{U}/step4_16"; POOLS=f"{U}/unique16/pools"
BOWTIE="/usr/bin/bowtie"; BUILD="/usr/bin/bowtie-build"
WORK="/tmp/snpbt"; os.makedirs(WORK, exist_ok=True)
IDX=f"{WORK}/all_cand_idx"
ALL=["C57BL_6NJ","BALB_cJ","A_J","FVB_NJ","C3H_HeJ","LP_J","129S1_SvImJ","DBA_2J","AKR_J","CBA_J",
     "NZO_HlLtJ","NOD_ShiLtJ","WSB_EiJ","CAST_EiJ","PWK_PhJ","SPRET_EiJ"]

# ---- all candidates -> globals + combined FASTA/index ----
SEQ=[]; HOME=[]; TPS=[]
for X in ALL:
    c=pd.read_csv(f"{S4}/{X}.candidates16.tsv.gz")
    SEQ+=c.sequence.tolist(); HOME+=[X]*len(c); TPS+=c.timepoints.astype(str).tolist()
N=len(SEQ); LEN=[len(s) for s in SEQ]


def run_Y(Y):
    """Stream Y's pool through bowtie -v3 -a; return per-candidate (min_mm, Y_allele) for candidates home!=Y."""
    minmm={}; allele={}
    p1=subprocess.Popen(["zcat",f"{POOLS}/{Y}.pool.txt.gz"],stdout=subprocess.PIPE)
    p2=subprocess.Popen([BOWTIE,"--norc","-v","3","-a","-r","-p","6",IDX,"-"],
                        stdin=p1.stdout,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,text=True,bufsize=1<<20)
    p1.stdout.close()
    for line in p2.stdout:
        F=line.rstrip("\n").split("\t")
        if len(F)<6: continue
        gi=int(F[2]); off=int(F[3]); rseq=F[4]
        if off!=0 or len(rseq)!=LEN[gi] or HOME[gi]==Y: continue
        mmdesc=F[7] if len(F)>7 else ""
        mm=0 if mmdesc=="" else mmdesc.count(",")+1
        if mm<minmm.get(gi,99): minmm[gi]=mm; allele[gi]=rseq
    p2.wait(); p1.wait()
    out=f"{WORK}/partial_{Y}.tsv"
    with open(out,"w") as f:
        for gi,m in minmm.items(): f.write(f"{gi}\t{m}\t{allele[gi]}\n")
    return Y, len(minmm)


if __name__=="__main__":
    if not os.path.exists(IDX+".1.ebwt"):
        fa=f"{WORK}/all_cand.fa"
        with open(fa,"w") as f:
            for i,s in enumerate(SEQ): f.write(f">{i}\n{s}\n")
        subprocess.run([BUILD,"--quiet","--threads","8",fa,IDX],check=True)
    print(f"index over {N:,} candidates ready", flush=True)
    with ProcessPoolExecutor(max_workers=16) as ex:
        for Y,n in ex.map(run_Y, ALL):
            print(f"[{Y}] {n:,} candidates matched (<=3, same-length, non-self)", flush=True)
    gmin=[99]*N; gY=[None]*N; gAl=[None]*N
    for Y in ALL:                                             # ALL order -> first strain at the global min wins ties
        with open(f"{WORK}/partial_{Y}.tsv") as f:
            for line in f:
                gi_s,m_s,al=line.rstrip("\n").split("\t"); gi=int(gi_s); m=int(m_s)
                if m<gmin[gi]: gmin[gi]=m; gY[gi]=Y; gAl[gi]=al
    rows=[]
    for gi in range(N):
        if 1<=gmin[gi]<=3:
            for tp in TPS[gi].replace(";",",").split(","):
                tp=tp.strip()
                if tp and tp!="nan":
                    rows.append((f"{HOME[gi]}|{tp}|{SEQ[gi]}",HOME[gi],gY[gi],SEQ[gi],gAl[gi],gmin[gi]))
    out=pd.DataFrame(rows,columns=["cand_id","home","variant_strain","home_seq","Y_allele","mm"])
    outcsv=f"{U}/unique16/snp_variant_refinement.bowtie.csv"; out.to_csv(outcsv,index=False)
    print(f"wrote {outcsv}: {len(out):,} rows, {out.cand_id.nunique():,} uniq cand_id; mm={out.mm.value_counts().sort_index().to_dict()}", flush=True)
    dp=f"{U}/unique16/snp_variant_refinement.delivered_orig.csv"
    dp=dp if os.path.exists(dp) else f"{U}/unique16/snp_variant_refinement.csv"
    D=set(pd.read_csv(dp,usecols=["cand_id"]).cand_id); R=set(out.cand_id)
    print(f"REPRODUCTION vs delivered: delivered={len(D):,} bowtie={len(R):,} overlap={len(R&D):,} "
          f"({100*len(R&D)/len(D):.1f}% of delivered); bowtie-extra={len(R-D):,} delivered-missed={len(D-R):,}", flush=True)
