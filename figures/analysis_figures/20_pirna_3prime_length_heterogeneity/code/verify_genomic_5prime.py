"""Is the sequence-prefix raggedness biologically real? Re-anchor on the 5' GENOMIC POSITION (cand_self16):
a genuine 3'-trim isoform shares the SAME 5' genomic coordinate (same locus) but a different length.
Compares (a) sequence-prefix raggedness, (b) genomic-5'-anchored raggedness, (c) of the sequence-prefix isoforms,
the fraction that ALSO share the 5' genomic position (concordance = how often the prefix match is a true same-locus trim)."""
import subprocess, pandas as pd, numpy as np
U="/mnt/home3/miska/nm667/scratch/inProgress/mice_PiRNA/analysis/claude_biomni_analysis/unique_pirna"
ST="/mnt/home3/miska/nm667/miniconda3/envs/ccTE/bin/samtools"
clean=pd.read_csv(f"{U}/unique16/final_classified_clean_2read.csv.gz",usecols=["sequence","strain","timepoint"])
clean["L"]=clean.sequence.str.len(); clean["cand_id"]=clean.strain+"|"+clean.timepoint+"|"+clean.sequence
WIN={"16.5dpc":[27],"12.5dpp":[27,30],"20.5dpp":[30]}; TPN={"16.5dpc":"E16.5","12.5dpp":"P12.5","20.5dpp":"P20.5"}
need=dict(zip(clean.cand_id,clean.sequence)); rows=[]
for X in sorted(clean.strain.unique()):
    out=subprocess.run([ST,"view",f"{U}/cand_self16/{X}.cand_self16.bam"],capture_output=True,text=True,timeout=900).stdout
    for ln in out.splitlines():
        f=ln.split("\t")
        if len(f)<6 or int(f[1]) not in (0,16) or f[0] not in need: continue
        seq=f[0].split("|")[2]; L=len(seq); pos=int(f[3]); strand="+" if int(f[1])==0 else "-"
        five=pos if strand=="+" else pos+L-1
        rows.append((f[0],f[0].split("|")[0],f[0].split("|")[1],seq,L,f[2].split("#")[-1],five,strand))
c=pd.DataFrame(rows,columns=["cand_id","strain","tp","seq","L","chrom","five","strand"]).drop_duplicates("cand_id")
print(f"coords for {len(c):,} of {len(need):,} clean candidates ({100*len(c)/len(need):.0f}%)")
for tp in ["16.5dpc","12.5dpp","20.5dpp"]:
    sub=c[c.tp==tp]; nlen=sub.groupby(["strain","chrom","five","strand"]).L.nunique().rename("nlen")
    bylen={L:set(sub[sub.L==L].seq) for L in range(23,35)}
    for L in WIN[tp]:
        win=sub[sub.L==L].merge(nlen,on=["strain","chrom","five","strand"])
        gen=100*(win.nlen>1).mean()                      # (b) genomic 5'-anchored raggedness
        pre={k:{y[:L] for y in bylen[L+k]} for k in (1,2,3)}
        win["seqrag"]=win.seq.apply(lambda x: any(x[:L-k] in bylen[L-k] for k in (1,2,3)) or any(x in pre[k] for k in (1,2,3)))
        seqr=100*win.seqrag.mean()                        # (a) sequence-prefix raggedness
        # (c) concordance: among seq-prefix-ragged window piRNAs, do they sit at a multi-length 5' position?
        conc=100*(win[win.seqrag].nlen>1).mean()
        print(f"{TPN[tp]} {L}nt (n={len(win):,}): seq-prefix={seqr:.0f}%  genomic-5'-anchored={gen:.0f}%  concordance(prefix-hits at multi-length 5' locus)={conc:.0f}%")
