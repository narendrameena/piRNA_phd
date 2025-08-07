import pysam
import numpy as np

bam = pysam.AlignmentFile(snakemake.input["bam"],'rb')
log = snakemake.input["log"]

with open(log) as fh:
    for r in fh:
        if "Uniquely mapped reads number" in r:
            N = int(r.split()[-1])

counts = []
lengths = []

for i,ref in enumerate(bam.references):
    counts.append(bam.count(ref))
    lengths.append(bam.get_reference_length(ref))

C = np.array(counts)
L = np.array(lengths)

rpms = (10**6 * C) / N
rpkms = (10**9 * C) / (N * L)

with open(snakemake.output[0],"w") as ofh:
    ofh.write("ref,count,length,rpm,rpkm\n")
    for a,b,c,d,e in zip(bam.references,counts,lengths,rpms,rpkms):
        ofh.write(",".join([str(s) for s in [a,b,c,d,e]]) + "\n")

bam.close()
