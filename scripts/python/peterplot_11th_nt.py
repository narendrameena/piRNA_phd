def parsearguments():
    import argparse,sys

    description = """=========================== Peterplot generator ============================"""
    epilog = """============================== Alea jacta est =============================="""
    parser = argparse.ArgumentParser(description=description, epilog=epilog,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('SEQFILE',help="Path to the fastq/collapsed fasta/bam file containing the reads")
    parser.add_argument('OUTFNAME',help="Output file name")
    parser.add_argument('--outraw',metavar='filename',help="Raw counts output file name")
    parser.add_argument('--minlength',type=int,default=18,metavar='N',help="Minimum length to report on the x axis")
    parser.add_argument('--maxlength',type=int,default=33,metavar='N',help="Maximum length to report on the x axis")
    parser.add_argument('--ylims',metavar='N,N',help="Limits for the y axis, comma separated (lower,upper)")
    parser.add_argument('--title',default="",help="Title for the plot")
    parser.add_argument('--uncollapse',action="store_true",help="Uncollapse the reads (requires fasta headers formatted as 'id-count'")
    parser.add_argument('--normfw',type=float,help="Normalisation factor for the forward strand",default=None)
    parser.add_argument('--normrv',type=float,help="Normalisation factor for the reverse strand",default=None)
    parser.add_argument('--normnreads',action="store_true",help="Normalise reads to library size (or mapped reads in case of BAM)")
    parser.add_argument('--stranded', help="Generate a stranded plot", action="store_true")
    parser.add_argument('--fwonly', help="Plot only forward strand", action="store_true")
    parser.add_argument('--rvonly', help="Plot only reverse strand", action="store_true")
    parser.add_argument('--noN', help="Don't display the counts for the 'N' nucleotide", action="store_true")
    parser.add_argument('--3p', help="Count the 3' nucleotides instead of the 5'", action="store_true")
    parser.add_argument('--chrnames', help="Count reads for the specified 'chromosome' only (BAM exclusive)", nargs="*")
    parser.add_argument('--nmaps', help="Maximum number of multi maps a read can have and still be counted (BAM exclusive)", type=int, default=1)

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    return(vars(parser.parse_args()))

def writeplot(OUTFNAME,fig):
    fig.savefig(OUTFNAME)

def plot(counts,minlength=18,maxlength=33,title=None,ylims=None,stranded=False,noN=True,fwonly=False,rvonly=False,colours=None,isrna=True,legend=True):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    import copy

    if not title:
        title = ""

    counts = copy.deepcopy(counts)

    if fwonly or rvonly:
        s = "rv" if fwonly else "fw"
        for l in counts[s]:
            for n in counts[s][l]:
                counts[s][l][n] = 0

    lengths = np.arange(minlength,maxlength+1)
    if stranded:
        lengths_rv = np.arange(minlength,maxlength+1)

    df = pd.DataFrame(counts["fw"])
    df.fillna(0,inplace=True)
    if stranded:
        df_rv = pd.DataFrame(counts["rv"]) * -1
        df_rv.fillna(0,inplace=True)

    fig, ax = plt.subplots()
    ax.set_title(title)
    if colours:
        nucs = [["A",colours["A"]],["C",colours["C"]],["G",colours["G"]],["T",colours["T"]]]
    else:
        nucs = [["A","red"],["C","blue"],["T","green"],["G","purple"]]
    if not noN:
        nucs += ["N","cyan"]

    bottom = [0] * len(df.columns)
    bottom_rv = [0] * len(df.columns)
    idx = np.arange(len(df.columns))
    width = 0.8

    for n in nucs:
        if n[0] in df.index:
            label = n[0] if not isrna or n[0] != "T" else "U"
            ax.bar(idx, df.loc[n[0]], color=n[1], bottom=bottom, width=width,label=label)
            bottom += df.loc[n[0]]
            if stranded:
                ax.bar(idx, df_rv.loc[n[0]], color=n[1], bottom=bottom_rv, width=width)
                bottom_rv += df_rv.loc[n[0]]

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(axis='both', direction='out')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    plt.xticks(idx+width/2.,lengths)
    plt.xlim([0,len(lengths)])

    if ylims:
        plt.ylim(ylims)
    if legend:
        plt.legend()

    return(fig)

def getcounts(SEQFILE,minlength=18,maxlength=33,uncollapse=False,normfw=None,normrv=None,normnreads=False,stranded=False,noN=False,c3p=False,chrnames=None,nmaps=1,t2u=False,rvonly=False,fwonly=False):

    from tstk.io import openfastx
    from Bio import SeqIO
    from Bio.Seq import Seq
    import numpy as np
    import sys,os,re

    fpath = SEQFILE
    ext = os.path.splitext(fpath)[1]

    lengths = np.arange(minlength,maxlength+1)
    if stranded:
        lengths_rv = np.arange(minlength,maxlength+1)
    counts = {}

    if noN:
        counts = {key: {"A":0,"C":0,"T":0,"G":0} for key in lengths}
        counts_rv = {key: {"A":0,"C":0,"T":0,"G":0} for key in lengths}
    else:
        counts = {key: {"A":0,"C":0,"T":0,"G":0,"N":0} for key in lengths}
        counts_rv = {key: {"A":0,"C":0,"T":0,"G":0,"N":0} for key in lengths}

    if c3p:
        nucpos = -1
    else:
        nucpos = 10 #0 index based also add 8 and 10

    nrecords = 0
    if ext == ".bam":
        import pysam
        bamfile = pysam.AlignmentFile(fpath)
        if chrnames:
            reads = []
            for chrname in chrnames:
                reads += bamfile.fetch(chrname)
        else:
            reads = bamfile.fetch()

        for read in reads:
            tags = dict(read.get_tags())
            if "NH" not in tags or tags["NH"] <= nmaps:
                if read.is_reverse:
                    seq = Seq(read.seq)
                    seq = seq.reverse_complement()
                else:
                    seq = read.seq

                firstnuc = seq[nucpos]
                if not noN or firstnuc != "N":
                    if len(seq) in lengths:
                        if stranded and read.is_reverse:
                                counts_rv[len(seq)][firstnuc] += 1
                        else:
                                counts[len(seq)][firstnuc] += 1

        if normnreads:
            if chrnames:
                nrecords = 0 
                for chrname in chrnames:
                    nrecords += sum(1 for r in pysam.AlignmentFile(fpath).fetch(chrname) if dict(r.get_tags())["HI"] == 1)
            else:
                nrecords = sum(1 for r in pysam.AlignmentFile(fpath) if dict(r.get_tags())["HI"] == 1)
    else:
        checkheader = True

        fpathnoext,fext,ftype,fh = openfastx(fpath)

        for record in SeqIO.parse(fh,ftype):
            if uncollapse:
                if checkheader:
                    reid = re.compile("^\d+-\d+$")
                    if not reid.match(record.id):
                        print("WARNING: the fasta header doesn't seem to match the 'id-count' format. Ignoring uncollapse request.")
                        uncollapse = False
                    checkheader = False

                if ftype == "fastq":
                    print("WARNING: file format is fastq, ignoring uncollapse request")
                    uncollapse = False
                    count = 1
                else:
                    count = int(record.id.split("-")[1])
            else:
                count = 1

            nrecords += count

            firstnuc = record.seq[nucpos]
            if not noN or firstnuc != "N":
                if len(record.seq) in lengths:
                    counts[len(record.seq)][firstnuc] += count

    normfactor_fw = None
    normfactor_rv = None

    if normnreads:
        normfactor_fw = nrecords
        normfactor_rv = nrecords

    if normfw:
        if normnreads:
            print("WARNING: normalisation to library size specified together with specific normalisation factor. Using the specific factor for FW reads")
        normfactor_fw = normfw

    if normrv:
        if normnreads:
            print("WARNING: normalisation to library size specified together with specific normalisation factor. Using the specific factor for RV reads")
        normfactor_rv = normrv

    if normfactor_fw: 
        for l in counts:
            for n in counts[l]:
                counts[l][n] = counts[l][n] / normfactor_fw

    if stranded and normfactor_rv:
        for l in counts_rv:
            for n in counts_rv[l]:
                counts_rv[l][n] = counts_rv[l][n] / normfactor_rv

    try:
        fh.close()
    except NameError:
        pass

    if t2u:
        for l in counts:
            counts[l]["U"] = counts[l].pop("T")
        for l in counts_rv:
            counts_rv[l]["U"] = counts_rv[l].pop("T")

    counts = {"fw":counts,"rv":counts_rv}

    return counts

def writeraw(fname, counts, stranded):
    if stranded:
        strands=["fw","rv"]
        pos_template="{l}_{s}\t"
    else:
        strands=["fw"]
        pos_template="{l}\t"

    some_length = next(iter(counts["fw"]))
    nucs = counts["fw"][some_length].keys()

    with open(fname, "w") as ofh:
        ofh.write("\t" + "\t".join(sorted(nucs)) + "\n")
        for s in strands:
            for l in sorted(counts[s]):
                ofh.write(pos_template.format(l=l,s=s) + "\t".join([str(counts[s][l][n]) for n in sorted(counts[s][l])]) + "\n")

def main():
    import os
    args = parsearguments()

    if args["uncollapse"]:
        args["OUTFNAME"] = args["OUTFNAME"].replace(".collapsed","")

    outformats = ['svg','png','jpg','eps','raw','ps','pdf','tif','pgf']

    outdir,outfname = os.path.split(args["OUTFNAME"])
    outfnamenoext,outfext = os.path.splitext(outfname)
    outformat = outfext[1:]
    if outformat not in outformats:
        raise ValueError("Output format not recognised form file extension '{}'. Please use one of {}".format(outformat,outformats))

    title = args["title"] if args["title"] else os.path.basename(args["OUTFNAME"])
    ylims = [float(lim) for lim in args["ylims"].split(",")] if args["ylims"] else None

    counts = getcounts(args["SEQFILE"],minlength=args["minlength"],maxlength=args["maxlength"],uncollapse=args["uncollapse"],normfw=args["normfw"],normrv=args["normrv"],normnreads=args["normnreads"],stranded=args["stranded"],noN=args["noN"],c3p=args["3p"],chrnames=args["chrnames"],nmaps=args["nmaps"],fwonly=args["fwonly"],rvonly=args["rvonly"])
    f = plot(counts,args["minlength"],args["maxlength"],title=title,ylims=ylims,stranded=args["stranded"],fwonly=args["fwonly"],rvonly=args["rvonly"])
    writeplot(args["OUTFNAME"],f)
    if args["outraw"]:
        writeraw(args["outraw"],counts,stranded=args["stranded"])

if __name__ == "__main__":
    import os

    main()