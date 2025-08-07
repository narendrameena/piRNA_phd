####### 
#@author narumeena
#@description This script converts a BAM file, which is a binary file format that stores nucleotide 
# sequence alignment data, into a bigWig file, which is a compact, indexed, and binary version of a bedGraph file.
#running the script
#python bigWigsFileFromBamFile.py --input_file input.bam --output_file output.bw

######


import argparse
import logging
import os
import pysam
import numpy as np
from collections import defaultdict

def bam2bw(bamfile, outfile):
    ext_len = 150

    # open BAM file
    try:
        bam = pysam.AlignmentFile(bamfile, "rb")
    except IOError as e:
        logging.error("Error opening BAM file: {}".format(e))
        return

    # create a defaultdict to store coverage data
    coverage = defaultdict(int)

    # iterate through reads in the BAM file
    try:
        for read in bam:
            # get the start and end positions of the read
            start = read.reference_start
            end = read.reference_end

            # extend the read by ext_len bases on both sides
            start = max(0, start - ext_len)
            end += ext_len

            # increment the coverage for each base in the extended read
            for i in range(start, end):
                coverage[i] += 1
    except Exception as e:
        logging.error("Error reading BAM file: {}".format(e))
        return

    # get the total number of reads
    total_reads = sum(coverage.values())

    # convert coverage to RPM
    rpm = {pos: (count / total_reads) * 1e6 for pos, count in coverage.items()}

    # write RPM data to a bigWig file
    try:
        with pysam.WigFile(outfile, "w", header={"description": "Coverage data"},reference=bam.references[0]) as wig:
            for pos, count in rpm.items():
                wig.write("chr1", pos, pos + 1, count)
    except Exception as e:
        logging.error("Error writing bigWig file: {}".format(e))
        return


# set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--input_bam", help="name of the input BAM file")
parser.add_argument("--output_bw", help="name of the output bigWig file")
args = parser.parse_args()

# convert BAM file to bigWig
try:
    bam2bw(args.input_bam, args.output_bw)
    logging.info("Successfully converted BAM file to bigWig")
except Exception as e:
    logging.error("Error converting BAM file to bigWig: {}".format(e))
