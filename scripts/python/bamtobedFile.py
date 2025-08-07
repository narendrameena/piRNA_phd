#!/usr/bin/env python

####### 
#@author narumeena
#@description This script will read the BAM file input.bam and write a BED12 file output.bed containing the reads in the BAM file.
#python bamtobedFile.py --input_file input.bam --output_file output.bw

######

import argparse
import logging
import pysam
from tqdm import tqdm

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--input_bam", help="input BAM file")
parser.add_argument("--output_bed", help="output BED12 file")
args = parser.parse_args()

# Open the BAM file
try:
    bam_file = pysam.AlignmentFile(args.input_bam, "rb")
except Exception as e:
    logging.error(f"Error opening BAM file: {e}")
    exit(1)

# Open a BED12 file for writing
try:
    bed_file = open(args.output_bed, "w")
except Exception as e:
    logging.error(f"Error opening BED12 file: {e}")
    bam_file.close()
    exit(1)

# Iterate over the reads in the BAM file
for read in tqdm(bam_file, desc="Converting BAM to BED12"):
    # Check if the read is mapped
    if not read.is_unmapped:
        # Get the chromosome, start, and end positions of the read
        chrom = bam_file.getrname(read.reference_id)
        start = read.reference_start
        end = read.reference_end
        # Get the name, score, and strand of the read
        name = read.query_name
        score = read.mapping_quality
        strand = "+" if read.is_reverse else "-"
        # Get the blocks (exons) of the read
        blocks = read.get_blocks()
        # Calculate the block sizes and starts
        block_sizes = ",".join([str(b[1] - b[0]) for b in blocks])
        block_starts = ",".join([str(b[0] - start) for b in blocks])
        # Write the read to the BED12 file
        bed_file.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{start}\t{end}\t0\t{len(blocks)}\t{block_sizes}\t{block_starts}\n")

# Close the BAM and BED12 files
bam_file.close()
bed_file.close()
