#!/usr/bin/env python

#######
#@author narumeena
#@description Python script to process bam file (mapped to mouse genome ) with a window size of 20kb 
# and steps 1kb on all chromosomes. Filtered out windows with fewer than 200 distinct counts read and 
# with more than 100  RPM counts. Only considered sequences with 24–32 lengths. Adjust the script to 
# handle cases where the window size or step size is not exact multiple of the length of the. And write 
# the filtred windows to a tab-limited file with TPM count.
#running the script
#python slidingWindowOnBam.py --input_file bamFile.bam --output_file outputFile.bed

######

import pysam
import os
import sys
import argparse
from tqdm import tqdm
import numpy as np
import logging
from collections import defaultdict


# Set up argument parsing
parser = argparse.ArgumentParser(description='Process BAM file to find TPM values.')
parser.add_argument('--input_file', help='input BAM file', required=True)
parser.add_argument('--output_file', help='output file', required=True)


# Parse command-line arguments
args = parser.parse_args()
input_file = args.input_file
output_file = args.output_file


# Constants
WIN_SIZE = 500000
STEP_SIZE = 20000
MIN_READ_COUNT = 20
MAX_TPM = 20
MIN_SEQ_LENGTH = 24
MAX_SEQ_LENGTH = 32



# Define a read callback function
def is_mapped(read):
    return read.is_unmapped == False

if WIN_SIZE <= 0:
    logging.error("Window size must be greater than zero")
    sys.exit(1)

if STEP_SIZE <= 0:
    logging.error("Step size must be greater than zero")
    sys.exit(1)

if not os.path.isfile(input_file):
    logging.error(f"Error: input file '{input_file}' does not exist")
    sys.exit(1)

# Open the BAM file
try:
    bamfile = pysam.AlignmentFile(input_file , "rb")
except pysam.SamtoolsError as e:
    logging.error(f'Error opening BAM file: {e}')
    sys.exit(1)

# Create a list of chromosomes in the BAM file
try:
    chroms = bamfile.references
except pysam.SamtoolsError as e:
    logging.error(f'Error getting chromosome list: {e}')
    sys.exit(1)


# Initialize a counter for the total number of mapping reads
total_mapped_reads = 0

# Iterate over the reads in the BAM file
for read in bamfile:
    # Check if the read is mapped (i.e., if it has a valid reference position)
    if read.reference_id != -1:
        # If the read is mapped, increment the counter
        total_mapped_reads += 1
print(total_mapped_reads)


# Open the output file
if os.path.isfile(output_file):
    # If the output file already exists, create a new empty file
    with open(output_file, "w") as outfile:
        pass
else:
    # If the output file doesn't exist, create it
    with open(output_file, "a+") as outfile:
        pass



# Open the output file
try:
    with open(output_file, "a+") as outfile:
        # Write the header line
        #outfile.write("chrom\twin_start\twin_end\tTPM\tdistinct_count\n")

        # Iterate over chromosomes
        for chrom in chroms:            
            # Get the length of the chromosome
            chrom_length = bamfile.get_reference_length(chrom)

            # Skip the chromosome if the window size is larger than the chromosome length
            if chrom_length < WIN_SIZE:
              logging.warning(f"Warning: window size ({WIN_SIZE}) is larger than chromosome length ({chrom_length}), skipping chromosome")
              continue

            # Calculate the number of windows
            num_wins = (chrom_length - WIN_SIZE) // STEP_SIZE + 1
            if num_wins * STEP_SIZE + WIN_SIZE > chrom_length:
              num_wins -= 1

            if num_wins < 0:
                logging.error("Number of windows must be greater than or equal to zero")
                sys.exit(1)

            # Iterate over windows
            for i in tqdm(range(num_wins), desc=f'Processing chromosome {chrom}'):
                # Calculate the start and end positions of the window
                win_start = i * STEP_SIZE
                win_end = win_start + WIN_SIZE

                # Skip the window if it extends beyond the end of the chromosome
                if win_end > chrom_length:
                    continue

                
                # Iterate over reads in the window
                reads_count = 0
                read_uniq = set()
                for read in bamfile.fetch(chrom, win_start, win_end):
                    # Skip the read if it doesn't meet the criteria for distinct counts, TPM counts, and sequence length
                    #if read.mapping_quality == 0 or read.query_length < MIN_SEQ_LENGTH or read.query_length > MAX_SEQ_LENGTH:
                    #    continue
                    reads_count += 1
                    read_uniq.add(read.query_sequence)
                
                # Skip the window if it has fewer than the minimum number of distinct counts
                #if len(read_uniq) < MIN_READ_COUNT:
                #    continue

                # Calculate the RPM (read per million) value for the window
                rpm = (reads_count / total_mapped_reads)
                tpm = rpm * 1,000,000 
                # Skip the window if the rpm value is too high
                #if tpm < MAX_TPM:
                #    continue
                
                # Write the results to the output file
                # outfile.write(f"{chrom}\t{win_start}\t{win_end}\t{tpm:.2f}\t{len(reads)}\n")
                # Write the window data to the output file in BED12 format
                #print(f"{chrom}\t{win_start}\t{win_end}\t{len(read_uniq)}\t{tpm[0]:.4f}\t+\t{win_start}\t{win_end}\t0,0,0\t1\t{WIN_SIZE}\t0\n")
                outfile.write(f"{chrom}\t{win_start}\t{win_end}\t{len(read_uniq)}\t{tpm[0]:.4f}\t+\t{win_start}\t{win_end}\t0,0,0\t{total_mapped_reads}\t{reads_count}\t0\n")


    # Close the file
    outfile.close()

except Exception as e:
  logging.error(f"Error processing BAM file: {e}")
  sys.exit(1)

# Close the BAM file
bamfile.close()

# Close the file
outfile.close()
