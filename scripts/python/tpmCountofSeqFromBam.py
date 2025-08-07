####### 
#@author narumeena
#@description This script calculates the Transcripts Per Million (TPM) values for each mapping
# sequence in a BAM file. The BAM file is a binary version of a SAM file, which is a text file 
# that contains sequence alignment data.
#running the script
#python tpmCountofSeqFromBam.py --bam_file bamFile.bam --output_file outputFile.tsv

######


import pysam
import sys
import os
import argparse
import logging
from tqdm import tqdm

# Set up logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Set up argparse to handle command line arguments
parser = argparse.ArgumentParser(description='Calculate TPM values for each mapping sequence in a BAM file')
parser.add_argument('--bam_file', type=str, help='Input BAM file')
parser.add_argument('--output_file', type=str, help='Output file')
args = parser.parse_args()

# Open bam file using pysam
try:
    bamfile = pysam.AlignmentFile(args.bam_file, "rb")
    logging.info(f'Successfully opened {args.bam_file}')
except Exception as e:
    logging.error(f'Error opening {args.bam_file}: {e}')
    sys.exit(1)

# Create empty dictionary for storing tpm counts
tpm_counts = {}

# Iterate through each mapping read in bam file
total_reads = bamfile.count()
logging.info(f'Calculating TPM values for {total_reads} reads')
for i, read in enumerate(tqdm(bamfile.fetch(), total=total_reads, desc="Processing reads")):
    if not read.is_unmapped:
        # Get the read's sequence
        sequence = read.query_sequence
        # If the sequence is not already in the dictionary, add it with a count of 1
        if sequence not in tpm_counts:
            tpm_counts[sequence] = 1
        # If the sequence is already in the dictionary, increment its count by 1
        else:
            tpm_counts[sequence] += 1

# Calculate the total number of mapping reads
total_mapped_reads = sum(tpm_counts.values())

# Calculate the tpm for each sequence
for sequence, count in tpm_counts.items():
    tpm = (count/total_mapped_reads) * 1000000
    tpm_counts[sequence] = tpm

# Open output file for writing
try:
    with open(args.output_file, "w") as out:
        # Write tpm counts to output file
        for sequence, tpm in tpm_counts.items():
            out.write(f"{sequence}\t{tpm}\n")
    logging.info(f'Successfully wrote TPM values to {args.output_file}')
except Exception as e:
    logging.error(f'Error writing to {args.output_file}: {e}')
    sys.exit(1)

# Close bam file
bamfile.close()
