

####### 
#@author narumeena
#@description This script looks correct and should generate a pseudo genome from a reference genome and selected regions 
# in a BED12 file. It will add "N" characters for non-selected regions when generating the pseudo genome, and it will 
# keep the same chromosome coordinates as the original genome when generating the pseudo genome.
#python pseudognome.py --ref_genome_file ref_genome.fasta --bed_file bed_file.bed --output_file pseudogenome.fasta

######


import pybedtools
import pyfaidx
import sys
import os 
import tqdm
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate a pseudo genome from a reference genome and selected regions in a BED12 file')
parser.add_argument('--ref_genome_file', help='Path to the reference genome FASTA file')
parser.add_argument('--bed_file', help='Path to the BED12 file')
parser.add_argument('--output_file', help='Path to the output pseudo genome FASTA file')
args = parser.parse_args()

ref_genome_file = args.ref_genome_file
bed_file = args.bed_file
output_file = args.output_file

try:
    # Load the BED12 file
    bed12 = pybedtools.BedTool(bed_file)

    # Load the reference FASTA file
    fasta = pyfaidx.Fasta(ref_genome_file)

    # Extract the regions of interest from the FASTA file using the BED12 coordinates
    extracted_regions = []
    for record in bed12:
        extracted_regions.append(fasta[record.chrom][record.start:record.end])

    # Convert the non-selected regions to N
    fasta_with_Ns = []
    for chrom in tqdm.tqdm(fasta.keys()):
        seq = str(fasta[chrom])
        for record in bed12:
            if record.chrom == chrom:
                seq = seq[:record.start] + 'N'*(record.end-record.start) + seq[record.end:]
        fasta_with_Ns.append(SeqRecord(Seq(seq), id=chrom))

    # Save the modified FASTA file
    with open(output_file, 'w') as f:
        SeqIO.write(fasta_with_Ns, f, 'fasta')

except Exception as e:
    print(f'Error: {e}')
