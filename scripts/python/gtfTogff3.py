

#!/usr/bin/env python

####### 
#@author narumeena
#@description This script will convert gtf to gff3 format
#python gtf_to_gff3.py --gtf_file input.gtf --gff3_file output.gff3

######


import argparse
import logging
import sys

# Set up the command line arguments
parser = argparse.ArgumentParser(description='Convert GTF file to GFF3 format')
parser.add_argument('--gtf_file', type=argparse.FileType('r'), help='the GTF file to convert')
parser.add_argument('--gff3_file', type=argparse.FileType('w'), help='the GFF3 file to write')

# Parse the command line arguments
args = parser.parse_args()

# Set up logging
logging.basicConfig(level=logging.INFO)

# Convert the GTF file to GFF3 format
with args.gtf_file, args.gff3_file:
  for line in args.gtf_file:
    cols = line.strip().split('\t')
    seqid = cols[0]
    source = cols[1]
    feature = cols[2]
    start = cols[3]
    end = cols[4]
    score = cols[5]
    strand = cols[6]
    phase = cols[7]
    attributes = {}
    for attr in cols[8].split(';'):
      attr = attr.strip()
      if not attr:
        continue
      key, value = attr.split()
      attributes[key] = value.strip('"')
    args.gff3_file.write('\t'.join([seqid, source, feature, start, end, score, strand, phase, attributes['gene_id']]) + '\n')

logging.info('Finished converting GTF file to GFF3 format')
