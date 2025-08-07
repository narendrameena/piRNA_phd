
####### 
#@author narumeena
#@description This script reads in two GFF3 (General Feature Format version 3) files as command line arguments, 
# one with annotations to borrow (borrow_lines) and one to add the borrowed annotations to (lines).
#python script.py input1.gff3 input2.gff3

######


import sys
import logging
import argparse

def overlap(start1, end1, start2, end2):
    """Calculate the overlap between two intervals [start1, end1] and [start2, end2]."""
    return max(0, min(end1, end2) - max(start1, start2))

def percent_overlap(start1, end1, start2, end2):
    """Calculate the percent overlap between two intervals [start1, end1] and [start2, end2]."""
    overlap_length = overlap(start1, end1, start2, end2)
    total_length = max(end1, end2) - min(start1, start2)
    return overlap_length / total_length

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Set up argument parser
parser = argparse.ArgumentParser()
parser.add_argument('input1', help='GFF3 file with annotations to borrow')
parser.add_argument('input2', help='GFF3 file to add borrowed annotations to')
parser.add_argument('-o', '--output', default='output.gff3', help='output file name (default: output.gff3)')
parser.add_argument('-p', '--percent', type=float, default=0.6, help='overlap percentage threshold (default: 0.6)')
args = parser.parse_args()

# Read the GFF3 file with the annotations we want to borrow
try:
    with open(args.input1, 'r') as f:
        borrow_lines = f.readlines()
except FileNotFoundError:
    logging.error(f'Error: Could not find file {args.input1}')
    sys.exit(1)

# Read the GFF3 file we want to add the annotations to
try:
    with open(args.input2, 'r') as f:
        lines = f.readlines()
except FileNotFoundError:
    logging.error(f'Error: Could not find file {args.input2}')
    sys.exit(1)

# Initialize a counter for the borrow_lines
borrow_line_counter = 0

# Initialize a list to store the modified lines
modified_lines = []

# Iterate over the lines in the GFF3 file we want to add annotations to
for i in range(len(lines)):
    line = lines[i]
    # If the line is a feature line (e.g. exon, CDS, etc.)
    if not line.startswith('#'):
        # Split the line into fields
        fields = line.strip().split('\t')
        # Get the type of the feature
        feature_type = fields[2]
        # Get the start and end positions of the feature
        start = int(fields[3])
        end = int(fields[4])
        # Get the next line in the borrow_lines
        borrow_line = borrow_lines[borrow_line_counter]
        # Split the borrow_line into fields
        borrow_fields = borrow_line.strip().split('\t')
        # Get the type of the borrow_line feature
        borrow_feature_type = borrow_fields[2]
        # Get the start and end positions of the borrow_line feature
        borrow_start = int(borrow_fields[3])
        borrow_end = int(borrow_fields[4])
        # If the feature types are the same and the borrow_line feature overlaps with the input feature by at least the specified percentage
        if feature_type == borrow_feature_type and percent_overlap(start, end, borrow_start, borrow_end) >= args.percent:
            # Append the line with the borrowed annotation to the modified_lines list
            modified_lines.append(line.strip() + '\t' + borrow_fields[8])
            # Increment the borrow_line_counter
            borrow_line_counter += 1
        else:
            # Append the original line to the modified_lines list
            modified_lines.append(line.strip())
    else:
        # Append the original line to the modified_lines list
        modified_lines.append(line.strip())

# Write the modified lines to the output file
with open(args.output, 'w') as f:
    f.write('\n'.join(modified_lines))

logging.info('Done!')

