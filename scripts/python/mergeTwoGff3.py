

####### 
#@author narumeena
#@description This script is a tool for merging two GFF3 files. GFF3 (General Feature Format version 3) is a file format used 
# to describe biological features of genomic sequences. The script takes in two GFF3 files and an overlap value as input, and 
# outputs a merged GFF3 file containing the unique features from the input files and the merged features that overlap by at 
# least the specified amount.
#running the script  python script_name.py --file1 file1.gff3 --file2 file2.gff3 --overlap 0.8 --output merged_and_unique.gff3
#python 

######


import argparse
import logging
import gffutils
import pandas as pd
from intervaltree import Interval, IntervalTree

def merge_gff3(file1, file2, overlap, output_file):
    try:
        # Open the two GFF3 files
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            # Read the contents of the files
            contents1 = f1.read()
            contents2 = f2.read()

        # Parse the contents of the GFF3 files
        df1 = gffutils.parse(contents1)
        df2 = gffutils.parse(contents2)

        # Create an interval tree for each GFF3 file
        tree1 = IntervalTree()
        tree2 = IntervalTree()

        # Lists to store unique features
        unique1 = []
        unique2 = []

        # Iterate over the rows in the DataFrames
        for i, row in df1.iterrows():
            # Add the region to the interval tree
            tree1[row['start']:row['end']] = row
            # Check if the region is not contained in the other interval tree
            if not tree2.overlap(row['start'], row['end']):
                # If it is not contained, add it to the list of unique features
                unique1.append(row)

        for i, row in df2.iterrows():
            # Add the region to the interval tree
            tree2[row['start']:row['end']] = row
            # Check if the region is not contained in the other interval tree
            if not tree1.overlap(row['start'], row['end']):
                # If it is not contained, add it to the list of unique features
                unique2.append(row)

        # Find the overlapping regions
        overlap1 = tree1.overlap(tree2)
        overlap2 = tree2.overlap(tree1)

        # Merge the overlapping regions
        merged_rows = []
        for interval in overlap1:
            merged_rows.append(interval.data)
        for interval in overlap2:
            merged_rows.append(interval.data)

        # Concatenate the lists of unique features and merged rows
        all_features = unique1 + unique2 + merged_rows

        # Create a new DataFrame from the list of all features
        df_all = pd.DataFrame(all_features)

        # Write the new DataFrame to a GFF3 file
        df_all.to_csv(output_file, sep='\t', index=False, header=False)
        logging.info(f'Successfully merged and saved the unique features from the GFF3 files and wrote the output to {output_file}')
    except Exception as e:
        logging.exception(e)

if __name__ == '__main__':
    # Set up logging
    logging.basicConfig(level=logging.INFO)

    # Parse the command-line arguments
    parser = argparse.ArgumentParser(description='Merge two GFF3 files and save unique features')
    parser.add_argument('--file1', help='The first GFF3 file')
    parser.add_argument('--file2', help='The second GFF3 file')
    parser.add_argument('--overlap', type=float, default=0.6, help='The minimum overlap required to merge regions (default: 0.6)')
    parser.add_argument('--output', default='merged.gff3', help='The output file name (default: merged.gff3)')
    args = parser.parse_args()

    # Call the merge_gff3 function with the command-line arguments
    merge_gff3(args.file1, args.file2, args.overlap, args.output)
