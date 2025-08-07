
####### 
#@author narumeena
#@description This script is a Python script that takes in two BED12 files as input and merges them, 
# keeping merged and  unique features.
# running the script
#python script_name.py --file1 file1.bed12 --file2 file2.bed12 --overlap 0.8 --output merged_and_unique.bed12
######



def merge_bed12(file1, file2, overlap, output_file):
    # Open the two BED12 files
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        # Read the contents of the files
        contents1 = f1.read()
        contents2 = f2.read()

    # Parse the contents of the BED12 files into DataFrames
    df1 = pd.read_csv(file1, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'])
    df2 = pd.read_csv(file2, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts'])

    # Create an interval tree for each BED12 file
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

    # Write the new DataFrame to a BED12 file
    df_all.to_csv(output_file, sep='\t', index=False, header=False)
    logging.info(f'Successfully merged and saved the unique features from the BED12 files and wrote the output to {output_file}')

if __name__ == '__main__':
    # Set up logging
    logging.basicConfig(level=logging.INFO)

    # Parse the command-line arguments
    parser = argparse.ArgumentParser(description='Merge two BED12 files and save unique features')
    parser.add_argument('--file1', help='The first BED12 file')
    parser.add_argument('--file2', help='The second BED12 file')
    parser.add_argument('--overlap', type=float, default=0.6, help='The minimum overlap required to merge regions (default: 0.6)')
    parser.add_argument('--output', default='merged.bed12', help='The output file name (default: merged.bed12)')
    args = parser.parse_args()

    # Call the merge_bed12 function with the command-line arguments
    merge_bed12(args.file1, args.file2, args.overlap, args.output)


