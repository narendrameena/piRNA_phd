
####### 
#@author narumeena
#@description This script takes in a list of filenames as command line arguments and combines 
# their contents into a single output file. It does this by reading each file, storing the 
# contents in a dictionary, and then writing the contents of the dictionary to the output file.
#running the script
#python python script.py --filenames input1.tsv input2.tsv -o output.tsv


######


import argparse
import csv

# Create the parser
parser = argparse.ArgumentParser()

# Add an argument for the list of filenames
parser.add_argument('--filenames', nargs='+', help='List of filenames to be combined')

# Add an argument for the output filename
parser.add_argument('-o', '--output', help='Name of the output file')

# Parse the arguments
args = parser.parse_args()

# Dictionary to store the data
data = {}

# Read the data from each file and store it in the dictionary
for filename in args.filenames:
    with open(filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        # Skip the header row
        next(reader)
        for row in reader:
            key = row[0]
            if key in data:
                data[key].append(row[1:])
            else:
                data[key] = [row[1:]]

# Write the combined data to a new file
with open(args.output, 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    for key, rows in data.items():
        for row in rows:
            writer.writerow([key] + row)
