

####### 
#@author narumeena
#@description This Python script reads a CSV file specified as the first command line argument, 
# transposes each row of the input file into a column, and writes the resulting columns to a 
# new CSV file specified as the second command line argument.
#running the script
#python rawsToColumnsACsvFile.py --input_file rawFile.csv --output_file columnFile.csv


#Import necessary libraries
import csv
import sys
import os
import argparse
import logging
from tqdm import tqdm

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

# Set up argparse
parser = argparse.ArgumentParser(description='Transpose rows to columns in a CSV file')
parser.add_argument('--input_file', type=str, help='input CSV file')
parser.add_argument('--output_file', type=str, help='output CSV file')
args = parser.parse_args()

# Set input and output files
input_file = args.input_file
output_file = args.output_file

# Open the input csv file
try:
    with open(input_file, 'r') as input_file:
        # Read the file as a csv
        reader = csv.reader(input_file)

        # Open the output csv file
        with open(output_file, 'w') as output:
            # Create a csv writer
            writer = csv.writer(output)

            # Iterate over the rows in the input file
            for row in tqdm(reader, desc='Transposing rows to columns'):
                # Transpose the row to create a column
                column = [row]

                # Write the column to the output file
                writer.writerows(column)

except FileNotFoundError:
    logging.error(f'{input_file} not found')
    sys.exit(1)
except PermissionError:
    logging.error(f'Permission denied for {input_file}')
    sys.exit(1)
