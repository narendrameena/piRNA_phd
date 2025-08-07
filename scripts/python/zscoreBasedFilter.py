
####### 
#@author narumeena
#@description python  script for markers slection  for each condition with respect to other conditions based 
# on z score from a tab limited file where columns are conditions and rows attributes and each cell contain 
# TPM values and save it to a file along with z-score 
#running the script
#python zscoreBasedFilter.py --input_file input_file --output_file output_file --z_score_threshold  z_score_threshold

######

import sys
import logging
import argparse
import pandas as pd

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Set up command line argument parsing
parser = argparse.ArgumentParser(description="Filter markers based on z-scores from a tab-separated file.")
parser.add_argument("--input_file", help="Input file (TSV format)")
parser.add_argument("--output_file", help="Output file (TSV format)")
parser.add_argument("--z_score_threshold", type=float, default=2.0, help="Z-score threshold for selecting markers (default: 2.0)")
args = parser.parse_args()

# Set input and output files
input_file = args.input_file
output_file = args.output_file
z_score_threshold = args.z_score_threshold

try:
    # read in the TPM values from the file
    df = pd.read_csv(input_file, sep='\t')

except FileNotFoundError:
    logging.exception(f"Error: Input file '{input_file}' not found.")
    sys.exit(1)

# compute the z-scores for each attribute
logging.info("Computing z-scores...")
z_scores = df.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

# select the markers for each condition based on the z-scores
logging.info("Selecting markers...")
markers = {}
for col in df.columns:
    markers[col] = z_scores[z_scores[col] > z_score_threshold].index.tolist()

# save the markers and z-scores to a file
logging.info("Saving results to file...")
try:
    with open(output_file , 'w') as f:
        for col, m in markers.items():
            f.write(f"{col}\t{','.join(m)}\n")
            z_scores[m].to_csv(f, sep='\t', header=False)
except OSError:
    logging.exception(f"Error: Could not write to output file '{output_file}'.")
    sys.exit(1)

logging.info("Process completed successfully.")
