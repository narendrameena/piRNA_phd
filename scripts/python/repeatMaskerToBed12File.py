

####### 
#@author narumeena
#@description This script converts a file in the RepeatMasker format to the BED12 format.
#running the script
#python repeatMaskerToBed12File.py --repeatmasker_file repeatmasker_file --bed12_file output.bed12

######

import sys
import pandas as pd 
import argparse
import logging
from tqdm import tqdm

def convert_to_bed12(repeatmasker_file, bed12_file):
    try:
        #Read in repeatmasker file
        repeatmasker_df = pd.read_csv(repeatmasker_file, sep='\s+', skiprows=3, header=None)

        #Define columns
        repeatmasker_df.columns = ['SW_score', 'perc_div', 'perc_del', 'perc_ins', 'query_name', 'query_start', 'query_end', 'query_left', 'strand', 'repeat_name', 'repeat_class', 'repeat_start', 'repeat_end', 'repeat_left', 'id']

        #Convert to bed12 format
        bed12_df = repeatmasker_df[['query_name', 'query_start', 'query_end', 'repeat_name', 'SW_score', 'strand', 'query_start', 'query_end', '0', '1', 'query_end-query_start', 'strand+","+repeat_name']]

        #Rename columns to match bed12 format
        bed12_df.columns = ['chrom', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']

        #Save bed12 file
        bed12_df.to_csv(bed12_file, sep='\t', index=False, header=False)
    except Exception as e:
        logging.exception("An error occurred during conversion")
        raise e

if __name__ == "__main__":
    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Convert a RepeatMasker file to BED12 format')
    parser.add_argument('--repeatmasker_file', help='Input file in RepeatMasker format')
    parser.add_argument('--bed12_file', help='Output file in BED12 format')
    args = parser.parse_args()

    # Convert file
    logging.info("Starting conversion")
    convert_to_bed12(args.repeatmasker_file, args.bed12_file)
    logging.info("Conversion complete")


