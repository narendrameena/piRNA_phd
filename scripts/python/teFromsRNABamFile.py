

####### 
#@author narumeena
#@description python script pipeline for analyzing paired-end sRNA (small ribonucleic acid) (STAR alinged)
# data and repeatmasker files to capture transposon element expression.
#running the script
#python script.py --repeatmasker /path/to/repeatmasker.txt --srna /path/to/srna.bam --output /path/to/output.tsv

######

import argparse
import logging
import pybedtools
import pandas as pd

# Set up argument parser
parser = argparse.ArgumentParser(description="Analyze paired-end sRNA data and repeatmasker files to capture transposon element expression")
parser.add_argument("--repeatmasker", required=True, help="Repeatmasker file in txt format")
parser.add_argument("--srna", required=True, help="sRNA alignment file in bam format")
parser.add_argument("--output", required=True, help="Output file in tsv format")

# Parse command line arguments
args = parser.parse_args()

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

try:
    # Load the repeatmasker file as a pandas dataframe
    df_repeatmasker = pd.read_csv(args.repeatmasker, sep="\t", skiprows=3, header=None, names=["chrom", "start", "end", "strand", "repeat_name", "repeat_class", "repeat_family", "score", "align_pos"])

    # Filter the dataframe to keep only the columns of interest
    df_repeatmasker = df_repeatmasker[["chrom", "start", "end", "repeat_name", "strand"]]

    # Convert the dataframe to a BedTool object
    repeatmasker = pybedtools.BedTool.from_dataframe(df_repeatmasker)

    # Load the sRNA alignment file as a BedTool object
    srna_alignment = pybedtools.BedTool(args.srna)

    # Intersect the sRNA alignment and repeatmasker file to count the number of reads
    # that map to each transposon element
    counts = srna_alignment.intersect(repeatmasker, c=True)

    # Convert the counts to a pandas dataframe
    df = counts.to_dataframe(names=["chrom", "start", "end", "transposon", "strand", "count"])

    # Calculate the total number of reads mapping to the genome
    total_reads = df["count"].sum()

    # Calculate the expression level of each transposon element in RPM
    df["rpm"] = df["count"] / total_reads * 1e6

    # Sort the dataframe by RPM in descending order
    df = df.sort_values("rpm", ascending=False)

except Exception as e:
    logging.error("An error occurred while processing the input files: %s", e)

try:
    # Save the results to a TSV file
    df.to_csv(args.output, sep="\t", index=False)

except Exception as e:
    logging.error("An error occurred while writing the output file: %s", e)
