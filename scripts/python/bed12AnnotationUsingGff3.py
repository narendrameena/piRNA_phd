
####### 
#@author narumeena
#@description The script then uses the pybedtools module to open the BED12 and GFF3 files and create BedTool 
# objects from them. It then uses the intersect method of the BedTool object to annotate the BED12 regions 
# with the GFF3 annotations. 
#running the script
#python script.py --bed12_file path/to/bed12_file.bed --gff3_file path/to/gff3_file.gff3 --output_file path/to/output_file.csv
######


import argparse
import csv
import logging
import pybedtools

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(message)s")

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--bed12_file", help="path to BED12 file")
parser.add_argument("--gff3_file", help="path to GFF3 file")
parser.add_argument("--output_file", help="path to output CSV file")
args = parser.parse_args()

try:
    # Open the BED12 file and create a BedTool object
    with open(args.bed12_file, "r") as bed12_file:
        bed12 = pybedtools.BedTool(bed12_file)

    # Open the GFF3 file and create a BedTool object
    with open(args.gff3_file, "r") as gff3_file:
        gff3 = pybedtools.BedTool(gff3_file)

    # Annotate the BED12 regions with the GFF3 annotations
    annotated = bed12.intersect(gff3, wa=True, wb=True)

    # Open the output CSV file
    with open(args.output_file, "w", newline="") as csv_file:
        # Create a CSV writer object
        writer = csv.writer(csv_file)

        # Write the header row
        writer.writerow(["chrom", "start", "end", "name", "score", "strand", "thick_start", "thick_end", "item_rgb", "block_count", "block_sizes", "block_starts", "annotation_chrom", "annotation_source", "annotation_type", "annotation_start", "annotation_end", "annotation_score", "annotation_strand", "annotation_phase", "annotation_attributes"])

        # Iterate over the annotated regions
        for region in annotated:
            # Write the region data to the CSV file
            writer.writerow(region)

except Exception as e:
    # Log an error message and exit
    logging.error(f"An error occurred: {e}")
    exit(1)

# Log a message indicating that the script has finished
logging.info("Finished annotating regions and saving to CSV file")
