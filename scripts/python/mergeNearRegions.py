


####### 
#@author narumeena
#@description This is a script that takes in a BED12 file (a file in the BED format that includes 12 tab-separated 
# fields) and merges regions that are within 1kb of each other.

#running the script
#python merge_regions.py --input_bed12_file <input_bed12_file> --output_bed12_file <output_bed12_file> --distance <distance>

######



import argparse
import logging
import pybedtools
import sys

def merge_regions(bed12_file, output_file, distance):
    # Create a dictionary to store the merged regions
    merged_regions = {}

    # Use pybedtools to parse the BED12 file
    try:
        bed12 = pybedtools.BedTool(bed12_file)
    except pybedtools.helpers.BEDToolsError as e:
        logger.error("Error parsing BED12 file: " + str(e))
        sys.exit(1)

    # Iterate through each region in the BED12 file
    for region in bed12:
        chrom = region.chrom
        start = region.start
        end = region.end
        strand = region.strand

        # Check if the chromosome is already in the dictionary
        if chrom not in merged_regions:
            # If it's not, create a new entry for the chromosome in the dictionary
            merged_regions[chrom] = []
            # Add the current region to the list of merged regions for the chromosome
            merged_regions[chrom].append((start, end, strand))
        else:
            # Get the last region in the list of merged regions for the chromosome
            last_region = merged_regions[chrom][-1]
            mr_start, mr_end, mr_strand = last_region
            # Check if the current region is within the specified distance of the last region
            # and if the strands are the same
            if (start - mr_end <= distance or mr_start - end <= distance) and strand == mr_strand:
                # If it is, update the last region with the new start and end coordinates
                mr_start = min(start, mr_start)
                mr_end = max(end, mr_end)
                merged_regions[chrom][-1] = (mr_start, mr_end, strand)
            else:
                # If the region is not within the specified distance of the last region, add it to the list
                merged_regions[chrom].append((start, end, strand))

    # Open the output file for writing
    try:
        with open(output_file, 'w') as f:
            # Iterate through the merged regions dictionary
            for chrom, regions in merged_regions.items():
                # Iterate through the list of merged regions for each chromosome
                for start, end, strand in regions:
                    # Write the merged region to the output file in BED12 format
                    f.write(f"{chrom}\t{start}\t{end}\t.\t.\t{strand}\t.\t.\t.\t.\t.\t.\n")
    except IOError as e:
        logger.error("Error writing to output file: " + str(e))
        sys.exit(1)



def main():
    # Set up logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # Parse the command line arguments
    parser = argparse.ArgumentParser(description='Merge regions in a BED12 file that are within a specified distance of each other')
    parser.add_argument('--input_bed12_file', help='input BED12 file')
    parser.add_argument('--output_bed12_file', help='output BED12 file')
    parser.add_argument('--distance', type=int, default=1000, help='maximum overlap distance default is 1000')
    args = parser.parse_args()

    # Check if the correct number of command line arguments were provided
    if len(sys.argv) < 3:
        logger.error("Usage: python merge_regions.py --input_bed12_file <input_bed12_file> --output_bed12_file <output_bed12_file> --distance <distance>")
        sys.exit(1)

    # Call the merge_regions function with the input and output file names
    try:
        merge_regions(args.input_bed12_file, args.output_bed12_file, args.distance)
        logger.info("Successfully merged regions.")
    except Exception as e:
        # Log the error message and traceback
        logger.exception("Error while merging regions:")
        # Exit with a non-zero exit code to indicate an error occurred
        sys.exit(1)

# Call the main function
if __name__ == '__main__':
    main()





    
