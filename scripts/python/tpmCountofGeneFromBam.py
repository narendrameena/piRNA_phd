
####### 
#@author narumeena
#@description This script calculates Transcripts Per Million (TPM) counts for transcripts 
# in a GFF3 file using mapped read counts from a BAM file.
#python tpmCountofGeneFromBam.py --bam_file bamFile.bam --gff3_file annotation.gff3 --output_file outputFile.tsv

######


import argparse
import csv
import logging
import pysam

def calculate_tpm_counts(bam_file, gff3_file):
    """Calculates TPM (Transcripts Per Million) counts for transcripts in a GFF3 file using mapped read counts from a BAM file"""
    tpm_counts = {}
    total_mapping_reads = 0

    # Open and parse the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Count the total number of mapping reads in the BAM file
    for read in bam:
        if not read.is_unmapped:
            total_mapping_reads += 1

    # Parse the GFF3 file
    with open(gff3_file, "r") as gff3:
        for line in gff3:
            # Skip comment lines
            if line.startswith("#"):
                continue
            # Parse the fields of each line
            fields = line.strip().split("\t")
            # Check if the line is a transcript
            if fields[2] == "transcript":
                # Extract the transcript identifier
                transcript_id = fields[8].split(";")[0].split("=")[1]
                # Count the number of reads that map to the transcript
                read_count = 0
                for read in bam.fetch(fields[0], int(fields[3]), int(fields[4])):
                    read_count += 1
                # Normalize the read count by the total number of mapping reads and the length of the transcript
                tpm = (read_count / total_mapping_reads) * (10**6) / (int(fields[4]) - int(fields[3]))
                tpm_counts[transcript_id] = tpm

    # Close the BAM file
    bam.close()

    return tpm_counts

def main(bam_file, gff3_file, output_file):
    """Calculates TPM counts and writes the results to a TSV file"""
    try:
        # Calculate TPM counts
        tpm_counts = calculate_tpm_counts(bam_file, gff3_file)

        # Write the TPM counts to a TSV file
        with open(output_file, "w", newline="") as tsv_file:
            writer = csv.writer(tsv_file, delimiter="\t")
            writer.writerow(["Transcript ID", "TPM Count"])
            for transcript_id, tpm in tpm_counts.items():
                writer.writerow([transcript_id, tpm])

    except IOError:
        # Log the error message and exit the program
        logging.exception("An error occurred while reading/writing a file")
        sys.exit(1)

    except ValueError:
        # Log the error message and exit the program
        logging.exception("An error occurred due to invalid input values")
        sys.exit(1)

    except Exception as e:
        logging.exception("An error occurred while calculating TPM counts")
        sys.exit(1)
        
    

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam_file", help="Path to the BAM file")
    parser.add_argument("--gff3_file", help="Path to the GFF3 file")
    parser.add_argument("--output_file", help="Path to the output TSV file")
    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    # Run the main function
    main(args.bam_file, args.gff3_file, args.output_file)



