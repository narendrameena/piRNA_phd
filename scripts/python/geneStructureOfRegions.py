
######

#@author narumeena
#@description This script processes RNA paired-end sequencing data and a reference genome to 
# generate a GFF3 file containing the gene expression levels in specified regions of interest.
#@tools required conda install samtools bedtools hisat2 stringtie gff3
#running script for paired end 
#python script.py --sequencing_data sequencing_data.fq --sequencing_data_2 read2.fq --sequencing_type paired --reference_genome reference_genome.fa --regions regions.bed --output output.gff3
#running the script for single end 
#python script.py --sequencing_data sequencing_data.fq --sequencing_type single --reference_genome reference_genome.fa --regions regions.bed --output output.gff3


#####


import argparse
import tempfile
import subprocess
import os
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser()
parser.add_argument("--sequencing_data", required=True, help="Path to the RNA sequencing data file")
parser.add_argument("--sequencing_data_2", help="Path to the second read file of the RNA paired-end sequencing data")
parser.add_argument("--sequencing_type", required=True, choices=["single", "paired"], help="Type of RNA sequencing data (single or paired)")
parser.add_argument("--reference_genome", required=True, help="Path to the reference genome")
parser.add_argument("--regions", required=True, help="Path to the BED file containing the regions of interest")
parser.add_argument("--output", required=True, help="Path to the output file")

args = parser.parse_args()
sequencing_data = args.sequencing_data
sequencing_data_2 = args.sequencing_data_2
sequencing_type = args.sequencing_type
reference_genome = args.reference_genome
regions = args.regions
output = args.output

with tempfile.TemporaryDirectory() as tmp_dir:
    # Check if the required tools are installed
    try:
        subprocess.run(["hisat2", "--version"], check=True)
        subprocess.run(["bedtools", "--version"], check=True)
        subprocess.run(["stringtie", "--version"], check=True)
        subprocess.run(["gff3", "--version"], check=True)
    except subprocess.CalledProcessError:
        print("One or more required tools are not installed or not in the PATH")
        exit(1)

    # Check if the input files exist
    if not os.path.exists(sequencing_data):
        print(f"Error: The file '{sequencing_data}' does not exist")
        exit(1)
    if sequencing_type == "paired" and not os.path.exists(sequencing_data_2):
        print(f"Error: The file '{sequencing_data_2}' does not exist")
        exit(1)
    if not os.path.exists(reference_genome):
        print(f"Error: The file '{reference_genome}' does not exist")
        exit(1)
    if not os.path.exists(regions):
        print(f"Error: The file '{regions}' does not exist")
        exit(1)

    # Generate the HISAT2 index in the temporary directory
    logger.info("Generating the HISAT2 index")
    subprocess.run(["hisat2-build", "-p", "8", reference_genome, f"{tmp_dir}/reference_genome_index"])

    # Align the reads to the reference genome using HISAT2
    logger.info("Aligning the reads to the reference genome using HISAT2")
    if sequencing_type == "single":
        subprocess.run(["hisat2", "-p", "8", "--dta", "-x", f"{tmp_dir}/reference_genome_index", "--single-end", sequencing_data, "-S", f"{tmp_dir}/aligned_reads.sam"])
    else:
        subprocess.run(["hisat2", "-p", "8", "--dta", "-x", f"{tmp_dir}/reference_genome_index", "-1", sequencing_data, "-2", sequencing_data_2, "-S", f"{tmp_dir}/aligned_reads.sam"])

    # Extract the relevant data from the reference genome
    logger.info("Extracting the relevant data from the reference genome")
    subprocess.run(["bedtools", "getfasta", "-fi", reference_genome, "-bed", regions, "-fo", f"{tmp_dir}/extracted_genome.fa"])

    # Generate a list of transcripts from the aligned reads and the extracted genome data using stringtie
    logger.info("Generating a list of transcripts using stringtie")
    subprocess.run(["stringtie", f"{tmp_dir}/aligned_reads.sam", "-G", f"{tmp_dir}/extracted_genome.fa", "-o", f"{tmp_dir}/transcripts.gtf"])

    # Convert the transcript data into a GFF3 file using gff3
    logger.info("Converting the transcript data into a GFF3 file")
    subprocess.run(["gff3", "-i", f"{tmp_dir}/transcripts.gtf", "-o", output])

