#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <fasta_file> <star_mapped_bam> <reference_genome>"
    exit 1
fi

# Assign arguments to variables
FASTA_FILE=$1
STAR_MAPPED_BAM=$2
REFERENCE_GENOME=$3

# Create a directory for STAR genome data
mkdir star_genome_dir

# Generate the STAR genome directory
STAR --runMode genomeGenerate --runThreadN 4 --genomeDir star_genome_dir --genomeFastaFiles $REFERENCE_GENOME

# Align your fasta sequences
STAR --runMode alignReads --runThreadN 4 --genomeDir star_genome_dir --readFilesIn $FASTA_FILE --outFileNamePrefix aligned_sequences_

# Convert SAM to BAM, sort, and index
samtools view -bS aligned_sequences_Aligned.out.sam > aligned_sequences.bam
samtools sort aligned_sequences.bam -o aligned_sequences.sorted.bam
samtools index aligned_sequences.sorted.bam

# Extract coordinates and retrieve reads from the STAR-mapped BAM file
samtools idxstats aligned_sequences.sorted.bam | awk 'BEGIN {FS="\t"}; $3>0 {print $1}' > coordinates.txt
while read coordinate; do
   samtools view -b $STAR_MAPPED_BAM $coordinate >> extracted_sequences.bam
done < coordinates.txt

# Count reads (You'd need a GTF annotation file for this. Replace 'your_annotation_file.gtf' with your file)
# featureCounts -a your_annotation_file.gtf -o counts.txt extracted_sequences.bam

# Cleanup
rm -r star_genome_dir
rm aligned_sequences_Aligned.out.sam
rm aligned_sequences.bam
rm aligned_sequences.sorted.bam.bai
rm coordinates.txt

echo "Extraction completed. The extracted BAM file is named extracted_sequences.bam"
