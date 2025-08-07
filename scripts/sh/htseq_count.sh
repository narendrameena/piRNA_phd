#!/bin/bash

# Check if correct number of arguments are passed
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <aligned_bam> <annotation_gtf> <output_folder>"
    exit 1
fi

# Arguments
aligned_bam=$1
annotation_gtf=$2
output_folder=$3

# Create a temporary folder for intermediate files
tmp_folder="$output_folder/tmp"
mkdir -p "$tmp_folder"

# Step 1: Use htseq-count to get counts for each transcript
htseq-count -f bam -r pos -s no -t transcript -i transcript_id "$aligned_bam" "$annotation_gtf" > "$tmp_folder/counts.txt"

# Step 2: Calculate total mapped reads in the BAM file
total_reads=$(samtools view -c -F 260 "$aligned_bam")

# Step 3: Calculate RPM values for each transcript
awk -v total="$total_reads" 'NR > 1 {print $1, $2, $2 / total * 1000000}' "$tmp_folder/counts.txt" > "$output_folder/rpm_values.txt"

# Step 4: Calculate RPKM values for each transcript
awk -v total="$total_reads" 'NR > 1 {print $1, $2, $2 / ($3 * (length($2)/1000)) * 1000000}' "$tmp_folder/counts.txt" > "$output_folder/rpkm_values.txt"

# Clean up temporary files
rm -rf "$tmp_folder"
