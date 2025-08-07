#!/bin/bash

# Check if correct number of arguments are passed
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_aligned_bam> <input_annotation_gtf> <input_count> <output_file>"
    exit 1
fi

aligned_bam=$1
annotation_gtf=$2
input_count=$3
output_file=$4

# Calculate total reads
total_reads=$(samtools view -c -F 260 "$aligned_bam")

# Calculate RPM and RPKM values and append as extended columns
awk -v total="$total_reads" 'BEGIN{OFS="\t"} NR > 2 {rpm = ($7 * 1000000) / total; rpkm = ($7 / (($6 / 1000) * (total / 1000000))); print $0, rpm, rpkm}' "$input_count" > "$output_file"
