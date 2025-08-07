#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <gene_annotation_file> <piRNA_input_file> <output_file>"
    exit 1
fi

# Input file paths
gene_annotation_file="$1"
piRNA_input_file="$2"
output_file="$3"

# Create a temporary directory for intermediate files
tmp_dir=$(mktemp -d)

# Intermediate file paths
genes_with_names="$tmp_dir/genes_with_names.bed"
converted="$tmp_dir/converted.bed"
genome_file="$tmp_dir/genome.txt"
genes_with_names_sorted="$tmp_dir/genes_with_names_sorted.bed"
sorted_converted="$tmp_dir/sorted_converted.bed"

# Command to process files
awk -F'\t' 'BEGIN{OFS="\t"} $3 == "gene" {split($9, a, ";"); for (i in a) {if (match(a[i], "Name=")) {split(a[i], name, "=")}}; print $1, $4-1, $5, name[2], ".", $7}' "$gene_annotation_file" > "$genes_with_names"

awk -F'\t' 'BEGIN {OFS="\t"} $3=="piRNA" {print $1, $4-1, $5, $9, ".", $7}' "$piRNA_input_file"| grep -v MT > "$converted"

awk '$1 == "##sequence-region" { print $2 "\t" $4 }' "$gene_annotation_file" | sort -k1,1 -V > "$genome_file"

bedtools sort -i "$genes_with_names" -g "$genome_file" > "$genes_with_names_sorted"

bedtools sort -i "$converted" -g "$genome_file" | grep -v MT > "$sorted_converted"

bedtools closest -a "$sorted_converted" -b "$genes_with_names_sorted" -g "$genome_file" -wo -s -sorted -k 1 -D a > "$output_file"

# Clean up temporary directory
rm -rf "$tmp_dir"

echo "Processing complete. Output saved in $output_file"
