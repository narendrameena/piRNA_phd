#!/bin/bash

# Description: This script extracts intergenic regions from a provided GFF3 input file and outputs them in GFF3 format.
# Requirements: awk, bedtools, gff2bed

# Check for required tools
command -v awk >/dev/null 2>&1 || { echo "awk is required but it's not installed. Aborting." >&2; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "bedtools is required but it's not installed. Aborting." >&2; exit 1; }
command -v gff2bed >/dev/null 2>&1 || { echo "gff2bed is required but it's not installed. Aborting." >&2; exit 1; }

# Check for the right number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.gff3> <output.intergenic.gff3>"
    exit 1
fi

input_gff="$1"
output_gff3="$2"

# Extract the directory of the output file
output_dir=$(dirname "$output_gff3")

# Create a temporary directory within the directory of the output file
tmp_dir=$(mktemp -d -p "$output_dir")

# Extract the base name of the input GFF3 file
input_base=$(basename "$input_gff" .gff3)

# Generate genome size file from GFF3 by finding the max end coordinate for each chromosome or sequence
awk '{if($1 !~ /^#/ && $5+0>0) arr[$1]=($5>arr[$1]?$5:arr[$1])} END{for (i in arr) print i "\t" arr[i]}' "$input_gff" | sort -k1,1 > "$tmp_dir/genome_size.txt"

# Convert GFF3 to BED using BEDOPS with memory constraint and sorting in the temporary directory
gff2bed --max-mem 2G --sort-tmpdir "$tmp_dir" < "$input_gff" > "$tmp_dir/${input_base}.unsorted.bed"

# Sort the BED file based on the order in the genome size file
bedtools sort -i "$tmp_dir/${input_base}.unsorted.bed" -g "$tmp_dir/genome_size.txt" > "$tmp_dir/${input_base}.sorted.bed"

# Get intergenic regions
bedtools complement  -L  -i "$tmp_dir/${input_base}.sorted.bed" -g "$tmp_dir/genome_size.txt" > "$tmp_dir/intergenic.bed"

# Convert intergenic BED to GFF3
awk 'BEGIN{OFS="\t"} {print $1, "intergenic", "region", $2+1, $3, ".", ".", ".", "Parent=intergenic_" NR ";"}' "$tmp_dir/intergenic.bed" > "$output_gff3"

# Clean up the temporary directory
rm -r "$tmp_dir"

echo "Intergenic GFF3 has been saved to $output_gff3"
