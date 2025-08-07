#!/bin/bash

# Check if required tools are installed
command -v bedtools >/dev/null 2>&1 || { echo >&2 "bedtools is required but it's not installed. Aborting."; exit 1; }

# Check if correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <feature_count_file> <gff3_file> <output_gene_names_file>"
    exit 1
fi

feature_count_file=$1
gff3_file=$2
output_gene_names_file=$3

# Create a BED file from feature count file with counts > 0
awk '$7 > 0 && !/;/' $feature_count_file | grep 'gene_gene'|grep 'Name'| awk -F'\\TBL' '{if (NF == 7){ match($1, /Name=[^:]+/);  namePart = substr($1, RSTART+5, RLENGTH-5); print  $2, $3, $4, namePart, $7, $5 }}'  |awk '$5 > 0 && !/;/'| tr " " "\t" > features_with_counts.bed

# Filter 3' UTR entries from GFF3 file and intersect with the BED file
# -s option is used to enforce strandness
# -f 1.0 ensures 100% of the gene overlaps with the gene region
awk '$3 == "gene"' $gff3_file | bedtools intersect -a stdin -b features_with_counts.bed -wa -wb -s -f 1.0 | awk '{print $9}' | sed 's/.*Name=\([^;]*\).*/\1/'| grep -v 'ID' | sort | uniq > $output_gene_names_file

# Cleanup
#rm features_with_counts.bed

echo "Gene names extracted to $output_gene_names_file"
