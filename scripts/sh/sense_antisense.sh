#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <path_to_bam> <path_to_gff3> <reference_fasta> <output_directory>"
    exit 1
fi

bam_file="$1"
gff3_file="$2"
reference_fasta="$3"
output_dir="$4"

mkdir -p "$output_dir"
tmp_dir=$(mktemp -d -p "$output_dir" sense_antisense-XXXXXX)

# Create genome size file
samtools faidx "$reference_fasta"
cut -f1,2 "$reference_fasta.fai" > "$tmp_dir/genome_size.txt"

# Extract sense and antisense regions
awk '$3=="similarity" && $7=="+"' "$gff3_file" > "$tmp_dir/positive_strand_genes.gff3"
awk '$3=="similarity" && $7=="-"' "$gff3_file" > "$tmp_dir/negative_strand_genes.gff3"

# Extract reads
bedtools intersect -abam "$bam_file" -b "$tmp_dir/positive_strand_genes.gff3" -s -u > "$tmp_dir/sense_positive.bam"
bedtools intersect -abam "$bam_file" -b "$tmp_dir/negative_strand_genes.gff3" -S -u > "$tmp_dir/sense_negative.bam"
bedtools intersect -abam "$bam_file" -b "$tmp_dir/positive_strand_genes.gff3" -S -u > "$tmp_dir/antisense_positive.bam"
bedtools intersect -abam "$bam_file" -b "$tmp_dir/negative_strand_genes.gff3" -s -u > "$tmp_dir/antisense_negative.bam"

# Merge BAM files
samtools merge "$tmp_dir/all_sense.bam" "$tmp_dir/sense_positive.bam" "$tmp_dir/sense_negative.bam"
samtools merge "$tmp_dir/all_antisense.bam" "$tmp_dir/antisense_positive.bam" "$tmp_dir/antisense_negative.bam"

# Copy desired BAM files to output directory
cp "$tmp_dir/all_sense.bam" "$output_dir/all_sense.bam"
cp "$tmp_dir/all_antisense.bam" "$output_dir/all_antisense.bam"

# Convert BAM to BED
bedtools bamtobed -i "$tmp_dir/all_sense.bam" > "$output_dir/sense.bed"
bedtools bamtobed -i "$tmp_dir/all_antisense.bam" > "$output_dir/antisense.bed"

# Convert BAM to FASTA
samtools fasta "$tmp_dir/all_sense.bam" > "$output_dir/sense.fasta"
samtools fasta "$tmp_dir/all_antisense.bam" > "$output_dir/antisense.fasta"

# Generate BEDGraph for sense
bedtools genomecov -ibam "$tmp_dir/all_sense.bam" -bg -split > "$tmp_dir/sense.bedgraph"
LC_COLLATE=C sort -k1,1 -k2,2n "$tmp_dir/sense.bedgraph" > "$tmp_dir/sorted_sense.bedgraph"
bedGraphToBigWig "$tmp_dir/sorted_sense.bedgraph" "$tmp_dir/genome_size.txt" "$output_dir/sense.bw"

# Generate BEDGraph for antisense
bedtools genomecov -ibam "$tmp_dir/all_antisense.bam" -bg -split > "$tmp_dir/antisense.bedgraph"
LC_COLLATE=C sort -k1,1 -k2,2n "$tmp_dir/antisense.bedgraph" > "$tmp_dir/sorted_antisense.bedgraph"
bedGraphToBigWig "$tmp_dir/sorted_antisense.bedgraph" "$tmp_dir/genome_size.txt" "$output_dir/antisense.bw"

# Count the reads
sense_count=$(( $(samtools view -c "$tmp_dir/sense_positive.bam") + $(samtools view -c "$tmp_dir/sense_negative.bam") ))
antisense_count=$(( $(samtools view -c "$tmp_dir/antisense_positive.bam") + $(samtools view -c "$tmp_dir/antisense_negative.bam") ))

# Calculate ratio
if [ $antisense_count -eq 0 ]; then
    ratio="undefined"
else
    ratio=$(echo "scale=2; $sense_count/$antisense_count" | bc)
fi



echo "sense_negative.bam and antisense_negative.bam have been saved to $output_dir"

# Write results
echo -e "Type\tCount" > "$output_dir/results.tsv"
echo -e "Sense\t$sense_count" >> "$output_dir/results.tsv"
echo -e "Antisense\t$antisense_count" >> "$output_dir/results.tsv"
echo -e "Sense/Antisense Ratio\t$ratio" >> "$output_dir/results.tsv"

echo "Results and BED files written to $output_dir"

# Cleanup
#rm "$reference_fasta.fai"
rm -rf "$tmp_dir"
