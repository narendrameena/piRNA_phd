import os
import sys
import subprocess
import pandas as pd
import pybedtools

def get_chromosomes_from_bed(bed_file):
    """Extract chromosome names from a BED file, excluding those containing 'MT'."""
    chromosomes = set()
    with open(bed_file, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith('#'):
                chrom = line.split('\t')[0]
                if 'MT' not in chrom:
                    chromosomes.add(chrom)
    return chromosomes

def get_chromosomes_from_gff3(gff3_file):
    """Extract chromosome names from a GFF3 file, excluding those containing 'MT'."""
    chromosomes = set()
    with open(gff3_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            if 'MT' not in chrom:
                chromosomes.add(chrom)
    return chromosomes


def get_common_chromosomes(gff3_file, bed_file):
    """Get chromosomes that are common between GFF3 and BED files."""
    gff3_chromosomes = get_chromosomes_from_gff3(gff3_file)
    bed_chromosomes = get_chromosomes_from_bed(bed_file)
    return gff3_chromosomes.intersection(bed_chromosomes)

def create_genome_file_from_gff3(gff3_file, common_chromosomes, genome_file):
    """Create a genome file from a GFF3 file, including only common chromosomes."""
    try:
        chrom_lengths = {}
        with open(gff3_file, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                chrom, end = parts[0], int(parts[4])
                if chrom in common_chromosomes:
                    if chrom not in chrom_lengths or chrom_lengths[chrom] < end:
                        chrom_lengths[chrom] = end

        with open(genome_file, 'w') as out:
            for chrom, length in sorted(chrom_lengths.items()):
                out.write(f"{chrom}\t{length}\n")

        print(f"Genome file created for common chromosomes: {genome_file}")
    except Exception as e:
        print(f"Error occurred while creating genome file: {e}", file=sys.stderr)
        sys.exit(1)



def filter_bed12_file(bed_file, common_chromosomes, output_dir):
    """Filter the BED12 file to include only common chromosomes."""
    filtered_bed_file = os.path.join(output_dir, os.path.basename(bed_file) + '.filtered.bed12')
    try:
        with open(bed_file, 'r') as infile, open(filtered_bed_file, 'w') as outfile:
            for line in infile:
                if line.strip() and not line.startswith('#'):
                    chrom = line.split('\t')[0]
                    if chrom in common_chromosomes:
                        outfile.write(line)
        return filtered_bed_file
    except Exception as e:
        print(f"Error occurred while filtering BED12 file: {e}", file=sys.stderr)
        sys.exit(1)



def gff3_to_bed(gff3_file, common_chromosomes, output_dir):
    """Convert GFF3 to BED format, keeping only genes with 'Name' attribute and common chromosomes."""
    bed_output = os.path.join(output_dir, 'genes.bed')
    try:
        with open(gff3_file, 'r') as gff3, open(bed_output, 'w') as bed_out:
            for line in gff3:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if parts[0] in common_chromosomes and parts[2] == 'gene':
                    gene_name = get_name_attribute(parts[8])
                    if gene_name:
                        start = int(parts[3]) - 1  # Convert to 0-based start for BED
                        bed_out.write(f"{parts[0]}\t{start}\t{parts[4]}\t{gene_name}\n")
        print(f"Conversion completed: GFF3 to BED with filtering for common chromosomes and 'Name' attribute.")
        return pybedtools.BedTool(bed_output)
    except Exception as e:
        print(f"Error occurred during conversion: {e}", file=sys.stderr)
        sys.exit(1)



def sort_bed(bed_file, genome_file, output_dir):
    """Sort a BED file using a genome file."""
    sorted_bed_file = os.path.join(output_dir, os.path.basename(bed_file) + '.sorted.bed')
    try:
        bed = pybedtools.BedTool(bed_file)
        bed.sort(g=genome_file).saveas(sorted_bed_file)
        return sorted_bed_file
    except Exception as e:
        print(f"Error occurred while sorting BED file: {e}", file=sys.stderr)
        sys.exit(1)

def rename_regions(bed_file, genes_bed, output_dir):
    """Rename regions based on overlapping or nearest genes, keeping locus info and gene 'Name', handling malformed BED lines."""
    try:
        bed = pybedtools.BedTool(bed_file)
        genes_bed = pybedtools.BedTool(genes_bed)

        # Process each region to determine the relationship with genes
        df_renamed = pd.DataFrame()
        for region in bed:
            try:
                chrom, start, end, name, score, strand = region.fields[0:6]

                # Find overlapping genes
                overlapping_genes = genes_bed.intersect(pybedtools.BedTool(region), wa=True, wb=True)
                if len(overlapping_genes) > 0:
                    # Overlapping gene found
                    gene_name = get_name_attribute(overlapping_genes[0].fields[8])
                    relationship = "overlap"
                else:
                    # Find the closest gene for non-overlapping regions
                    closest_gene = genes_bed.closest(pybedtools.BedTool(region), d=True, t='first')[0]
                    gene_name = get_name_attribute(closest_gene.fields[8])
                    distance = int(closest_gene.fields[-1])
                    relationship = "upstream" if distance >= 0 else "downstream"

                # Construct new name with relationship info
                new_name = f"{name} ({relationship} with {gene_name})"

                # Append to the dataframe
                df_renamed = df_renamed.append({'chrom': chrom, 'start': start, 'end': end, 'name': new_name, 'score': score, 'strand': strand}, ignore_index=True)

            except pybedtools.cbedtools.MalformedBedLineError as e:
                print(f"Malformed BED line encountered: {e}")
                continue  # Skip this line and continue processing

        # Sort and save the renamed BED file
        df_renamed.sort_values(by=['chrom', 'start', 'end'], inplace=True)
        renamed_bed = pybedtools.BedTool.from_dataframe(df_renamed)
        renamed_bed.saveas(os.path.join(output_dir, 'renamed_regions.bed'))
        print(f"Renamed BED file with {len(df_renamed)} entries created as {os.path.join(output_dir, 'renamed_regions.bed')}")

    except Exception as e:
        print(f"Unexpected error occurred during region renaming: {e}", file=sys.stderr)

def get_name_attribute(gff3_attributes):
    """Extract the 'Name' attribute from GFF3 attributes."""
    attributes = dict(item.split('=') for item in gff3_attributes.split(';') if '=' in item)
    return attributes.get('Name', 'NA')


def main(gff3_file, bed_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print("Identifying common chromosomes...")
    common_chromosomes = get_common_chromosomes(gff3_file, bed_file)
    
    genome_file = os.path.join(output_dir, "genome.txt")
    
    # Create genome file for common chromosomes with their lengths
    create_genome_file_from_gff3(gff3_file, common_chromosomes, genome_file)
    
    print("Converting GFF3 to BED...")
    genes_bed = gff3_to_bed(gff3_file, common_chromosomes, output_dir)
    
    print("Filtering input BED12 file...")
    filtered_bed12_file = filter_bed12_file(bed_file, common_chromosomes, output_dir)

    print("Sorting filtered BED12 file...")
    sorted_bed12_file = sort_bed(filtered_bed12_file, genome_file, output_dir)
    
    print("Renaming BED regions...")
    rename_regions(sorted_bed12_file, genes_bed, output_dir)
    
    print("Process completed.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <genes.gff3> <regions.bed> <output_directory>")
        sys.exit(1)
    
    gff3_file_path = sys.argv[1]
    bed_file_path = sys.argv[2]
    output_directory = sys.argv[3]

    main(gff3_file_path, bed_file_path, output_directory)