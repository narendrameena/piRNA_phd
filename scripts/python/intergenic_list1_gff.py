import argparse

def extract_entities_from_gff3(input_gff3, output_gff3):
    chromosomes = {}
    output_entries = []

    with open(input_gff3, 'r') as file:
        for line in file:
            if not line.startswith("#"):
                parts = line.strip().split("\t")
                chrom, source, feature_type, start, end, score, strand, phase, attributes = parts

                # Parse the ID or Parent attribute, depending on feature type
                attributes_dict = {attr.split("=")[0]: attr.split("=")[1] for attr in attributes.split(";")}
                if feature_type == "gene":
                    identifier = attributes_dict.get("ID", None)
                    gene_name = attributes_dict.get("Name", None)  # Get the "Name" attribute for genes
                elif feature_type == "lnc_RNA":
                    identifier = attributes_dict.get("Parent", None)
                    gene_name = None
                else:
                    continue

                # Adjusting the identifier by prepending the type and including both ID and Parent attributes
                identifier = f"{feature_type}:{identifier}"
                attributes_dict["ID"] = f"{feature_type}_{attributes_dict.get('ID', '')}_{attributes_dict.get('Parent', '')}"

                # Add the "Name" attribute for genes if available
                if gene_name:
                    attributes_dict["Name"] = gene_name

                # Adjusted attributes line
                adjusted_attributes = ":".join([f"{key}={value}" for key, value in attributes_dict.items()])

                # Modify this line to include the "Name" attribute after the feature type
                output_entries.append(f"{chrom}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{adjusted_attributes}")
                #
                if chrom not in chromosomes:
                    chromosomes[chrom] = []
                chromosomes[chrom].append((int(start), int(end), identifier, strand))

    # Now, compute the intergenic regions
    for chrom, entities in chromosomes.items():
        sorted_entities = sorted(entities, key=lambda x: x[0])

        for i in range(1, len(sorted_entities)):
            # Ensure entities are on the same strand
            if sorted_entities[i-1][3] != sorted_entities[i][3]:
                continue

            intergenic_start = sorted_entities[i-1][1] + 1
            intergenic_end = sorted_entities[i][0] - 1
            
            # Ensure the intergenic region has positive length
            if intergenic_start + 10 < intergenic_end:
                identifier = f"intergenic:{sorted_entities[i][2]}"  # Using the identifier of the upcoming gene/lncRNA
                output_entries.append(f"{chrom}\tcustom\tintergenic\t{intergenic_start}\t{intergenic_end}\t.\t{sorted_entities[i][3]}\t.\tID={identifier}")
            else:
                print(f"Warning: Detected zero or negative-length intergenic region on {chrom} between {intergenic_start} and {intergenic_end}. Skipping...")

    # Write to the output GFF3
    with open(output_gff3, 'w') as output:
        output.write("##gff-version 3\n")
        for entry in output_entries:
            output.write(entry + "\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate GFF3 file with gene, lnc_RNA, and intergenic regions.')
    parser.add_argument('input_gff3', help='Input GFF3 file')
    parser.add_argument('output_gff3', help='Output GFF3 file')
    args = parser.parse_args()
    
    extract_entities_from_gff3(args.input_gff3, args.output_gff3)
