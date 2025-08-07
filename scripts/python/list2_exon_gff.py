import argparse

def extract_features_and_introns_from_gff3(input_gff3, output_gff3):
    features = {}
    gene_and_lncRNA_regions = []

    with open(input_gff3, 'r') as file:
        for line in file:
            if not line.startswith("#"):
                parts = line.strip().split("\t")
                chrom, _, feature_type, start, end, _, strand, _, attributes = parts
                
                # Recording gene and lncRNA regions
                if feature_type in ["gene"]:
                    gene_and_lncRNA_regions.append((chrom, int(start), int(end)))
                
                if feature_type in ["exon"]:
                    attributes_dict = {attr.split("=")[0]: attr.split("=")[1] for attr in attributes.split(";")}
                    feature_id = attributes_dict.get("ID", None)
                    parent_id = attributes_dict.get("Parent", None)
                    
                    # Combine the "ID", "Parent" attributes, and locus coordinates under the "ID" field
                    combined_id = f"{feature_id}_{parent_id}_{start}_{end}" if parent_id else f"{feature_id}_{start}_{end}"
                    attributes_dict["ID"] = f"{feature_type}:{combined_id}"
                    
                    modified_attributes = ":".join([f"{k}={v}" for k, v in attributes_dict.items() if k == "ID"])
                    
                    if feature_type not in features:
                        features[feature_type] = []
                    features[feature_type].append((chrom, start, end, strand, modified_attributes))

    # Inferring introns from exons
    if "exon" in features:
        exons = sorted(features["exon"], key=lambda x: x[1])
        for i in range(1, len(exons)):
            chrom, _, _, strand, attributes = exons[i-1]
            _, prev_end, next_start, _, _ = exons[i]
            intron_start = int(prev_end) + 1
            intron_end = int(next_start) - 1
            
            # Check for valid intron coordinates
            if intron_start +10 < intron_end:  # Ensure the intron has positive length
                combined_id = attributes.split(';')[0].split('=')[1].split(':')[1]  # Extract the combined ID
                # Create a combined ID for introns with locus coordinates
                intron_attributes = f"ID=intron:{combined_id}__{chrom}_{intron_start}_{intron_end}"
                if "intron" not in features:
                    features["intron"] = []
                features["intron"].append((chrom, str(intron_start), str(intron_end), strand, intron_attributes))
            else:
                print(f"Warning: Invalid intron coordinates detected between exons ending at {prev_end} and starting at {next_start}.")

    with open(output_gff3, 'w') as output:
        output.write("##gff-version 3\n")
        for feature_type, entries in features.items():
            for entry in entries:
                chrom, start, end, strand, attributes = entry
                start, end = int(start), int(end)
                
                # Check if the feature lies within gene or lncRNA regions
                in_region = any(s <= start <= e and s <= end <= e for c, s, e in gene_and_lncRNA_regions if c == chrom)
                
                if in_region:
                    output.write(f"{chrom}\tcustom\t{feature_type}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract specific features and introns from a GFF3 file.')
    parser.add_argument('input_gff3', help='Input GFF3 file')
    parser.add_argument('output_gff3', help='Output GFF3 file with specified features and introns')
    args = parser.parse_args()
    
    extract_features_and_introns_from_gff3(args.input_gff3, args.output_gff3)
