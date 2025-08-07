import argparse

def extract_gene_features_from_gff3(input_gff3, output_gff3):
    features = {}
    gene_info = {}
    transcript_to_gene = {}
    gene_regions = {}

    with open(input_gff3, 'r') as file:
        for line in file:
            if not line.startswith("#"):
                parts = line.strip().split("\t")
                chrom, _, feature_type, start, end, _, strand, _, attributes = parts
                
                # Extracting attributes
                attributes_dict = {attr.split("=")[0]: attr.split("=")[1] for attr in attributes.split(";") if "=" in attr}
                feature_id = attributes_dict.get("ID", None)
                parent_id = attributes_dict.get("Parent", None)

                # Store gene information, transcript-to-gene mapping, and gene regions
                if feature_type == "gene":
                    gene_id = feature_id.split(":")[1]
                    gene_info[gene_id] = attributes_dict
                    gene_regions[(chrom, gene_id)] = (int(start), int(end))
                elif feature_type == "transcript":
                    gene_id = parent_id.split(":")[1]
                    transcript_id = feature_id.split(":")[1]
                    transcript_to_gene[transcript_id] = gene_id

                # Process exons, UTRs, and other features linked to transcripts
                if feature_type in ["exon", "three_prime_UTR", "five_prime_UTR"] and parent_id:
                    transcript_id = parent_id.split(":")[1]
                    gene_id = transcript_to_gene.get(transcript_id, None)
                    if gene_id and gene_id in gene_info:
                        # Add gene 'Name' to feature attributes if available
                        if "Name" in gene_info[gene_id]:
                            attributes_dict["Name"] = gene_info[gene_id]["Name"]

                    combined_id = f"{feature_id}_{start}_{end}"
                    attributes_dict["ID"] = f"{feature_type}:{combined_id}"

                    modified_attributes = ";".join([f"{k}={v}" for k, v in attributes_dict.items()])
                        
                    if feature_type not in features:
                        features[feature_type] = []
                    features[feature_type].append((chrom, start, end, strand, modified_attributes))

    # Inferring introns from exons
    if "exon" in features:
        exons = sorted(features["exon"], key=lambda x: (x[0], int(x[1])))
        for i in range(1, len(exons)):
            chrom, _, _, strand, attributes = exons[i-1]
            _, prev_end, next_start, _, _ = exons[i]
            intron_start = int(prev_end) + 1
            intron_end = int(next_start) - 1

            # Check for valid intron coordinates
            if intron_start < intron_end:
                combined_id = attributes.split(';')[0].split('=')[1].split(':')[1]  # Extract the combined ID
                intron_attributes = f"ID=intron:{combined_id}__{chrom}_{intron_start}_{intron_end}"
                if "intron" not in features:
                    features["intron"] = []
                features["intron"].append((chrom, str(intron_start), str(intron_end), strand, intron_attributes))

    # Writing the output with gene region check
    with open(output_gff3, 'w') as output:
        output.write("##gff-version 3\n")
        for feature_type, entries in features.items():
            for entry in entries:
                chrom, start, end, strand, attributes = entry
                start, end = int(start), int(end)

                # Check if the feature lies within gene regions
                in_region = any(start >= s and end <= e for (c, _), (s, e) in gene_regions.items() if c == chrom)
                
                if in_region:
                    output.write(f"{chrom}\tcustom\t{feature_type}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract exons, UTRs, and introns exclusively for genes from a GFF3 file.')
    parser.add_argument('input_gff3', help='Input GFF3 file')
    parser.add_argument('output_gff3', help='Output GFF3 file with specified features for genes only')
    args = parser.parse_args()
    
    extract_gene_features_from_gff3(args.input_gff3, args.output_gff3)
