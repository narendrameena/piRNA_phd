import argparse
import os
from collections import defaultdict

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        fasta_dict = {}
        header = None
        sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    fasta_dict[header] = sequence
                header = line
                sequence = ""
            else:
                sequence += line
        if header:
            fasta_dict[header] = sequence
    return fasta_dict

def main():
    parser = argparse.ArgumentParser(description='Process FASTA file and group sequences by name prefix.')
    parser.add_argument('fasta_file', help='Input FASTA file')
    parser.add_argument('output_dir', help='Output directory for grouped FASTA files')
    args = parser.parse_args()

    # Ensure output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    sequences = read_fasta(args.fasta_file)
    
    # Group sequences by the prefix before "_markers_"
    grouped_sequences = defaultdict(list)
    for header, sequence in sequences.items():
        prefix = header.split("_markers_")[0][1:]  # Extracting prefix and removing ">"
        grouped_sequences[prefix].append((header, sequence))

    # Write to separate multifasta files in the specified directory based on prefix
    for prefix, seq_list in grouped_sequences.items():
        output_file_name = os.path.join(args.output_dir, f"{prefix}.fasta")
        with open(output_file_name, 'w') as out_file:
            for header, sequence in seq_list:
                out_file.write(header + "\n" + sequence + "\n")

if __name__ == "__main__":
    main()
