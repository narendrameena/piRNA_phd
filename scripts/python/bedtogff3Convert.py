import argparse

def bed_to_gff(bed_file, gff_file):
    with open(bed_file, 'r') as bed, open(gff_file, 'w') as gff:
        for line in bed:
            if line.strip():
                parts = line.strip().split('\t')

                chrom = parts[0]
                chrom_start = int(parts[1])
                chrom_end = int(parts[2])
                name = parts[3]
                score = parts[4]
                strand = parts[5]

                # Calculate the length of the entry
                entry_length = chrom_end - chrom_start

                # Check if the length is at least 27 bp, not zero or negative, and both chrom_start and chrom_end are not zero
                if entry_length >= 27 and entry_length > 0 and chrom_start != 0 and chrom_end != 0:
                    locus_name = f"{chrom}:{chrom_start}-{chrom_end}:{name}"
                    gff_cols = [
                        chrom,
                        'TRINITY',
                        'piRNA',
                        str(chrom_start),
                        str(chrom_end),
                        score,
                        strand,
                        '.',
                        f"ID={locus_name};Name={name}"
                    ]
                    gff.write('\t'.join(gff_cols) + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert BED file to GFF format.')
    parser.add_argument('input', help='Input BED file')
    parser.add_argument('output', help='Output GFF file')

    args = parser.parse_args()

    bed_to_gff(args.input, args.output)
