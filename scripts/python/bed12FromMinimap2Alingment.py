

####### 
#@author narumeena
#@description This script is a command-line tool that aligns two FASTA files using the minimap2 alignment 
# tool, and converts the alignment output to BED12 format.
#running the script
#python align_fasta.py --fasta1 path/to/fasta1.fa --fasta2 path/to/fasta2.fa --output path/to/output.bed
######


import argparse
import logging
import subprocess
import tempfile

def main():
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description='Align two FASTA files using minimap2 and save the output in BED12 format')
    parser.add_argument('--fasta1', required=True, help='Path to the refrance FASTA file')
    parser.add_argument('--fasta2', required=True, help='Path to the second FASTA file')
    parser.add_argument('--output', required=True, help='Path to the output BED12 file')
    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(level=logging.INFO)

    # Create a temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            # Run minimap2 to align the FASTA files
            minimap2_cmd = ['minimap2', '-x', 'asm5', '-c', args.fasta1, args.fasta2]
            result = subprocess.run(minimap2_cmd, capture_output=True)

            # Parse the minimap2 output and convert it to BED12 format
            bed12_lines = []
            for line in result.stdout.decode().splitlines():
                # Skip lines that don't contain alignments
                if line[0] == '@':
                    continue

                # Split the line into fields
                fields = line.split('\t')

                # Extract the relevant fields for BED12 format
                chrom = fields[5]
                start = int(fields[7])
                end = int(fields[8])
                name = fields[0]
                score = 1000
                strand = fields[4]
                thick_start = start
                thick_end = end
                item_rgb = '0,0,0'
                block_count = 1
                block_sizes = [end - start]
                block_starts = [0]

                # Convert the fields to a BED12 line
                bed12_line = '\t'.join([chrom, str(start), str(end), name, str(score), strand,
                                        str(thick_start), str(thick_end), item_rgb,
                                        str(block_count), ','.join(map(str, block_sizes)), ','.join(map(str, block_starts))])

                # Add the BED12 line to the list
                bed12_lines.append(bed12_line)

            # Write the BED12 lines to the output file
            with open(args.output, 'w') as f:
                f.write('\n'.join(bed12_lines))

        except FileNotFoundError as e:
            # Log the exception and exit with an            
            # 
            # Write the BED12 lines to the output file
            with open(args.output, 'w') as f:
                f.write('\n'.join(bed12_lines))

        except FileNotFoundError as e:
            # Log the exception and exit with an error code
            logging.error(e)
            sys.exit(1)
        except Exception as e:
            # Log the exception and exit with an error code
            logging.error(e)
            sys.exit(1)

if __name__ == '__main__':
    main()

