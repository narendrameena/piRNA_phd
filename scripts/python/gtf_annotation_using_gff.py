import pybedtools
import argparse
import os
import tempfile

def annotate_piRNAs(gff3_file, piRNA_gtf_file, output_file):
    # Load GFF3 and piRNA GTF files as BedTools objects
    gff3_bed = pybedtools.BedTool(gff3_file)
    
    # Create a temporary tabix index file for the piRNA GTF file
    with tempfile.NamedTemporaryFile(suffix=".gtf.gz", delete=False) as tmp_index:
        piRNA_bed = pybedtools.BedTool(piRNA_gtf_file).saveas(tmp_index.name).tabix()

    # Annotate piRNAs based on overlapping features in GFF3
    annotated_overlapping = piRNA_bed.intersect(gff3_bed, wo=True, stream=True)

    # Filter piRNAs that were not annotated by overlap
    unannotated_piRNAs = piRNA_bed.intersect(gff3_bed, v=True, stream=True)

    # Annotate remaining piRNAs based on closest features in GFF3
    annotated_closest = unannotated_piRNAs.closest(gff3_bed, d=True, stream=True)

    # Merge annotations
    annotated_piRNAs = annotated_overlapping.cat(annotated_closest, postmerge=False)

    # Generate a new output file name if the specified file already exists
    output_path, output_ext = os.path.splitext(output_file)
    new_output_file = output_file
    count = 1
    while os.path.exists(new_output_file):
        new_output_file = f"{output_path}_annotated_{count}{output_ext}"
        count += 1

    annotated_piRNAs.saveas(new_output_file)

def main():
    parser = argparse.ArgumentParser(description='piRNA annotation based on GFF3 file')
    parser.add_argument('gff3_file', help='Input GFF3 file for annotation')
    parser.add_argument('piRNA_gtf_file', help='Input GTF-like piRNA file')
    parser.add_argument('output_file', help='Output annotated piRNA file')
    args = parser.parse_args()

    annotate_piRNAs(args.gff3_file, args.piRNA_gtf_file, args.output_file)

if __name__ == "__main__":
    main()
