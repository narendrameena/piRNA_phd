import argparse
import pybedtools

def process_bed12_file(input_bed, output_bed):
    # Load the BED12 file
    bed = pybedtools.BedTool(input_bed)

    # Find overlapping intervals, preserving BED12 format
    overlaps = bed.cluster()

    # Create a dictionary to map each original interval to its cluster ID
    cluster_map = {}
    for interval in overlaps:
        cluster_id = f"pi-{interval[-1]}"  # The cluster ID is in the last column
        key = (interval.chrom, interval.start, interval.end)
        cluster_map[key] = cluster_id

    # Write original intervals with their cluster ID
    with open(output_bed, 'w') as f_out:
        for interval in bed:
            key = (interval.chrom, interval.start, interval.end)
            pi_name = cluster_map.get(key, "no_overlap")
            # Concatenate the interval fields and the pi_name as a single line
            line = '\t'.join(str(x) for x in interval.fields) + '\t' + pi_name + '\n'
            f_out.write(line)

def main():
    # Argument parser
    parser = argparse.ArgumentParser(description='Process BED12 file to assign unique identifiers to overlapping intervals while keeping original rows.')
    parser.add_argument('input_bed', type=str, help='Input BED12 file path')
    parser.add_argument('output_bed', type=str, help='Output BED12 file path with additional column')

    # Parse arguments
    args = parser.parse_args()

    # Process BED12 file
    process_bed12_file(args.input_bed, args.output_bed)

if __name__ == "__main__":
    main()
