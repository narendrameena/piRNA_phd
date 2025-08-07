import sys
import os
import csv
import pysam
import numpy as np
from scipy.stats import zscore
from multiprocessing import Pool, Manager

def check_and_reindex_bam(bam_path):
    bam_file_date = os.path.getmtime(bam_path)
    bai_file_path = bam_path + ".bai"
    if os.path.exists(bai_file_path):
        bai_file_date = os.path.getmtime(bai_file_path)
        if bai_file_date < bam_file_date:
            print("Re-indexing the BAM file as the current index is older...")
            pysam.index(bam_path)
    else:
        print("Index not found for the BAM file. Creating a new one...")
        pysam.index(bam_path)

def process_reads_chunk(start, end, bam_path, min_length, max_length, return_dict):
    bam = pysam.AlignmentFile(bam_path, "rb")
    overlap_count = np.zeros(25)

    for read in bam.fetch(start=start, stop=end):
        if min_length <= len(read.seq) <= max_length:
            for read2 in bam.fetch(start=start, stop=end):
                if read != read2:
                    overlap = get_overlap(read.seq, read2.seq)
                    if overlap > 0:
                        overlap_count[overlap - 1] += 1

    bam.close()
    return_dict[start] = overlap_count

def get_overlap(seq1, seq2):
    reversed_complement_seq2 = complement(seq2[::-1])
    for i in range(25, 0, -1):
        if seq1[:i] == reversed_complement_seq2[-i:]:
            return i
    return 0

def complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join([comp[base] for base in seq])

def get_zscores_from_overlaps(overlap_counts):
    return zscore(overlap_counts)

def main(bam_path, output_filename, min_length, max_length, num_threads):
    try:
        check_and_reindex_bam(bam_path)
        bam = pysam.AlignmentFile(bam_path, "rb")
        total_reads = bam.count()
        bam.close()
    except Exception as e:
        print(f"Error reading or indexing BAM file: {e}")
        sys.exit(1)
    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
        total_reads = bam.count()
        bam.close()
    except Exception as e:
        print(f"Error reading BAM file: {e}")
        sys.exit(1)

    if total_reads == 0:
        print("Error: The BAM file appears to have no readable alignments.")
        sys.exit(1)

    reads_per_chunk = total_reads // num_threads

    with Manager() as manager:
        return_dict = manager.dict()
        with Pool(num_threads) as pool:
            pool.starmap(process_reads_chunk, [(start, start + reads_per_chunk, bam_path, min_length, max_length, return_dict) for start in range(0, total_reads, reads_per_chunk)])

    total_overlap_counts = np.zeros(25)
    for value in return_dict.values():
        total_overlap_counts += value

    overlap_zscores = get_zscores_from_overlaps(total_overlap_counts)

    with open(output_filename, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Overlap_Length", "Count", "Z-score"])
        for i, (count, zscore_val) in enumerate(zip(total_overlap_counts, overlap_zscores)):
            csvwriter.writerow([i + 1, count, zscore_val])

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: script_name.py <BAM_path> <output_filename.csv>  <num_threads>")
        sys.exit(1)
    
    bam_path = sys.argv[1]
    output_filename = sys.argv[2]
    num_threads = int(sys.argv[3])
    min_length = 25
    max_length = 33

    main(bam_path, output_filename, min_length, max_length, num_threads)
