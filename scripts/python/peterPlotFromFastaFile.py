#!/usr/bin/env python

####### 
#@author narumeena
#@description This modified script will read a FASTA file, extract the sequences and their corresponding lengths, 
# calculate the nucleotide proportions for each sequence, filter the sequences by the minimum number of reads, 
# save the statistics to a tab-delimited file, and generate a bar plot showing the proportions of each nucleotide 
# at each position in the aligned sequences.
#running the script
#python script.py --fasta_file path/to/fasta_file --output_file path/to/output_file --statistics_file path/to/statistics_file [-m min_reads]

######


#!/usr/bin/env python

import argparse
import logging
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import Counter

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def peterPlotFromFASTA():
    # Parse the command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_file", required=True, help="Path to the FASTA file")
    parser.add_argument("--output_file", required=True, help="Path to the output file")
    args = parser.parse_args()

    fasta_file_path = args.fasta_file
    output_file = args.output_file

    # Open the FASTA file
    logging.info("Parsing FASTA file")
    fasta_file = open(fasta_file_path, "r")

    # Use the SeqIO.parse function to parse the FASTA file
    sequences = []
    lengths = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Extract the sequence and its corresponding length
        sequence = read.query_sequence
        length = len(sequence)
        sequences.append(sequence)
        lengths.append(length)

    # Calculate the nucleotide proportions for each sequence
    proportions = []
    for sequence in sequences:
        # Use a Counter object to count the nucleotides
        nucleotide_counts = Counter(sequence)

        # Add all the keys if they are not present in the Counter object
        for nucleotide in ['A', 'T', 'G', 'C', 'N']:
            if nucleotide not in nucleotide_counts:
                nucleotide_counts[nucleotide] = 0.0

        # Replace 'T' with 'U'
        if 'T' in nucleotide_counts:
            nucleotide_counts['U'] = nucleotide_counts.pop('T')

        total_nucleotides = sum(nucleotide_counts.values())
        nucleotide_proportions = {nucleotide: count / total_nucleotides for nucleotide, count in nucleotide_counts.items()}
        proportions.append(nucleotide_proportions)

    # Calculate the mean nucleotide proportions
    mean_proportions = {}
    for i in range(len(sequences[0])):
        # Use a list comprehension to extract the i-th nucleotide from each sequence
        nucleotides = [sequence[i] for sequence in sequences]

        # Use a Counter object to count the nucleotides
        nucleotide_counts = Counter(nucleotides)

        # Add all the keys if they are not present in the Counter object
        for nucleotide in ['A', 'T', 'G', 'C', 'N']:
            if nucleotide not in nucleotide_counts:
                nucleotide_counts[nucleotide] = 0.0

        # Replace 'T' with 'U'
        if 'T' in nucleotide_counts:
            nucleotide_counts['U'] = nucleotide_counts.pop('T')

        total_nucleotides = sum(nucleotide_counts.values())
        nucleotide_proportions = {nucleotide: count / total_nucleotides for nucleotide, count in nucleotide_counts.items()}
        mean_proportions[i] = nucleotide_proportions

    # Generate the bar plot
    logging.info("Generating bar plot")
    fig, ax = plt.subplots()

    # Set the bar width
    bar_width = 0.8

    # Set the positions of the bars on the x-axis
    positions = range(len(mean_proportions))

    # Iterate over the nucleotides
    for nucleotide in ['A', 'C', 'G', 'U', 'N']:
        # Use a list comprehension to extract the nucleotide proportions for each position
        nucleotide_proportions = [mean_proportions[position][nucleotide] for position in positions]

        # Plot the bar for the current nucleotide
        ax.bar(positions, nucleotide_proportions, bar_width, label=nucleotide)

    # Add the x-tick labels
    plt.xticks(positions, positions)

    # Add the y-tick labels
    plt.yticks(np.arange(0, 1.1, 0.1))

