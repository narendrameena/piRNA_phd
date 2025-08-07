#!/usr/bin/env python


####### 
#@author narumeena
#@description This script is a command-line tool that analyzes a BAM/SAM file and generates a bar 
# plot showing the proportions of each nucleotide (A, C, G, T, and N) at each position in the 
# aligned sequences. 
#python script.py --bam_file path/to/bam_file --output_file path/to/output_file --statistics_file path/to/statistics_file [-m min_reads]

######



import sys
import pysam
import matplotlib.pyplot as plt

# Parse the input BAM/SAM file
bam_file = pysam.AlignmentFile(sys.argv[1], "rb")

# Initialize dictionaries to store the counts of each nucleotide at each position
counts = {}
totals = {}

# Iterate through the reads in the BAM/SAM file
for read in bam_file.fetch():
    # Iterate through the nucleotides in the read
    for i in range(read.query_length):
        # Get the nucleotide at the current position
        nucleotide = read.query_sequence[i]
        # Increment the count for this nucleotide at this position
        if i in counts:
            if nucleotide in counts[i]:
                counts[i][nucleotide] += 1
            else:
                counts[i][nucleotide] = 1
        else:
            counts[i] = {nucleotide: 1}
        # Increment the total count for this position
        if i in totals:
            totals[i] += 1
        else:
            totals[i] = 1

# Calculate the average proportions of each nucleotide at each position
proportions = {}
for i in counts:
    proportions[i] = {}
    for nucleotide in counts[i]:
        proportions[i][nucleotide] = counts[i][nucleotide] / totals[i]

# Get the list of positions
positions = list(proportions.keys())

# Get the list of nucleotides
nucleotides = list(proportions[positions[0]].keys())

# Initialize the plot
fig, ax = plt.subplots()

# Set the bar width
bar_width = 0.1

# Set the colors for the bars
colors = ['b', 'g', 'r', 'c', 'm']

# Set the X-axis tick marks to be the positions
plt.xticks(positions)

# Set the X-axis label
plt.xlabel("Position")

# Set the Y-axis label
plt.ylabel("Proportion")

# Set the title
plt.title("Average Proportions of Nucleotides")

# Initialize the X-axis offset
offset = 0

# Iterate through the nucleotides
for i in range(len(nucleotides)):
    # Get the current nucleotide
    nucleotide = nucleotides[i]
    # Get the proportions for this nucleotide
    values = [proportions[pos][nucleotide] for pos in positions]
    # Create a bar plot for this nucleotide
    plt.bar(positions, values, bar_width, bottom=offset, color=colors[i])
    # Update the X-axis offset
    offset += values

# Create a list of patches for the legend
patches = []
for i in range(len(nucleotides)):
    patches.append(mpatches.Patch(color=colors[i], label=nucleotides[i]))

# Add the legend to the plot
plt.legend(handles=patches)

# Save the plot to the specified file
plt.savefig(sys.argv[2])

# Open the output file for writing
with open(sys.argv[3], 'w') as f:
    # Write the header line
    f.write("position,A,C,G,U,N\n")
    # Iterate through the positions
    for pos in positions:
        # Write the position and the count for each nucleotide
        f.write(str(pos))
        for nucleotide in nucleotides:
            f.write("," + str(counts[pos][nucleotide]))
        f.write("\n")

# Close the BAM/SAM file
bam_file.close()
