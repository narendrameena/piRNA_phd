
####### 
#@author narumeena
#@description python script to plot circos plot visualisation of of many  bamFiles
#running the script
#python plot_circos.py  --bam_files file1.bam file2.bam file3.bam -o output.png

######

import argparse
import pysam
import pycircos

# Parse the command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--bam_files', nargs='+', help='BAM files to plot')
parser.add_argument('--output', '-o', help='Output plot file')
args = parser.parse_args()

# Read in the BAM files and extract the necessary information
data = []
for bam_file in args.bam_files:
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for read in bam:
            data.append({
                'start': read.reference_start,
                'end': read.reference_end,
                'ref': read.reference_name,
                'read': read.query_sequence
            })

# Create the Circos plot
plot = pycircos.CircosPlot(data, stroke_color='black')

# Add labels to the plot
plot.label(r0='1.1r', r1='1.2r', position='center', color='black', text=args.bam_files)

# Save the plot to the output file
if args.output:
    plot.save(args.output)
else:
    plot.show()





