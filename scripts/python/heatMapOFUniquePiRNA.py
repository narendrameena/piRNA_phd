
####### 
#@author narumeena
#@description Python script that  take two tab limited files as argpharse 1)  first file: a) columns are 
# condition b) raws are empty c) each cell value is a piRNA seq 2) the second file: b) columns are condition  
# b) raw are piRNA seq  c) each cell value is  rpm value. Combine both of them based on piRNA seq.  Plot 
# heatmap with wide cell heatmap with the option of changing colour and width and height of the cell where 
# conditions is column and piRNA seq is raw. Also, divide the heatmap vertically based criterion on condition. 
# also, do dendrogram on rows and columns of the heatmap. Save the output from the output of argpharse
#running the script
#python heatmap.py file1.csv file2.csv file3.csv --transcript_ids transcript1 transcript2 --file_name_criteria criteria1 criteria2 --output_file heatmap.png
######



import argparse
import pandas as pd
import seaborn as sns
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

# Parse command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('file1', help='First tab-delimited file')
parser.add_argument('file2', help='Second tab-delimited file')
parser.add_argument('output_file', help='Output file name')
args = parser.parse_args()

# Parse first file and store data in dictionary
file1_dict = {}
with open(args.file1, 'r') as f:
    for line in f:
        condition, seq = line.strip().split('\t')
        file1_dict[seq] = condition

# Parse second file and store data in dictionary
file2_dict = {}
with open(args.file2, 'r') as f:
    for line in f:
        condition, seq, rpm = line.strip().split('\t')
        file2_dict[seq] = rpm

# Combine dictionaries
combined_dict = {}
for seq, condition in file1_dict.items():
    combined_dict[seq] = {'condition': condition, 'rpm': file2_dict[seq]}

# Create DataFrame from combined dictionary
df = pd.DataFrame.from_dict(combined_dict, orient='index')

# Create heatmap
sns.heatmap(df, cmap='YlGnBu', linewidths=.5, annot=True)

# Compute dendrogram for rows and columns
row_dendrogram = hierarchy.dendrogram(hierarchy.linkage(df, method='ward'))
col_dendrogram = hierarchy.dendrogram(hierarchy.linkage(df.T, method='ward'))

# Save figure to file
plt.savefig(args.output_file)
