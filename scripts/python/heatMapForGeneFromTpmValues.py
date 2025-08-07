
####### 
#@author narumeena
#@description This script takes in one or more input CSV files, filters the rows based on the specified 
# transcript IDs, divides the data into separate dataframes based on the file name criteria, pivots the 
# data, calculates the correlation matrix and p-values, and creates a heatmap using the correlation data 
# and p-values.
#running the script
#python heatmap.py file1.csv file2.csv file3.csv --transcript_ids transcript1 transcript2 --file_name_criteria criteria1 criteria2 --output_file heatmap.png
######



import argparse
import pandas as pd
import matplotlib.pyplot as plt
import logging
import seaborn as sns
import numpy as np
import scipy.stats as stats

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', help='input files')
parser.add_argument('--transcript_ids_1', nargs='+', help='transcript IDs to plot in heatmap 1')
parser.add_argument('--transcript_ids_2', nargs='+', help='transcript IDs to plot in heatmap 2')
parser.add_argument('--file_name_criteria', nargs='+', help='file name criteria for dividing the heatmap vertically')
parser.add_argument('--output_file', help='output file name')
parser.add_argument('--transcript_id_col', default='Transcript ID', help='column name for transcript IDs')
parser.add_argument('--tpm_count_col', default='TPM Count', help='column name for TPM counts')
args = parser.parse_args()

try:
    logging.info('Started processing input files')

    df_list = []
    for file in args.files:
        df = pd.read_csv(file)
        df_sub = df[[args.transcript_id_col, args.tpm_count_col]]
        df_list.append(df_sub)

    df_all = pd.concat(df_list)

    # Select only the rows with the specified transcript IDs
    df_filtered_1 = df_all[df_all[args.transcript_id_col].isin(args.transcript_ids_1)]
    df_filtered_2 = df_all[df_all[args.transcript_id_col].isin(args.transcript_ids_2)]

    # Divide the data into separate dataframes based on the file name criteria
    df_list_1 = []
    df_list_2 = []
    for criteria in args.file_name_criteria:
        df_sub_1 = df_filtered_1[df_filtered_1.index.str.contains(criteria)]
        df_sub_2 = df_filtered_2[df_filtered_2.index.str.contains(criteria)]
        df_list_1.append(df_sub_1)
        df_list_2.append(df_sub_2)

    df_pivot_list_1 = []
    df_pivot_list_2 = []
    for df in df_list_1:
        df_pivot = df.pivot(index=args.transcript_id_col, columns=df.index, values=args.tpm_count_col)
        df_pivot_list_1.append(df_pivot)
    for df in df_list_2:
        df_pivot = df.pivot(index=args.transcript_id_col, columns=df.index, values=args.tpm_count_col)
        df_pivot_list_2.append(df_pivot)

        # Create the figure and the subplots
    fig, axs = plt.subplots(2, len(df_pivot_list_1), figsize=(20, 10))

    # Loop through each pivoted dataframe and plot the heatmap in each subplot
    for i, df_pivot in enumerate(df_pivot_list_1):
        # Calculate the correlation matrix and p-values for the pivoted data
        corr = df_pivot.corr()
        p_values = np.empty_like(corr)
        for i in range(len(corr)):
            for j in range(len(corr)):
                p_values[i, j] = stats.pearsonr(corr.iloc[i], corr.iloc[j])[1]

        # Create the heatmap
        sns.clustermap(corr, annot=True, mask=p_values < 0.01, cmap='coolwarm', linewidths=0.5, square=True, ax=axs[0, i])
        axs[0, i].set_title(args.file_name_criteria[i])
    for i, df_pivot in enumerate(df_pivot_list_2):
        # Calculate the correlation matrix and p-values for the pivoted data
        corr = df_pivot.corr()
        p_values = np.empty_like(corr)
        for i in range(len(corr)):
            for j in range(len(corr)):
                p_values[i, j] = stats.pearsonr(corr.iloc[i], corr.iloc[j])[1]

        # Create the heatmap
        sns.clustermap(corr, annot=True, mask=p_values < 0.01, cmap='coolwarm', linewidths=0.5, square=True, ax=axs[1, i])
        axs[1, i].set_title(args.file_name_criteria[i])

    # Adjust the layout and save the figure
    plt.tight_layout()
    plt.savefig(args.output_file)

    logging.info('Completed creating correlation plot')

except Exception as e:
    logging.exception('Error occurred while processing the input files')







