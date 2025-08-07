import pandas as pd
import sys

# Get input and output paths from Snakemake
excel_file = sys.argv[1]
bed_file = sys.argv[2]

# Load the third sheet (index 2, 0-based indexing)
sheet_name = 2  # Third sheet
df = pd.read_excel(excel_file, sheet_name=sheet_name)

# Define BED6 columns
bed6_columns = ['seqnames', 'start', 'end', 'type', 'width', 'strand']

# Identify extra columns to include as attributes
extra_columns = [col for col in df.columns if col not in bed6_columns]

# Combine extra columns into a single attributes string
df['attributes'] = df[extra_columns].apply(
    lambda row: ';'.join([f"{col}={row[col]}" for col in extra_columns]), axis=1
)

# Select only BED6 columns and the attributes column
bed_with_attributes = df[bed6_columns + ['attributes']]

# Save the DataFrame to a BED file
bed_with_attributes.to_csv(bed_file, sep='\t', index=False, header=False)