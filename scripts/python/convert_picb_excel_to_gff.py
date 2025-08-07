import pandas as pd
import sys

# Get input and output paths from Snakemake
excel_file = sys.argv[1]
gff_file = sys.argv[2]

# Load the third sheet (index 2, 0-based indexing)
sheet_name = 2  # Third sheet
df = pd.read_excel(excel_file, sheet_name=sheet_name)

# Define mandatory GFF3 columns
gff_columns = ['seqnames', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase']

# Ensure the GFF3 mandatory fields are present or create defaults
df['source'] = 'custom_script'  # Example source, modify as needed
df['score'] = '.'  # Placeholder for score if not provided
df['phase'] = '.'  # Placeholder for phase if not provided

# Identify extra columns to include as attributes
extra_columns = [col for col in df.columns if col not in gff_columns]

# Combine extra columns into a single attributes string
df['attributes'] = df[extra_columns].apply(
    lambda row: ';'.join([f"{col}={row[col]}" for col in extra_columns]), axis=1
)

# Reorder the DataFrame to match GFF3 format
gff_df = df[['seqnames', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']]

# Save the DataFrame to a GFF3 file
gff_df.to_csv(gff_file, sep='\t', index=False, header=False)

print(f"GFF3 file saved to: {gff_file}")