
#!/usr/bin/env python
####### 
#@author narumeena
#@description This script reads a GFF3 file and converts it into a BED12 file. GFF3 
# (General Feature Format version 3) is a standard file format used to store genomic 
# feature annotations, while BED12 (Browser Extensible Data version 12) is another 
# standard file format used to represent genomic intervals with additional annotation fields.
#python gffToBedFile.py in.gff3 > output.bed12


#!/usr/bin/env python3

import logging

# Set the logging level
logging.basicConfig(level=logging.INFO)

# Open the input GFF3 file for reading
with open(sys.argv[1], 'r') as gff_file:
    # Iterate over each line in the file
    for line in gff_file:
        # Skip comments and blank lines
        if line.startswith('#') or line.strip() == '':
            continue

        # Split the line on tabs to get the GFF3 fields
        fields = line.split('\t')

        # Extract the relevant fields: seqid, start, end, strand, and attributes
        seqid = fields[0]
        start = fields[3]
        end = fields[4]
        strand = fields[6]
        attributes = fields[8]

        # Parse the attributes field to extract the ID and Name attributes
        attributes_dict = {}
        for attribute in attributes.split(';'):
            if attribute.strip() == '':
                continue
            key, value = attribute.split('=')
            attributes_dict[key.strip()] = value.strip()

        # Build the BED12 line by concatenating the fields
        bed12_line = '\t'.join([seqid, start, end, attributes_dict['ID'], '0', strand, start, end, '0', '1', '1', '1'])

        # Print the BED12 line to standard output
        print(bed12_line)
        logging.info(f"Converted line: {bed12_line}")
