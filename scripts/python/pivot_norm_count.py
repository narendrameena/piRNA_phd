import csv
import sys

def convert_table_to_csv(input_file, output_file):
    try:
        # Read data from the input file
        with open(input_file, mode='r') as input_csv:
            reader = csv.reader(input_csv)
            data = [row for row in reader]

        # Transpose the data
        transposed_data = list(map(list, zip(*data)))

        # Write to CSV
        with open(output_file, mode='w', newline='') as output_csv:
            writer = csv.writer(output_csv)
            writer.writerow(["Gene_Name", "Condition", "Norm_count"])
            for row in transposed_data[1:]:
                for i in range(1, len(row)):
                    writer.writerow([transposed_data[0][i], row[0], row[i]])

        print(f"Conversion successful. Data written to {output_file}")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file.csv> <output_file.csv>")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_path = sys.argv[2]

    convert_table_to_csv(input_file_path, output_file_path)
