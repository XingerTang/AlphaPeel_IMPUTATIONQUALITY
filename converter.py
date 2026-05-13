import pandas as pd
import argparse
import subprocess

def convert_file_format(input_filename, output_filename):
    # Read the input file
    data = pd.read_csv(input_filename, sep=r'\s+', header=None)

    # Get the number of columns (loci) and the number of individuals
    num_loci = len(data.columns) - 1
    num_individuals = len(data) // 3

    # Create a list to store the DataFrames
    dfs = []

    # Convert the data
    for i in range(num_loci):
        locus_values = []
        for j in range(num_individuals):
            individual_values = data.iloc[j*3:(j+1)*3, i+1].tolist()
            locus_values.extend(individual_values)
        dfs.append(pd.DataFrame([[f'loci_{i+1}', f'loci_{i+1}', i+1] + locus_values], columns=['loci_name', 'loci_name_2', 'loci_id'] + [f'individual_{k}_value_{m}' for k in range(num_individuals) for m in range(3)]))

    # Concatenate the DataFrames
    converted_data = pd.concat(dfs, ignore_index=True)

    # Write the converted data to the output file
    converted_data.to_csv(output_filename, sep=' ', header=False, index=False)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-input_file", type=str, required=True)
    parser.add_argument("-output_file", type=str, required=True)
    parser.add_argument("-params_file", type=str, required=True)
    parser.add_argument("-program_path", type=str, required=True)

    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    params_file = args.params_file
    program_path = args.program_path

    convert_file_format(input_file, output_file)
    subprocess.run([
        "python",
        program_path,
        "-f",
        params_file
    ], input=b"\n", stdout=subprocess.PIPE)


if __name__ == "__main__":
    main()