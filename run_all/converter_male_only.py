import pandas as pd
import argparse
import subprocess

output_filename = "converted_multi_x_chr_male_switched.txt"

def convert_file_format(input_filename, output_filename, ped_filename, switched=True):
    data = pd.read_csv(input_filename, sep=r'\s+', header=None)
    num_loci = len(data.columns) - 1
    pedigree = pd.read_csv(ped_filename, sep=r'\s+', header=None)

    males = pedigree[pedigree.iloc[:, 3] == 0].iloc[:, 0].tolist()

    num_females = len(males)

    dfs = []

    if switched: 
        for i in range(num_loci):
            locus_values = []
            for j in range(num_females):
                ind_id = males[j]
                aa_prob = data.iloc[(ind_id - 1)*3, i+1]
                hetero_prob = data.iloc[(ind_id - 1)*3 + 1, i+1]
                AA_prob = data.iloc[(ind_id - 1)*3 + 2, i+1]
                individual_values = [aa_prob, AA_prob, hetero_prob]
                locus_values.extend(individual_values)
            dfs.append(pd.DataFrame([[f'loci_{i+1}', f'loci_{i+1}', i+1] + locus_values], columns=['loci_name', 'loci_name_2', 'loci_id'] + [f'individual_{k}_value_{m}' for k in males for m in range(3)]))
    else:
        for i in range(num_loci):
            locus_values = []
            for j in range(num_females):
                ind_id = males[j]
                individual_values = data.iloc[(ind_id - 1)*3:ind_id*3, i+1].tolist()
                locus_values.extend(individual_values)
            dfs.append(pd.DataFrame([[f'loci_{i+1}', f'loci_{i+1}', i+1] + locus_values], columns=['loci_name', 'loci_name_2', 'loci_id'] + [f'individual_{k}_value_{m}' for k in males for m in range(3)]))

    converted_data = pd.concat(dfs, ignore_index=True)

    converted_data.to_csv(output_filename, sep=' ', header=False, index=False)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-input_file", type=str, required=True)
    parser.add_argument("-output_file", type=str, required=True)
    parser.add_argument("-ped_file", type=str, required=True)
    parser.add_argument("-params_file", type=str, required=True)
    parser.add_argument("-program_path", type=str, required=True)
    parser.add_argument("-switched", type=int, required=True)

    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    ped_file = args.ped_file
    params_file = args.params_file
    program_path = args.program_path
    switched = bool(args.switched)

    convert_file_format(input_file, output_file, ped_file, switched)
    subprocess.run([
        "python",
        program_path,
        "-f",
        params_file
    ], input=b"\n", stdout=subprocess.PIPE)


if __name__ == "__main__":
    main()