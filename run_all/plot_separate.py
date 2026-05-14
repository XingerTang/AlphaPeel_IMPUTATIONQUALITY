import pandas as pd
import matplotlib.pyplot as plt


def plot_figure(filename):

    # Define the column names
    column_names = ['SNP_no', 'SNP', 'position', 'N', 'MAF', 'Iam_chance', 'Iam_hwe', 'hiQ', 'accuracy', 'info_IMPUTE', 'r2_MACH', 'r2_BEAGLE']

    # Read the file into a DataFrame, skipping the first two rows (creation time and column names if they're not directly parsed correctly)
    df = pd.read_csv(filename, 
                    sep=r'\s+',  # To handle multiple spaces as a single delimiter
                    skiprows=2,  # Skip the first two rows (creation time and possibly the column names row if it causes parsing issues)
                    names=column_names,  # Use the manually defined column names
                    index_col=False)  # Don't set any column as the index

    # Convert the 'MAF' column from percentage to float
    df['MAF'] = df['MAF'].apply(lambda x: float(x.strip('%')) / 100)

    # Create a figure with multiple subplots
    fig, axs = plt.subplots(8, 1, figsize=(10, 20))

    # Define the columns to plot
    columns_to_plot = ['MAF', 'Iam_chance', 'Iam_hwe', 'hiQ', 'accuracy', 'info_IMPUTE', 'r2_MACH', 'r2_BEAGLE']

    # Iterate over the columns and create a line plot for each
    for i, column in enumerate(columns_to_plot):
        if column == 'accuracy':
            # Create a categorical plot for accuracy
            accuracy_levels = ['cold', 'tepid', 'hot', 'very hot']
            accuracy_map = {level: j for j, level in enumerate(accuracy_levels)}
            df['accuracy_level'] = df['accuracy'].map(accuracy_map)
            axs[i].bar(df['SNP_no'], df['accuracy_level'])
            axs[i].set_yticks(range(len(accuracy_levels)))
            axs[i].set_yticklabels(accuracy_levels)
        else:
            axs[i].plot(df['SNP_no'], df[column])
            axs[i].set_ylabel(column)
        
        axs[i].set_title(f'{column} against Loci')
        axs[i].set_xlabel('Loci')
        axs[i].grid(True)

    # Layout so plots do not overlap
    fig.tight_layout()

    fig.suptitle(f"{filename[:-9]}")

    plt.savefig(f"figures/{filename[:-9]}.png", dpi=300)
    

def main():

    with open("accuracy_list.txt") as fh:
        files = fh.readlines()
    
    for file in files:
        file = file[:-1]
        plot_figure(file)


if __name__ == "__main__":
    main()