import pandas as pd
import matplotlib.pyplot as plt
import plotly.subplots as sp
import plotly.graph_objects as go
import plotly.offline as pyo
import dash
from dash import dcc, html
from dash.dependencies import Input, Output


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

    plt.savefig(f"figures_1/{filename[:-9]}.png", dpi=300)
    

    return df


def plot_multi_file(files, dfs):

    # Define the columns to plot
    columns_to_plot = ['MAF', 'Iam_chance', 'Iam_hwe', 'hiQ', 'accuracy', 'info_IMPUTE', 'r2_MACH', 'r2_BEAGLE']

    # Create a Dash app
    app = dash.Dash(__name__)

    # Define the layout
    app.layout = html.Div([
        html.H1('Interactive Plot'),
        dcc.Dropdown(
            id='column-dropdown',
            options=[{'label': column, 'value': i} for i, column in enumerate(columns_to_plot)],
            value=0,  # default value
            style={'width': '200px'}
        ),
        dcc.Checklist(
            id='file-checklist',
            options=[{'label': file[:-10], 'value': i} for i, file in enumerate(files)],
            value=[i for i in range(len(files))],  # default values
            labelStyle={'display': 'inline-block'}
        ),
        dcc.Graph(id='plot')
    ])

    # Define the callback function
    @app.callback(
        Output('plot', 'figure'),
        [Input('column-dropdown', 'value'),
        Input('file-checklist', 'value')]
    )
    def update_plot(selected_column, selected_files):
        fig = sp.make_subplots(rows=1, cols=1)
        for i, df in enumerate(dfs):
            if i in selected_files:
                filename = files[i][:-10]
                if columns_to_plot[selected_column] == 'accuracy':
                    accuracy_levels = ['cold', 'tepid', 'hot', 'very hot']
                    accuracy_map = {level: j for j, level in enumerate(accuracy_levels)}
                    df['accuracy_level'] = df['accuracy'].map(accuracy_map)
                    fig.add_trace(go.Bar(x=df['SNP_no'], y=df['accuracy_level'], name=filename), row=1, col=1)
                    fig.update_yaxes(ticktext=accuracy_levels, tickvals=[0, 1, 2, 3], row=1, col=1)
                else:
                    fig.add_trace(go.Scatter(x=df['SNP_no'], y=df[columns_to_plot[selected_column]], name=filename), row=1, col=1)
                    fig.update_yaxes(title=columns_to_plot[selected_column], row=1, col=1)
                fig.update_xaxes(title='Loci', row=1, col=1)
                fig.update_layout(title=f'{columns_to_plot[selected_column]} against Loci', showlegend=True)
        return fig

    app.run()

    # Create a figure
    fig = update_plot(0, [i for i in range(len(files))])

    # Save the figure as an HTML file
    pyo.plot(fig, filename='figures_1/all.html', auto_open=False)


def main():

    with open("accuracy_list.txt") as fh:
        files = fh.readlines()

    dfs = []
    
    for file in files:
        file = file[:-1]
        dfs.append(plot_figure(file))

    plot_multi_file(files, dfs)


if __name__ == "__main__":
    main()