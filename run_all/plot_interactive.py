import pandas as pd
import plotly.subplots as sp
import plotly.graph_objects as go
import plotly.offline as pyo
import dash
from dash import dcc, html
from dash.dependencies import Input, Output


with open("accuracy_list.txt") as fh:
    filename = fh.readlines()


dfs = []


for file in filename:

    file = file[:-1]
    # Define the column names
    column_names = ['SNP_no', 'SNP', 'position', 'N', 'MAF', 'Iam_chance', 'Iam_hwe', 'hiQ', 'accuracy', 'info_IMPUTE', 'r2_MACH', 'r2_BEAGLE']

    # Read the file into a DataFrame, skipping the first two rows (creation time and column names if they're not directly parsed correctly)
    df = pd.read_csv(file, 
                    sep=r'\s+',  # To handle multiple spaces as a single delimiter
                    skiprows=2,  # Skip the first two rows (creation time and possibly the column names row if it causes parsing issues)
                    names=column_names,  # Use the manually defined column names
                    index_col=False)  # Don't set any column as the index

    # Convert the 'MAF' column from percentage to float
    df['MAF'] = df['MAF'].apply(lambda x: float(x.strip('%')) / 100)

    dfs.append(df)


for i, file in enumerate(filename):
    filename[i] = file[10:-10]



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
        options=[{'label': file, 'value': i} for i, file in enumerate(filename)],
        value=[i for i in range(len(filename))],  # default values
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
    global filename

    fig = sp.make_subplots(rows=1, cols=1)
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'black', 'gray', 'yellow']
    for i, df in enumerate(dfs):
        if i in selected_files:
            if columns_to_plot[selected_column] == 'accuracy':
                accuracy_levels = ['cold', 'tepid', 'hot', 'very hot']
                accuracy_map = {level: j for j, level in enumerate(accuracy_levels)}
                df['accuracy_level'] = df['accuracy'].map(accuracy_map)
                fig.add_trace(go.Bar(x=df['SNP_no'], y=df['accuracy_level'], name=filename[i]), row=1, col=1)
                fig.update_yaxes(ticktext=accuracy_levels, tickvals=[0, 1, 2, 3], row=1, col=1)
            else:
                fig.add_trace(go.Scatter(mode='markers', x=df['SNP_no'], y=df[columns_to_plot[selected_column]], marker=dict(color=colors[i % len(colors)]), name=filename[i]), row=1, col=1)
                fig.update_yaxes(title=columns_to_plot[selected_column], row=1, col=1)
            fig.update_xaxes(title='Loci', row=1, col=1)
            fig.update_layout(title=f'{columns_to_plot[selected_column]} against Loci', showlegend=True)
    return fig


# Create a figure
fig = update_plot(0, [i for i in range(len(filename))])

# Save the figure as an HTML file
pyo.plot(fig, filename='figures/all.html', auto_open=False)


if __name__ == "__main__":
    app.run(debug=True)
    