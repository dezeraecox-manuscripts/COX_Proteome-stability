import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib

from loguru import logger
from GEN_Utils import FileHandling


smoothing_path = f'results/lysate_denaturation/smoothing/processed_data.xlsx'
cluster_path = f'results/lysate_denaturation/clustering/clustered.xlsx'
output_folder = f'results/lysate_denaturation/plot_per_protein/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 14 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

# Choose data type, number of clusters, set colours
cluster_colors = {4: 'royalblue', 2: 'firebrick', 3: 'rebeccapurple', 1: 'darkorange'}
cluster_cmaps = {4: 'Blues', 2: 'Reds', 3: 'Purples', 1: 'Oranges'}


def visualise_transformation(test_dfs, sequences, range_vals=None, output_folder=output_folder):
    """Visualise what each transformation does to the data
    test_df: dict mapping str(data_type): (df, cols)"""
    for sequence in sequences:
        fig, ax = plt.subplots()
        for data_type, (df, cols) in test_dfs.items():
            data = df[df['Sequence'] == sequence]
            if data.shape[0] > 0:
                sns.lineplot(x='channel', y='ratio', data=pd.melt(data, id_vars=None, value_vars=cols, var_name='channel', value_name='ratio'), ci='sd', label=data_type)
        if range_vals:
            plt.ylim(*range_vals)
        plt.legend()
        plt.title(sequence)
        ax.axhline(1, c='r', linestyle='--')
        ax.axhline(0, c='grey', linestyle='--')
        plt.tight_layout()
        plt.savefig(f'{output_folder}{sequence[0:5]}.png')
        plt.savefig(f'{output_folder}{sequence[0:5]}.svg')
        plt.show()
        plt.clf()

def population_plot(df, cols, title, yrange):
    df[cols].T.plot(legend=None)
    plt.ylim(*yrange)
    plt.title(title)
    plt.savefig(f'{output_folder}{title}.png')
    plt.show()


# ------------------------------Plot normalisation for individual peptides------------------------------
# Read in processed data
processed_data = pd.read_excel(f'{smoothing_path}', sheet_name=None)

# Visualise the transformation on individual sequences as examples
test_sequences = ['CTPACVSFGPK', 'AITIAGVPQSVTECVK', 'YTVQDESHSEWVSCVR', 'DVQIGDIVTVGECRPLSK']
dfs_to_plot = ['raw', 'pval_smooth', 'pval_scaled', 'pval_loess']

for_plotting = {}
for key, df in processed_data.items():
    if key in dfs_to_plot:
        for_plotting[key] = (df, [col for col in df.columns.tolist() if type(col) != str])

visualise_transformation(
    test_dfs=for_plotting,
    sequences=test_sequences,
    range_vals=(-2, 2)
    )

# -------------------------Plot example peptides for HSP70 and DNAJB1 peptides-------------------------
# Read in processed data
processed_data = pd.read_excel(f'{cluster_path}', sheet_name=None)

proteins_of_interest = ['P63017', 'Q9QYJ3', 'P14152', 'P08249']
sequences = {sequence: protein for sequence, protein in processed_data['clustered'][['Sequence', 'Proteins']].values if protein in proteins_of_interest}

quant_cols = [col for col in processed_data['sigmoid_summary'].columns.tolist() if type(col) != str]

# plot test sequences
for sequence, df in processed_data['sigmoid_summary'].set_index('Sequence').iterrows():
    if sequence in sequences.keys():
        # generate fitted values
        seq_vals = dict(zip(df.index.tolist(), df.values.tolist()))
        
        # add plot elements
        fig, ax = plt.subplots(figsize=(4, 3))
        # sns.scatterplot(comparison, list(df[comparison_cols]), marker = 'o', label=None, color='lightgrey')
        sns.lineplot(quant_cols, list(df[quant_cols]), label=None, color=cluster_colors[seq_vals['cluster']], linewidth=5)
        plt.legend('')
        plt.xlabel('Urea Concentration (M)')
        plt.ylabel(r'$\Delta$ Corrected cys ratio')
        plt.title(sequence)
        # plt.suptitle(seq_vals['Proteins'])
        plt.ylim(-1.5, 1.5)
        plt.savefig(f'{output_folder}{sequence}.png')
        plt.savefig(f'{output_folder}{sequence}.svg')
        plt.show()

# Generate per-protein plot
for protein in proteins_of_interest:
    # add plot elements
    fig, ax = plt.subplots(figsize=(4, 3))
    for sequence, df in processed_data['sigmoid_summary'].set_index('Sequence').iterrows():
        if sequence in sequences.keys():
            if sequences[sequence] == protein:
                # generate fitted values
                seq_vals = dict(zip(df.index.tolist(), df.values.tolist()))
                
                # sns.scatterplot(comparison, list(df[comparison_cols]), marker = 'o', label=None, color='lightgrey')
                sns.lineplot(quant_cols, list(df[quant_cols]), color=cluster_colors[seq_vals['cluster']], linewidth=5, label=sequence)
    # plt.legend('')
    plt.xlabel('Urea Concentration (M)')
    plt.ylabel(r'$\Delta$ Corrected cys ratio')
    plt.title(protein)
    # plt.suptitle(seq_vals['Proteins'])
    plt.ylim(-1.5, 1.5)
    plt.savefig(f'{output_folder}{protein}.png')
    plt.savefig(f'{output_folder}{protein}.svg')
    plt.show()