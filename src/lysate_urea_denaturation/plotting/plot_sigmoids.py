import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

from GEN_Utils import FileHandling
from loguru import logger

logger.info('Import ok')

def denaturant(urea, top, bottom, cM, m):
    # adapted from https://en.wikipedia.org/wiki/Equilibrium_unfolding, by keeping terms for bottom as in Savitski, subbing deltaG into standard equation, and reintroducing bottom term as per boltzmann

    temp_constant = 298.15
    gas_constant = 8.31446261815324
    constant = temp_constant * gas_constant

    y = bottom + ((top - bottom) / (1 + np.exp((m*(cM-urea)/constant))))

    # deltaG can then be calculated as m(cM-urea) - generally calculated at 0M urea therefore m(cM)

    return y


input_path = f'results/lysate_denaturation/sigmoid_fitting/sigmoid_fits.xlsx'
output_folder = f'results/lysate_denaturation/plot_sigmoids/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Choose data type, number of clusters, set colours
cluster_colors = {4: 'royalblue', 2: 'firebrick', 3: 'rebeccapurple', 1: 'darkorange'}
cluster_cmaps = {4: 'Blues', 2: 'Reds', 3: 'Purples', 1: 'Oranges'}
test_sequences = ['CTPACVSFGPK', 'AITIAGVPQSVTECVK', 'YTVQDESHSEWVSCVR', 'DVQIGDIVTVGECRPLSK',  'CNEIISWLDK', 'GPAVGIDLGTTYSCVGVFQHGK', 'ILDKCNEIISWLDK', 'VCNPIITK',]

font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 14 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'


# read in raw data
raw_data = pd.read_excel(f'{input_path}', sheet_name=None)

summary = raw_data['summary'].copy()
summary.drop([col for col in summary.columns.tolist() if 'Unnamed: ' in str(col)], axis=1, inplace=True)
quant_cols = [col for col in summary.columns.tolist() if type(col) != str]

if not test_sequences:
    test_sequences = summary['Sequence'][0:20]

# visualise r_squared - shows high bias toward good fits in cluster as anticipated
for cluster, df in summary.groupby('cluster'):
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.distplot(df['r_squared'], color=cluster_colors[cluster], bins=10)
    plt.xlabel('Goodness of fit (R2)')
    plt.ylabel('Density')
    plt.savefig(f'{output_folder}{cluster}_sigmoid_rsquared.png')
    plt.savefig(f'{output_folder}{cluster}_sigmoid_rsquared.svg')
    plt.ylim(0, 4.5)
    plt.show()

# visualise overall fitted results for each cluster
for cluster, data in summary.groupby('cluster'):
    fig, ax = plt.subplots(figsize=(4, 3))
    for sequence, df in data.iterrows():
        # generate fitted values
        (bottom, top, cM, m, r_squared, protein, sequence) = tuple(df[['bottom_value', 'top_value', 'cM_value', 'm_value', 'r_squared', 'Proteins', 'Sequence']])
        y_vals = denaturant(np.array(quant_cols), top, bottom, cM, m)
        sns.lineplot(quant_cols, y_vals, color=cluster_colors[cluster], linewidth=5, alpha=0.25)
    plt.legend('')
    plt.xlabel('Urea conc (M)')
    plt.ylabel(r'$\Delta Corrected cys ratio')
    plt.ylim(-1, 1)
    plt.savefig(f'{output_folder}sigmoid_{cluster}.png')
    plt.savefig(f'{output_folder}sigmoid_{cluster}.svg')
    plt.show()

# plot test sequences
for sequence, df in summary.set_index('Sequence').iterrows():
    if sequence in test_sequences:
        # generate fitted values
        seq_vals = dict(zip(df.index.tolist(), df.values.tolist()))
        y_vals = denaturant(np.array(quant_cols), seq_vals['top_value'], seq_vals['bottom_value'], seq_vals['cM_value'], seq_vals['m_value'])
        
        # add plot elements
        fig, ax = plt.subplots(figsize=(4, 3))
        # sns.scatterplot(comparison, list(df[comparison_cols]), marker = 'o', label=None, color='lightgrey')
        sns.lineplot(quant_cols, list(df[quant_cols]), label=None, color='lightgrey', linewidth=5, alpha=0.5)
        sns.lineplot(quant_cols, y_vals, color=cluster_colors[seq_vals['cluster']], linewidth=5)
        plt.annotate(f'R2={round(seq_vals["r_squared"], 2)}', (0.8, 0.85), xycoords = 'figure fraction')
        plt.legend('')
        plt.xlabel('Urea conc (M)')
        plt.ylabel(r'$\Delta Corrected cys ratio')
        plt.title(sequence)
        plt.ylim(-1, 1)
        plt.savefig(f'{output_folder}{sequence}.png')
        plt.savefig(f'{output_folder}{sequence}.svg')
        plt.show()
