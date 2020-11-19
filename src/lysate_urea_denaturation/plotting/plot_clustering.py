
import os, functools
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import skfuzzy as fuzz
from GEN_Utils import FileHandling
from loguru import logger
from pykalman import KalmanFilter
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

logger.info("Import ok")

input_folder = f'results/lysate_denaturation/clustering/clustered.xlsx'
output_folder = f'results/lysate_denaturation/plot_clustering/'

cluster_colors = {4: 'royalblue', 2: 'firebrick', 3: 'rebeccapurple', 1: 'darkorange'}
cluster_cmaps = {4: 'Blues', 2: 'Reds', 3: 'Purples', 1: 'Oranges'}

if not os.path.exists(output_folder):
    os.makedirs(output_folder)    

font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 14 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

# Read in clustered data
clustered = pd.read_excel(f'{input_folder}', sheet_name=None)
clustered.update({key: df.drop([col for col in df.columns.tolist() if 'Unnamed: ' in str(col)], axis=1) for key, df in clustered.items()})

clustered_data = pd.merge(clustered['summary'], clustered['clustered'], on=['Proteins', 'Sequence'])
quant_cols = [col for col in clustered_data.columns.tolist() if type(col) != str]

# visualise number of peptides in each cluster
fig, ax = plt.subplots(figsize=(4, 3))
sns.barplot(data=clustered_data.groupby(['cluster']).count().reset_index(), x='cluster', y='Sequence', color='cluster', palette=cluster_colors)
plt.ylabel('Number of Peptides')
plt.xlabel('Cluster')
plt.savefig(f'{output_folder}peptides_per_cluster.png')
plt.savefig(f'{output_folder}peptides_per_cluster.svg')
plt.show()

# plot PCA
fig, ax = plt.subplots(figsize=(4, 3))
sns.scatterplot(data=clustered_data, x='PC1', y='PC2', hue='cluster', palette=cluster_colors, size='score', legend=None)
plt.xlim(-10, 10)
plt.ylim(-10, 10)
plt.savefig(f'{output_folder}PCA.png')
plt.savefig(f'{output_folder}PCA.svg')
plt.show()


# plot clustered peptides
for_plotting = pd.melt(clustered_data, id_vars=['Sequence', 'cluster', 'score'], value_vars=quant_cols, var_name='channel', value_name='ratio')

for cluster_number, df in for_plotting.groupby(f'cluster'):
    color = cluster_cmaps[cluster_number]
    fig, axes = plt.subplots(figsize=(4, 3))
    sns.lineplot(data=df, x='channel', y='ratio', hue='score', palette=color, legend=None)
    plt.ylim(-1.5, 1.5)
    plt.xlabel(f'Urea Concentration (M)')
    plt.ylabel(r'$\Delta$ Corrected cys ratio')
    plt.savefig(f'{output_folder}cluster_{cluster_number}.png')
    plt.savefig(f'{output_folder}cluster_{cluster_number}.svg')
    plt.show()

# Generate simple schematic plot from means

fig, axes = plt.subplots(figsize=(4, 3))
sns.lineplot(data=for_plotting.groupby(['cluster', 'channel']).mean().reset_index(), x='channel', y='ratio', hue='cluster', palette=cluster_colors, legend=None, linewidth=10, alpha=0.75)
plt.ylim(-1, 1)
plt.xlabel(None)
plt.xticks([], [])
plt.ylabel(None)
plt.yticks([], [])
plt.savefig(f'{output_folder}cluster_schematic.png')
plt.savefig(f'{output_folder}cluster_schematic.svg')
plt.show()