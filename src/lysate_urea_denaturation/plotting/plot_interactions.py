import os, re
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

input_folder = 'results/lysate_denaturation/protein_interactions/'
output_folder = 'results/lysate_denaturation/plot_interactions/'

cluster_colors = { '1': 'darkorange', '2': 'firebrick', '3': 'rebeccapurple', '4': 'royalblue', 1: 'darkorange', 2: 'firebrick', 3: 'rebeccapurple', 4: 'royalblue', 'multiple': 'grey', 'clusters': 'black'}

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 14 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

def color_significant(val):
    """
    Takes a scalar and returns a string with
    the css property `'color: red'` for values
    less than threshold, black otherwise.
    """
    color = '#ba2d2d' if float(val) < 0.05 else None
    return 'background: %s' % color

def color_index(sequence):
    cluster = sequence.split('_')[-1]
    color = cluster_colors[cluster]
    return f'background: {color}'


# Test better boxplot for plotting
def better_boxplot(dataframe, col_to_plot, cluster_col='cluster', cluster_colors=None, order=None, normed=False, title=False, output_folder=output_folder):

    df = dataframe.copy()

    if normed:
        df[col_to_plot] = (df[col_to_plot] - df[col_to_plot].mean()) / df[col_to_plot].std()

    pop_mean = df[col_to_plot].mean()

    error_df = df.groupby(cluster_col).mean().reset_index()
    error_df['error'] = df.groupby(cluster_col).std().reset_index()[col_to_plot]
    if order:
        error_df = error_df.set_index([cluster_col]).loc[order]
    error_df['x_pos'] = np.arange(0, len(error_df))


    fig, ax = plt.subplots(figsize=(4, 3))
    sns.stripplot(x=df[cluster_col], y=df[col_to_plot], palette=cluster_colors, color=df[cluster_col], alpha=0.5, order=order)
    # sns.boxplot(x=cluster_col, y=col_to_plot, data=df.groupby(cluster_col).mean().reset_index(), palette=cluster_colors, color=df[cluster_col], linewidth=2, ax=ax, order=order)
    plt.errorbar(error_df['x_pos'], error_df[col_to_plot], yerr=error_df['error'], ecolor='black', elinewidth=2.0, capsize=10, barsabove=True, capthick=1.0, marker='_', linewidth=0, markersize=20, markeredgecolor='black')
    ax.axhline(pop_mean, linewidth=1, color='grey', linestyle='--')

    plt.xlabel(cluster_col)
    plt.ylabel(col_to_plot)
    if title:
        plt.title(title)
        plt.savefig(f'{output_folder}{title}{col_to_plot}.png')
        plt.savefig(f'{output_folder}{title}{col_to_plot}.svg')
    plt.show()


# ---------------------------Read in precalculated datasets---------------------------
interactions  = pd.read_excel(f'{input_folder}interaction_summary.xlsx', sheet_name=None)
interactions.update({key: value.drop([col for col in value.columns.tolist() if 'Unnamed: ' in col], axis=1) for key, value in interactions.items()})

interaction_summary = interactions['interaction_summary'].copy()
interaction_summary[['group', 'cluster']] = interaction_summary['cluster_filter_type'].str.split('_', expand=True)

firstshell_interactors = interactions['firstshell_interactors'].copy()

manual_node_degree = interactions['manual_node_degrees'].copy()

# -------------------------------visualise interactions-------------------------------
for intersection_type in ['inside_v_outside', 'inside_cluster', 'outside_cluster', 'average_node_degree', 'node_degree_diff']:
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.barplot(data=interaction_summary, x='group', y=intersection_type, hue='cluster', palette=cluster_colors, hue_order=['1', '2', '3', '4', 'multiple'])
    plt.legend(bbox_to_anchor=[1.0, 1.0])
    plt.title(intersection_type)
    plt.savefig(f'{output_folder}{intersection_type}_interactions.png')
    plt.savefig(f'{output_folder}{intersection_type}_interactions.svg')


for group, data in firstshell_interactors.groupby('group'):
    fig, ax = plt.subplots(figsize=(4, 3))
    sns.violinplot(x='cluster', y='interaction_count', data=data, palette=cluster_colors, order=[1, 2, 3, 4, 'multiple'])
    plt.setp(ax.collections, alpha=.3)
    sns.swarmplot(x='cluster', y='interaction_count', data=data, palette=cluster_colors, alpha=0.5, order=[1, 2, 3, 4, 'multiple'])
    plt.ylabel('External interactors')
    plt.xlabel('Cluster')
    plt.title(group)
    plt.savefig(f'{output_folder}{group}_external_interactors.png')
    plt.savefig(f'{output_folder}{group}_external_interactors.svg')
    plt.show()

# plot better bxplot for unique # interactors
better_boxplot(firstshell_interactors[firstshell_interactors['group'] == 'unique'], col_to_plot='interaction_count', cluster_col='cluster', cluster_colors=cluster_colors, order=[1, 2, 3, 4, 'multiple'], normed=False, title='num_interactors', output_folder=output_folder)

# plot better barplot for manually calulcated average node degree (where score > 0.7) 
better_boxplot(manual_node_degree[manual_node_degree['group'] == 'unique'], col_to_plot='interaction_degree', cluster_col='cluster', cluster_colors=cluster_colors, order=[1, 2, 3, 4, 'multiple'], normed=False, title='Node interaction degree', output_folder=output_folder)