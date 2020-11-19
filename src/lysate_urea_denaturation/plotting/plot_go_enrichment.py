
# --------------------------------from enrichment barplots--------------------------------

import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 14 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'


def filter_levels(summary_df, level=None, level_min=0.0, level_max=None, proteins_min=2):
    # calculate log2 fold enrichment
    summary_df['log2_fold_enrichment'] = np.log2(summary_df['fold_enrichment'])

    # assign filtering criteria to determine what level (and optionally how many proteins must be in that level) for plotting
    filtered = summary_df.copy()

    if level:
        logger.info(f'filtering to only level {level}')
        filtered = filtered[filtered['level'] == level]
    if proteins_min:
        logger.info(f'filtering min number of proteins {proteins_min}')
        filtered = filtered[filtered['number_in_list'] >= proteins_min]
    if level_min:
        logger.info(f'filtering min level {level_min}')
        filtered = filtered[filtered['level'] >= level_min]
    if level_max:
        logger.info(f'filtering max level {level_max}')
        filtered = filtered[filtered['level'] <= level_max]
    selected = filtered.groupby(['search_type', 'column', 'family'])['level'].min().reset_index()
    selected['selected'] = 'yes'
    filtered = pd.merge(filtered, selected, how='outer', on=['search_type', 'column', 'family', 'level'])
    filtered = filtered[filtered['selected'] == 'yes']
    return filtered

def plot_enrichment(filtered, colour_dict=None, filter_top=5, output_folder=False):
    # generate bar plots

    # Initialize the matplotlib figure
    for (category, column), df in filtered.groupby(['search_type', 'column']):
        source = df.copy()
        source = source.sort_values(['log2_fold_enrichment'], ascending=False)
        if filter_top:
            source = source.iloc[0:filter_top]
        fig, ax = plt.subplots(figsize=(8, source.shape[0]*0.45))
        if colour_dict:
            sns.barplot(x="log2_fold_enrichment", y="term_label", data=source, color=colour_dict[column])
        else:
            sns.barplot(x="log2_fold_enrichment", y="term_label", data=source)
        # Add term_ids
        for position, term_id in enumerate(source['term_id']):
            plt.annotate(term_id, (2, position), color='white')
        plt.xlim(0, 7)
        plt.title(f'GO enrichment for {category}', )
        plt.ylabel(None)
        plt.xlabel("Log2 Fold Enrichment")
        sns.despine(right=True, top=True)
        if output_folder:
            plt.savefig(f'{output_folder}/{column}_{category.replace(" ", "_")}.svg')
            plt.tight_layout()
            plt.savefig(f'{output_folder}/{column}_{category.replace(" ", "_")}.png')
        plt.show()


def plot_compiled_enrichment(filtered, colour_dict=None, filter_top=5, output_folder=False):
    # generate bar plots

    # Initialize the matplotlib figure
    for column, df in filtered.groupby(['column']):
        source = df.copy()
        source = source.sort_values(['log2_fold_enrichment'], ascending=False)
        if filter_top:
            source = source.iloc[0:filter_top]
        fig, ax = plt.subplots(figsize=(8, source.shape[0]*0.45))
        if colour_dict:
            sns.barplot(x="log2_fold_enrichment", y="term_label", data=source, color=colour_dict[column])
        else:
            sns.barplot(x="log2_fold_enrichment", y="term_label", data=source)
        # Add term_ids
        for position, term_id in enumerate(source['term_id']):
            plt.annotate(term_id, (2, position), color='white')
        plt.xlim(0, 7)
        plt.title(f'GO enrichment' )
        plt.ylabel(None)
        plt.xlabel("Log2 Fold Enrichment")
        sns.despine(right=True, top=True)
        if output_folder:
            plt.savefig(f'{output_folder}/{column}.svg')
            plt.tight_layout()
            plt.savefig(f'{output_folder}/{column}.png')
        plt.show()

if __name__ == "__main__":

    input_path = 'results/lysate_denaturation/go_enrichment/cluster_enrichment.xlsx'
    output_folder = f'results/lysate_denaturation/plot_go_enrichment/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # read in summary statistics
    protein_enrichment = pd.read_excel(f'{input_path}', sheet_name=None)

    filtered_dfs = {}
    for clustered_type, df in protein_enrichment.items():

        if not os.path.exists(f'{output_folder}{clustered_type}/'):
            os.makedirs(f'{output_folder}{clustered_type}/')

        filtered_df = filter_levels(df, level=None, level_min=0.0, level_max=None, proteins_min=1)
        filtered_dfs[clustered_type] = filtered_df
        cluster_colors = {1: 'darkorange', 2: 'firebrick', 3: 'rebeccapurple', 4: 'royalblue', '1_clusters': '#9cc5a1', '2_clusters': '#49a078', '3_clusters': '#216869', '4_clusters': '#1f2421', 'multiple': 'grey'}

        plot_enrichment(filtered_df, colour_dict=cluster_colors, filter_top=5, output_folder=f'{output_folder}{clustered_type}/')
        plot_compiled_enrichment(filtered_df, colour_dict=cluster_colors, filter_top=False, output_folder=f'{output_folder}{clustered_type}/')
        


