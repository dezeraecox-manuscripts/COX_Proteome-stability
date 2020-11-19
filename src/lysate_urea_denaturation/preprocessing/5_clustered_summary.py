import os, functools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from GEN_Utils import FileHandling
from loguru import logger

logger.info('Import OK')


input_path = f'results/lysate_denaturation/clustering/cluster_summary.xlsx'
output_folder = f'results/lysate_denaturation/clustering/'

cluster_type, cluster_number = ('pval_loess', 4)


if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# import cluster data
clusters = pd.read_excel(input_path, sheet_name=None)
clustered = clusters[cluster_type].drop([col for col in clusters[cluster_type].columns.tolist() if 'Unnamed: ' in str(col)], axis=1).copy()

# to generate cluster map
# In order to ease figure plotting for paper, map cluster numbers into desired order
cols = [col for col in clustered.columns.tolist() if type(col) != str]
df = pd.melt(clustered.rename(columns={f'member_{cluster_number}': 'cluster', f'score_{cluster_number}': 'score'}), id_vars=['Sequence', 'cluster'], value_vars=cols, var_name='conc', value_name='ratio')
sns.lineplot(data=df.groupby(['cluster', 'conc']).mean().reset_index(), x='conc', y='ratio', hue='cluster', palette='muted')
cluster_map = {0: 2, 1: 4, 2: 1, 3: 3}

cluster_df = clusters[cluster_type][['Sequence', 'Proteins', f'member_{cluster_number}', f'score_{cluster_number}']].rename(columns={f'member_{cluster_number}': 'cluster', f'score_{cluster_number}': 'score'})

if cluster_map:
    clustered[f'member_{cluster_number}'] = clustered[f'member_{cluster_number}'].map(cluster_map)
    cluster_df['cluster'] = cluster_df['cluster'].map(cluster_map)

cluster_counts = cluster_df.groupby(['Proteins', f'cluster']).count().groupby('Proteins').count()

# collect all proteins with at least one peptide assigned to that cluster
proteins_mixed = pd.pivot(cluster_df, index=None, columns='cluster', values='Proteins')
proteins_mixed = pd.DataFrame([list(proteins_mixed[cluster].dropna().unique()) for cluster in proteins_mixed.columns.tolist()], index=proteins_mixed.columns.tolist()).T

# collect all proteins only found in one cluster
unique = cluster_counts[cluster_counts['Sequence'] == 1].index.tolist()
proteins_unique = pd.pivot(cluster_df[cluster_df['Proteins'].isin(unique)].groupby('Proteins').mean()['cluster'].reset_index(), index=None, columns='cluster', values='Proteins')
proteins_unique = pd.DataFrame([list(proteins_unique[cluster].dropna().unique()) for cluster in proteins_unique.columns.tolist()], index=proteins_unique.columns.tolist()).T

# collect proteins found in multiple clusters
proteins_multiple = pd.DataFrame(cluster_counts[cluster_counts['Sequence'] > 1].index.tolist())
proteins_multiple.columns = ['multiple']

# collect proteins according to how many clusters they are found in
cluster_count = pd.pivot(cluster_counts.reset_index(), index=None, columns='Sequence', values='Proteins')
cols = cluster_count.columns.tolist()
cluster_count = pd.DataFrame([list(cluster_count[cluster].dropna().unique()) for cluster in cluster_count.columns.tolist()], index=cluster_count.columns.tolist()).T
cluster_count.columns = [f'{cluster}_clusters' for cluster in cols]

# generate summary df
summary = cluster_df.copy()
summary['mixed'] = summary['cluster'].copy()
summary['unique'] = ['multiple' if protein in proteins_multiple['multiple'].tolist() else cluster for protein, cluster in summary[['Proteins', 'cluster']].values]
summary['count'] = summary['Proteins'].map(dict(zip(cluster_counts.reset_index()['Proteins'], cluster_counts.reset_index()['Sequence'])))

FileHandling.df_to_excel(
    output_path=f'{output_folder}clustered.xlsx',
    sheetnames= ['clustered', 'proteins_mixed', 'proteins_unique', 'proteins_multiple', 'clusters_count', 'summary'],
    data_frames= [clustered, proteins_mixed, proteins_unique, proteins_multiple, cluster_count, summary]
    )
