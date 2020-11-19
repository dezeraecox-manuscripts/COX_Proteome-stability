import re
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import statsmodels.api as sa
import statsmodels.formula.api as sfa
import scikit_posthocs as sp
import networkx as nx


from loguru import logger
from GEN_Utils import FileHandling

from utilities.database_collection import network_interactions, all_interactions, interaction_enrichment

logger.info('Import OK')

input_path = f'results/lysate_denaturation/clustering/clustered.xlsx'
output_folder = 'results/lysate_denaturation/protein_interactions/'

confidence_threshold = 0.7

if not os.path.exists(output_folder):
    os.makedirs(output_folder)



# ------------------------------Read in clustered data------------------------------
# Read in standard components - hits & background
proteins = pd.read_excel(f'{input_path}', sheet_name='summary')
proteins = proteins.drop([col for col in proteins.columns.tolist() if 'Unnamed: ' in col], axis=1)[['Proteins', 'mixed', 'unique', 'count']]
proteins = pd.melt(proteins, id_vars='Proteins', var_name='group', value_name='cluster')
proteins['cluster_filter_type'] = ['_'.join([var, str(val)]) for var, val in proteins[['group', 'cluster']].values]

cluster_summary = proteins.groupby('cluster_filter_type').count()['Proteins'].reset_index()

# Test 1: Get intra-cluster interactions (i.e. interactions within a cluster)
intra_cluster_interactions = {}
for cluster_type, df in proteins.groupby('cluster_filter_type'):
    gene_ids = df['Proteins'].unique()
    intra_cluster_interactions[cluster_type] = network_interactions(gene_ids, tax_id=10090, id_type='uniprot')

# calculate number of interactions for which evidence is > 0.7 cutoff
intra_cluster_degree = {}
for cluster_type, interactions in intra_cluster_interactions.items():
    filtered_ints = interactions[interactions['score'].astype(float) > confidence_threshold]
    intra_cluster_degree[cluster_type] = len(filtered_ints)

cluster_summary['number_within_cluster'] = cluster_summary['cluster_filter_type'].map(intra_cluster_degree)
cluster_summary['normalised_within_cluster'] = cluster_summary['number_within_cluster'] / cluster_summary['Proteins']

# Test 2: Get intra-cluster interactions within whole interaction dataset vs inter-cluster interactions
gene_ids = proteins['Proteins'].unique()
interactions = network_interactions(gene_ids, tax_id=10090, id_type='uniprot')
interactions = interactions[interactions['score'].astype(float) > confidence_threshold] # less than half remain!

# calculate number of interactions for which evidence is > 0.7 cutoff
inter_vs_intra = {}
for cluster_type, df in proteins.groupby('cluster_filter_type'):
    gene_ids = df['Proteins'].unique()
    cluster_ints = interactions.copy()
    cluster_ints['int_A'] = [1 if protein in gene_ids else 0 for protein in cluster_ints['originalId_A']]
    cluster_ints['int_B'] = [1 if protein in gene_ids else 0 for protein in cluster_ints['originalId_B']]
    cluster_ints['int_type'] = cluster_ints['int_A'] + cluster_ints['int_B']
    inter_vs_intra[cluster_type] = cluster_ints['int_type'].value_counts()
inter_vs_intra = pd.DataFrame(inter_vs_intra).T.reset_index()
inter_vs_intra.columns = ['cluster_filter_type', 'not_in_cluster', 'outside_cluster', 'inside_cluster']

cluster_summary = pd.merge(cluster_summary, inter_vs_intra, on='cluster_filter_type')
cluster_summary['inside_v_outside'] = cluster_summary['inside_cluster'] / cluster_summary['outside_cluster']


# Test 3: total interactions for all genes (up to 1000 per gene)
gene_ids = proteins['Proteins'].unique()
interaction_results = all_interactions(genes=gene_ids, tax_id=10090, max_partners=1000, id_type='uniprot', confidence_threshold=0.7)

interaction_count = {}
for protein in gene_ids:
    interaction_genes = interaction_results[(interaction_results['uniprot_A'] == protein) | (interaction_results['uniprot_B'] == protein)][['uniprot_A', 'uniprot_B']].dropna().drop_duplicates().values
    interaction_genes = set([item for sublist in interaction_genes for item in sublist if item != protein])
    interaction_count[protein] = len(interaction_genes)
interaction_count = pd.DataFrame(interaction_count.values(), index=interaction_count.keys()).reset_index()
interaction_count.columns = ['Proteins', 'interaction_count']

proteins = pd.merge(proteins, interaction_count, on='Proteins', how='left')

# Test 4: interaction enrichment provided by STRING
enrichment = {}
for cluster_type, df in proteins.groupby('cluster_filter_type'):
    gene_ids = df['Proteins'].unique()
    enrichment[cluster_type] = interaction_enrichment(genes=gene_ids, tax_id=10090, id_type='uniprot')
enrichment['all_proteins'] = interaction_enrichment(genes=proteins['Proteins'].unique(), tax_id=10090, id_type='uniprot')
enrichment_summary = pd.concat(list(enrichment.values()))
enrichment_summary['cluster_filter_type'] = list(enrichment.keys())

summary = enrichment_summary.reset_index()[['cluster_filter_type', 'number_of_nodes',  'number_of_edges',  'average_node_degree',  'local_clustering_coefficient',  'expected_number_of_edges',  'p_value']].copy()

cluster_summary = pd.merge(cluster_summary, enrichment_summary, on='cluster_filter_type')

cluster_summary['node_degree_diff'] = (cluster_summary['number_of_edges'].astype(float) / cluster_summary['number_of_nodes'].astype(float)) / (cluster_summary['expected_number_of_edges'].astype(float) / cluster_summary['number_of_nodes'].astype(float))

# Test 5: calculate average node interaction manually for each unique cluster
node_degrees = []
for cluster_type, df in proteins.groupby('cluster_filter_type'):
    gene_ids = df['Proteins'].unique()
    interactions = network_interactions(gene_ids, tax_id=10090, id_type='uniprot')
    interactions = interactions[interactions['score'].astype(float) > confidence_threshold]
    # build a network
    interaction_network = nx.from_pandas_edgelist(interactions, 'originalId_A', 'originalId_B')
    interaction_degrees = {protein: val for (protein, val) in interaction_network.degree()}
    df['interaction_degree'] = df['Proteins'].map(interaction_degrees)
    node_degrees.append(df)
node_degrees = pd.concat(node_degrees)
node_degrees = node_degrees.drop_duplicates()

FileHandling.df_to_excel(
    output_path=f'{output_folder}interaction_summary.xlsx',
    data_frames=[proteins, cluster_summary, node_degrees],
    sheetnames=['firstshell_interactors', 'interaction_summary', 'manual_node_degrees'])

# ----------------------compelete statistical test for # external interactors----------------------

for group, df in proteins.groupby('group'):
    logger.info(group)
    lm = sfa.ols('interaction_count ~ C(cluster)', data=df).fit()
    anova = sa.stats.anova_lm(lm)
    logger.info(anova)
    ph_test = sp.posthoc_ttest(df, val_col='interaction_count', group_col='cluster', p_adjust='bonf')
    logger.info(ph_test)

# compare single vs multiple proteins
interaction_count_unique = proteins[proteins['group'] == 'unique'].copy()
interaction_count_unique['cluster_num'] = ['single' if clusters == 1 else 'multiple' for clusters in interaction_count_unique['cluster']]
interaction_count_unique = interaction_count_unique[['Proteins', 'interaction_count', 'cluster_num']].drop_duplicates()


lm = sfa.ols('interaction_count ~ C(cluster_num)', data=interaction_count_unique).fit()
anova = sa.stats.anova_lm(lm)
logger.info(anova)
ph_test = sp.posthoc_ttest(interaction_count_unique, val_col='interaction_count', group_col='cluster_num', p_adjust='bonf')
logger.info(ph_test)