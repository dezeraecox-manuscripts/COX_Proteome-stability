"""Test the enrichment of the entire dataset, or specific clusters against gene ontologies associated with complexes"""

import re
import os
import pandas as pd
import numpy as np

from collections import defaultdict
from scipy import stats
from utilities.database_map_and_filter import ortholog_map, uniprot_go_genes
from utilities.statistical_tests import fischers_test, apply_oneway_anova

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

cluster_path = f'results/lysate_denaturation/clustering/clustered.xlsx'
background_path = f'results/lysate_denaturation/normalised/normalised_summary.xlsx'
ontology_path = 'results/lysate_denaturation/gene_ontology_datasets/go_term_summary.xlsx'
size_path = 'results/lysate_denaturation/gene_ontology_datasets/size_summary.xlsx'
output_folder = 'results/lysate_denaturation/go_enrichment/'
resource_folder = 'resources/bioinformatics_databases/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


# -----------------Read in standard components-----------------
# clustered data
clustered_data = pd.read_excel(f'{cluster_path}', sheet_name='summary')
clustered_data.drop([col for col in clustered_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# raw data as background
background_genes = pd.read_excel(f'{background_path}', sheet_name='raw')['Proteins'].unique()

# read in ontologies
ontology_genes = pd.read_excel(f'{ontology_path}')
ontology_genes.drop([col for col in ontology_genes.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
ontology_genes.drop('homologue_id', inplace=True, axis=1)


# -----------Test ontologies against each filter type with Fischers-----------
filter_types = ['mixed', 'unique', 'count']
results = []
for col in ontology_genes.columns.tolist():
    go_genes = ontology_genes[col].dropna().tolist()
    contingencies = {}
    test_results = {}
    for filter_name in filter_types:
        for group, df in clustered_data.groupby(filter_name):
            hit_genes = df['Proteins'].dropna()
            oddsratio, pval, test_result = fischers_test(go_genes, hit_genes, background_genes, output_folder=None)
            test_results[f'{filter_name}_{group}'] = [oddsratio, pval]
            contingencies[f'{filter_name}_{group}'] = test_result
    test_results = pd.DataFrame(test_results, index=['oddsratio', 'pval']).T
    test_results['ontology'] = col
    results.append(test_results)   

    FileHandling.df_to_excel(output_path=f'{output_folder}{col}_fisher.xlsx', data_frames=[test_results]+list(contingencies.values()), sheetnames=['test_results']+list(contingencies.keys()))

go_summary_results = pd.concat(results)
FileHandling.df_to_excel(output_path=f'{output_folder}clustered_goterm_enrichment.xlsx', data_frames=[go_summary_results], sheetnames=['test_results'])

# ----------------Perform Panther enrichment test----------------
from utilities.statistical_tests import apply_enrichment
# read in clustered data, and background data - this uses all the proteins that contained identified cys peptides
combinations = pd.read_excel(f'{cluster_path}', sheet_name=['proteins_mixed', 'proteins_unique', 'proteins_multiple'])
ref_list = pd.read_excel(f'{background_path}', sheet_name='raw')['Proteins'].unique()

# perform enrichment test
protein_enrichment = {}
for combination, df in combinations.items():
    df.drop([col for col in df if 'Unnamed: ' in str(col)], axis=1, inplace=True)
    protein_enrichment[combination] = apply_enrichment(df, searches=None, obo_path=f'{resource_folder}PANTHERGOslim.obo', organism='10090', refOrganism='10090', enrichmentTestType='FISHER', correction='FDR', min_proteins=2, reference=ref_list)

# Save all to excel
FileHandling.df_to_excel(
    output_path=f'{output_folder}cluster_enrichment.xlsx',
    sheetnames= list(protein_enrichment.keys()),
    data_frames= list(protein_enrichment.values())
    )


