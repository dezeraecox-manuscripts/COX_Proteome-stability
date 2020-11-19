import os, re, string
import functools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from statsmodels.stats.outliers_influence import summary_table
from scipy.stats import pearsonr, spearmanr

from utilities.statistical_tests import correlation_df

from GEN_Utils import FileHandling
from loguru import logger

logger.info('Import OK')

input_path = 'results/published_stability_correlations/published_stability_datasets.xlsx'
output_folder = f'results/published_stability_correlations/correlation/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# -----------------------------------read in cleaned data-----------------------------------

cleaned_stability = pd.read_excel(input_path, sheet_name=None)
dataset_details = cleaned_stability.pop('index')
dataset_details.drop([col for col in dataset_details.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
dataset_categories = cleaned_stability.pop('dataset_categories')
dataset_categories.drop([col for col in dataset_categories.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
cleaned_stability.update({dataset: df.drop([col for col in df.columns.tolist() if 'Unnamed: ' in col], axis=1) for dataset, df in cleaned_stability.items()})

# -------------------------------calculate dataset statistics-------------------------------
# proportion mapped
proportions = {}
for dataset in dataset_details['dataset_id']:
    df = cleaned_stability[dataset]
    # - From original protein name to UniProt ID
    protein_proportion = len(df['Proteins'].dropna()) / len(df) * 100

    # proportion with fitted stability
    stability_col = [col for col in df.columns.tolist() if 'Tm' in col][0]
    stability_proportion = len(df[stability_col].dropna()) / len(df) * 100

    # drop all non-mapped proteins to only look at the mapped versions
    df = df.dropna(subset=['Proteins'])
    # - To KO, OrthoDB, MIG
    ko_proportion = len(df['KO'].dropna()) / len(df) * 100
    ortho_proportion = len(df['OrthoDB'].dropna()) / len(df) * 100
    mig_proportion = len(df['MIG_ID'].dropna()) / len(df) * 100

    mapped_stability_proportion = len(df[stability_col].dropna()) / len(df) * 100

    proportions[dataset] = [protein_proportion, ko_proportion, ortho_proportion, mig_proportion, stability_proportion, mapped_stability_proportion]
proportions = pd.DataFrame(proportions).T.reset_index()
proportions.columns = ['dataset_id', 'protein_proportion', 'ko_proportion', 'ortho_proportion', 'mig_proportion', 'stability_overall_proportion', 'stability_mapped_proportion']

dataset_details = pd.merge(dataset_details, proportions, on='dataset_id')
dataset_details['index_key'] = dataset_details['dataset_id'].str.rstrip(string.ascii_uppercase)
dataset_details = pd.merge(dataset_details, dataset_categories[['index_key', 'technique', 'species', 'fit_available']], on='index_key', how='inner')

FileHandling.df_to_excel(
    output_path=f'{output_folder}dataset_details.xlsx',
    data_frames=[dataset_details],
    sheetnames=['dataset_details']
)

# ------------------------merge dfs according to each mapping system------------------------
dfs = []
for dataset, df in cleaned_stability.items():
    stability_col = [col for col in df if 'Tm' in col][0]
    df.rename(columns={stability_col:'stability'}, inplace=True)
    df['dataset_id'] = dataset
    dfs.append(df)
compiled = pd.concat(dfs).dropna(subset=['Proteins'])

compiled_ko = pd.pivot_table(compiled, index='KO', columns='dataset_id', values='stability')
compiled_ortho = pd.pivot_table(compiled, index='OrthoDB', columns='dataset_id', values='stability')
compiled_mig = pd.pivot_table(compiled, index='MIG_ID', columns='dataset_id', values='stability')

FileHandling.df_to_excel(
    output_path=f'{output_folder}compiled_dfs.xlsx',
    data_frames=[compiled, dataset_details, compiled_ko, compiled_ortho, compiled_mig],
    sheetnames=['compiled', 'dataset_details', 'compiled_ko', 'compiled_ortho', 'compiled_mig']
)


# ------------------------------calculate pairwise correlation------------------------------

# KO
# Correlation calculation for all published datasets using limit of minimum 5 overlapping proteins
correlation_summary = {}
for mapping_type, mapping_df in zip(['correlation_KO', 'correlation_ortho', 'correlation_mig'],[compiled_ko, compiled_ortho, compiled_mig]):
    correlation_data = {}
    for resource_1 in dataset_details['dataset_id']:
        for resource_2 in dataset_details['dataset_id']:
            correlation_data[f'{resource_1}-{resource_2}'] = correlation_df(mapping_df[f'{resource_1}'], mapping_df[f'{resource_2}'])
    correlation = pd.concat(list(correlation_data.values()), axis=1)
    correlation.columns = list(correlation_data.keys())
    correlation = correlation.T.reset_index().rename(columns={'index':'comparison'})
    correlation[['dataset_1', 'dataset_2']] = correlation['comparison'].str.split('-', expand=True)
    correlation_summary[mapping_type] = correlation
    # Visualise table with significance markers        
    correlation.style.bar(subset=['pearsons_r'], align='mid', color=['#de251f', '#1f78de']).background_gradient(subset=['pearsons_pval'], cmap='Reds_r')


# save compiled dfs
FileHandling.df_to_excel(
    output_path=f'{output_folder}correlations.xlsx',
    data_frames=list(correlation_summary.values()),
    sheetnames=list(correlation_summary.keys())
)

