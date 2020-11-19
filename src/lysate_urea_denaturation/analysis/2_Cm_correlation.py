import os, functools
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.preprocessing import minmax_scale

from utilities.statistical_tests import stability_summary_calculator, correlation_df
from utilities.database_map_and_filter import create_uniprot_xref

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

resource_folder = 'resources/bioinformatics_databases/'
published_dataset = f'resources/published_datasets_compiled.xlsx'
measured_dataset = f'results/lysate_denaturation/clustering/clustered.xlsx'

output_folder = 'results/lysate_denaturation/cM_correlation/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


# Generate standard databases for human, mouse and ortholog proteins

mouse_ortho = create_uniprot_xref(input_path=resource_folder, tax_id='10090', gene_ids=[], id_type='OrthoDB')
ortho = dict(zip(mouse_ortho['UniProtKB-AC'], mouse_ortho['ID']))

mouse_keggo = create_uniprot_xref(input_path=resource_folder, tax_id='10090', gene_ids=[], id_type='KO')
# note some non-unique keys here...
keggo = dict(zip(mouse_keggo['UniProtKB-AC'], mouse_keggo['ID']))

homology_db = pd.read_table(f'{resource_folder}HOM_MouseHumanSequence.txt')
homology_db.dropna(subset=['SWISS_PROT IDs', 'HomoloGene ID'], inplace=True)
sp_to_homoid = dict(zip(homology_db['SWISS_PROT IDs'], homology_db['HomoloGene ID']))

    
# -----------------------------Prepare measured datasets-----------------------------

measured = pd.read_excel(measured_dataset, sheet_name='sigmoid_summary')
measured = measured.drop([col for col in measured.columns.tolist() if 'Unnamed: ' in str(col)], axis=1)
sample_cols = [col for col in measured.columns.tolist() if type(col) == float]

# Generate summary cM at peptide and protein levels
measured = measured[['Sequence', 'Proteins', 'cM_value', 'cluster', 'filtered']]
measured.columns = ['pep_sequence', 'Proteins', 'peptide_stability', 'cluster', 'category']
# add homology information from orthoDB, KeggO, MIG
measured['KO'] = measured['Proteins'].map(keggo)
measured['OrthoDB'] = measured['Proteins'].map(ortho)
measured['MIG_ID'] = measured['Proteins'].map(sp_to_homoid)
measured['resource'] = 'measured'

# calculate summary for each protein in each cluster with no filtering
measured['protein_stability'] = np.nan
cluster_data = {}
for cluster, df in measured.groupby('cluster'):
    df['resource'] = f'measured_{cluster}'
    cluster_data[f'{cluster}'] = stability_summary_calculator(df, id_col='KO')
cluster_data['all'] = stability_summary_calculator(measured, id_col='KO')

# calculate summary with filtering
for cluster, df in measured[measured['category'] == 1].groupby('cluster'):
    df['resource'] = f'filtered_{cluster}'
    cluster_data[f'{cluster}_filtered'] = stability_summary_calculator(df, id_col='KO')
filtered = measured[measured['category'] == 1].copy()
filtered['resource'] = 'filtered'
cluster_data['all_filtered'] = stability_summary_calculator(filtered, id_col='KO')

measured = pd.concat(list(cluster_data.values()))

measured_summary = [df[['KO', 'mean_stability', 'median_stability']].rename(columns={'mean_stability': f'{resource}_mean_stability', 'median_stability': f'{resource}_median_stability'}).drop_duplicates() for resource, df in measured.groupby('resource')]

measured_summary = functools.reduce(lambda left, right: pd.merge(left, right, on=['KO'], how='outer'), measured_summary)

FileHandling.df_to_excel(
    output_path=f'{output_folder}measured_summary.xlsx',
    data_frames=[measured_summary, measured],
    sheetnames=['measured_summary', 'measured']
)

#---------------------- Prepare published datasets----------------------

stability_summary = pd.read_excel(f'{published_dataset}', sheet_name='compiled_ko')
stability_summary.drop([col for col in stability_summary.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

resources = stability_summary.set_index('KO').columns.tolist()


# --------------------------Calculate correlations--------------------------

# generate scaled dataframes
measured_scaled = measured.copy()
# measured_scaled['mean_stability'] = minmax_scale(measured_scaled['mean_stability'])
# measured_scaled['median_stability'] = minmax_scale(measured_scaled['median_stability'])
# measured_scaled['peptide_stability'] = minmax_scale(measured_scaled['peptide_stability'])
stability_scaled = pd.DataFrame(minmax_scale(stability_summary.set_index('KO')), index=stability_summary['KO'], columns=stability_summary.set_index('KO').columns)


# Correlation of measured data against each resource, after min-max scaling
measured_comparison = pd.merge(measured_scaled[measured_scaled['resource'] == 'measured'], stability_scaled, on='KO')
measured_correlation_data = {}
for resource in resources:
    measured_correlation_data[f'{resource}-all'] = correlation_df(measured_comparison['peptide_stability'], measured_comparison[f'{resource}'])
    for cluster, df in measured_comparison.groupby('cluster'):
        measured_correlation_data[f'{resource}-{cluster}'] = correlation_df(df['peptide_stability'], df[f'{resource}'])
measured_cluster_correlation = pd.concat(list(measured_correlation_data.values()), axis=1)
measured_cluster_correlation.columns = list(measured_correlation_data.keys())
measured_cluster_correlation = measured_cluster_correlation.T

# Visualise table with significance markers  
measured_cluster_correlation.style.bar(subset=['pearsons_r'], align='mid', color=['#de251f', '#1f78de']).background_gradient(subset=['pearsons_pval'], cmap='Reds_r')
# --> Shows no significant correlation with any previous datasets

# Correlation of filtered data (those that fit sigmoid criteria) against each resource, after min-max scaling
filtered_comparison = pd.merge(measured_scaled[measured_scaled['resource'] == 'filtered'], stability_scaled, on='KO')
filtered_correlation_data = {}
for resource in resources:
    filtered_correlation_data[f'{resource}-all'] = correlation_df(filtered_comparison['peptide_stability'], filtered_comparison[f'{resource}'])
    for cluster, df in filtered_comparison.groupby('cluster'):
        filtered_correlation_data[f'{resource}-{cluster}'] = correlation_df(df['peptide_stability'], df[f'{resource}'])
filtered_cluster_correlation = pd.concat(list(filtered_correlation_data.values()), axis=1)
filtered_cluster_correlation.columns = list(filtered_correlation_data.keys())
filtered_cluster_correlation = filtered_cluster_correlation.T

# Visualise table with significance markers  
filtered_cluster_correlation.style.bar(subset=['pearsons_r'], align='mid', color=['#de251f', '#1f78de']).background_gradient(subset=['pearsons_pval'], cmap='Reds_r')
# --> Shows no significant correlation with any previous datasets

# save output to excel
FileHandling.df_to_excel(output_path=f'{output_folder}correlations.xlsx', sheetnames=['measured_scaled', 'stability_scaled', 'measured_comparison', 'measured_cluster_correlation', 'filtered_comparison', 'filtered_cluster_correlation'], data_frames=[measured_scaled, stability_scaled, measured_comparison, measured_cluster_correlation.reset_index(), filtered_comparison, filtered_cluster_correlation.reset_index()])

