"""Process and clean predicted and calculated features for proteins"""
import os, re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from GEN_Utils import FileHandling
from IPython.display import Image
from pymol import cmd
from scipy.stats import normaltest, kruskal, fisher_exact, ttest_ind
from statsmodels.stats.multitest import multipletests
from upsetplot import generate_counts, from_contents
from collections import defaultdict
import functools
import statsmodels.api as sa
import statsmodels.formula.api as sfa
import scikit_posthocs as sp
from itertools import combinations


from loguru import logger

from utilities.database_collection import uniprot_features, asa_calculator
from utilities.database_map_and_filter import create_uniprot_db,  uniprot_summary
from utilities.decorators import ProgressBar
from utilities.statistical_tests import apply_oneway_anova, apply_kruskal, fischers_test, apply_chisquare, apply_fisher_test

logger.info('Import OK')

cluster_path = f'results/lysate_denaturation/clustering/clustered.xlsx'
background_path = f'results/lysate_denaturation/normalised/normalised_summary.xlsx'
predicted_folder = f'results/lysate_denaturation/predicted_features/'
calculated_folder = f'results/lysate_denaturation/uniprot_features/'
complex_path = 'results/lysate_denaturation/protein_complexes/'
output_folder = 'results/lysate_denaturation/protein_features/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


def result_filter(result_type, entries, property_name=False, input_folder=output_folder):
    result = pd.read_table(f'{input_folder}{result_type}.tsv').set_index('#')
    if property_name:
        property_cols = [col for col in result.columns.tolist() if property_name in col]
    else:
        property_cols = result.columns.tolist()
    filtered = result[property_cols].reset_index()
    filtered = filtered[filtered['#'].isin(entries)]
    return filtered

def set_intersection_without_replacement(sets):
    set_list = sets.values()
    subsets = []
    names = []
    for i in range(len(sets), 0, -1):
        for combo in combinations(set_list, i):
            overlap = set.intersection(*combo)
            subsets.append(overlap)
            set_list = [subset.difference(overlap) for subset in set_list]
        for combo in combinations(sets.keys(), i):
            names.append(''.join(combo))

    return subsets, names

# ------------------------------Read in clustered data------------------------------
clustered_data = pd.read_excel(f'{cluster_path}', sheet_name='summary')
clustered_data.drop([col for col in clustered_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

col='unique'

# -----------------Determine overlapping proteins for venn diagram-----------------

sets = {str(cluster): set(df['Proteins']) for cluster, df in clustered_data.groupby('mixed')}
# subsets = [set.intersection(*combo) for i in range(len(sets), 0, -1) for combo in combinations(sets.values(), i)]
names = [''.join(combo) for i in range(len(sets), 0, -1) for combo in combinations(sets.keys(), i)]
subsets, names = set_intersection_without_replacement(sets)
counts = [len(intersects) for intersects in reversed(subsets)]
names = [f'Cluster {name}' for name in reversed(names)]
protein_venn = pd.DataFrame(counts, index=names).reset_index().rename(columns={'index':'name', 0:'count'})


# Determine proteins with multiple peptides across clusters
multiple_peps = clustered_data.groupby('Proteins').count().reset_index()
multiple_peps = multiple_peps[multiple_peps['Sequence'] > 1]['Proteins'].tolist()
# based on distribution, roughly half (251/545) proteins have multiple peptides 
# and of those ~75% of proteins with multiple cysteines in multiple clusters
distribution = clustered_data[clustered_data['Proteins'].isin(multiple_peps)][['Proteins', 'unique']].drop_duplicates().groupby('unique').count()
multiple_peps_counts = clustered_data[clustered_data['Proteins'].isin(multiple_peps)][['Proteins', 'unique']].groupby('Proteins').count().reset_index()
multiple_peps = pd.merge(clustered_data[clustered_data['Proteins'].isin(multiple_peps)][['Proteins', 'unique']].drop_duplicates(), multiple_peps_counts.rename(columns={'unique': 'peptide_count'}), on='Proteins', how='outer')

# check for correlation between number of peptides identified and number of clusters
multiple_pep_clustered = pd.merge(multiple_peps, clustered_data[['Proteins', 'mixed', 'count']], on='Proteins', how='outer')

FileHandling.df_to_excel(output_path=f'{output_folder}protein_venn.xlsx', sheetnames=['protein_venn', 'multiple_peptide_distribution', 'multiple_peptides_clustered'], data_frames=[protein_venn, multiple_peps, multiple_pep_clustered])


# -----------------Processing iFeature predicted protein info-----------------
# Collect proteins of interest into summary dataframe --> using unique here as stats test require independent samples
#  and we care primarily about whether there is a difference between proteins found in one vs multiple clusters
proteins = clustered_data[['Proteins', col]].drop_duplicates()

# amino acid composition: frequency of all natural amino acids --> heatmap
prot_amino_acids = result_filter(result_type='AAC', entries=proteins['Proteins'].unique().tolist(), property_name=False, input_folder=f'{predicted_folder}proteins/')

prot_ctdc = result_filter(result_type='CTDC', entries=proteins['Proteins'].unique().tolist(), input_folder=f'{predicted_folder}proteins/')
column_dict = {
    'polarity.G1': 'Polar', 'polarity.G2': 'Neutral', 'polarity.G3': 'Hydrophobic',
    'hydrophobicity_PRAM900101.G1': 'Low', 'hydrophobicity_PRAM900101.G2': 'Medium', 'hydrophobicity_PRAM900101.G3': 'High',
    'charge.G1': 'charge_positive', 'charge.G2': 'charge_neutral', 'charge.G3': 'charge_negative',
    'secondarystruct.G1': 'Helix', 'secondarystruct.G2': 'Strand', 'secondarystruct.G3': 'Coil',
    'solventaccess.G1': 'Buried', 'solventaccess.G2': 'Exposed', 'solventaccess.G3': 'Intermediate',
    }
prot_ctdc.rename(columns=column_dict, inplace=True)

predicted = pd.merge(prot_amino_acids, prot_ctdc, on='#').rename(columns={'#': 'Proteins'})

predicted_features = pd.merge(proteins, predicted.drop_duplicates().reset_index(drop=True), on=['Proteins'], how='inner')
predicted_features['polar_sum'] = predicted_features[['S', 'T', 'Y', 'N', 'Q']].sum(axis=1) # mimics figure from meltome dataset
predicted_features['hydrophobic_sum'] = predicted_features[['A', 'V', 'I', 'L', 'M', 'F', 'W']].sum(axis=1) # mimics figure from meltome dataset

features = list(column_dict.values())
features.append('polar_sum')
features.append('hydrophobic_sum')


# Due to non-parametric nature of the data, perform Kruskal–Wallis test with Dunn’s post hoc test, multiple testing correction according to Benjamini–Hochberg
predicted_stats = {}
predicted_stats['compiled_data'] = predicted_features
for feature in features:
    logger.info(f'Processing {feature}')
    test_val, pval, posthoc = apply_kruskal(predicted_features, val_col=feature, group_col=col, p_adjust='bonferroni', sig_cutoff=0.05)
    predicted_stats[f'{feature}_kruskal'] = pd.DataFrame([test_val, pval], index=['test_statistic', 'pvalue'])
    predicted_stats[f'{feature}_posthoc'] = posthoc
FileHandling.df_to_excel(output_path=f'{output_folder}predicted_features_kruskal.xlsx', sheetnames=list(predicted_stats.keys()), data_frames=list(predicted_stats.values()))


# Perform one-way anova with Bonferonni post-hoc tests
predicted_anovas = {}
predicted_anovas['compiled_data'] = predicted_features
for feature in column_dict.values():
    logger.info(f'Processing {feature}')
    predicted_features[[col, feature]]
    anova, pairwise = apply_oneway_anova(predicted_features[[feature, col]], xcol=feature, group_col=col)
    predicted_anovas[f'{feature}_unique_anova'] = anova 
    predicted_anovas[f'{feature}_unique_pairwise'] = pairwise
FileHandling.df_to_excel(output_path=f'{output_folder}predicted_features_anova.xlsx', sheetnames=list(predicted_anovas.keys()), data_frames=list(predicted_anovas.values()))

# test normality
# statistic, pval = normaltest(a, axis=0, nan_policy='omit')

# ----------------------------Processing IUPred Disorder prediction----------------------------

disorder_prediction = pd.read_excel(f'{predicted_folder}disorder_prediction.xlsx')
disorder_prediction.drop([col for col in disorder_prediction.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

disorder_proportion = (disorder_prediction.groupby(['Proteins', 'disordered']).count() / disorder_prediction.groupby(['Proteins']).count() * 100)['disorder_probability'].reset_index().rename(columns={'disorder_probability': 'proportion'})
disorder_proportion = pd.merge(clustered_data[['Proteins', col]].copy().drop_duplicates(), disorder_proportion, on='Proteins', how='left')

# Due to non-parametric nature of the data, perform Kruskal–Wallis test with Dunn’s post hoc test, multiple testing correction according to Benjamini–Hochberg
disorder_stats = {}
disorder_stats['compiled_data'] = disorder_proportion
disorder_stats['disorder_probability'] = disorder_prediction
test_val, pval, posthoc = apply_kruskal(disorder_proportion, val_col='proportion', group_col=col, p_adjust='bonferroni', sig_cutoff=0.05)
disorder_stats[f'kruskal'] = pd.DataFrame([test_val, pval], index=['test_statistic', 'pvalue'])
disorder_stats[f'posthoc'] = posthoc
FileHandling.df_to_excel(output_path=f'{output_folder}disorder_prediction_kruskal.xlsx', sheetnames=list(disorder_stats.keys()), data_frames=list(disorder_stats.values()))


# --------------------------Process cys residue features from UniProt--------------------------
uniprot_stats = {}
peptide_features = pd.read_excel(f'{calculated_folder}uniprot_feature_summary.xlsx', sheet_name='peptides')

# Molecular weight
# cleanup to remove non-cys peptides, and collect unique peptides
peptide_features = peptide_features[~peptide_features['cys_in_peptide'].isnull()]
peptide_features = peptide_features.drop([col for col in peptide_features.columns.tolist() if 'Unnamed: ' in col], axis=1).rename(columns={'pep_sequence': 'Sequence'})
peptide_features = pd.merge(peptide_features, clustered_data, on=['Sequence', 'Proteins'], how='right')
peptide_features['length'] = [len(seq) for seq in peptide_features['full_sequence']]
peptide_features['cluster_num'] = ['single' if clusters == 1 else 'multiple' for clusters in peptide_features['count']]
peptide_features = peptide_features[['Proteins', 'MW', 'length', 'cluster_num', col]].drop_duplicates()

# Perform one-way anova with Bonferonni post-hoc tests
mw_anovas, mw_anova_posthocs = apply_oneway_anova(peptide_features, xcol='MW', group_col=col)
# compare single-cluster to multiple cluster proteins
mw_multiple_ttest = ttest_ind(
    peptide_features[peptide_features['cluster_num'] == 'single']['MW'].tolist(),
    peptide_features[peptide_features['cluster_num'] == 'multiple']['MW'].tolist())

uniprot_stats['mw_anovas'] = mw_anovas
uniprot_stats['mw_anova_ph'] = mw_anova_posthocs
uniprot_stats['mw_multiple_ttest'] = pd.DataFrame(mw_multiple_ttest, index=['ttest_val', 'p-val'])

# Length
# Perform one-way anova with Bonferonni post-hoc tests
len_anovas, len_anova_posthocs = apply_oneway_anova(peptide_features, xcol='length', group_col=col)
# compare single-cluster to multiple cluster proteins
len_multiple_ttest = ttest_ind(
    peptide_features[peptide_features['cluster_num'] == 'single']['length'].tolist(),
    peptide_features[peptide_features['cluster_num'] == 'multiple']['length'].tolist())

uniprot_stats['len_anovas'] = len_anovas
uniprot_stats['len_anova_ph'] = len_anova_posthocs
uniprot_stats['len_multiple_anovas'] = pd.DataFrame(len_multiple_ttest, index=['ttest_val', 'p-val'])

# Assess representation of domains in uni vs multidomain proteins
pfam_domains = pd.read_csv(f'{calculated_folder}whole_protein_pfam_domains.csv')
pfam_domains.drop([col for col in pfam_domains.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
proteins = clustered_data[['Proteins', col]].drop_duplicates()
pfam_domains = pd.merge(proteins, pfam_domains, how='left', on=['Proteins'])
pfam_domains = pfam_domains.groupby('Proteins').count()['pfam_id'].reset_index() # calculate number of domains per protein
peptide_features = pd.merge(peptide_features, pfam_domains, how='left', on=['Proteins'])

# Perform one-way anova with Bonferonni post-hoc tests
pfam_anovas, pfam_anova_posthocs = apply_oneway_anova(peptide_features, xcol='pfam_id', group_col=col)
# compare single-cluster to multiple cluster proteins
pfam_multiple_ttest = ttest_ind(
    peptide_features[peptide_features['cluster_num'] == 'single']['pfam_id'].tolist(),
    peptide_features[peptide_features['cluster_num'] == 'multiple']['pfam_id'].tolist())

uniprot_stats['pfam_anovas'] = pfam_anovas
uniprot_stats['pfam_anova_ph'] = pfam_anova_posthocs
uniprot_stats['pfam_multiple_anovas'] = pd.DataFrame(pfam_multiple_ttest, index=['ttest_val', 'p-val'])

uniprot_stats['len_mw_summary'] = peptide_features


# read in uniprot cys peptide features
cys_peptides = pd.read_excel(f'{calculated_folder}cys_peptide_features.xlsx')
cys_peptides.drop([col for col in cys_peptides.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# drop non-mapped rows
calculated_features = cys_peptides.copy().dropna(subset=['pdb_aa']).rename(columns={'proteins': 'Proteins', 'pep_sequence': 'Sequence'})
# add cluster info
calculated_features = pd.merge(calculated_features, clustered_data, on=['Proteins', 'Sequence'], how='right')
calculated_features['mapped'] = [1 if cys_pos == 'C' else 0 for cys_pos in calculated_features['pdb_aa']]

# Number of sequences with mapped structures - chi-squared test
proportions = calculated_features[['Proteins', 'mapped', 'unique', 'count']].copy().drop_duplicates()

proportion = (proportions.groupby([col, 'mapped']).count())['Proteins'].reset_index().rename(columns={'Proteins': 'proportion', col: 'cluster'})
proportion = pd.pivot_table(proportion, values='proportion', index='mapped', columns='cluster')
test_val, pval, dof, contingency, chi_test = apply_chisquare(proportion, p_adjust='bonferroni', sig_cutoff=0.05)
chi_test['group'] = col

uniprot_stats['calculated_features'] = calculated_features
uniprot_stats['mapped_chi_proportions'] = proportion
uniprot_stats['mapped_chi_enrichment'] = chi_test


FileHandling.df_to_excel(
    output_path=f'{output_folder}calculated_protein_features_summary.xlsx',
    sheetnames=list(uniprot_stats.keys()),
    data_frames=list(uniprot_stats.values())
)


# -------------------------------------Chaperone Proteins-------------------------------------
chaperone_stats = {}
chaperone_data = pd.read_excel(f'{calculated_folder}chaperone_proteins.xlsx', sheet_name=None)

# Read in background data to generate list of protein sequences to map
raw_data = pd.read_excel(background_path, sheet_name='raw')
raw_data = raw_data.drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in str(col)], axis=1)
background_genes = raw_data['Proteins'].copy().unique()

clustered = clustered_data.copy()
clustered['chaperone'] = [1 if protein in chaperone_data['chaperone_proteins']['Proteins'].tolist() else 0 for protein in clustered['Proteins']]

# investigate chaperones in each cluster
chaperones = clustered[clustered['chaperone'] == 1].copy()
chaperones = pd.merge(chaperones, chaperone_data['chaperone_details'][['Entry', 'Entry name', 'Protein names', 'Gene names']], left_on='Proteins', right_on='Entry', how='left')
chaperone_membership = {cluster: df['Proteins'].unique().tolist() for cluster, df in chaperones.groupby('mixed')}
counts = from_contents(chaperone_membership)
cluster_counts = counts.reset_index().set_index('id')
cluster_counts['sum'] = cluster_counts.sum(axis=1)
cluster_counts.sort_values('sum')
chaperone_stats['cluster_counts'] = cluster_counts
chaperone_stats['clustered_chaperones'] = clustered
chaperone_stats['chaperone_details'] = chaperones


# Determine overlapping chaperones for venn diagram
sets = {str(cluster): set(df['Proteins']) for cluster, df in chaperones.groupby('mixed')}
# subsets = [set.intersection(*combo) for i in range(len(sets), 0, -1) for combo in combinations(sets.values(), i)]
names = [''.join(combo) for i in range(len(sets), 0, -1) for combo in combinations(sets.keys(), i)]
subsets, names = set_intersection_without_replacement(sets)
counts = [len(intersects) for intersects in reversed(subsets)]
names = [name for name in reversed(names)]
chaperone_stats['overlap'] = pd.DataFrame(counts, index=[f'Cluster {name}' for name in names]).reset_index().rename(columns={'index':'name', 0:'count'})


# save summaries to excel
FileHandling.df_to_excel(output_path=f'{output_folder}chaperone_enrichment_summary.xlsx', data_frames=list(chaperone_stats.values()), sheetnames=list(chaperone_stats.keys()))