"""Process and clean predicted and calculated features for peptides"""
import os, re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from GEN_Utils import FileHandling
from IPython.display import Image
from pymol import cmd
from scipy.stats import normaltest, kruskal, fisher_exact
from statsmodels.stats.multitest import multipletests
from collections import defaultdict
import functools

from loguru import logger

from utilities.database_collection import uniprot_features, asa_calculator
from utilities.database_map_and_filter import create_uniprot_db,  uniprot_summary
from utilities.decorators import ProgressBar
from utilities.statistical_tests import apply_oneway_anova, apply_kruskal, fischers_test, apply_chisquare

logger.info('Import OK')

cluster_path = f'results/lysate_denaturation/clustering/clustered.xlsx'
predicted_folder = f'results/lysate_denaturation/predicted_features/'
calculated_folder = f'results/lysate_denaturation/uniprot_features/'
output_folder = 'results/lysate_denaturation/peptide_features/'

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

# ------------------------------Read in clustered data------------------------------
clustered_data = pd.read_excel(f'{cluster_path}', sheet_name='summary')
clustered_data.drop([col for col in clustered_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# --------------------------Process predicted peptide info--------------------------

# Create peptide number map to add to sequences of interest
peptides = pd.read_table(f'{predicted_folder}peptides.fasta', header=None)
sequence_map = dict(zip(peptides[0].str.strip('>')[np.arange(0, len(peptides), 2)], peptides[0][np.arange(1, len(peptides), 2)]))

# amino acid composition: frequency of all natural amino acids
pep_amino_acids = result_filter(result_type='AAC', entries=list(sequence_map.keys()), property_name=False, input_folder=f'{predicted_folder}peptides/')

pep_ctdc = result_filter(result_type='CTDC', entries=list(sequence_map.keys()), input_folder=f'{predicted_folder}peptides/')
column_dict = {
    'polarity.G1': 'Polar', 'polarity.G2': 'Neutral', 'polarity.G3': 'Hydrophobic',
    'hydrophobicity_PRAM900101.G1': 'Low', 'hydrophobicity_PRAM900101.G2': 'Medium', 'hydrophobicity_PRAM900101.G3': 'High',
    'charge.G1': 'charge_positive', 'charge.G2': 'charge_neutral', 'charge.G3': 'charge_negative',
    'secondarystruct.G1': 'Helix', 'secondarystruct.G2': 'Strand', 'secondarystruct.G3': 'Coil',
    'solventaccess.G1': 'Buried', 'solventaccess.G2': 'Exposed', 'solventaccess.G3': 'Intermediate',
    }
pep_ctdc.rename(columns=column_dict, inplace=True)

predicted = pd.merge(pep_amino_acids, pep_ctdc[list(column_dict.values())+['#']], on='#').rename(columns={'#': 'pep_id'})
predicted['Sequence'] = predicted['pep_id'].map(sequence_map)
predicted['polar_sum'] = predicted[['S', 'T', 'Y', 'N', 'Q']].sum(axis=1) # mimics figure from meltome dataset
predicted['hydrophobic_sum'] = predicted[['A', 'V', 'I', 'L', 'M', 'F', 'W']].sum(axis=1) # mimics figure from meltome dataset
predicted_features = pd.merge(clustered_data, predicted.drop('pep_id', axis=1).drop_duplicates().reset_index(drop=True), on=['Sequence'], how='inner')

features = list(column_dict.values())
features.append('polar_sum')
features.append('hydrophobic_sum')

# Due to non-parametric nature of the data, perform Kruskal–Wallis test with Dunn’s post hoc test, multiple testing correction according to Benjamini–Hochberg
predicted_stats = {}
predicted_stats['compiled_data'] = predicted_features
for feature in features:
    logger.info(f'Processing {feature}')
    for col in ['mixed', 'unique', 'count']:
        test_val, pval, posthoc = apply_kruskal(predicted_features, val_col=feature, group_col=col, p_adjust='bonferroni', sig_cutoff=0.05)
        predicted_stats[f'{feature}_{col}_kruskal'] = pd.DataFrame([test_val, pval], index=['test_statistic', 'pvalue'])
        predicted_stats[f'{feature}_{col}_posthoc'] = posthoc

FileHandling.df_to_excel(output_path=f'{output_folder}predicted_features_summary.xlsx', sheetnames=list(predicted_stats.keys()), data_frames=list(predicted_stats.values()))

# --------------------------Process calculated peptide info--------------------------

# read in uniprot cys peptide features
cys_peptides = pd.read_excel(f'{calculated_folder}cys_peptide_features.xlsx')
cys_peptides.drop([col for col in cys_peptides.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# drop non-mapped rows
uniprot_cols = [col for col in cys_peptides.columns.tolist() if col not in ['proteins', 'pep_sequence', 'cys_in_peptide']]
calculated_features = cys_peptides.copy().replace('None', np.nan).dropna(how='all', subset=uniprot_cols).dropna(subset=['cys_in_peptide']).rename(columns={'proteins': 'Proteins', 'pep_sequence': 'Sequence'})
calculated_features = calculated_features[(calculated_features['pdb_aa'].isin(['C', np.nan])) & (calculated_features['dssp_aa'].isin(['C', np.nan]))]


# calculate summary info - mean of numerical entries, highest frequency for non-numerical
calculated_summary = calculated_features.copy().groupby(['Proteins', 'Sequence', 'cys_in_peptide']).mean().reset_index()

calculated_summary = functools.reduce(lambda left, right: pd.merge(left, right, on=['Proteins', 'Sequence', 'cys_in_peptide'], how='outer'), [calculated_summary, pd.DataFrame(calculated_features.groupby(['Proteins', 'Sequence', 'cys_in_peptide'])['pdb_structure'].apply(list)), pd.DataFrame(calculated_features.groupby(['Proteins', 'Sequence', 'cys_in_peptide'])['pdb_disorder'].apply(list)), pd.DataFrame(calculated_features.groupby(['Proteins', 'Sequence', 'cys_in_peptide'])['dssp_structure'].apply(list)), pd.DataFrame(calculated_features.groupby(['Proteins', 'Sequence', 'cys_in_peptide'])['pfam_domains'].apply(list))])

def max_freq(elements):
    elements = [element for element in elements if type(element) == str]
    if len(elements) == 0:
        return np.nan
    return max(set(elements), key = elements.count)

calculated_summary['pdb_structure_freq'] = calculated_summary['pdb_structure'].map(max_freq)
calculated_summary['pdb_disorder_freq'] = calculated_summary['pdb_disorder'].map(max_freq)
calculated_summary['dssp_structure_freq'] = calculated_summary['dssp_structure'].map(max_freq)
calculated_summary['pfam_domains_freq'] = [[item for item in elements if type(item) == str] for elements in calculated_summary['pfam_domains']]
calculated_summary['pfam_domains_freq'] = calculated_summary['pfam_domains_freq'].map(max_freq)
calculated_summary['pfam_domains_freq'] = calculated_summary['pfam_domains_freq'].replace('[]', np.nan)

# add cluster info
calculated_summary = pd.merge(calculated_summary, clustered_data, on=['Proteins', 'Sequence'], how='outer')
# Using mixed, as these are per-residue features
col = 'mixed' 
stats = {}
stats = {}

# Number of sequences with mapped structures - chi-squared test
proportions = pd.merge(cys_peptides[['proteins', 'pep_sequence', 'pdb_aa']].drop_duplicates().copy().rename(columns={'proteins': 'Proteins', 'pep_sequence': 'Sequence'}).replace('None', np.nan), clustered_data, on=['Proteins', 'Sequence'], how='right') # only clustered peptides
proportions['feature_legend'] = [1 if structure == 'C' else 0 for structure in proportions['pdb_aa'].fillna('None')]
proportions = proportions.groupby(['Proteins', 'Sequence', 'mixed', 'unique', 'count']).max()['feature_legend'].reset_index()
proportion = (proportions.groupby([col, 'feature_legend']).count())['Proteins'].reset_index().rename(columns={'Proteins': 'proportion', col: 'cluster'})
proportion = pd.pivot_table(proportion, values='proportion', index='feature_legend', columns='cluster')
test_val, pval, dof, contingency, chi_test = apply_chisquare(proportion, p_adjust='bonferroni', sig_cutoff=0.05)
chi_test['group'] = col
stats['mapped_chi_ct'] = pd.DataFrame(contingency)
stats['mapped_chi_observed'] = proportion
stats['mapped_chi'] = chi_test

# limit to clustered residues with mapped sequences for remaining tests
calculated_clustered = calculated_summary[(~calculated_summary['cys_in_peptide'].isnull()) &(~calculated_summary['cluster'].isnull())].copy()

# proportion of residues with each structural feature - chi-squared test
proportions = calculated_clustered.copy()
simple_structure = {'H' : 'helix' , 'B' : 'strand' , 'E' : 'strand' , 'G' : 'helix' , 'I' : 'helix' , 'T' : 'turn' , 'S' : 'bend' , '-': '-'}

proportions['pdb_simple_structure'] = proportions['pdb_structure_freq'].map(simple_structure)
proportion_pdb = (proportions.groupby([col, 'pdb_simple_structure']).count())['Proteins'].reset_index().rename(columns={'Proteins': 'proportion', col: 'cluster'})
proportion_pdb = pd.pivot_table(proportion_pdb, values='proportion', index='pdb_simple_structure', columns='cluster')
test_val, pval, dof, contingency, chi_test = apply_chisquare(proportion_pdb, p_adjust='bonferroni', sig_cutoff=0.05)
chi_test['group'] = col
stats['pdb_structure_chi_ct'] = pd.DataFrame(contingency)
stats['pdb_structure_chi_observed'] = proportion_pdb
stats['pdb_structure_chi'] = chi_test

proportions['dssp_simple_structure'] = proportions['pdb_structure_freq'].map(simple_structure)
proportion_dssp = (proportions.groupby([col, 'dssp_simple_structure']).count())['Proteins'].reset_index().rename(columns={'Proteins': 'proportion', col: 'cluster'})
proportion_dssp = pd.pivot_table(proportion_dssp, values='proportion', index='dssp_simple_structure', columns='cluster')
test_val, pval, dof, contingency, chi_test = apply_chisquare(proportion_dssp, p_adjust='bonferroni', sig_cutoff=0.05)
chi_test['group'] = col
stats['dssp_structure_chi_ct'] = pd.DataFrame(contingency)
stats['dssp_structure_chi_observed'] = proportion_dssp
stats['dssp_structure_chi'] = chi_test


# proportion of residues with disorder - chi-squared test
proportions = calculated_clustered.copy().dropna(subset=['pdb_disorder_freq'])

proportion = (proportions.groupby([col, 'pdb_disorder_freq']).count())['Proteins'].reset_index().rename(columns={'Proteins': 'proportion', col: 'cluster'})
proportion = pd.pivot_table(proportion, values='proportion', index='pdb_disorder_freq', columns='cluster')
test_val, pval, dof, contingency, chi_test = apply_chisquare(proportion, p_adjust='bonferroni', sig_cutoff=0.05)
chi_test['group'] = col
stats['disorder_chi_ct'] = pd.DataFrame(contingency)
stats['disorder_chi_observed'] = proportion
stats['disorder_chi'] = chi_test


# proportion of residues within structured domains - chi-squared test
proportions = calculated_clustered.copy()
proportions['feature_legend'] = [0 if structure == 'None' else 1 for structure in proportions['pfam_domains_freq'].fillna('None')]
proportion = (proportions.groupby([col, 'feature_legend']).count())['Proteins'].reset_index().rename(columns={'Proteins': 'proportion', col: 'cluster'})
proportion = pd.pivot_table(proportion, values='proportion', index='feature_legend', columns='cluster')
test_val, pval, dof, contingency, chi_test = apply_chisquare(proportion, p_adjust='bonferroni', sig_cutoff=0.05)
chi_test['group'] = col
stats['domains_chi_ct'] = pd.DataFrame(contingency)
stats['domains_chi_observed'] = proportion
stats['domains_chi'] = chi_test


# mean solvent exposure (DSSP/PDB) - kruskal-Wallis/ANOVA w Bonferroni correction
accessibility = calculated_clustered.copy()

test_val, pval, posthoc = apply_kruskal(accessibility, val_col='pdb_asa', group_col=col, p_adjust='bonferroni', sig_cutoff=0.05)
accessibility_pdb = pd.DataFrame([test_val, pval], index=['test_statistic', 'pvalue']).T
accessibility_pdb.index = [f'accessibility_pdb']
accessibility_pdb['group'] = col
stats['pdb_asa_kruskal'] = accessibility_pdb

test_val, pval, posthoc = apply_kruskal(accessibility, val_col='dssp_asa', group_col=col, p_adjust='bonferroni', sig_cutoff=0.05)
accessibility_dssp = pd.DataFrame([test_val, pval], index=['test_statistic', 'pvalue']).T
accessibility_dssp.index = [f'accessibility_dssp']
accessibility_dssp['group'] = col
stats['dssp_asa_kruskal'] = accessibility_dssp

# add summary dfs
stats['compiled_data'] = calculated_features
stats['compiled_summary'] = calculated_summary
stats['compiled_clustered'] = calculated_clustered



# --------------------------Process cys residue features from UniProt/PFAM--------------------------
uniprot_features = pd.read_excel(f'{calculated_folder}uniprot_feature_summary.xlsx', sheet_name=None)

# cleanup to remove non-cys peptides, and collect unique peptides
residue_features = uniprot_features['residues'].copy()
residue_features = residue_features[~residue_features['cys_in_peptide'].isnull()]
residue_features = residue_features.drop([col for col in residue_features.columns.tolist() if 'Unnamed: ' in col], axis=1).rename(columns={'pep_sequence': 'Sequence'})

feature_legend = {'None': 0, 'active site': 1, 'metal ion-binding site': 1, 'modified residue': 1, 'sequence conflict': 0, 'disulfide bond': 1, 'lipid moiety-binding region': 1, 'mutagenesis site': 0, 'binding site': 1}
residue_features['feature_legend'] = residue_features['feature_type'].map(feature_legend)
residue_features = residue_features[['Proteins', 'Sequence', 'cys_in_full', 'feature_legend']].drop_duplicates().groupby(['Proteins', 'Sequence', 'cys_in_full']).max().reset_index() # remove duplicates, collect max
clustered_features = pd.merge(clustered_data, residue_features, on=['Sequence', 'Proteins'], how='left')

# Proportion on cys residues with each feature type - chi-squared
proportions = clustered_features.copy().groupby([col, 'feature_legend']).count()['Proteins'].reset_index().rename(columns={'Proteins': 'proportion', col: 'cluster'})
proportion = pd.pivot_table(proportions, values='proportion', index='feature_legend', columns='cluster')
test_val, pval, dof, contingency, chi_test = apply_chisquare(proportion, p_adjust='bonferroni', sig_cutoff=0.05)
chi_test['group'] = col
stats['residue_features_chi_ct'] = pd.DataFrame(contingency)
stats['residue_features_chi_observed'] = proportion
stats['residue_features_chi'] = chi_test

stats['residue_features'] = residue_features

# Complete domain-level analysis
# cleanup to remove non-cys peptides, and collect unique peptides
domain_features = uniprot_features['domains'].copy()
domain_features = domain_features[~domain_features['cys_in_peptide'].isnull()]
domain_features = domain_features.drop([col for col in domain_features.columns.tolist() if 'Unnamed: ' in col], axis=1).rename(columns={'pep_sequence': 'Sequence'})

feature_legend = {'None': 0, 'chain': 1, 'strand': 1, 'domain': 1, 'region of interest': 0,
       'compositionally biased region': 2, 'helix': 1, 'topological domain': 1,
       'intramembrane region': 1, 'repeat': 0, 'coiled-coil region': 1,
       'zinc finger region': 1, 'splice variant': 0, 'short sequence motif': 0,
       'disulfide bond': 2, 'nucleotide phosphate-binding region': 2,
       'mutagenesis site': 0, 'peptide': 0, 'turn': 1, 'transmembrane region': 1,
       'sequence conflict': 0, 'DNA-binding region': 2, 'propeptide': 0}
domain_features['feature_legend'] = domain_features['feature_type'].map(feature_legend)
domain_features = domain_features[['Proteins', 'Sequence', 'cys_in_full', 'feature_legend']].drop_duplicates().groupby(['Proteins', 'Sequence', 'cys_in_full']).max().reset_index() # remove duplicates, collect max
clustered_features = pd.merge(clustered_data, domain_features, on=['Sequence', 'Proteins'], how='left')

# Proportion on cys residues with each feature type - chi-squared
proportions = clustered_features.copy().groupby([col, 'feature_legend']).count()['Proteins'].reset_index().rename(columns={'Proteins': 'proportion', col: 'cluster'})
proportion = pd.pivot_table(proportions, values='proportion', index='feature_legend', columns='cluster')
test_val, pval, dof, contingency, chi_test = apply_chisquare(proportion, p_adjust='bonferroni', sig_cutoff=0.05)
chi_test['group'] = col
stats['domain_features_chi_ct'] = pd.DataFrame(contingency)
stats['domain_features_chi_observed'] = proportion
stats['domain_features_chi'] = chi_test

stats['domain_features'] = domain_features


# Complete domain-level analysis
# cleanup to remove non-cys peptides, and collect unique peptides
pfam = uniprot_features['pfam'].copy()
pfam = pfam[~pfam['cys_in_peptide'].isnull()]
pfam = pfam.drop([col for col in pfam.columns.tolist() if 'Unnamed: ' in col], axis=1).rename(columns={'pep_sequence': 'Sequence', 'pfam_in_domain': 'feature_legend'})

pfam = pfam[['Proteins', 'Sequence', 'cys_in_full', 'feature_legend']].drop_duplicates().groupby(['Proteins', 'Sequence', 'cys_in_full']).max().reset_index() # remove duplicates, collect max
clustered_features = pd.merge(clustered_data, pfam, on=['Sequence', 'Proteins'], how='left')

# Proportion on cys residues with each feature type - chi-squared
proportions = clustered_features.copy().groupby([col, 'feature_legend']).count()['Proteins'].reset_index().rename(columns={'Proteins': 'proportion', col: 'cluster'})
proportion = pd.pivot_table(proportions, values='proportion', index='feature_legend', columns='cluster')
test_val, pval, dof, contingency, chi_test = apply_chisquare(proportion, p_adjust='bonferroni', sig_cutoff=0.05)
chi_test['group'] = col
stats['pfam_chi_ct'] = pd.DataFrame(contingency)
stats['pfam_chi_observed'] = proportion
stats['pfam_chi'] = chi_test

stats['pfam'] = pfam

# ----------------------------Save to excel----------------------------

FileHandling.df_to_excel(
    output_path=f'{output_folder}calculated_peptide_features_summary.xlsx',
    sheetnames=list(stats.keys()),
    data_frames=list(stats.values())
)

