
import pandas as pd
import numpy as np
import os
import re

from ProteomicsUtils import FileHandling, StatUtils, PlotUtils, DataWrangling

from loguru import logger

logger.info('Import OK')

pd.set_option('display.max_columns', 500)


def tmt_peptides(input_path, sample_names=None, pep_cols=[]):

    peptides = pd.read_table(input_path, sep='\t')
    logger.info(
        f'Imported peptides from {input_path}. {peptides.shape[0]} entries found.')
    cleaned_dfs = {}
    standard_cols = ['Sequence', 'Proteins', 'Gene names',
                     'Protein names', 'Unique (Groups)', 'Unique (Proteins)', ] + pep_cols

    # filter out non-unique, reverse and contaminant peptides
    filtered_pep = peptides[(peptides['Reverse'] != '+') & (
        peptides['Potential contaminant'] != '+') & (peptides['Unique (Proteins)'] == 'yes')]
    # add cys rank
    filtered_pep['cys_rank'] = [
        1 if 'C' in pep else 0 for pep in filtered_pep['Sequence']]

    if sample_names is None:
        logger.info(f'Sample names not set. Collecting all samples.')
        logger.debug(f'Columns found: {peptides.columns.tolist()}')
        sample_names = [x.replace('Experiment ', '').split(
            '_')[:-1] for x in peptides.columns.tolist() if 'Experiment ' in x]
        sample_names = list(set([('_').join(x) for x in sample_names]))
        logger.info(f'Samples detected: {sample_names}')

    cleaned_dfs = {}
    standard_cols = ['Sequence', 'Proteins', 'Gene names',
                     'Protein names', 'Unique (Groups)', 'Unique (Proteins)', ] + pep_cols
    for sample in sample_names:
        sample_df = filtered_pep[standard_cols +
                                 [x for x in filtered_pep.columns.tolist() if sample in x]]
        sample_cols = [col for col in sample_df.columns.tolist(
        ) if 'Reporter intensity corrected' in col]
        # filter those without any values for value in variable:
        # In theory (and in previous MQ cleanup), drop any peptides with any missing values here? To start with better to only drop ones that are all missing
        sample_df = sample_df.replace(0, np.nan).dropna(
            axis=0, how='all', subset=sample_cols)
        cleaned_dfs[sample] = sample_df[standard_cols + sample_cols]
    logger.info(f'Successfully cleaned peptide dataframe.')

    return cleaned_dfs


def tmt_proteins(input_path, sample_names=None, prot_cols=[]):

    logger.info(f'Collecting proteins')
    proteins = pd.read_table(input_path, sep='\t')
    logger.info(
        f'Imported proteins from {input_path}. {proteins.shape[0]} entries found.')

    # remove contaminant and reverse proteins
    proteins = proteins[(proteins['Reverse'] != '+') &
                        (proteins['Potential contaminant'] != '+')]
    logger.info(
        f'Removed contaminant and reverse proteins: {proteins.shape[0]} entries remain.')

    cleaned_prot_dfs = {}
    standard_cols = ['Protein IDs', 'Gene names',
                     'Protein names', 'Number of proteins'] + prot_cols

    for sample in sample_names:
        sample_cols = standard_cols + \
            [x for x in proteins.columns.tolist() if sample in x]
        sample_df = proteins[sample_cols]
        ratio_cols = [col for col in sample_df.columns.tolist(
        ) if 'Reporter intensity corrected' in col]
        logger.debug(f'Ratio cols: {ratio_cols}')
        #collect columns of interest
        sample_vals = proteins[sample_cols + ratio_cols]
        #collect only proteins with at least one quantification in that sample
        sample_reps = sample_df.replace(0, np.nan).dropna(
            axis=0, how='all', subset=ratio_cols)
        logger.debug(f'Sample reps: {sample_reps.head(10)}')
        # collect only proteins which are master proteins
        master_proteins = sample_reps[sample_reps['Number of proteins'] == 1]
        logger.debug(f'Master proteins: {master_proteins.head(10)}')
        cleaned_prot_dfs[sample] = master_proteins

    logger.info(f'Successfully cleaned proteins dataframe.')

    return cleaned_prot_dfs

def tmt_cleaner(input_folder, output_path, sample_names=None, proteins_file='proteinGroups.txt', peptides_file='peptides.txt', prot_cols=[], pep_cols=[]):

    cleaned_peptides = tmt_peptides(input_path=f'{input_folder}{peptides_file}', sample_names=sample_names, pep_cols=pep_cols)
    cleaned_proteins = tmt_proteins(input_path=f'{input_folder}{proteins_file}', sample_names=sample_names, prot_cols=prot_cols)

    logger.info(f'Sorting cleaned data per sample...')
    ## Collecting specific results for each set of samples for further processing
    for sample in sample_names:
        #collect peptide dataframe, rename relevant columns
        pep_dataframe = cleaned_peptides[sample]
        prot_dataframe = cleaned_proteins[sample]
        #save to individual excel spreadsheets
        FileHandling.df_to_excel(f'{output_folder}{sample}_Compiled.xlsx', sheetnames=[
                                 'Proteins', 'Peptides'], data_frames=[prot_dataframe, pep_dataframe])

    logger.info(f'Proteins and peptides successfully cleaned. Dataframes save to {output_folder}.')

    return cleaned_peptides, cleaned_proteins

if __name__ == "__main__":

    input_folder = 'raw_data/inhibitor_urea_denaturation/txt/'
    output_folder = 'results/inhibitor_urea_denaturation/initial_cleanup/'

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    tmt_cleaner(input_folder, output_folder, sample_names=['Control', 'VER'], proteins_file='proteinGroups.txt', peptides_file='peptides.txt', prot_cols=[], pep_cols=[])
