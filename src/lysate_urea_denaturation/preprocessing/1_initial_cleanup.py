
import pandas as pd
import numpy as np
import os, re, collections
import functools

from GEN_Utils.CalcUtils import sorted_nicely
from GEN_Utils import FileHandling

from loguru import logger

logger.info('Import OK')

input_folder = 'raw_data/lysate_urea_denaturation/'
output_folder = 'results/lysate_denaturation/initial_cleanup/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def dict_tree():
    return collections.defaultdict(dict_tree)



def pd_peptides(input_path, sample_names=None, pep_cols=[], peptides_file='_Peptides.xlsx'):

    peptides = pd.read_excel(f'{input_path}{peptides_file}')
    logger.info(
        f'Imported peptides from {input_path}. {peptides.shape[0]} entries found.')
    
    # filter out non-unique, missed cleavage, reverse and contaminant peptides
    filtered_pep = peptides[(peptides['Quan Info'] == 'Unique')].rename(
        columns={'Master Protein Accessions': 'Proteins'})
    # filtered_pep = filtered_pep[filtered_pep['# Missed Cleavages'] == 0]

    # Clean sequence column
    filtered_pep['Sequence'] = filtered_pep['Annotated Sequence'].str.split('.').str[1]
    
    
    # add cys rank
    filtered_pep['cys_rank'] = [
        1 if 'C' in pep else 0 for pep in filtered_pep['Sequence']]

    # Collect sample columns and rename
    logger.info(f'Sample names not set. Collecting all samples.')
    logger.debug(f'Columns found: {filtered_pep.columns.tolist()}')
    sample_cols = [x for x in filtered_pep.columns.tolist() if 'Abundance Ratio: (' in x]
    sample_names = [x.replace('Abundance Ratio: (', '').split(
        ',')[0] for x in sample_cols]
    filtered_pep.rename(columns=dict(zip(sample_cols, sample_names)), inplace=True)
    logger.info(f'Samples detected: {sample_names}')

    # filter those without any values for value in variable:
    # In theory (and in previous MQ cleanup), drop any peptides with any missing values here? To start with better to only drop ones that are all missing
    filtered_pep = filtered_pep.replace(0, np.nan).dropna(
        axis=0, how='all', subset=sample_names)

    # Calculate average ratio for peptides which have multiple different modifications
    # Note this may not be the best solution - probably more appropriate to repeat search and report only the summed modified abundance

    filtered_pep = filtered_pep.groupby(['Sequence', 'Proteins']).mean().reset_index()

    cleaned_df = filtered_pep[['Sequence', 'Proteins'] + sample_names]

    logger.info(f'Successfully cleaned peptide dataframe.')

    return cleaned_df


def pd_proteins(input_path, sample_names=None, prot_cols=[], proteins_file='_Proteins'):

    logger.info(f'Collecting proteins')
    proteins = pd.read_excel(f'{input_path}{proteins_file}')
    logger.info(
        f'Imported proteins from {input_path}. {proteins.shape[0]} entries found.')

    # remove non-unique protein groups
    filtered_prot = proteins[(proteins['# Protein Groups'] == 1)].rename(
        columns={'Accession': 'Proteins'})
    logger.info(
        f'Removed non-unique proteins: {proteins.shape[0]} entries remain.')

    # Collect sample columns and rename
    logger.info(f'Sample names not set. Collecting all samples.')
    logger.debug(f'Columns found: {filtered_prot.columns.tolist()}')
    sample_cols = [x for x in filtered_prot.columns.tolist(
    ) if 'Abundance Ratio: (' in x]
    sample_names = [x.replace('Abundance Ratio: (', '').split(
        ',')[0] for x in sample_cols]
    filtered_prot.rename(columns=dict(
        zip(sample_cols, sample_names)), inplace=True)
    logger.info(f'Samples detected: {sample_names}')

    # filter those without any values for value in variable:
    # In theory (and in previous MQ cleanup), drop any peptides with any missing values here? To start with better to only drop ones that are all missing
    filtered_prot = filtered_prot.replace(0, np.nan).dropna(
        axis=0, how='all', subset=sample_names)

    cleaned_df = filtered_prot[['Proteins', 'Description'] + sample_names]

    logger.info(f'Successfully cleaned proteins dataframe.')

    return cleaned_df

def pd_cleaner(input_folder, output_path, proteins_file='_Proteins.xlsx', peptides_file='_Peptides.xlsx', prot_cols=[], pep_cols=[]):


    sample_list = [filename.split(peptides_file)[0] for filename in os.listdir(input_folder) if peptides_file in filename]

    cleaned_dfs = dict_tree()
    for sample in sample_list:
        sample_name = sample.split('_')[-1]
        cleaned_dfs[sample_name]['Peptides'] = pd_peptides(input_path=f'{input_folder}{sample}', pep_cols=pep_cols, peptides_file='_Peptides.xlsx')
        cleaned_dfs[sample_name]['Proteins'] = pd_proteins(input_path=f'{input_folder}{sample}', prot_cols=prot_cols, proteins_file='_Proteins.xlsx')

    logger.info(f'Sorting cleaned data per sample...')

    ## Collecting specific results for each set of samples for compiling
    prots = []
    peps = []
    replicate = 1
    for sample_name in cleaned_dfs.keys():

        #collect dataframes, rename relevant columns to account for replicates
        pep_dataframe = cleaned_dfs[sample_name]['Peptides']
        prot_dataframe = cleaned_dfs[sample_name]['Proteins']

        old_cols = sorted_nicely([col for col in pep_dataframe.columns.tolist() if 'F' in col])
        new_cols = [f'{x}_{replicate}' for x in range(1, len(old_cols)+1)]
        peps.append(pep_dataframe.rename(columns=dict(zip(old_cols, new_cols))))
        prots.append(prot_dataframe.rename(columns=dict(zip(old_cols, new_cols))))

        replicate += 1

    # compile replicates into single file
    proteins = functools.reduce(lambda left, right: pd.merge(
        left, right, on=['Proteins', 'Description'], how='outer'), prots)

    peptides = functools.reduce(lambda left, right: pd.merge(left, right, on=['Sequence', 'Proteins'], how='outer'), peps)

        #save to individual excel spreadsheets
    FileHandling.df_to_excel(f'{output_folder}Compiled.xlsx', sheetnames=[
                                'Proteins', 'Peptides'], data_frames=[proteins, peptides])

    logger.info(f'Proteins and peptides successfully cleaned. Dataframes save to {output_folder}.')

    return cleaned_dfs

if __name__ == "__main__":


    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    pd_cleaner(input_folder, output_folder, proteins_file='_Proteins.xlsx',
                peptides_file='_Peptides.xlsx', prot_cols=[], pep_cols=[])
