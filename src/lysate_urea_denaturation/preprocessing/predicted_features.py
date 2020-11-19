"""
Using iFeature online version, as the command-line version looks hard to wrangle
In this case, collect fasta sequences and peptide sequences for entire 'background' dataset
and search using online interface.
"""

import gzip, shutil, os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from distutils.dir_util import copy_tree

from Bio import SeqIO
import phylopandas as ph

from utilities.database_map_and_filter import tar_file_to_folder
from utilities.database_collection import disorder_prediction_iupred
from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')


background_path = f'results/lysate_denaturation/normalised/normalised_summary.xlsx'
input_folder = f'raw_data/lysate_urea_denaturation/predicted_features/'
output_folder = 'results/lysate_denaturation/predicted_features/'

resource_folder = 'resources/bioinformatics_databases/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def build_peptide_fasta(peptides, output_path):
    peptide_file = []
    for protein, sequence_df in peptides.groupby('Proteins'):
        sequences = sequence_df['Sequence'].tolist()
        for x, sequence in enumerate(sequences):
            peptide_file.append(f'>{protein}_{x}')
            peptide_file.append(f'{sequence}')
    with open(f'{output_path}', 'w') as filehandle:
        filehandle.writelines(f"{line}\n" for line in peptide_file)


def build_protein_fasta(genes, tax_id, output_path, resource_folder=resource_folder):
    full_fasta = ph.read_fasta(f'{resource_folder}{tax_id}.fasta')
    full_fasta['Entry'] = full_fasta['id'].str.split('|').str[1]
    search_genes = '|'.join(genes)
    fasta = full_fasta[full_fasta['Entry'].str.contains(search_genes)]
    protein_file = []
    for protein, df in fasta.groupby('Entry'):
        sequence = df['sequence'].tolist()[0]
        protein_file.append(f'>{protein}')
        protein_file.append(f'{sequence}')
    with open(f'{output_path}', 'w') as filehandle:
        filehandle.writelines(f"{line}\n" for line in protein_file)



def result_filter(result_type, entries, property_name=False, input_folder=output_folder):
    result = pd.read_table(f'{input_folder}{result_type}.tsv').set_index('#')
    if property_name:
        property_cols = [col for col in result.columns.tolist() if property_name in col]
    else:
        property_cols = result.columns.tolist()
    filtered = result[property_cols].reset_index()
    filtered = filtered[filtered['#'].isin(entries)]
    return filtered


# ----------------Generate FASTA files to search with iFeature----------------

# # Read in background proteins/peptides
# background = pd.read_excel(background_path, sheet_name='control_norm_cys_ratio')
# background = background.copy().drop([col for col in background.columns.tolist() if 'Unnamed: ' in col], axis=1)

# # generate fasta file for all peptides
# peptides = background[['Proteins', 'Sequence']].drop_duplicates()
# build_peptide_fasta(peptides, output_path=f'{output_folder}peptides.fasta')

# # generate fasta file for all proteins
# proteins = background['Proteins'].unique()
# build_protein_fasta(proteins, tax_id='10090', output_path=f'{output_folder}proteins.fasta')

# process peptide and protein fasta files using iFeature web server: http://ifeature.erc.monash.edu/center.php?page=descriptors
# Selected all options that do not require symmetrical sequences, including K-means clustering and DBSCAN clustering
# Protein job ID is: 20200517100852
# Peptide job ID is: 20200517103850
# collected resulting folder and placed into output folder

# Open tar file into output folder
# tar_file_to_folder(input_path=f'{output_folder}proteins.tar', output_path=f'{input_folder}proteins')
# tar_file_to_folder(input_path=f'{output_folder}peptides.tar', output_path=f'{input_folder}peptides')

# Copy into results folder
copy_tree(f'{input_folder}peptides/', f'{output_folder}peptides/')
copy_tree(f'{input_folder}proteins/', f'{output_folder}proteins/')

# --------------------------------------Disorder Prediction from IUPred--------------------------------------
# using pre-imported protein list from raw data 

disorder_results = disorder_prediction_iupred(proteins, output_folder)
FileHandling.df_to_excel(output_path=f'{output_folder}disorder_prediction.xlsx', data_frames=[disorder_results], sheetnames=['IUPRed'])
