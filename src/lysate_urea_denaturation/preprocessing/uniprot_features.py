import gzip
import json
import os
import re
import shutil
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import phylopandas as ph
import seaborn as sns
from Bio import pairwise2
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from fuzzysearch import find_near_matches
from IPython.display import Image
from pymol import cmd

from loguru import logger
from GEN_Utils import FileHandling

from utilities.database_map_and_filter import create_uniprot_db,  uniprot_summary
from utilities.database_collection import (uniprot_features, uniprot_sequences, uniprot_function, uniprot_download, asa_calculator, pfam_domains)
from utilities.decorators import ProgressBar

logger.info('Import OK')

background_path = f'results/lysate_denaturation/smoothing/processed_data.xlsx'
resource_folder = 'resources/bioinformatics_databases/'

output_folder = 'results/lysate_denaturation/uniprot_features/'
tax_id = '10090'


if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def mol_weight(df, seq_col):
    """Calculate molecular weight from sequence"""
    weights = []
    for sequence in df[seq_col].tolist():
        try:
            analysed_seq = ProteinAnalysis(sequence)
            weights.append(analysed_seq.molecular_weight())
        except:
            weights.append(np.nan)
    df['MW'] = weights
    return df

# ------------------------------Read in background proteins------------------------------
raw_data = pd.read_excel(background_path, sheet_name='raw')
raw_data = raw_data.drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in str(col)], axis=1)
proteins = raw_data['Proteins'].unique().tolist()

# --------------------Generate UniProt summary with mapped structure--------------------
uniprot_info = uniprot_summary(tax_id=tax_id, resource_folder=resource_folder, genes=proteins, reviewed=True)
uniprot_info = pd.merge(raw_data[['Proteins', 'Sequence']].drop_duplicates(), uniprot_info, left_on='Proteins', right_on='Entry')
pdb_info = uniprot_info[['Entry', 'PDB']].copy().dropna().drop_duplicates()
pdb_info['structures'] = pdb_info['PDB'].str.split('; ')
pdb_info = pdb_info.explode('structures')
pdb_info[['pdb_id', 'chain']] = pdb_info['structures'].str.split(':', expand=True)

# ----------------------------Collect generic structural features----------------------------

# Collect complete sequences
full_fasta = ph.read_fasta(f'{resource_folder}{tax_id}.fasta')
full_fasta['Entry'] = full_fasta['id'].str.split('|').str[1]
full_sequences = dict(zip(full_fasta['Entry'], full_fasta['sequence']))

uniprot_info['full_sequence'] = uniprot_info['Entry'].map(full_sequences)
uniprot_info = mol_weight(uniprot_info, 'full_sequence')

uniprot_info.to_csv(f'{output_folder}uniprot_summary.csv')

# --------------------Option 1: Collect ASA and Stucture from DSSP--------------------
# Generate file list to grab from dssp
structures = '\n'.join([f'{entry}.dssp' for entry in pdb_info['pdb_id'].str.lower().unique().tolist()])
# with open('C:/Users/Dezerae/Desktop/dssp_files.txt', 'w') as f:
#     f.write(structures)
# # in Ubuntu (WSL):
# rsync -avz rsync://rsync.cmbi.umcn.nl/dssp/ </path/to/dssp/folder/> --files-from=files.txt
# move into resource folder

# Retrieve dssp info for structures of interest from pre-computed files
structures = pdb_info['pdb_id'].str.lower().unique().tolist()
dssp_dfs = []
header_columns = [0, 6, 13, 15, 26, 30, 35, 39, 51, 62, 73, 84, 92, 98, 104, 110, 116, 123, 130, 137, 154, 165]
header = ['#', 'RESIDUE', 'AA', 'STRUCTURE', 'BP1', 'BP2', 'ACC', 'N-H-->O', 'O-->H-N', 'N-H-->O', 'O-->H-N', 'TCO', 'KAPPA', 'ALPHA', 'PHI', 'PSI', 'X-CA', 'Y-CA', 'Z-CA', 'CHAIN', 'AUTHCHAIN']

for structure in structures:
    try:
        dssp_raw = pd.read_table(f'{resource_folder}dssp/{structure}.dssp', skiprows=27, header=None)
        dssp_raw.drop(0, inplace=True) # remove header row
        dssp_raw[~dssp_raw[0].str.contains("!")]# remove extra rows indicating chain swaps
        dssp_raw = dssp_raw[0].str.split('', expand=True)
        dssp = pd.DataFrame(columns=header)
        for x, entry in enumerate(header):
            start, stop = (header_columns[x], header_columns[x+1])
            dssp[entry] = dssp_raw.iloc[:, start:stop].astype(str).apply(lambda x: ''.join(x), axis=1)
        dssp['pdb_id'] = structure

        dssp_dfs.append(dssp)
    except:
        logger.info(f'No file located for {structure}.')

dssp_summary = pd.concat(dssp_dfs)
dssp_summary.rename(columns={'#': 'dssp_#',
                    'AA': 'residue_aa',
                    'STRUCTURE': 'residue_structure',
                    'ACC': 'residue_asa',
                    'CHAIN': 'chain'}, inplace=True)
dssp_summary['structure_type'] = dssp_summary['residue_structure'].str.split('').str[3].str.replace(' ', '-')
dssp_summary['chain'] = dssp_summary['chain'].str.strip()
dssp_summary['pdb_id'] = dssp_summary['pdb_id'].str.upper()
dssp_summary['residue_position'] = dssp_summary['RESIDUE'].str.extract('(\d+)').astype(float)
dssp_summary.to_csv(f'{resource_folder}dssp/dssp_summary.csv')

# -----------------------Option 2A: Collect ASA from pymol-----------------------

# Collect ASA values calculated by PyMol - if no asa values found, calculate for all otherwise those
# that are missing
asa_folder = f'{resource_folder}pymol_asa/'
if not os.path.exists(asa_folder):
    asa_calculator(pdb_info, resource_folder=resource_folder)
file_list = [filename.strip('.csv') for filename in os.listdir(f'{resource_folder}pymol_asa/') if '.csv' in filename]
for structure in pdb_info['pdb_id'].unique():
    if structure not in file_list:
        asa_calculator(pdb_info[pdb_info['pdb_id'] == structure], resource_folder=asa_folder)

# compile calculations for all structures
asa_summary = pd.concat([pd.read_csv(f'{asa_folder}{filename}', dtype={'pdb': 'str'}) for filename in os.listdir(f'{asa_folder}') if '.csv' in filename])
asa_summary.drop([col for col in asa_summary.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# Collect sequence information
sequences = ph.read_fasta(f'{resource_folder}pdb_seqres.txt')
sequences[['pdb_id', 'chain']] = sequences['id'].str.split('_', expand=True)
sequences['pdb_id'] = sequences['pdb_id'].str.upper()

# match sequences to structures
dfs = []
sequence_failed = []
for (structure, chain), df in asa_summary.groupby(['pdb', 'chain']):
    try:
        sequence = sequences[sequences['id'].str.upper() == f'{structure}_{chain}']['sequence'].values[0]
        asa = df.copy()
        asa = asa[asa['aa'] == chain]
        asa['aa'] = [aa for aa in sequence]
        dfs.append(asa)
    except:
        sequence_failed.append((structure, chain)) # 14/756 failed
asa_summary = pd.concat(dfs)
asa_summary.to_csv(f'{output_folder}asa_summary.csv')

#----------Option 2B: Collect structural information from PDB secondary structure---------

description = []
sequence = []
linetext = []
with open(f'{resource_folder}ss_dis.txt.') as f:    
    for line in f:
        if '>' in line:
            sequence.append(linetext)
            description.append(line.strip('>').strip('\n'))
            linetext = []

        else:
            linetext.append(line.strip('\n').replace(' ', '-'))
    sequence.append(linetext) # as the last sequence is not appended unless it goes through the loop

seqstructure = pd.DataFrame([''.join(entry) for entry in sequence[1:]], columns=['sequence'])
seqstructure['description'] = description
seqstructure[['pdb_id', 'chain', 'parameter']] = seqstructure['description'].str.split(':', expand=True)

seq_structure = seqstructure.set_index(['pdb_id', 'chain', 'parameter'])['sequence'].unstack(level=-1).reset_index()
seq_structure.rename(columns={'sequence': 'chain_sequence'}, inplace=True)

# Add full sequence info
seq_map = uniprot_info[['PDB', 'full_sequence']].copy()
seq_map['structures'] = seq_map['PDB'].str.split('; ')
seq_map = seq_map.explode('structures')
seq_map['structures'] = seq_map['structures'].str.split(':').str[0]
seq_map = dict(zip(seq_map['structures'], seq_map['full_sequence']))
seq_structure['full_sequence'] = seq_structure['pdb_id'].map(seq_map)
seq_structure.dropna(subset=['full_sequence'], inplace=True)
seq_structure.to_csv(f'{output_folder}secondary_structure_summary.csv')

# combine asa and structure info
secstr = seq_structure.set_index(['pdb_id', 'chain']).to_dict(orient='index')
pdb_summary = []
for (pdb_id, chain), df in asa_summary.groupby(['pdb', 'chain']):
    try:
        df['disorder'] = [entry for entry in secstr[(pdb_id, chain)]['disorder']]
        df['secstr'] = [entry for entry in secstr[(pdb_id, chain)]['secstr']]
        pdb_summary.append(df)
    except:
        df['disorder'] = np.nan
        df['secstr'] = np.nan
pdb_summary = pd.concat(pdb_summary)
pdb_summary.to_csv(f'{output_folder}pdb_summary.csv')

# ------------------Collect PFAM domain information from pdb or uniprot------------------

pfam_chains = pd.read_table(f'{resource_folder}hmmer_pdb_all.txt')
pfam_chains.rename(columns={'PDB_ID': 'pdb_id', 'CHAIN_ID': 'chain', 'PdbResNumStart': 'pfam_start', 'PdbResNumEnd': 'pfam_end', 'PFAM_ACC': 'pfam_id', 'PFAM_Name': 'pfam_name', 'PFAM_desc': 'pfam_description', 'eValue': 'evalue'}, inplace=True)
pfam_chains['chain'] = pfam_chains['chain'].astype(str)

pfam_chains.to_csv(f'{output_folder}pfam_structure.csv')




# -------------------------------Prepare info for peptides-------------------------------

# Collect per-peptide information
peptides = uniprot_info[['Entry', 'name', 'full_sequence', 'Length', 'MW', 'PDB', 'Sequence']].copy().rename(columns={'Sequence': 'pep_sequence', 'Entry': 'proteins', 'Length': 'full_length', 'PDB': 'pdb_id'}).drop_duplicates()
peptides['pep_start'] = [full.find(sequence) + 1 for full, sequence in zip(peptides['full_sequence'], peptides['pep_sequence'])]
peptides['pep_stop'] = [start + len(sequence) - 1 for start, sequence in zip(peptides['pep_start'], peptides['pep_sequence'])]

# Find cys positions in individual peptides
cys_peptide = []
for (index, peptide_info) in peptides.iterrows():
    peptide = peptide_info['pep_sequence']
    pep_start, pep_stop = (peptide_info['pep_start'], peptide_info['pep_stop'])
    cys_peptide.append([pos for pos, char in enumerate(peptide) if char == 'C'])
peptides['cys_in_peptide'] = cys_peptide

# Split PDB structures and explode peptides df
# peptides.dropna(subset=['pdb_id'], inplace=True)
info_cols = ['proteins', 'name', 'full_sequence', 'full_length', 'full_mw', 'pep_sequence', 'pep_start', 'pep_stop', 'cys_in_peptide']
peptides['pdb_id'] = peptides['pdb_id'].str.split('; ')
peptides = peptides.explode('pdb_id')#.explode('cys_in_peptide')
peptides[['pdb_id', 'pdb_chain']] = peptides['pdb_id'].str.split(':', expand=True)
# peptides['cys_in_full'] = peptides['cys_in_peptide'] + peptides['pep_start']

# Add chain info from pdb
pdb_pivot = []
for (structure, chain), df in pdb_summary.groupby(['pdb', 'chain']):
    pdb_pivot.append([structure, chain, ''.join(df['aa'].str.strip(' ').tolist()), df['secstr'].tolist(),  df['disorder'].tolist(), df['asa'].tolist()])
pdb_pivot = pd.DataFrame(pdb_pivot)
pdb_pivot.columns = ['pdb_id','pdb_chain','pdb_chain_sequence','pdb_chain_secstr', 'pdb_chain_disorder', 'pdb_chain_asa']

peptides = pd.merge(peptides, pdb_pivot, on=['pdb_id', 'pdb_chain'], how='left')

# Add chain info from dssp
dssp_pivot = []
for (structure, chain), df in dssp_summary.groupby(['pdb_id', 'chain']):
    dssp_pivot.append([structure, chain, ''.join(df['residue_aa'].str.strip(' ').tolist()), df['structure_type'].tolist(), df['residue_asa'].tolist()])
dssp_pivot = pd.DataFrame(dssp_pivot)
dssp_pivot.columns = ['pdb_id','pdb_chain','dssp_chain_sequence','dssp_chain_secstr','dssp_chain_asa']
dssp_pivot = dssp_pivot[~dssp_pivot['dssp_chain_sequence'].str.contains('!')]

peptides = pd.merge(peptides, dssp_pivot, on=['pdb_id', 'pdb_chain'], how='left')

# replace nan sequences with None
peptides.fillna('None', inplace=True)

# locate peptides within each chain sequence
peptides['pep_in_pdb_chain'] = [chain.find(sequence) for chain, sequence in zip(peptides['pdb_chain_sequence'], peptides['pep_sequence'])]
peptides['pep_in_pdb_match'] = [find_near_matches(sequence, chain, max_l_dist=2)[0].start if len(find_near_matches(sequence, chain, max_l_dist=2)) > 0 else -1 for chain, sequence in zip(peptides['pdb_chain_sequence'], peptides['pep_sequence'])]

peptides['pep_in_dssp_chain'] = [chain.find(sequence) for chain, sequence in zip(peptides['dssp_chain_sequence'], peptides['pep_sequence'])]
peptides['pep_in_dssp_match'] = [find_near_matches(sequence, chain, max_l_dist=2)[0].start if len(find_near_matches(sequence, chain, max_l_dist=2)) > 0 else -1 for chain, sequence in zip(peptides['dssp_chain_sequence'], peptides['pep_sequence'])]

position_cols = ['pep_in_pdb_chain','pep_in_pdb_match','pep_in_dssp_chain','pep_in_dssp_match']
for col in position_cols:
    peptides[col] = peptides[col].replace(-1, np.nan)

# # Find chain sequence within full protein
# peptides['pdb_chain_in_full'] = [full.find(chain) for full, chain in zip(peptides['full_sequence'], peptides['pdb_chain_sequence'])]
# peptides['pdb_chain_in_full_match'] = [find_near_matches(chain, full, max_l_dist=2)[0].start if len(find_near_matches(chain, full, max_l_dist=2)) > 0 else -1 for full, chain in zip(peptides['full_sequence'], peptides['pdb_chain_sequence'])]
# peptides['pdb_chain_in_full_align'] = [pairwise2.align.globalxx(full, chain) for full, chain in zip(peptides['full_sequence'], peptides['pdb_chain_sequence'])]

# peptides['dssp_chain_in_full'] = [full.find(chain) for full, chain in zip(peptides['full_sequence'], peptides['dssp_chain_sequence'])]
# peptides['dssp_chain_in_full_match'] = [find_near_matches(chain, full, max_l_dist=2)[0].start if len(find_near_matches(chain, full, max_l_dist=2)) > 0 else -1 for full, chain in zip(peptides['full_sequence'], peptides['dssp_chain_sequence'])]
# peptides['dssp_chain_in_full_align'] = [pairwise2.align.globalxx(full, chain) for full, chain in zip(peptides['full_sequence'], peptides['dssp_chain_sequence'])]

# locate cysteines within each peptide sequence
peptides = peptides.explode('cys_in_peptide')
peptides['cys_in_pdb_chain'] = peptides['cys_in_peptide'] + peptides['pep_in_pdb_match']
peptides['cys_in_dssp_chain'] = peptides['cys_in_peptide'] + peptides['pep_in_dssp_match']
peptides['cys_in_full'] = peptides['cys_in_peptide'] + peptides['pep_start']

peptides.to_csv(f'{output_folder}peptides_summary.csv')

# -------------Prepare specific residue information for each cys position-------------
pfam_summary = pfam_chains.copy().set_index(['pdb_id', 'chain'])
cys_info = []
for (index, peptide_info) in peptides.iterrows():

    protein, peptide, cys_in_peptide, cys_in_full, pdb_id, pdb_chain = (peptide_info['proteins'], peptide_info['pep_sequence'], peptide_info['cys_in_peptide'], peptide_info['cys_in_full'], peptide_info['pdb_id'], peptide_info['pdb_chain'])
    try:
        pdb_cys_position = int(peptide_info['cys_in_pdb_chain'])
        pdb_aa, pdb_structure, pdb_disorder, pdb_asa = (peptide_info['pdb_chain_sequence'][pdb_cys_position], peptide_info['pdb_chain_secstr'][pdb_cys_position], peptide_info['pdb_chain_disorder'][pdb_cys_position], peptide_info['pdb_chain_asa'][pdb_cys_position])
    except ValueError:
        pdb_aa, pdb_structure, pdb_disorder, pdb_asa = (np.nan, np.nan, np.nan, np.nan)

    try:
        dssp_cys_position = int(peptide_info['cys_in_dssp_chain'])
        dssp_aa, dssp_structure, dssp_asa = (peptide_info['dssp_chain_sequence'][dssp_cys_position], peptide_info['dssp_chain_secstr'][dssp_cys_position], peptide_info['dssp_chain_asa'][dssp_cys_position])
    except ValueError:
        dssp_aa, dssp_structure, dssp_asa = (np.nan, np.nan, np.nan)

    # PFAM domains based on position of cys in total protein
    try:
        domains = pfam_summary.loc[(pdb_id, pdb_chain)]
        domains = domains[(domains['pfam_start'].astype(int) < cys_in_full) & (cys_in_full < domains['pfam_end'].astype(int))]['pfam_id'].tolist()
    except KeyError:
        f'{pdb_structure, pdb_chain, cys_in_full} structure not found.'
        domains = np.nan

    cys_info.append([protein, peptide, cys_in_peptide, cys_in_full, pdb_id, pdb_chain, pdb_aa, pdb_structure, pdb_disorder, pdb_asa, dssp_aa, dssp_structure, dssp_asa, domains])
cys_info = pd.DataFrame(cys_info, columns=['proteins', 'pep_sequence', 'cys_in_peptide', 'cys_in_full', 'pdb_id', 'pdb_chain', 'pdb_aa', 'pdb_structure', 'pdb_disorder', 'pdb_asa', 'dssp_aa', 'dssp_structure', 'dssp_asa', 'pfam_domains'])

cys_info.to_csv(f'{output_folder}cys_info_summary.csv')

FileHandling.df_to_excel(output_path=f'{output_folder}cys_peptide_features.xlsx', data_frames=[cys_info], sheetnames=['cys_info'])

# ---------Collect residue-specific and domain-specific features from uniprot, domains from PFAM---------

# Read in background data to generate list of protein sequences to map
peptides = raw_data[['Proteins', 'Sequence']].rename(columns={'Sequence': 'pep_sequence'})
peptides = peptides.drop_duplicates()

# add full sequence info
genes = peptides['Proteins'].unique()
full_sequences = uniprot_sequences(genes, resource_folder)
peptides['full_sequence'] = peptides['Proteins'].map(full_sequences)
peptides = mol_weight(peptides, seq_col='full_sequence')


peptides['pep_start'] = [full.find(sequence) + 1 for full, sequence in zip(peptides['full_sequence'], peptides['pep_sequence'])]
peptides['pep_stop'] = [start + len(sequence) - 1 for start, sequence in zip(peptides['pep_start'], peptides['pep_sequence'])]

# Find cys positions in individual peptides
cys_peptide = []
for (index, peptide_info) in peptides.iterrows():
    peptide = peptide_info['pep_sequence']
    pep_start, pep_stop = (peptide_info['pep_start'], peptide_info['pep_stop'])
    cys_peptide.append([pos for pos, char in enumerate(peptide) if char == 'C'])
peptides['cys_in_peptide'] = cys_peptide

# Split PDB structures and explode peptides df
peptides = peptides.explode('cys_in_peptide')
peptides['cys_in_full'] = peptides['cys_in_peptide'] + peptides['pep_start']

# Collect protein active site features from uniprot
genes = peptides['Proteins'].unique()
protein_features, genes_unprocessed = uniprot_features(genes, resource_folder=f'{resource_folder}')

# collect PFAM domains
pfam_domain_info = pfam_domains(proteins, resource_folder=resource_folder)
pfam_domain_info.to_csv(f'{output_folder}whole_protein_pfam_domains.csv')

residue_features = []
region_features = []
pfam_features = []
for gene, position in peptides[['Proteins', 'cys_in_full']].values:
    if gene in genes_unprocessed:
        continue
    features = protein_features[gene]
    if 'position' in features.columns.tolist():
        positions = features[~features['position'].isnull()]
        positions = positions[positions['position'] == float(position)]
        residue_features.append(positions[['feature_type', 'feature_description']].astype(str).agg('_'.join, axis=1).tolist())
    else:
        residue_features.append(['None'])

    if 'begin' in features.columns.tolist():
        ranges = features[~features['begin'].isnull()]
        ranges = ranges[(ranges['begin'] <= position) & (ranges['end'] >= position)]
        region_features.append(ranges[['feature_type', 'feature_description']].astype(str).agg('_'.join, axis=1).tolist())
    
    else:
        region_features.append(['None'])

    if gene in pfam_domain_info['Proteins'].unique():
        domains = pfam_domain_info[(pfam_domain_info['Proteins'] == gene) &(pfam_domain_info['pfam_start'].astype(int) <= position) & (pfam_domain_info['pfam_end'].astype(int) >= position)]
        pfam_features.append(domains['pfam_acc'].tolist())
    else:
        pfam_features.append([np.nan])

# handle residue-specific features
residues = peptides[~peptides['Proteins'].isin(genes_unprocessed)].copy()
residues['residue_features'] = residue_features
residues = residues.explode('residue_features')
residues[['feature_type', 'feature_description']] = residues['residue_features'].str.split('_', expand=True)
residues['feature_type'] = residues['feature_type'].fillna('None')

# handle domain features
domains = peptides[~peptides['Proteins'].isin(genes_unprocessed)].copy()
domains['domain_features'] = region_features
domains = domains.explode('domain_features')
domains[['feature_type', 'feature_description']] = domains['domain_features'].str.split('_', expand=True)
domains['feature_type'] = domains['feature_type'].fillna('None')

# handle PFAM domains
pfam = peptides[~peptides['Proteins'].isin(genes_unprocessed)].copy()
pfam['pfam_domains'] = pfam_features
pfam['pfam_in_domain'] = [1 if len(domain) > 0 else 0 for domain in pfam['pfam_domains']]

# Save to excel
FileHandling.df_to_excel(output_path=f'{output_folder}uniprot_feature_summary.xlsx', data_frames=[residues, domains, pfam, peptides], sheetnames=['residues', 'domains', 'pfam', 'peptides'])


# ---------------------------------------------Chaperone proteins---------------------------------------------

# Read in background data to generate list of protein sequences to map
raw_data = pd.read_excel(f'{background_path}', sheet_name='raw')
raw_data = raw_data.drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in str(col)], axis=1)
peptides = raw_data[['Proteins', 'Sequence']].copy().drop_duplicates()

# Collect protein keywords
genes = peptides['Proteins'].unique()
protein_keywords = uniprot_function(genes, resource_folder)
peptides['keywords'] = peptides['Proteins'].map(protein_keywords)
peptides['keywords'] = [['_'.join([k, v]) for k, v in words.items()] if type(words) == dict else np.nan for words in peptides['keywords']]
peptides = peptides.explode('keywords')
chaperones = peptides[peptides['keywords'] == 'KW-0143_Chaperone']['Proteins'].drop_duplicates()

chaperone_details = uniprot_summary(tax_id='10090', resource_folder=resource_folder, genes=chaperones, reviewed=True)

FileHandling.df_to_excel(
    output_path=f'{output_folder}chaperone_proteins.xlsx',
    sheetnames=['chaperone_proteins', 'protein_keywords', 'chaperone_details'],
    data_frames=[chaperones, peptides, chaperone_details]
)