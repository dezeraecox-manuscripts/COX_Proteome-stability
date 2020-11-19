import os, re, string
import pandas as pd
import numpy as np
import functools

from GEN_Utils import FileHandling
from loguru import logger

from utilities.database_map_and_filter import uniprot_summary, create_uniprot_xref

logger.info('Import OK')

resource_folder = 'resources/bioinformatics_databases/'
input_folder = 'raw_data/published_stability_correlations/'
output_folder = f'results/published_stability_correlations/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# -----------------------------------clean-up each resource-----------------------------------

cleaned_stability = {}
dataset_details = {}

dataset_categories = pd.read_excel(f'{input_folder}supp_info_investigation_annotated.xlsx', sheet_name='index')

# -----------------------------------------Jarzab_2020-----------------------------------------
dataset = 'Jarzab_2020'
input_path = f'{input_folder}/41592_2020_801_MOESM4_ESM.xlsx'
sheet = None
identifier_type = 'Protein ID'

# read raw data, collect relevant columns
raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet)
index = pd.read_excel(f'{input_folder}/41592_2020_801_MOESM3_ESM.xlsx', sheet_name='Meltome data set')
index['Dataset ID'] = index['Dataset ID'].str.split(':').str.join('_')
index['species'] = index['Organism'].map({'Homo sapiens': 'human', 'Homo sapiens ': 'human','Mus musculus': 'mouse', 'Mus musculus ': 'mouse'})
index_details = dict(zip(index['Dataset ID'], [tuple(entry) for entry in index[['Cells/lysate', 'Strain/tissue', 'species']].values]))

dfs = []
for dataset_id, df in raw_data.items():
    if dataset_id == 'README':
        continue
    data = df[['Protein ID', 'Melting point [°C]']]
    data['details'] = [index_details[dataset_id]] * len(data)
    data['Dataset ID'] = dataset_id
    dfs.append(data)
dfs = pd.concat(dfs)
mean = dfs.groupby(['details', 'Protein ID']).mean().reset_index()
mean[['Cells/lysate', 'Strain/tissue', 'species']] = pd.DataFrame(mean['details'].tolist(), index=mean.index)

# drop species other than human or mouse
mean = mean.dropna(subset=['species'])

# map identifiers to Swissprot ID's
# all of the datasets use different identifiers - some are UniProt names, some UniProt Acc. and some a combination of both..
# try to capture these, but no guarantees!
# There are also some examples where info is assigned to unreviewed UniProt ID's - these have been left as unmapped
mouse_id = uniprot_summary(tax_id='10090', resource_folder=resource_folder, genes=[], reviewed=True)
mouse_id_map = dict(zip(mouse_id['Entry name'].str.split('_').str[0].tolist(), mouse_id['Entry']))
mouse_id_map.update(dict(zip(mouse_id['Entry name'], mouse_id['Entry'])))
mouse_id_map.update(dict(zip(mouse_id['Entry'], mouse_id['Entry'])))
mouse_genes = mouse_id[['Gene names', 'Entry']].copy()
mouse_genes['Gene names'] = mouse_genes['Gene names'].str.upper().str.split(' ')
mouse_genes = mouse_genes.explode('Gene names')
mouse_id_map.update(dict(zip(mouse_genes['Gene names'].str.split('_').str[0].tolist(), mouse_genes['Entry'])))

human_id = uniprot_summary(tax_id='9606', resource_folder=resource_folder, genes=[], reviewed=True)
human_id_map = dict(zip(human_id['Entry name'].str.split('_').str[0].tolist(), human_id['Entry']))
human_id_map.update(dict(zip(human_id['Entry name'], human_id['Entry'])))
human_id_map.update(dict(zip(human_id['Entry'], human_id['Entry'])))
human_genes = human_id[['Gene names', 'Entry']].copy()
human_genes['Gene names'] = human_genes['Gene names'].str.upper().str.split(' ')
human_genes = human_genes.explode('Gene names')
human_id_map.update(dict(zip(human_genes['Gene names'].str.split('_').str[0].tolist(), human_genes['Entry'])))

mapped = []
for species, df in mean.groupby('species'):
    df['identifier'] = df['Protein ID'].str.split('_').str[0].str.split('-').str[0]
    if species == 'mouse':
        df['Proteins'] = df['identifier'].map(mouse_id_map)
    if species == 'human':
        df['Proteins'] = df['identifier'].map(human_id_map)
    mapped.append(df)
mean = pd.concat(mapped)

# to check proportion of unmapped protein IDs
for group, df in mean.groupby('details'):
    logger.info(f'Unmapped IDs for {group}: {len(df[df["Proteins"].isnull()])} / {len(df)}')


# grab columns of interest for each dataset and add to cleaned dictionary with A, B, ... appended
# add details to details
dataset_labels = list(string.ascii_uppercase)
for i, (details, df) in enumerate(mean.groupby('details')):
    stability = df[['Proteins', 'Melting point [°C]']].copy()
    stability.columns = ['Proteins', f'{dataset}{dataset_labels[i]}_Tm']
    cleaned_stability[f'{dataset}{dataset_labels[i]}'] = stability
    dataset_details[f'{dataset}{dataset_labels[i]}'] = details

# ------------------------------------------Ball_2020------------------------------------------
dataset = 'Ball_2020'
input_path = f'{input_folder}/42003_2020_795_MOESM2_ESM.xlsx'
sheet = 'Staurosporine_TPP_data'
identifier_type = 'UniProtKB-AC'
# read raw dta, collect relevant columns
raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet)
stability = raw_data[['Uniprot.ID', 'Average_Filtered_meltPoint_DMSO_Webb_TPP']].copy()
stability.columns = ['Proteins', f'{dataset}_Tm']
cleaned_stability[f'{dataset}'] = stability

dataset_details[f'{dataset}'] = ('lysate', 'K562', 'human')

# ---------------------------------------Sridharan_2019----------------------------------------
dataset = 'Sridharan_2019'
input_path = f'{input_folder}/41467_2019_9107_MOESM4_ESM.xlsx'
sheet = 'Table S1'
identifier_type = 'UniProtKB-AC'
tax_id = '9606'

# read raw data
raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet)

# Collect relevant columns
raw_data['mean_meltPoint'] = raw_data[['meltPoint_crude_1', 'meltPoint_crude_2']].mean(axis=1)
raw_data['mean_meltPoint'] = raw_data[['meltPoint_crude_1', 'meltPoint_crude_2']].mean(axis=1)
stability = raw_data[['protein_ID', 'mean_meltPoint', ]].copy()
stability.columns = ['Proteins', f'{dataset}_Tm']
cleaned_stability[f'{dataset}'] = stability

dataset_details[f'{dataset}'] = ('lysate', 'Jurkat', 'human')

# -----------------------------------------Walker_2019-----------------------------------------
dataset = 'Walker_2019'
input_path = f'{input_folder}/pnas.1819851116.sd04.xlsx'
sheet = '-TMAO'
identifier_type = 'UniProtKB-AC'
tax_id = '9606'

# read raw data
raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet)

# Collect relevant columns
stability = raw_data[['Protein Names', '[GdmCl]1/2']].copy()
stability.columns = ['Proteins', f'{dataset}_Tm']
cleaned_stability[f'{dataset}'] = stability

dataset_details[f'{dataset}'] = ('lysate', 'HCA2-hTert', 'human')

# ---------------------------------------Miettinen_2018----------------------------------------
dataset = 'Miettinen_2018'
input_path = f'{input_folder}/embj201798359-sup-0003-tableev2.xlsx'
sheet = 'Table EV2'
identifier_type = 'UniProt'
tax_id = '10090'

# read raw data
raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet, skiprows=3)

# map identifiers to Swissprot ID's
raw_data['Protein'] = raw_data['Uniprot'].str.upper().str.split('-').str[0]

# Collect relevant columns
stability = raw_data[['Protein', 'ctrl.Tm']].copy()
stability.columns = ['Proteins', f'{dataset}_Tm']
cleaned_stability[f'{dataset}'] = stability

dataset_details[f'{dataset}'] = ('cells', 'MCF-7', 'human')

# ----------------------------------------Savitski_2018----------------------------------------
dataset = 'Savitski_2018'
input_path = f'{input_folder}Supplementary Dataset 6_reference_melting curves.xlsx'
sheet = 'Jurkat_processed'
identifier_type = 'Accession No.' # This appears to be International Protein Identifier (IPI) --> very difficult to map to UniProt as it seems like all the databases have been depreciated! Attempt to map using protein names?
tax_id = '9606'


# read raw data
raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet, skiprows=1)

# map identifiers to Swissprot ID's
ids = uniprot_summary(tax_id=tax_id, resource_folder=resource_folder, genes=[], reviewed=True)
id_map = dict(zip(ids['Entry name'].str.split('_').str[0].tolist(), ids['Entry']))
ids['Genes'] = ids['Gene names'].str.split(' ')
ids = ids.explode('Genes')
id_map.update(zip(ids['Genes'], ids['Entry']))
raw_data['Protein'] = raw_data['Protein Name'].map(id_map)
logger.info(f'Unmapped IDs for {dataset}A: {len(raw_data[raw_data["Protein"].isnull()])} / {len(raw_data)}')

# Collect relevant columns
stability = raw_data[['Protein', 'median_meltPoint']].copy()
stability.columns = ['Proteins', f'{dataset}A_Tm']
cleaned_stability[f'{dataset}A'] = stability


input_path = f'{input_folder}Supplementary Dataset 6_reference_melting curves.xlsx'
sheet = 'T-cell_processed'
# This appears to be International Protein Identifier (IPI) --> very difficult to map to UniProt as it seems like all the databases have been depreciated! Attempt to map using protein names?
identifier_type = 'Accession No.'
tax_id = '9606'

# read raw data
raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet, skiprows=1)

# map identifiers to Swissprot ID's
raw_data['Protein'] = raw_data['Protein Name'].map(id_map)
logger.info(f'Unmapped IDs for {dataset}B: {len(raw_data[raw_data["Protein"].isnull()])} / {len(raw_data)}')

# Collect relevant columns
stability = raw_data[['Protein', 'median_meltPoint']].copy()
stability.columns = ['Proteins', f'{dataset}B_Tm']
cleaned_stability[f'{dataset}B'] = stability

dataset_details[f'{dataset}A'] = ('cells', 'Jurkat', 'human')
dataset_details[f'{dataset}B'] = ('cells', 'T-cells', 'human')

# -----------------------------------------Ogburn_2017-----------------------------------------
dataset = 'Ogburn_2017'
input_path = f'{input_folder}/pr7b00442_si_003.xlsx'
sheet = 'Sheet1'
identifier_type = 'Protein Number'
tax_id = '9606'

# read raw data
raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet, skiprows=2, usecols=[1, 15])

# map identifiers to Swissprot ID's
raw_data['Protein'] = raw_data['Protein Number'].str.split('-').str[0]
logger.info(f'Unmapped IDs for {dataset}A: {len(raw_data[raw_data["Protein"].isnull()])} / {len(raw_data)}')

# Collect relevant columns
stability = raw_data[['Protein', 'C 1/2']].copy()
stability.columns = ['Proteins', f'{dataset}A_Tm']
cleaned_stability[f'{dataset}A'] = stability

input_path = f'{input_folder}/pr7b00442_si_004.xlsx'
sheet = 'Sheet1'
identifier_type = 'Protein Number'
tax_id = '9606'

# read raw data
raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet, skiprows=2, usecols=[1, 15])

# map identifiers to Swissprot ID's
raw_data['Protein'] = raw_data['Protein Number'].str.split('-').str[0]
logger.info(f'Unmapped IDs for {dataset}B: {len(raw_data[raw_data["Protein"].isnull()])} / {len(raw_data)}')

# Collect relevant columns
stability = raw_data[['Protein', 'C 1/2']].copy()
stability.columns = ['Proteins', f'{dataset}B_Tm']
cleaned_stability[f'{dataset}B'] = stability

dataset_details[f'{dataset}A'] = ('lysate', 'MCF-7', 'human')
dataset_details[f'{dataset}B'] = ('lysate', 'MCF-7', 'human')

# --------------------------------------Leuenberger_2017---------------------------------------
dataset = 'Leuenberger_2017'
input_path = f'{input_folder}/aai7825_Leuenberger_Table-S3.xlsx'
sheet = 'Human HeLa Cells'
identifier_type = 'Protein_ID'
tax_id = '9606'

# read raw data
raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet)

# map identifiers to Swissprot ID's
raw_data['Protein'] = raw_data['Protein_ID'].str.split(';').str[-1]
logger.info(f'Unmapped IDs for {dataset}: {len(raw_data[raw_data["Protein"].isnull()])} / {len(raw_data)}')

# Collect relevant columns
stability = raw_data[['Protein', 'Tm Protein']].copy()
stability.columns = ['Proteins', f'{dataset}_Tm']
cleaned_stability[f'{dataset}'] = stability

dataset_details[f'{dataset}'] = ('lysate', 'HeLa', 'human')

# ----------------------------------------Roberts_2016-----------------------------------------
dataset = 'Roberts_2016'
input_path = f'{input_folder}/pr6b00927_si_003.xlsx'
sheet = None
identifier_type = 'UniProtKB-AC'
tax_id = '10090'

# read raw data
raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet, skiprows=1)

# data is provided per-mouse, with no filtering for Cm indicated, 
# however, within the methods they indicate filtering R2 values > 0.8
# apply this here, leaving peptides with <0.8 as NaN
mouse_dfs = []
for mouse, df in raw_data.items():
    mouse_number = mouse.strip('Mouse') 
    mouse_data = df[['Sequence', 'Protein Accessions', 'C 1/2', 'Adjusted R Squared', 'Isolation Interference']]
    mouse_data['mouse'] = mouse.strip('Mouse')
    mouse_dfs.append(mouse_data)
mouse_dfs = pd.concat(mouse_dfs)
mouse_dfs['C 1/2'] = [val if fit >= 0.8 else np.nan for val, fit in mouse_dfs[['C 1/2', 'Adjusted R Squared']].values]
mouse_dfs['C 1/2'] = [val if val >= 0 else np.nan for val in mouse_dfs['C 1/2']] # there is one aweful outlier with R2 0.95 but C1/2 -4.8e6
mouse_dfs['C 1/2'] = [val if fit <= 30 else np.nan for val, fit in mouse_dfs[['C 1/2', 'Isolation Interference']].values] # second filtering criteria applied in study
mouse_dfs['age'] = ['old' if int(mouse) > 8 else 'young' for mouse in mouse_dfs['mouse']]
mouse_mean = mouse_dfs.groupby(['Sequence', 'Protein Accessions', 'age']).mean().reset_index()

# map identifiers to Swissprot ID's
mouse_mean['Protein'] = mouse_mean['Protein Accessions'].str.split('-').str[0]
logger.info(f'Unmapped IDs for {dataset}: {len(mouse_mean[mouse_mean["Protein"].isnull()])} / {len(mouse_mean)}')

# Collect relevant columns
stability_A = mouse_mean[mouse_mean['age'] == 'young'][['Protein', 'C 1/2']].copy()
stability_A.columns = ['Proteins', f'{dataset}A_Tm']
cleaned_stability[f'{dataset}A'] = stability_A

stability_B = mouse_mean[mouse_mean['age'] == 'old'][['Protein', 'C 1/2']].copy()
stability_B.columns = ['Proteins', f'{dataset}B_Tm']
cleaned_stability[f'{dataset}B'] = stability_B

dataset_details[f'{dataset}A'] = ('lysate', 'brain', 'mouse')
dataset_details[f'{dataset}B'] = ('lysate', 'brain', 'mouse')

# -----------------------------------------Becher_2016-----------------------------------------
dataset = 'Becher_2016'
input_path = f'{input_folder}/41589_2016_BFnchembio2185_MOESM256_ESM.xlsx'
sheet = 'SupplDataset3'
identifier_type = 'IPI'# This appears to be International Protein Identifier (IPI) --> very difficult to map to UniProt as it seems like all the databases have been depreciated! Attempt to map using protein names?
tax_id = '9606'

# read raw data
raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet, skiprows=2)

# map identifiers to Swissprot ID's
ids = uniprot_summary(tax_id=tax_id, resource_folder=resource_folder, genes=[], reviewed=True)
id_map = dict(zip(ids['Entry name'].str.split('_').str[0].tolist(), ids['Entry']))
ids['Genes'] = ids['Gene names'].str.split(' ')
ids = ids.explode('Genes')
id_map.update(zip(ids['Genes'], ids['Entry']))
raw_data['Protein'] = raw_data['clustername'].map(id_map)
logger.info(f'Unmapped IDs for {dataset}: {len(raw_data[raw_data["Protein"].isnull()])} / {len(raw_data)}')

raw_data['average_Tm'] = raw_data[['meltPoint_Vehicle_1', 'meltPoint_Vehicle_2']].mean(axis=1)

# Collect relevant columns
stability = raw_data[['Protein', 'average_Tm']].copy()
stability.columns = ['Proteins', f'{dataset}_Tm']
cleaned_stability[f'{dataset}'] = stability

dataset_details[f'{dataset}'] = ('cells', 'HepG2', 'human')

# ----------------------------------------Franken_2015-----------------------------------------
dataset = 'Franken_2015'
input_path = f'{input_folder}/41596_2015_BFnprot2015101_MOESM411_ESM.xlsx'
sheet = 'Results'
identifier_type = 'Protein Name'
tax_id = '9606'

# read raw data
raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet)

# map identifiers to Swissprot ID's
# appears that their dataset contains protein families (e.g. CDK11A|CDK11B)
# As these are challenging to map without duplicating, drop any proteins that are not unique
id_map = uniprot_summary(tax_id=tax_id, resource_folder=resource_folder, genes=[], reviewed=True)
id_map['Gene names'] = id_map['Gene names'].str.split(' ')
id_map = id_map.explode('Gene names')
id_map = dict(zip(id_map['Gene names'], id_map['Entry']))
raw_data = raw_data[raw_data['Protein_ID'].str.split('|').str.len() == 1]
raw_data['Protein'] = raw_data['Protein_ID'].map(id_map)
logger.info(f'Unmapped IDs for {dataset}: {len(raw_data[raw_data["Protein"].isnull()])} / {len(raw_data)}')

raw_data['average_Tm'] = raw_data[['meltPoint_Vehicle_1', 'meltPoint_Vehicle_2']].mean(axis=1)

# Collect relevant columns
stability = raw_data[['Protein', 'average_Tm']].copy()
stability.columns = ['Proteins', f'{dataset}_Tm']
cleaned_stability[f'{dataset}'] = stability

dataset_details[f'{dataset}'] = ('cells', 'K562', 'human')

# ----------------------------------------Savitski_2014----------------------------------------
dataset = 'Savitski_2014'
input_path = f'{input_folder}/Table_S4_Thermal_Profiling_Staurosporine_cell_extract.xlsx'
sheet = 'staurosporine_in_cell_extract'
identifier_type = 'IPI acc. no.'
tax_id = '9606'

# read raw data
raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet)

# map identifiers to Swissprot ID's
ids = uniprot_summary(tax_id=tax_id, resource_folder=resource_folder, genes=[], reviewed=True)
ids['Genes'] = ids['Gene names'].str.split(' ')
ids = ids.explode('Genes')
id_map = dict(zip(ids['Genes'], ids['Entry']))
id_map.update(zip(ids['Protein names'], ids['Entry']))
raw_data['Protein'] = raw_data['Protein NAME'].map(id_map)
logger.info(f'Unmapped IDs for {dataset}A: {len(raw_data[raw_data["Protein"].isnull()])} / {len(raw_data)}')

raw_data['average_Tm'] = raw_data[['meltP_Vehicle_Expt1', 'meltP_Vehicle_Expt2']].mean(axis=1)

# Collect relevant columns
stability = raw_data[['Protein', 'average_Tm']].copy()
stability.columns = ['Proteins', f'{dataset}A_Tm']
cleaned_stability[f'{dataset}A'] = stability


input_path = f'{input_folder}Table_S3_Thermal_Profiling_ATP_cell_extract.xlsx'
sheet = 'ATP'
# This appears to be International Protein Identifier (IPI) --> very difficult to map to UniProt as it seems like all the databases have been depreciated! Attempt to map using protein names?
identifier_type = 'Accession No.'
tax_id = '9606'

# read raw data
raw_data = pd.read_excel(f'{input_path}', sheet_name=sheet)

# map identifiers to Swissprot ID's
raw_data['Protein'] = raw_data['Protein NAME'].map(id_map)
logger.info(f'Unmapped IDs for {dataset}B: {len(raw_data[raw_data["Protein"].isnull()])} / {len(raw_data)}')

raw_data['average_Tm'] = raw_data[['meltP_Vehicle_1', 'meltP_Vehicle_2']].mean(axis=1)

# Collect relevant columns
stability = raw_data[['Protein', 'average_Tm']].copy()
stability.columns = ['Proteins', f'{dataset}B_Tm']
cleaned_stability[f'{dataset}B'] = stability

dataset_details[f'{dataset}A'] = ('lysate', 'K562', 'human')
dataset_details[f'{dataset}B'] = ('lysate', 'K562', 'human')


# --------------------------------Add homology info---------------------------------

# Generate standard databases for human, mouse and ortholog proteins

human_ortho = create_uniprot_xref(input_path=resource_folder, tax_id='9606', gene_ids=[], id_type='OrthoDB')
mouse_ortho = create_uniprot_xref(input_path=resource_folder, tax_id='10090', gene_ids=[], id_type='OrthoDB')
ortho = dict(zip(human_ortho['UniProtKB-AC'], human_ortho['ID']))
ortho.update(dict(zip(mouse_ortho['UniProtKB-AC'], mouse_ortho['ID'])))

human_keggo = create_uniprot_xref(input_path=resource_folder, tax_id='9606', gene_ids=[], id_type='KO')
mouse_keggo = create_uniprot_xref(input_path=resource_folder, tax_id='10090', gene_ids=[], id_type='KO')
# note some non-unique keys here...
keggo = dict(zip(human_keggo['UniProtKB-AC'], human_keggo['ID']))
keggo.update(dict(zip(mouse_keggo['UniProtKB-AC'], mouse_keggo['ID'])))

homology_db = pd.read_table(f'{resource_folder}HOM_MouseHumanSequence.txt')
homology_db.dropna(subset=['SWISS_PROT IDs', 'HomoloGene ID'], inplace=True)
sp_to_homoid = dict(zip(homology_db['SWISS_PROT IDs'], homology_db['HomoloGene ID']))

# add homology information from orthoDB, KeggO, MIG
for index, df in cleaned_stability.items():
    df['KO'] = df['Proteins'].map(keggo)
    df['OrthoDB'] = df['Proteins'].map(ortho)
    df['MIG_ID'] = df['Proteins'].map(sp_to_homoid)
    cleaned_stability[index].update(df)

# Generate summary details table
dataset_details = pd.DataFrame(dataset_details).T.reset_index()
dataset_details.columns = ['dataset_id', 'sample_type', 'sample_origin', 'sample_species']

cleaned_stability['index'] = dataset_details
cleaned_stability['dataset_categories'] = dataset_categories

# Save cleaned datasets to single excel file
FileHandling.df_to_excel(
    output_path=f'{output_folder}published_stability_datasets.xlsx',
    data_frames=list(cleaned_stability.values()),
    sheetnames=list(cleaned_stability.keys())
    )
