"""Collect specific gene ontologies, and additional background/complex information """
import os
import re
import functools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sa
import statsmodels.formula.api as sfa
from GEN_Utils import FileHandling
from loguru import logger
from scipy import stats

import scikit_posthocs as sp
from utilities.database_map_and_filter import (go_lineage_tracer,go_term_details, 
                                         ontology_wordfinder, uniprot_go_genes, create_uniprot_xref, ortholog_map)

logger.info('Import OK')

clustered_path = f'results/lysate_denaturation/clustering/clustered.xlsx'
background_path = f'results/lysate_denaturation/normalised/normalised_summary.xlsx'
resource_folder = f'resources/bioinformatics_databases/'
output_folder = 'results/lysate_denaturation/gene_ontology_datasets/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def add_homolog_id(datasets):
    for name, data in datasets.items():
        data['homologue_id'] = data['Proteins'].map(swiss_to_homologue_id)
        data.to_csv(f'{output_folder}{name}.csv')
        datasets.update({name: data})
    return datasets


# MIG database for homology
homology_db = pd.read_table(f'{resource_folder}HOM_MouseHumanSequence.txt')
homology_db.dropna(subset=['SWISS_PROT IDs', 'HomoloGene ID'], inplace=True)
swiss_to_homologue_id = dict(zip(homology_db['SWISS_PROT IDs'], homology_db['HomoloGene ID']))

# Dataset 1: Any proteins associated with generic "protein complex" GO term GO:0032991
# collect proteins from ontology 
complex_genes = uniprot_go_genes(tax_id='10090', go_term='GO:0032991', child_terms=True, direct=False, output='list')
complex_genes = pd.DataFrame(complex_genes).rename(columns={0: 'Proteins'})

# Dataset 2: Against proteins associated with specific complexes: proteasome
# find terms associated with proteasome
potential_proteasome_terms = ontology_wordfinder(['proteasome']) # decided on "GO:0000502: proteasome complex"
proteasome_genes = uniprot_go_genes(tax_id='10090', go_term='GO:0000502', child_terms=False, direct=True, output='list')
proteasome_genes = pd.DataFrame(proteasome_genes).rename(columns={0: 'Proteins'})

# Dataset 3: Against proteins associated with specific complexes: ribosome
# find terms associated with proteasome
potential_terms = ontology_wordfinder(['ribosome']) # decided on "GO:0003735 structural constituent of ribosome""
ribosome_genes = uniprot_go_genes(tax_id='10090', go_term='GO:0003735', child_terms=False, direct=True, output='list')
ribosome_genes = pd.DataFrame(ribosome_genes).rename(columns={0: 'Proteins'})

# Dataset 4: Against proteins associated with specific complexes: DNA repair complex
# find terms associated with proteasome
potential_terms = ontology_wordfinder(['DNA repair complex']) # decided on "GO:1990391 DNA repair complex""
dna_genes = uniprot_go_genes(tax_id='10090', go_term='GO:1990391', child_terms=False, direct=True, output='list')
dna_genes = pd.DataFrame(dna_genes).rename(columns={0: 'Proteins'})

# Dataset 5: Against proteins associated with specific complexes: nuclear pore
# find terms associated with proteasome
potential_terms = ontology_wordfinder(['nuclear pore']) # decided on "GO:0005643 Nuclear pore"
pore_genes = uniprot_go_genes(tax_id='10090', go_term='GO:0005643', child_terms=True, direct=False, output='list')
pore_genes = pd.DataFrame(pore_genes).rename(columns={0: 'Proteins'})

# Dataset 6: Against proteins associated with protein folding (chaperone)
# find terms associated with proteasome
potential_terms = ontology_wordfinder(['chaperone']) # decided on "GO:0061077 chaperone-mediated protein folding"
chaperone_genes = uniprot_go_genes(tax_id='10090', go_term='GO:0061077', child_terms=True, direct=False, output='list')
chaperone_genes = pd.DataFrame(chaperone_genes).rename(columns={0: 'Proteins'})

# Add homologue id's to each df
datasets = dict(zip(['complex_genes', 'proteasome_genes', 'ribosome_genes', 'dna_genes', 'pore_genes', 'chaperone_genes'], [complex_genes, proteasome_genes, ribosome_genes, dna_genes, pore_genes, chaperone_genes]))
datasets = add_homolog_id(datasets)

# save summaries to excel
file_list = [filename for filename in os.listdir(output_folder) if '.csv' in filename]
datasets = dict(zip([name.split('.csv')[0] for name in file_list], [pd.read_csv(f'{output_folder}{filename}') for filename in file_list]))
summary_df = functools.reduce(lambda left, right: pd.merge(left, right, on='homologue_id', how='outer'), [df[['Proteins', 'homologue_id']].rename(columns={'Proteins': name}) for name, df in datasets.items()])

FileHandling.df_to_excel(output_path=f'{output_folder}go_term_summary.xlsx', data_frames=[summary_df], sheetnames=['summary'])
