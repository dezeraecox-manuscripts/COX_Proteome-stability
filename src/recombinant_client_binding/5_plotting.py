import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')


font = {'family' : 'Arial',
'weight' : 'normal',
'size'   : 14 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

def main(sample_name, input_path):
    # Read in normalised data
    raw_data = pd.read_excel(input_path, sheet_name='native_norm')
    raw_data.drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

    # Select peptides of interest
    peptides = [
        'GPAVGIDLGTTYSCVGVFQHGK', #cys17
        'CNEIINWLDK', #cys574
        'VCNPIITK', #cys603
        'VSLEEIYSGCTK', #cys179
        'EALCGCTVNVPTLDGR', #cys267/269
        ]

    peptide_data = raw_data[raw_data['Sequence'].isin(peptides)]

    # Melt df for plotting
    peptide_data = pd.melt(peptide_data, id_vars=['Sequence', 'Proteins', 'replicate'], value_vars=['2', '4'], var_name='sample', value_name='ratio')

    # generate plots for individual proteins
    hue_order = ['4', '2']
    order = {
        'P11142': ['GPAVGIDLGTTYSCVGVFQHGK', 'CNEIINWLDK', 'VCNPIITK'],
        'P25685': ['VSLEEIYSGCTK', 'EALCGCTVNVPTLDGR']
    } 

    for protein, df in peptide_data.groupby('Proteins'):

        fig, ax = plt.subplots(figsize=(12, 3))

        sns.swarmplot(data=df, x='Sequence', y='ratio', hue='sample', color='#878787', dodge=True, order=order[protein], hue_order=hue_order, size=10)
        sns.boxplot(data=df.groupby(['Sequence', 'sample']).mean().reset_index(), x='Sequence', y='ratio', hue='sample', color='black', dodge=True, order=order[protein], hue_order=hue_order)
        sns.pointplot(data=df, x='Sequence', y='ratio', hue='sample', color='black', dodge=0.38, order=order[protein], hue_order=hue_order, join=False, errwidth=0.5, capsize=0.1, markers=[' ', ' '])
        ax.axhline(1.0, color='grey', linestyle='--')
        plt.xticks(rotation=20)
        # plt.ylim(0.9, 1.1)
        plt.legend('')
        plt.savefig(f'{output_folder}{sample_name}_{protein}_native_normed.png')
        plt.savefig(f'{output_folder}{sample_name}_{protein}_native_normed.svg')

if __name__ == '__main__':
    
    sample_names = ['Heat', 'Urea']
    output_folder = 'results/recombinant_client_assay/plot_native_norm/'

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    
    for sample_name in sample_names:
        input_path = f'results/recombinant_client_assay/native_normalisation/{sample_name}_native_normed.xlsx'
        main(sample_name, input_path)