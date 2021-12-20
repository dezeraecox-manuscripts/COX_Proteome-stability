import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 12 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

input_path = 'results/inhibitor_urea_denaturation/detect_outliers/outlier_summary.xlsx'
output_folder = 'results/inhibitor_urea_denaturation/plot_proteins/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# ----------------------Prepare raw data----------------------
compiled = pd.read_excel(input_path, sheet_name=None)
compiled.update({key: df.drop([col for col in df.columns.tolist(
) if 'Unnamed: ' in str(col)], axis=1) for key, df in compiled.items()})

# collect outliers
outliers = compiled['outliers'].copy()

# Collect all peptides for each outlier protein
all_peps = compiled['summary'].copy()
all_peps = all_peps[all_peps['count'] >= 7].copy()
all_peps['noncys_ratio_std'] = all_peps['noncys_ratio_std'].replace(
    0, np.nan) # remove single-noncys-peptide proteins
all_peps.dropna(subset=['noncys_ratio_std'], inplace=True)

for protein in outliers['Proteins'].unique().tolist():
    protein_df = all_peps[all_peps['Proteins'] == protein].copy()
    protein_df = protein_df[protein_df['Sequence'].str.contains('C')].copy()
    palette = {
        seq: 'black' if 'C' in seq else 'darkgrey' for seq in protein_df['Sequence'].tolist()}
    dashes = {
        seq: '' if (outlier == 1) else (2, 2) for seq, outlier in protein_df[['Sequence', 'outlier']].values}
    # dashes = {1: '', 0: (2, 2)}

    for_plotting = pd.melt(
        protein_df,
        id_vars=['Sequence', 'outlier'],
        value_vars=[col for col in protein_df if type(col) != str],
        var_name='Urea Conc. (M)',
        value_name='Corrected cys ratio'
        ).dropna()

    thresh_band = protein_df['noncys_ratio_std'].tolist()[0]

    fig, ax = plt.subplots()
    sns.lineplot(
        data=for_plotting,
        x='Urea Conc. (M)',
        y='Corrected cys ratio',
        hue='Sequence',
        palette=palette,
        style='Sequence',
        dashes=dashes,
        marker='o'
    )
    plt.fill_between(
        x=np.arange(0, 7),
        y1=-thresh_band,
        y2=thresh_band,
        color='lightgrey',
        alpha=0.5)
    plt.title(f'{protein}')
    plt.ylim(-1.25, 1.25)
    plt.legend(bbox_to_anchor=(1.0, 1.0))
    plt.savefig(f'{output_folder}{protein}.svg')
    plt.savefig(f'{output_folder}{protein}.png')
    plt.show()
    plt.clf()
