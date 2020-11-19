import os, re, string
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import bokeh.palettes as bp

from loguru import logger
from GEN_Utils import FileHandling
from ProteomicsUtils import StatUtils

logger.info('Import OK')

input_folder = 'results/recombinant_denaturation/initial_cleanup/'
output_path = 'results/recombinant_denaturation/plot_kinetics/'

if not os.path.exists(output_path):
    os.mkdir(output_path)

## Import and cleanup raw plate reader data for each sample
file_list = [filename for filename in os.listdir(input_folder)]

urea_subset = [0.5, 1.5, 2.0, 2.5, 3.5]

for filename in file_list:
    
    sample_name = filename.split('_')[0]
    cleaned_kinetic = pd.read_excel(f'{input_folder}{filename}', sheet_name=None)
    cleaned_kinetic.update({key: df.drop([col for col in df.columns.tolist() if 'Unnamed: ' in col], axis=1) for key, df in cleaned_kinetic.items()})

    # Generate line plots for original raw fluorescence data
    for sample_type, df in cleaned_kinetic['TPE'].groupby('samples'):

        # melt df for plotting
        plot = df.drop(['Well\nRow', 'Well\nCol'], axis=1)
        plotting = pd.melt(plot, id_vars=['samples', 'urea'], var_name='time', value_name='fluorescence',)
        plotting['time'] = plotting['time'].astype(int)
        if urea_subset:
            plotting = plotting[plotting['urea'].isin(urea_subset)]

        # Generate  plot
        fig, ax = plt.subplots()
        sns.lineplot(x='time', y='fluorescence', hue='urea', data=plotting, ci=None, marker='o', palette='magma_r', markeredgecolor="none", hue_norm=(0, 6))
        # plt.ylim(0, 120000)
        plt.xlabel('Time (min)')
        plt.ylabel('TPE Fluorescence (A.U.)')
        plt.legend(bbox_to_anchor=(1.25, 1.00))
        plt.title(sample_type)
        plt.autoscale()
        plt.tight_layout()
        plt.savefig(f'{output_path}{sample_name}_{sample_type}_raw.png')
        plt.show()