import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

from scipy.stats import ttest_ind
from GEN_Utils import FileHandling
from loguru import logger

font = {'family': 'normal',
        'weight': 'normal',
        'size': 10}
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'
logger.info('Import ok')

input_path = 'results/recombinant_denaturation/denaturant_fitting/cM_summary.xlsx'
output_folder = 'results/recombinant_denaturation/plot_Cm/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

def twosamp_test(df_1):
    """Applies ttest to each pair of columns within df"""
    results = {}

    for x, col in enumerate(df_1.columns.tolist()):
        vals_a = list(df_1[col])
        test_cols = df_1.columns.tolist()[x+1:]
        for test in test_cols:
            vals_b = list(df_1[test])
            ttest = ttest_ind(vals_a, vals_b, nan_policy='omit')
            results[(col, test)] = list(ttest)

    result = pd.DataFrame.from_dict(results).T
    result.columns = ['statistic', 'pvalue']

    return result

# Read in raw data
raw_data = pd.read_excel(input_path)
cM_summary = raw_data.copy().drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in col], axis=1)

# Add colour information
colour_dict = {'TRP': '#23b023', 'TPE': '#c90265'}
cM_summary['color'] = cM_summary['label_type'].map(colour_dict)
cM_summary['sample_name'] = cM_summary[['label_type', 'time']].agg('_'.join, axis=1)


# Generate plot for 4 h and 24 h separately
for time, df in cM_summary.groupby('time'):

    fig, ax = plt.subplots(figsize=(3, 4))
    plt.ylim(0, 6)
    plt.xlim(-0.5, 1.5)
    sns.swarmplot(x='label_type', y='cM', palette=colour_dict, data=df, size=10)
    sns.pointplot(data=df, x='label_type', y='cM', color='black', join=False, errwidth=1, capsize=0.075, markers='_')
    plt.ylabel('Fitted cM ([Urea] M)')
    plt.xlabel(None)
    plt.tight_layout()
    plt.savefig(f'{output_folder}{time}_cM.svg')
    plt.savefig(f'{output_folder}{time}_cM.png')
    plt.show()

# Generate plot for 4 h and 24 h together
colour_dict = {'TRP_24h': '#23b023', 'TPE_24h': '#c90265',
               'TRP_4h': '#23b023', 'TPE_4h': '#c90265'}

fig, ax = plt.subplots(figsize=(3, 4))
plt.ylim(0, 6)
plt.xlim(-0.5, 3.5)
sns.swarmplot(x='sample_name', y='cM', palette=colour_dict, data=cM_summary,
              size=10, order=['TPE_4h', 'TRP_4h',  'TPE_24h', 'TRP_24h'])
sns.pointplot(data=cM_summary, x='sample_name', y='cM', color='black',
              join=False, errwidth=1, capsize=0.075, markers='_', order=['TPE_4h', 'TRP_4h',  'TPE_24h', 'TRP_24h'])
plt.ylabel('Fitted cM ([Urea] M)')
plt.xlabel(None)
plt.tight_layout()
plt.savefig(f'{output_folder}summary_cM.svg')
plt.savefig(f'{output_folder}summary_cM.png')
plt.show()

