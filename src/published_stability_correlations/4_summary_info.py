
import os, functools
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib
from sklearn.preprocessing import minmax_scale


from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

input_path = f'results/published_stability_correlations/'
output_folder = 'results/published_stability_correlations/summary_details/'

# Datasets selected for use in comparison
datasets = ['Leuenberger_2017', 'Ogburn_2017A', 'Ogburn_2017B', 'Walker_2019', 'Roberts_2016A', 'Roberts_2016B', 'Jarzab_2020F', 'Jarzab_2020G', 'Jarzab_2020H', 'Becher_2016', 'Franken_2015', 'Miettinen_2018', 'Savitski_2018A', 'Savitski_2018B', 'Ball_2020', 'Jarzab_2020N', 'Jarzab_2020O', 'Savitski_2014A', 'Savitski_2014B', 'Sridharan_2019', 'Jarzab_2020M']


if not os.path.exists(output_folder):
    os.makedirs(output_folder)


# ---------------------------import datasets---------------------------
dataset_details = pd.read_excel(f'{input_path}correlation/dataset_details.xlsx')
dataset_details.drop([col for col in dataset_details.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

article_details = pd.read_excel(f'raw_data/published_stability_correlations/supp_info_investigation_annotated.xlsx')
article_details.drop([col for col in article_details.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# combine dataset details with publication details
details = pd.merge(dataset_details, article_details, on='index_key', how='left')


# Select only datasets of interest
details = details[details['dataset_id'].isin(datasets)]

# Assign dataset number from plot
details['dataset_number'] = details['dataset_id'].map(dict(zip(datasets, np.arange(1, len(datasets)+1))))
details = details.sort_values('dataset_number')

# collect columns of interest
details = details[['dataset_number', 'dataset_id', 'sample_type', 'sample_origin',  'sample_species', 'index_key', 'DOI', 'citation',]].copy().rename(columns={'index_key': 'Reference'})


FileHandling.df_to_excel(
    output_path=f'{output_folder}summary_details.xlsx',
    sheetnames=['summary'],
    data_frames=[details]
)