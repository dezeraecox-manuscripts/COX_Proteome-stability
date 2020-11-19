import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.optimize import curve_fit
from GEN_Utils import FileHandling
from loguru import logger

input_path_TPE = 'results/recombinant_denaturation/kinetic_fitting/'
input_path_TRP = 'results/recombinant_denaturation/initial_cleanup/'
output_folder = 'results/recombinant_denaturation/normalisation/'

urea_conc = list(np.arange(0, 6.5, 0.5))
urea_conc.pop(-2)

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# read in TPE data - in this case, normalise gradients to min/max
file_list = [filename for filename in os.listdir(input_path_TPE) if 'gradient_summary.xlsx' in filename]

tpe_dict = {}
for filename in file_list:
    sample_name = filename.split('_')[0]
    raw_data = pd.read_excel(f'{input_path_TPE}{filename}', sheet_name='TPE')
    raw_data.drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
    raw_data.set_index('urea', inplace=True)
    raw_data = raw_data - raw_data.loc[0.0]
    raw_data = raw_data / np.max(np.max(raw_data))
    tpe_dict[sample_name] =raw_data

# read in TRP data - in this case, take t0 as endpoint data and normalise to min/max
file_list = [filename for filename in os.listdir(input_path_TRP) if '.xlsx' in filename]

trp_dict = {}
for filename in file_list:
    sample_name = filename.split('_')[0]
    raw_data = pd.read_excel(f'{input_path_TRP}{filename}',sheet_name='TRP')
    raw_data.drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
    # collect t0 data and pivot
    endpoint_data = raw_data[['Well\nRow', 'Well\nCol', '0']].copy()
    endpoint_data = pd.pivot_table(endpoint_data, values='0', index='Well\nCol', columns='Well\nRow').rename(columns={'A': 'TRP_1', 'B': 'TRP_2','C': 'TRP_3', 'G': 'TRP_control'}).reset_index(drop=True)
    endpoint_data.columns = endpoint_data.columns.tolist()
    endpoint_data['urea'] = urea_conc
    endpoint_data.set_index('urea', inplace=True)
    # min, mac normalisation
    endpoint_data = endpoint_data - endpoint_data.loc[0.0]
    endpoint_data = endpoint_data / np.max(np.max(endpoint_data))
    trp_dict[sample_name] = endpoint_data

# Combine dictionaries into single output
dataframes = list(tpe_dict.values()) + list(trp_dict.values())
sheetnames = [f'TPE_{key}' for key in list(tpe_dict.keys())] + [f'TRP_{key}' for key in list(tpe_dict.keys())]

FileHandling.df_to_excel(
    data_frames=dataframes,
    sheetnames=sheetnames,
    output_path=f'{output_folder}normalised_summary.xlsx')
