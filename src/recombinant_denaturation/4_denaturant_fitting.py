import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.optimize import curve_fit
from GEN_Utils import FileHandling
from loguru import logger

logger.info('Import ok')

input_path = 'results/recombinant_denaturation/normalisation/normalised_summary.xlsx'
output_folder = 'results/recombinant_denaturation/denaturant_fitting/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)

# Define the fitting elements as per paper example
def denaturant_curve(urea, yf, mf, yu, mu, deltaG, m):
    temp = 298.15
    R_const = 8.31446261815324
    return (
        (yf + mf * urea)
        + (
            (yu + mu * urea)
            * np.exp(-((deltaG - m * urea) / (R_const * temp)))
        )
    ) / (1 + np.exp(-(deltaG - m * urea) / (R_const * temp)))

def linear(x, m, b):
    return m*x + b

# Read in data to dictionary
raw_data = pd.read_excel(f'{input_path}', sheet_name=None)
raw_data.keys()

sample_params = {}
for key, value in raw_data.items():
    test_data = value.copy()
    x_vals = np.array(test_data['urea'])[:-1]
    test_data.set_index(['urea'], inplace=True)

    # establish dictionary for this sample
    fit_params = {}
    cols = [col for col in test_data.columns.tolist() if 'control' not in col]
    for col in cols:
        y_vals = np.array(test_data[col])[:-1]

        # Attempt fitting and plot result
        logger.info(f'Attempting fit for {col}')
        try:
            popt, pcov = curve_fit(denaturant_curve, x_vals, y_vals, p0=[
                0, 0.1, 2, -0.2, -1, -1])
            fit_params[col] = popt
            x = np.linspace(0, np.max(x_vals), 50)
            y = denaturant_curve(x, *popt)       
            logger.info(f'Fit found for {col}.')

        except:
            logger.info(f'No fit found for {col}.')

    sample_params[key] = fit_params

# Save parameters to excel:
fitted_dfs = {}
for key, value in sample_params.items():
    fitted_df = pd.DataFrame.from_dict(value, orient='columns')
    fitted_df['parameter'] = ['yf', 'mf', 'yu', 'mu', 'deltaG', 'm']
    fitted_dfs[key] = fitted_df

FileHandling.df_to_excel(data_frames=list(fitted_dfs.values()), sheetnames=list(fitted_dfs.keys()), output_path=f'{output_folder}denat_fitting_params.xlsx')

# Calculate Cm values
# Read in data to dictionary
cM_vals = []
for timepoint, fit_df in fitted_dfs.items():
    fit_data = fit_df.copy().set_index('parameter').T
    fit_data['cM'] = fit_data['deltaG'] / fit_data['m']
    fit_data['time'] = timepoint.split('_')[1]
    fit_data['label_type'] = [label.split('_')[0] for label in fit_data.index.tolist()]
    cM_vals.append(fit_data.reset_index(drop=True))

cM_summary = pd.concat(cM_vals)

FileHandling.df_to_excel(data_frames=[cM_summary], sheetnames=['summary'], output_path=f'{output_folder}cM_summary.xlsx')