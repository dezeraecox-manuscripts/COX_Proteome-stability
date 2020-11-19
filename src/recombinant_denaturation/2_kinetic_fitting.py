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
output_folder = 'results/recombinant_denaturation/kinetic_fitting/'

urea_conc = list(np.arange(0, 6.5, 0.5))
urea_conc.pop(-2)
urea_pos = np.arange(1, len(urea_conc)+1, 1)
colours = list(reversed(bp.magma(12)))


if not os.path.exists(output_folder):
    os.mkdir(output_folder)

## Import and cleanup raw plate reader data for each sample
file_list = [filename for filename in os.listdir(input_folder)]

for filename in file_list:
    sample_name = filename.split('_')[0]
    cleaned_kinetic = pd.read_excel(f'{input_folder}{filename}', sheet_name=None)
    cleaned_kinetic.update({key: df.drop([col for col in df.columns.tolist() if 'Unnamed: ' in col], axis=1) for key, df in cleaned_kinetic.items()})

    # Regenerate single df with replicate labels
    clean_data = pd.concat((cleaned_kinetic.values()))
    samples = ['TRP_1', 'TRP_2', 'TRP_3', 'TPE_1', 'TPE_2', 'TPE_3', 'TRP_control', 'TPE_control']
    sample_pos = [x for x in string.ascii_uppercase[0:len(samples)]]
    sample_map = dict(zip(sample_pos, samples))
    clean_data['replicate'] = clean_data['Well\nRow'].map(sample_map)
    clean_data.rename(columns={'Well\nCol':'sample_col', 'Well\nRow':'sample_row'}, inplace=True)

    # For TPE at each concentration of urea, collect first 9 minutes and fit with linear regression
    data_for_fit = {}
    time_range = np.arange(0, 9, 1)
    for replicate, df in clean_data.groupby('replicate'):
        ## Cleaning dataframe for plotting
        test_sample = df.copy().drop(['samples', 'urea', 'replicate'] , axis=1).set_index(['sample_row', 'sample_col']).T
        test_sample.index = [np.arange(0, 60, 1)]
        data_for_fit[replicate] = test_sample.loc[time_range]

    fit_param_dict = {}
    for sample in samples:
        fit_params = pd.DataFrame(index=['gradient', 'intercept', 'r_squared'])
        for x, col in enumerate(data_for_fit[sample].columns.tolist()):
            # Collecting data
            x_data = time_range.astype(float)
            y_data = data_for_fit[sample][col].astype(float).tolist()
            norm_y_data = [value - y_data[0] for value in y_data]

            # Calculating fit parameters
            x_fit, y_fit, x_data_fitted, y_data_fitted, popt, pcov = StatUtils.fit_calculator(x_data, norm_y_data, StatUtils.linear)
            residuals = norm_y_data - StatUtils.linear(x_data, *popt)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((norm_y_data-np.mean(norm_y_data))**2)
            r_squared = 1 - (ss_res / ss_tot)
            fit_params[col] = list(popt) + [r_squared]
        fit_param_dict[sample] = fit_params

    sheetnames = [f'{x}' for x in list(fit_param_dict.keys())]
    dataframes = [df.reset_index() for df in list(fit_param_dict.values())]
    FileHandling.df_to_excel(
        f'{output_folder}{sample_name}_Norm_fit_params.xlsx',
        sheetnames=sheetnames, 
        data_frames=dataframes)

    # ## Collecting gradients for each concentration at each sample
    gradients_dict = {
        sample: fit_param_dict[sample].loc['gradient'].reset_index(drop=True)
        for sample in samples
    }

    gradients = pd.DataFrame.from_dict(gradients_dict, orient='columns')
    gradients['urea'] = urea_conc

    tpe_cols = ['urea'] + [col for col in gradients.columns.tolist() if 'TPE' in col]
    trp_cols = ['urea'] + [col for col in gradients.columns.tolist() if 'TRP' in col]

    data_frames = [gradients[tpe_cols], gradients[trp_cols]]

    FileHandling.df_to_excel(
        output_path=f'{output_folder}{sample_name}_gradient_summary.xlsx',
        data_frames=data_frames,
        sheetnames=['TPE', 'TRP'])

