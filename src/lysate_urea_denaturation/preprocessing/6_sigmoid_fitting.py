import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from lmfit import Model, Parameters, minimize, report_fit
from scipy.optimize import curve_fit
from scipy import stats

from utilities.statistical_tests import r_squared_calculator

from GEN_Utils import FileHandling
from loguru import logger

logger.info('Import ok')


# Define the fitting functions
def sigmoid(x, bottom, top, X50):
    return bottom + ((top - bottom) / (1 + np.exp((X50 - x))))


def boltzmann(x, bottom, top, V50, slope):
    return bottom + ((top - bottom) / (1 + np.exp((V50 - x)/slope)))


def denaturant(urea, top, bottom, cM, m):
    # adapted from https://en.wikipedia.org/wiki/Equilibrium_unfolding, by keeping terms for bottom as in Savitski, subbing deltaG into standard equation, and reintroducing bottom term as per boltzmann

    temp_constant = 298.15
    gas_constant = 8.31446261815324
    constant = temp_constant * gas_constant

    y = bottom + ((top - bottom) / (1 + np.exp((m*(cM-urea)/constant))))

    # deltaG can then be calculated as m(cM-urea) - generally calculated at 0M urea therefore m(cM)

    return y


def denaturant_fit(compiled, info_cols, quant_cols):
    """Attempts to fit a sigmoid to each row. Returns sample_params dict where keys are sequences"""

    fit_params = {}

    for info, quant_data in compiled.set_index(info_cols).iterrows():

        # extract x and y vals for fitting
        y_vals = np.array(list(quant_data[quant_cols]))
        x_vals = np.array([float(x) for x in quant_cols])

        # Attempt fitting
        try:
            model = Model(denaturant)
            params = model.make_params(
                bottom=-1, top=1, cM=3, m=-10000)
            result = model.fit(y_vals, params, urea=x_vals)
            r_squared = r_squared_calculator(
                x_vals, y_vals, denaturant, result.values.values())

            # Collect fitted parameters
            fit_stats = pd.DataFrame()
            for parameter, details in result.params.items():
                fit_stats[f'{parameter}_value'] = [details.value]
                fit_stats[f'{parameter}_stderr'] = [details.stderr]
                fit_stats[f'{parameter}_relerr'] = fit_stats[f'{parameter}_stderr'].values[0] / \
                    fit_stats[f'{parameter}_value'].values[0] * 100

            # add r-squared value, key info
            fit_stats['r_squared'] = r_squared
            fit_stats['key'] = [info]

            fit_params[info] = fit_stats
        except:
            logger.info(f'No fit found for {info}')
    return fit_params


def sigmoid_filter(summary, filter_R2=True, filter_range=True, filter_cM=True, filter_relerr=True, filter_direction=True):
    
    # apply filtering criteria
    filtered = summary.copy()

    if filter_R2:
        # Remove R2 < filter threshold
        filtered['filter_R2'] = [1 if R2 > 0.75 else 0 for R2 in filtered['r_squared']]
        logger.info(f"R2 filter: {filtered['filter_R2'].sum()}")

    if filter_range:
        # Remove top/bottom outside range - threshold = 10?
        filtered = filtered[(abs(filtered['top_value']) < 10) & (abs(filtered['bottom_value']) < 10)]
        filtered['filter_range'] = [1 if (abs(val_1) < 10) & (abs(val_2) < 10) else 0 for val_1, val_2 in filtered[['top_value', 'bottom_value']].values]
        logger.info(f"Range filter: {filtered['filter_range'].sum()}")

    if filter_cM:
        # Remove cM outside range tested
        filtered['filter_cM'] = [1 if (val < 6) & (val > 0) else 0 for val in filtered['cM_value']]
        logger.info(f"cM filter: {filtered['filter_cM'].sum()}")

    if filter_relerr:
        # Remove fits with > 50% uncertainty in cM fit
        filtered['filter_relerr'] = [1 if val < 50 else 0 for val in filtered['cM_relerr']]
        logger.info(f"Relative cM error: {filtered['filter_relerr'].sum()}")

    if filter_direction:
        # Remove sigmoids that trend upward
        filtered['filter_direction'] = [1 if val_0 > val_6 else 0 for val_0, val_6 in zip(filtered['0M_value'], filtered['6M_value'])]
        logger.info(f"Sigmoid direction: {filtered['filter_direction'].sum()}")

    filter_cols = [col for col in filtered.columns.tolist() if 'filter_' in str(col)]
    filtered['filter_count'] = filtered[filter_cols].sum(axis=1)
    filtered['filter_all'] = [1 if num == len(filter_cols) else 0 for num in filtered['filter_count']]
    logger.info(f"All filters: {filtered['filter_all'].sum()}")

    # add filtering info to original df
    summary['filtered'] = filtered['filter_all']

    return summary, filtered


if __name__ == '__main__':

    filter_cols = []

    input_path = f'results/lysate_denaturation/clustering/clustered.xlsx'
    output_folder = f'results/lysate_denaturation/sigmoid_fitting/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Read in cluster data
    clusters_summary = pd.read_excel(input_path, sheet_name=None)
    cluster_number = clusters_summary['summary']['cluster'].max()
    clusters = clusters_summary['clustered'].copy()
    clusters.drop([col for col in clusters.columns.tolist() if 'Unnamed: ' in str(col)], axis=1, inplace=True)

    info_cols = ['Sequence', 'Proteins', 'PC1', 'PC2', f'score_{cluster_number}', f'member_{cluster_number}']
    quant_cols = [col for col in clusters.columns.tolist() if type(col) == float]
    clusters = clusters[info_cols+quant_cols].rename(columns={f'member_{cluster_number}': 'cluster', f'score_{cluster_number}': 'score'})
    info_cols = ['Sequence', 'Proteins', 'PC1', 'PC2', 'cluster', 'score']

    # complete denaturant fit
    fit_params = denaturant_fit(clusters, info_cols=info_cols, quant_cols=quant_cols)
    fitting_parameters = pd.concat(fit_params.values()).reset_index(drop=True)

    # add back useful info
    fitting_parameters[info_cols] = pd.DataFrame(fitting_parameters['key'].tolist(), index=fitting_parameters.index)
    summary = pd.merge(clusters, fitting_parameters, on=info_cols, how='inner')

    # generate "fitted" results
    sigmoid_fitted_vals = {}
    for sequence, df in summary.iterrows():
        # generate fitted values
        (bottom, top, cM, m, r_squared, cluster, protein, sequence) = tuple(df[['bottom_value', 'top_value', 'cM_value', 'm_value', 'r_squared', 'cluster', 'Proteins', 'Sequence']])
        y_vals = denaturant(np.array(quant_cols), top, bottom, cM, m)
        sigmoid_fitted_vals[sequence] = y_vals
    sigmoid_fitted_vals = pd.DataFrame(sigmoid_fitted_vals).T
    sigmoid_fitted_vals.columns = quant_cols

    # Add upper and lower bound (to determine up vs down sigmoid)
    lower_bound = []
    upper_bound = []
    for sequence, df in summary.iterrows():
        # generate fitted values
        (bottom, top, cM, m, r_squared, cluster, protein, sequence) = tuple(df[['bottom_value', 'top_value', 'cM_value', 'm_value', 'r_squared', 'cluster', 'Proteins', 'Sequence']])
        lower_bound.append(denaturant(0, top, bottom, cM, m))
        upper_bound.append(denaturant(6, top, bottom, cM, m))
    summary['0M_value'] = lower_bound
    summary['6M_value'] = upper_bound

    summary, filtered = sigmoid_filter(summary, filter_R2=True, filter_range=True, filter_cM=True, filter_relerr=True, filter_direction=True)

    # save results to sigmoids excel
    FileHandling.df_to_excel(output_path=f'{output_folder}sigmoid_fits.xlsx', data_frames=[summary, sigmoid_fitted_vals, filtered], sheetnames=['summary', 'sigmoid_fitted_vals', 'filtering_steps'])

    # add results to the cluster summary excel doc
    clusters_summary['sigmoid_summary'] = summary
    clusters_summary['sigmoid_fitted_vals'] = sigmoid_fitted_vals
    clusters_summary['sigmoid_filtering_steps'] = filtered

    clusters_summary['summary'] = pd.merge(clusters_summary['summary'], summary[['Sequence', 'Proteins', 'filtered']].rename(columns={'filtered': 'sigmoid_filtered'}), on=['Sequence', 'Proteins'], how='outer')

    FileHandling.df_to_excel(output_path=input_path, data_frames=list(clusters_summary.values()), sheetnames=list(clusters_summary.keys()))