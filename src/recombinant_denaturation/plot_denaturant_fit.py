import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

from lmfit import Model
from scipy.optimize import curve_fit
from GEN_Utils import FileHandling
from loguru import logger

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 10 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

logger.info('Import ok')

input_path = 'results/recombinant_denaturation/normalisation/normalised_summary.xlsx'
output_folder = 'results/recombinant_denaturation/plot_denaturant_fit/'

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


def fitting_function(x_vals, y_vals, denaturant_curve):

    df_model = Model(denaturant_curve)
    # logger.info(df_model.param_names, df_model.independent_vars)
    params = df_model.make_params(yf=0, mf=0.2, yu=10, mu=-5, deltaG=-1, m=-1)
    result = df_model.fit(y_vals, params, urea=x_vals)
    logger.info(result.fit_report())

    return result


def fit_plotter(result, x_vals, y_vals, summary_plot, colour, label):
    # first collect relevant parameters
    fit_params = [result.params[parameter].value for parameter in list(
        result.params.keys())]
    x = np.linspace(0, np.max(summary_plot.index.tolist()), 50)
    y_pre = linear(
        x, result.params['mf'].value, result.params['yf'].value)
    y_post = linear(
        x, result.params['mu'].value, result.params['yu'].value)
    y = result.eval(urea=x)

    # Add plot with each component (pre, post, transition)
    plt.plot(summary_plot.index, summary_plot[(
        'fluorescence', 'mean')], 'o', color=colour)
    plt.errorbar(summary_plot.index, summary_plot[('fluorescence', 'mean')], yerr=summary_plot[(
        'fluorescence', 'std')], fmt='none', ecolor=colour, capsize=5)
    plt.plot(x, y, color=colour, marker=None, label=label)
    plt.plot(x, y_pre, marker=None, linestyle='--', c=colour, alpha=0.5)
    plt.plot(x, y_post, marker=None, linestyle='--',  c=colour, alpha=0.5)

    return plt.gcf()


# Read in data to dictionary
raw_data = pd.read_excel(input_path, sheet_name=None)

# establish dictionary to collect sample info
fit_params = {}
model_dict = {}
for key, value in raw_data.items():
    test_data = value.copy()
    cols = [col for col in test_data.set_index(
        'urea').columns.tolist() if 'control' not in col]
    x_vals = np.array([test_data['urea'][:-1] for col in cols]).flatten()
    y_vals = np.array([test_data[col][:-1] for col in cols]).flatten()

    for_plotting = test_data.melt(
        id_vars='urea', value_vars=cols, var_name='replicate', value_name='fluorescence')
    for_plotting[['urea', 'fluorescence']
                 ] = for_plotting[['urea', 'fluorescence']].astype('float')
    summary_plot = for_plotting.groupby('urea').agg(['mean', 'std'])
    if 'TRP' in key:
        colour = '#23b023'
        label = 'TRP'
    elif 'TPE' in key:
        colour = '#c90265'
        label = 'TPE'

    fit_params[key] = [x_vals, y_vals, summary_plot, colour, label]
    logger.info(f'Attempting fit for {key}')
    try:
        model_dict[key] = fitting_function(x_vals, y_vals, denaturant_curve)
    except:
        logger.info(f'No fit found for {key}.')
        
# Generate plot for 4 h and 24 h separately
for time in [4, 24]:
    samples = [key for key in raw_data.keys() if f'_{time}h' in key]

    fig_mixed, ax = plt.subplots(figsize=(6, 5))
    for key in samples:
        fit_plotter(model_dict[key], *fit_params[key])
    plt.legend()
    plt.ylabel('Fraction unfolded')
    plt.xlabel('Urea (M)')
    plt.ylim(-0.05, 1.05)
    # plt.autoscale()
    plt.tight_layout()
    plt.savefig(f'{output_folder}{time}h_mixedfit.png')
    plt.savefig(f'{output_folder}{time}h_mixedfit.svg')
    plt.show()
