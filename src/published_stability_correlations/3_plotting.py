
import os, functools
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
import statsmodels as sm
import matplotlib.pyplot as plt
from sklearn.preprocessing import minmax_scale
from matplotlib import colors
from matplotlib import cm
import matplotlib

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

input_path = f'results/published_stability_correlations/correlation/'
output_folder = 'results/published_stability_correlations/plot_correlation/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

font = {'family' : 'arial',
'weight' : 'normal',
'size'   : 14 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

# define plotting functions
def linear_regression(x, y, data, title=None):
    """
    Plots scatter plot with linear regression and 95% ci's
    Returns prediction dataframe

    Inputs:
    x, y: str column names within dataframe
    data: dataframe
    title (optional): str
    """

    y_data = data[y]
    x_data = data[x]
    X = sm.add_constant(x_data)
    res = sm.OLS(y_data, X).fit()

    st, data, ss2 = summary_table(res, alpha=0.05)
    preds = pd.DataFrame.from_records(data, columns=[s.replace('\n', ' ') for s in ss2])
    preds['x_data'] = list(x_data)
    preds['y_data'] = list(y_data)
    preds = preds.sort_values(by='x_data')

    fig, ax = plt.subplots()
    plt.scatter(x_data, y_data)
    plt.plot(preds['x_data'], preds['Predicted Value'], c='red', alpha=0.5, linestyle='--')
    ax.fill_between(preds['x_data'], preds['Mean ci 95% low'], preds['Mean ci 95% upp'], color='red', alpha=0.2)

    plt.ylabel(y)
    plt.xlabel(x)
    if title:
        plt.title(title)
    
    return preds


def hexbin(x, y, **kwargs):
    plt.hexbin(x, y, gridsize=20, linewidths=0, **kwargs)


def corrfunc(x, y, r_xy=(.1, .9), p_xy=(.5, .9), **kws):
    data = pd.DataFrame()
    data['y'] = y
    data['x'] = x

    data.dropna(inplace=True)

    (r, p) = pearsonr(data['x'], data['y'])
    ax = plt.gca()
    ax.annotate(f'r = {r:.2f}',
                xy=r_xy, xycoords=ax.transAxes)
    ax.annotate(f'p = {p:.3f}',
                xy=p_xy, xycoords=ax.transAxes)


def scatter_plotter(x_col, y_col, ax, data, color_df):
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)


def hex_plotter(x_col, y_col, ax, data, color_df):
    df = data[[x_col, y_col]].copy().dropna()
    pval = correlation[(correlation['dataset_1'] == x_col) & (correlation['dataset_2'] == y_col)]['spearmans_pval'].tolist()[0]
    color = colors.LinearSegmentedColormap.from_list('new_map', colors=['#ffffff', '#000000'])

    # ax.hexbin(df[x_col], df[y_col], gridsize=20, linewidths=0, cmap=colors.ListedColormap([(1, 1, 1)] + sns.color_palette(color, 512)[256:]))
    ax.hexbin(df[x_col], df[y_col], gridsize=10, linewidths=0, cmap=color)

    corr_color = correlation[(correlation['dataset_1'] == x_col) & (correlation['dataset_2'] == y_col)]['corr_color'].tolist()[0]
    sns.regplot(x_col, y_col, data, scatter_kws={'alpha': 0}, line_kws={'color': '#f00004', 'linewidth': 2}, ax=ax)

    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)


def kde_plotter(col, ax, data, color_df):
    sns.distplot(data[col], color='black', ax=ax, kde=True, hist=False)


def remove_labels(ax, i, j, xlabel, ylabel):
    keep_x, keep_y = (None, None)
    if j == 0:
        keep_y=True
    if i == (len(datasets)-1):
        keep_x=True
    if keep_y is None:
        keep_y=False
    if keep_x is None:
        keep_x=False

    if keep_x == False:
        ax.set_xlabel('')
    else:
        ax.set_xlabel(xlabel)

    if keep_y == False:
        ax.set_ylabel('')
    else:
        ax.set_ylabel(ylabel)

    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_yticklabels([])
    ax.set_yticks([])


def color_picker(val, val_min=-1, val_max=1):
    norm = colors.Normalize(vmin=val_min, vmax=val_max)
    rgb = cmap(norm(val))[:3] # will return rgba, we take only first 3 so we get rgb
    return colors.rgb2hex(rgb)


# --------------------import datasets--------------------

compiled_datasets = pd.read_excel(f'{input_path}compiled_dfs.xlsx', sheet_name=None)
dataset_details = compiled_datasets.pop('dataset_details')
dataset_details.drop([col for col in dataset_details.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
compiled_datasets.update({dataset_id: df.drop([col for col in df.columns.tolist() if 'Unnamed: ' in col], axis=1) for dataset_id, df in compiled_datasets.items()})

# define list of datasets of interest (not practical to include all of the Jarzab datasets)
# currently sorted by technique, species, type (cells/lysate)
datasets = ['Leuenberger_2017', 'Ogburn_2017A', 'Ogburn_2017B', 'Walker_2019', 'Roberts_2016A', 'Roberts_2016B', 'Jarzab_2020F', 'Jarzab_2020G', 'Jarzab_2020H', 'Becher_2016', 'Franken_2015', 'Miettinen_2018', 'Savitski_2018A', 'Savitski_2018B', 'Ball_2020', 'Jarzab_2020N', 'Jarzab_2020O', 'Savitski_2014A', 'Sridharan_2019', 'Jarzab_2020M']

# -------------plot custom hexbin correlation-------------
positions = dict(zip(datasets, np.arange(1, len(datasets)+1)))
# Generate custom cmap, add colour column for mapping correlation to cmap
cmap =  sns.diverging_palette(200, 340, s=100, l=20, n=9, center="light", as_cmap=True)

## Prepare scaled version of the dataframe, map colours to correlation
compiled_ko = compiled_datasets['compiled_ko'].set_index('KO')[datasets].copy()
scaled_ko = pd.DataFrame(minmax_scale(compiled_ko), columns=[col for col in compiled_ko], index=compiled_ko.index).reset_index()

correlation = pd.read_excel(f'{input_path}correlations.xlsx', sheet_name='correlation_KO')
correlation = correlation[(correlation['dataset_1'].isin(datasets)) & (correlation['dataset_2'].isin(datasets))]
correlation['corr_color'] = correlation['spearmans_r'].map(color_picker)

fig, axes = plt.subplots(len(datasets), len(datasets), figsize=(18, 18))
for i, dataset_1 in enumerate(datasets):
    for j, dataset_2 in enumerate(datasets):
        ax = axes[i][j]
        if i == j:
            # kde plot
            # kde_plotter(dataset_1, ax, data=scaled_ko, color_df=correlation)
            pass
        if i < j:
            # hex_plotter(x_col=dataset_2, y_col=dataset_1, ax=ax, data=scaled_ko)
            corr_color = correlation[(correlation['dataset_1'] == dataset_1) & (correlation['dataset_2'] == dataset_2)]['corr_color'].tolist()[0]
            pval = correlation[(correlation['dataset_1'] == dataset_1) & (correlation['dataset_2'] == dataset_2)]['spearmans_pval'].tolist()[0]
            ax.set_facecolor(corr_color)
            if pval < 0.01:
                ax.annotate('**', (0.25, 0.2), fontsize=40)
            elif pval < 0.05:
                ax.annotate('*', (0.35, 0.2), fontsize=40)
            else:
                pass
        if i > j:
            # hex plot
            hex_plotter(x_col=dataset_2, y_col=dataset_1, ax=ax, data=scaled_ko, color_df=correlation)
            pass

        # axes fix labels
        remove_labels(ax, i, j, f'{positions[dataset_2]}', f'{positions[dataset_1]}')

plt.tight_layout(pad=-0.3)
plt.savefig(f'{output_folder}custom_correlation_heatmap.png')
plt.savefig(f'{output_folder}custom_correlation_heatmap.svg')

# create custom colorbar
fig, ax = plt.subplots(figsize=(6, 1))
fig.subplots_adjust(bottom=0.5)
norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)
cb1 = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='horizontal')
cb1.set_label("Spearman's R")
plt.savefig(f'{output_folder}custom_colourbar_heatmap.svg')