import os, functools
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import ListedColormap
from sklearn.preprocessing import minmax_scale
from matplotlib import colors


from utilities.statistical_tests import stability_summary_calculator, correlation_df

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

published_stabilities = f'results\lysate_denaturation/published_datasets/dataset_details.xlsx'
measured_stabilitites = f'results/lysate_denaturation/cM_correlation/measured_summary.xlsx'
correlations = f'results/lysate_denaturation/cM_correlation/correlations.xlsx'

output_folder = 'results/lysate_denaturation/plot_correlation/'

cluster_colors = {4: 'royalblue', 2: 'firebrick', 3: 'rebeccapurple', 1: 'darkorange'}

font = {'family' : 'arial',
'weight' : 'normal',
'size'   : 14 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)



def remove_labels(ax, i, j, xlabel, ylabel):
    keep_x, keep_y = (None, None)
    if j == 0:
        keep_y=True
    if i == 4:
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


# define plotting function
def hexbin(x, y, **kwargs):
    plt.hexbin(x, y, gridsize=20, linewidths=0, **kwargs)


def corrfunc(x, y, r_xy=(.1, .9), p_xy=(.5, .9), **kws):
    data = pd.DataFrame()
    data['y'] = y
    data['x'] = x

    data.dropna(inplace=True)

    (r, p) = spearmanr(data['x'], data['y'])
    ax = plt.gca()
    ax.annotate(f'r = {r:.2f}',
                xy=r_xy, xycoords=ax.transAxes)
    ax.annotate(f'p = {p:.3f}',
                xy=p_xy, xycoords=ax.transAxes)


def visualise_correlations(comparison, correlations, filter_type, output_folder):

    if not os.path.exists(f'{output_folder}{filter_type}/'):
        os.makedirs(f'{output_folder}{filter_type}/')

    # 2. Compare dataset against all published
    comparison['color'] = comparison['cluster'].map(cluster_colors)
    for resource in resources:
        fig, ax = plt.subplots()
        sns.scatterplot(x='peptide_stability', y=f'{resource}', data=comparison, hue='cluster', palette=cluster_colors, s=100)
        sns.regplot(comparison['peptide_stability'], comparison[f'{resource}'], scatter_kws={'s': 0}, color='grey', robust=True)
        plt.legend(bbox_to_anchor=(1.0, 1.0))
        plt.ylabel('Normalised stability')
        plt.xlabel('Measured Cm (M)')
        plt.title(resource)
        plt.savefig(f'{output_folder}{filter_type}/{resource}.png')
        plt.show()

    # 3. Compare each cluster against published datasets 
    for resource in resources:
        for cluster, df in comparison[['peptide_stability', f'{resource}', 'cluster']].groupby('cluster'):
            if len(df) > 5: # prevent plotting clusters with only 2 or 3 values that look like false correlations
                fig, ax = plt.subplots()
                sns.scatterplot(x='peptide_stability', y=f'{resource}', data=df, color=cluster_colors[cluster], s=100, label=cluster)
                ## add reg lines for individual cluster_colors
                sns.regplot(df['peptide_stability'], df[f'{resource}'], scatter_kws={'s': 0}, line_kws={'linestyle': '--'},color=cluster_colors[cluster], truncate=False)
        sns.regplot(comparison['peptide_stability'], comparison[f'{resource}'], scatter_kws={'s': 0}, color='grey', truncate=False, robust=True)
        plt.legend(bbox_to_anchor=(1.0, 1.0))
        plt.ylabel('Normalised stability')
        plt.xlabel('Measured Cm (M)')
        plt.title(resource)
        plt.savefig(f'{output_folder}{filter_type}/clustered_{resource}.png')
        plt.savefig(f'{output_folder}{filter_type}/clustered_{resource}.svg')
        plt.show()

    # 4. Plot individual panels for each cluster against each resource
    for resource in resources:
        for cluster, df in comparison[['peptide_stability', f'{resource}', 'cluster']].groupby('cluster'):
            fig, ax = plt.subplots()
            sns.scatterplot(x='peptide_stability', y=f'{resource}', data=df, color=cluster_colors[cluster], s=100)
            # sns.regplot(comparison['peptide_stability'], comparison[f'{resource}'], scatter_kws={'s': 0}, color='grey', truncate=False, robust=True, label='All data')
            if len(df) > 5: # prevent plotting clusters with only 2 or 3 values that look like false correlations
                ## add reg lines for individual cluster_colors
                sns.regplot(df['peptide_stability'], df[f'{resource}'], scatter_kws={'s': 0}, line_kws={'linestyle': '--'},color=cluster_colors[cluster], truncate=False, label=f'Cluster {cluster}')
                corrfunc(df['peptide_stability'], df[f'{resource}'], r_xy=(0.82, 0.92), p_xy=(0.82, 0.85))
            plt.title(f'Cluster {cluster}')
            plt.ylim(-0.05, 1.05)
            # plt.xlim(-0.05, 6.05)
            plt.ylabel('Normalised stability')
            plt.xlabel('Measured Cm (M)')
            plt.title(f'{resource}')
            plt.savefig(f'{output_folder}{filter_type}/{cluster}_{resource}.png')
            plt.show()


def visualise_summary(comparison, correlations, filter_type, output_folder, labels='text'):

    cmap =  sns.diverging_palette(200, 340, s=100, l=30, n=9, center="light", as_cmap=True)

    # 5. Generate summary plot of published correlation (including R and significance)
    # generate scatterplot
    fig, ax = plt.subplots(figsize=(20, 5))
    sns.scatterplot(x='x_pos', y='y_pos', data=correlations, size='size', hue='spearmans_r', palette=cmap, sizes=(100, 1000), hue_norm=(-1, 1))
    # h, l = ax.get_legend_handles_labels()
    plt.legend(bbox_to_anchor=(1.0, 1.0), )
    if labels == 'text':
        plt.yticks(ticks=list(np.arange(0, len(correlations['y_pos'].unique()))), labels=reversed(['All', 'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4']))
        plt.xticks(ticks=list(positions.values()), labels=(list(positions.keys())), rotation=90)
    elif labels == 'numeric':
        plt.yticks(ticks=list(np.arange(0, len(correlations['y_pos'].unique()))), labels=reversed(['All', 'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4']))
        plt.xticks(ticks=np.arange(0, len(positions.values())), labels=(np.arange(1, len(positions.values())+1))) # to cope with reversed positions above
    plt.ylim(-0.5, 4.5)
    plt.xlabel(None)
    plt.ylabel(None)

    plt.savefig(f'{output_folder}{filter_type}correlation_summary.png')
    plt.savefig(f'{output_folder}{filter_type}correlation_summary.svg')

        
    # create custom colorbar
    fig, ax = plt.subplots(figsize=(6, 1))
    fig.subplots_adjust(bottom=0.5)
    norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)
    cb1 = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, orientation='horizontal')
    cb1.set_label("Spearman's R")
    plt.savefig(f'{output_folder}custom_colourbar.svg')


def visualise_heatmap(comparison, correlations, filter_type, output_folder, labels='text'):

    cmap =  sns.diverging_palette(200, 340, s=100, l=30, n=9, center="light", as_cmap=True)
    inverse_positions = {value: key for key, value in positions.items()}
    comparison_positions = {'4': 'Cluster 4', '3': 'Cluster 3',  '2': 'Cluster 2', '1': 'Cluster 1', 'all': 'All', }

    fig, axes = plt.subplots(len(correlations['dataset_2'].unique()), len(correlations['dataset_1'].unique()), figsize=(20, 5))#, sharex=True, sharey=True)
    for j, dataset_1 in enumerate(correlations['dataset_1'].unique()):
        for i, dataset_2 in enumerate(correlations['dataset_2'].unique()):
            # logger.info(f'{i}: {dataset_1}, {j}: {dataset_2}')
            j = positions[dataset_1]
            ax = axes[i][j]
            corr_color = correlations[(correlations['dataset_1'] == dataset_1) & (correlations['dataset_2'] == dataset_2)]['corr_color'].tolist()[0]
            pval = correlations[(correlations['dataset_1'] == dataset_1) & (correlations['dataset_2'] == dataset_2)]['spearmans_pval'].tolist()[0]
            ax.set_facecolor(corr_color)
            if pval < 0.01:
                ax.annotate('**', (0.3, 0.3), fontsize=28)
            if pval < 0.05:
                ax.annotate('*', (0.3, 0.3), fontsize=28)
            else:
                pass

            # axes fix labels
            remove_labels(ax, i, j, f'{positions[dataset_1]}', f'{comparison_positions[dataset_2]}')
    plt.tight_layout(pad=-0.5)
    plt.savefig(f'{output_folder}{filter_type}correlation_heatmap.png')
    plt.savefig(f'{output_folder}{filter_type}correlation_heatmap.svg')



# ------------------------------------------------Read in standard components------------------------------------------------

resource_details = pd.read_excel(f'{published_stabilities}')
resource_details.drop([col for col in resource_details.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
resources = resource_details['dataset_id'].tolist()

stability_summary = pd.read_excel(f'{correlations}', sheet_name=None)
stability_summary.update({key: value.drop([col for col in value.columns.tolist() if 'Unnamed: ' in str(col)], axis=1) for key, value in stability_summary.items()})

datasets = ['Leuenberger_2017', 'Ogburn_2017A', 'Ogburn_2017B', 'Walker_2019', 'Roberts_2016A', 'Roberts_2016B', 'Jarzab_2020F', 'Jarzab_2020G', 'Jarzab_2020H', 'Becher_2016', 'Franken_2015', 'Miettinen_2018', 'Savitski_2018A', 'Savitski_2018B', 'Ball_2020', 'Jarzab_2020N', 'Jarzab_2020O', 'Savitski_2014A', 'Sridharan_2019', 'Jarzab_2020M']

cmap =  sns.diverging_palette(200, 340, s=100, l=30, n=9, center="light", as_cmap=True)
def color_picker(val, val_min=-1, val_max=1):
    norm = colors.Normalize(vmin=val_min, vmax=val_max)
    rgb = cmap(norm(val))[:3] # will return rgba, we take only first 3 so we get rgb
    return colors.rgb2hex(rgb)
#---------------------- Prepare measured correlation datasets----------------------
comparison = stability_summary['measured_comparison'].copy().set_index('KO')[datasets].copy()

correlations = stability_summary['measured_cluster_correlation'].copy().rename(columns={'index': 'datasets'})
correlations[['dataset_1', 'dataset_2']] = correlations['datasets'].str.split('-', expand=True)
# generate positions and size calculations for plotting
correlations['size'] = pd.cut(correlations['spearmans_pval'], bins=[0.0, 0.01, 0.05, 1.0], labels=[0.01, 0.05, 1], include_lowest=True)
# correlations.dropna(subset=['size'], inplace=True) # remove 1:1 comparisons
# positions = resource_details.copy().sort_values(['technique', 'sample_species', 'sample_type'])
positions = dict(zip(datasets, np.arange(0, len(datasets))))
correlations['x_pos'] = correlations['dataset_1'].map(positions)  # to mimic order of the correlation facet plot
correlations['y_pos'] = [0 if cluster == 'all' else int(cluster) for cluster in correlations['dataset_2']]

# visualise_correlations(comparison, correlations, filter_type='measured', output_folder=output_folder)
visualise_summary(comparison, correlations, filter_type='measured', output_folder=output_folder)

#---------------------- Prepare filtered correlation datasets----------------------
comparison = stability_summary['filtered_comparison'].set_index('KO')[datasets].copy()


correlations = stability_summary['filtered_cluster_correlation'].copy().rename(columns={'index': 'datasets'})
correlations[['dataset_1', 'dataset_2']] = correlations['datasets'].str.split('-', expand=True)
correlations = correlations[(correlations['dataset_1'].isin(datasets))]
# generate positions and size calculations for plotting
# correlations['size'] = - np.log10(correlations['pearsons_pval']).replace([np.inf, -np.inf], np.nan)
correlations['size'] = pd.cut(correlations['spearmans_pval'], bins=[0.0, 0.01, 0.05, 1.0], labels=[0.01, 0.05, 1], include_lowest=True)
# correlations.dropna(subset=['size'], inplace=True) # remove 1:1 comparisons
# positions = resource_details.copy().sort_values(['technique', 'sample_species', 'sample_type'])
# positions = dict(zip(positions['dataset_id'], np.arange(0, len(positions['dataset_id']))))
positions = dict(zip(datasets, np.arange(0, len(datasets))))
correlations['x_pos'] = correlations['dataset_1'].map(positions)  # to mimic order of the correlation facet plot
correlations['y_pos'] = [len(correlations['dataset_2'].unique()) -1 - entry for entry in [0 if cluster == 'all' else int(cluster) for cluster in correlations['dataset_2']]]
correlations['corr_color'] = correlations['spearmans_r'].map(color_picker)


# visualise_correlations(comparison, correlations, filter_type='filtered', output_folder=output_folder)
visualise_summary(comparison, correlations, filter_type='filtered', output_folder=output_folder, labels='numeric')

# visualise_correlations(comparison, correlations, filter_type='filtered', output_folder=output_folder)
visualise_heatmap(comparison, correlations, filter_type='filtered', output_folder=output_folder, labels='numeric')

