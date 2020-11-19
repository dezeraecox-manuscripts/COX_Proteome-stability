import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
# test upset plot for chaperones in each cluster
from upsetplot import generate_counts, from_contents
from upsetplot import plot
from simple_venn import venn4
from matplotlib.colors import LinearSegmentedColormap


from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

peptide_folder = f'results/lysate_denaturation/peptide_features/'
protein_folder = f'results/lysate_denaturation/protein_features/'
output_folder = f'results/lysate_denaturation/plot_features/'


if not os.path.exists(output_folder):
    os.makedirs(output_folder)

cluster_colors = { '1': 'darkorange', '2': 'firebrick', '3': 'rebeccapurple', '4': 'royalblue', '0': 'grey', 1: 'darkorange', 2: 'firebrick', 3: 'rebeccapurple', 4: 'royalblue', 0: 'grey', 'multiple': 'grey', 'clusters': 'black'}

font = {'family' : 'normal',
'weight' : 'normal',
'size'   : 14 }
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

def generate_plot(df, x_col, y_col, group_col, xlabel, ylabel, title, yrange=False):
    for group, data in df.groupby(group_col):
        fig, ax = plt.subplots(figsize=(4, 3))
        sns.violinplot(x=x_col, y=y_col, data=data, palette=cluster_colors)
        plt.setp(ax.collections, alpha=.3)
        sns.swarmplot(x=x_col, y=y_col, data=data, palette=cluster_colors)
        if yrange:
            plt.ylim(*yrange)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.title(f'{group}_{title}')
        plt.savefig(f'{output_folder}{group}_{title}.png')
        plt.savefig(f'{output_folder}{group}_{title}.svg')
        plt.show()


# Test better boxplot for plotting
def better_boxplot(dataframe, col_to_plot, cluster_col='cluster', cluster_colors=None, order=None, normed=False, title=False, output_folder=output_folder):

    df = dataframe.copy()

    if normed:
        df[col_to_plot] = (df[col_to_plot] - df[col_to_plot].mean()) / df[col_to_plot].std()

    pop_mean = df[col_to_plot].mean()

    error_df = df.groupby(cluster_col).mean().reset_index()
    error_df['error'] = df.groupby(cluster_col).std().reset_index()[col_to_plot]
    if order:
        error_df = error_df.set_index([cluster_col]).loc[order]
    error_df['x_pos'] = np.arange(0, len(error_df))


    fig, ax = plt.subplots(figsize=(4, 3))
    sns.stripplot(x=df[cluster_col], y=df[col_to_plot], palette=cluster_colors, color=df[cluster_col], alpha=0.5, order=order)
    # sns.boxplot(x=cluster_col, y=col_to_plot, data=df.groupby(cluster_col).mean().reset_index(), palette=cluster_colors, color=df[cluster_col], linewidth=2, ax=ax, order=order)
    plt.errorbar(error_df['x_pos'], error_df[col_to_plot], yerr=error_df['error'], ecolor='black', elinewidth=2.0, capsize=10, barsabove=True, capthick=1.0, marker='_', linewidth=0, markersize=20, markeredgecolor='black')
    ax.axhline(pop_mean, linewidth=1, color='grey', linestyle='--')

    plt.xlabel(cluster_col)
    plt.ylabel(col_to_plot)
    if title:
        plt.title(title)
        plt.savefig(f'{output_folder}{title}{col_to_plot}.png')
        plt.savefig(f'{output_folder}{title}{col_to_plot}.svg')
    plt.show()


def layered_distplot(dataframe, col_to_plot, cluster_col='cluster', cluster_colors=None, order=None, normed=False, title=False, output_folder=output_folder):
    data = dataframe.copy()

    if normed:
        data[col_to_plot] = (data[col_to_plot] - data[col_to_plot].mean()) / data[col_to_plot].std()

    fig, ax = plt.subplots(len(data[cluster_col].unique()), sharex=True)
    for i, (group, df) in enumerate(data.groupby(cluster_col)):
        sns.kdeplot(df[col_to_plot], ax=ax[i], color=cluster_colors[group], shade=True)
        ax[i].axvline(df[col_to_plot].mean(), color='grey', linestyle='--')
        ax[i].get_legend().remove()

    plt.xlabel(col_to_plot)
    plt.ylabel('Density')
    if title:
        plt.title(title)
    # plt.savefig(f'{output_folder}{col_to_plot}.png')
    plt.show()


def generate_protein_heatmap(df, quant_cols, title, cmap='Greys', normed=False, output_folder=output_folder):
    
    dataframe = df.copy()

    if normed:
        dataframe[quant_cols] = (dataframe[quant_cols] - dataframe[quant_cols].mean()) / dataframe[quant_cols].std()

    for group, df in dataframe.groupby('group'):
        fig, ax = plt.subplots(figsize=(4, 3))
        sns.heatmap(df.groupby('cluster_filter_type').mean()[quant_cols], cmap=cmap)
        plt.title(title)
        plt.savefig(f'{output_folder}{title}_{group}.png')
        plt.savefig(f'{output_folder}{title}_{group}.svg')
        plt.show()


def generate_heatmap(df, quant_cols, title, cluster_col='cluster', cmap='Greys', center=None, vmin=None, vmax=None, normed=False, output_folder=output_folder):
    
    dataframe = df.copy()

    if normed:
        dataframe[quant_cols] = (dataframe[quant_cols] - dataframe[quant_cols].mean()) / dataframe[quant_cols].std()

    fig, ax = plt.subplots(figsize=(4, 3))
    sns.heatmap(dataframe.groupby(cluster_col).mean()[quant_cols], cmap=cmap, center=center, vmin=vmin, vmax=vmax)
    plt.title(title)
    plt.savefig(f'{output_folder}{title}.png')
    plt.savefig(f'{output_folder}{title}.svg')
    plt.show()


# ---------------predicted peptide features---------------
if not os.path.exists(f'{output_folder}peptides/'):
    os.makedirs(f'{output_folder}peptides/')

# read in features
peptide_predicted = pd.read_excel(f'{peptide_folder}predicted_features_summary.xlsx', sheet_name='compiled_data')
peptide_predicted.drop([col for col in peptide_predicted.copy().columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# Map colours - in the case of peptides, using the 'mixed' col
# as this relates to the specific peptides added to each cluster
cols_to_plot = ['Polar', 'Neutral', 'Hydrophobic', 'Low', 'Medium', 'High', 'charge_positive', 'charge_neutral', 'charge_negative', 'Helix', 'Strand', 'Coil', 'Buried', 'Exposed', 'Intermediate', 'polar_sum', 'hydrophobic_sum']


for title, group_cols in {'Amino acid class': ['Polar', 'Neutral', 'Hydrophobic'], 'Hydrophobicity': ['Low', 'Medium', 'High'], 'Charge': ['charge_positive', 'charge_neutral', 'charge_negative'], 'Secondary structure': ['Helix', 'Strand', 'Coil'], 'Solvent exposure': ['Buried', 'Exposed', 'Intermediate']}.items():
    # generate_heatmap(peptide_predicted, group_cols, title=title, cmap='Greys', normed=False, output_folder=output_folder)
    generate_heatmap(peptide_predicted, group_cols, title=f'Peptide {title}', cmap=sns.diverging_palette(180, 300, sep=80, n=7), center=0, vmin=-0.1, vmax=0.1, normed=True, output_folder=f'{output_folder}peptides/')

# --------------------calculated peptide features--------------------
peptide_calculated = pd.read_excel(f'{peptide_folder}calculated_peptide_features_summary.xlsx', sheet_name=None)
peptide_calculated.update({key: value.drop([col for col in value.copy().columns.tolist() if 'Unnamed: ' in str(col)], axis=1) for key, value in peptide_calculated.items()})


# Proportion of peptides in each cluster mapped to PDB residue
mapped = peptide_calculated['mapped_chi_observed'].copy().drop('feature_legend', axis=1).T.reset_index().rename(columns={'index':'cluster'})
mapped['proportion'] = mapped[1] / mapped[0] * 100

fig, ax = plt.subplots(figsize=(4, 3))
sns.barplot(x=mapped['cluster'], y=mapped['proportion'], color=mapped['proportion'], palette=cluster_colors)
plt.ylim(0, 100)
plt.ylabel('Proportion of peptides with PDB structures')
plt.xlabel('Cluster')
plt.savefig(f'{output_folder}mapped_peptides.png')
plt.savefig(f'{output_folder}mapped_peptides.svg')

# secondary structure
structure = peptide_calculated['dssp_structure_chi_observed'].copy().set_index('dssp_simple_structure')
structure_proportions = structure.copy() / structure.sum() * 100
structure_proportions = pd.melt(structure_proportions.reset_index(), id_vars='dssp_simple_structure', value_vars=structure_proportions.columns, var_name='cluster', value_name='percent').fillna(0)
# H = α-helix
# B = residue in isolated β-bridge
# E = extended strand, participates in β ladder
# G = 3-helix (310 helix)
# I = 5 helix (π-helix)
# T = hydrogen bonded turn
# S = bend

# Define some colours
colours = {'turn': '#46474a', 'strand': '#84868a', 'helix': '#5f6163', 'bend': '#a6a8ab', '-': '#232324'}
structure_proportions['baseline'] = 0
baselines = dict(zip(structure_proportions['cluster'], structure_proportions['baseline']))

fig, ax = plt.subplots(figsize=(4, 3))
for x, (structure, df) in enumerate(structure_proportions.groupby('dssp_simple_structure')):
    logger.info(x)
    position = (len(structure_proportions['dssp_simple_structure'].unique()) - x)
    df['baseline'] = df['cluster'].map(baselines)
    df['percent'] = df['baseline'] + df['percent']
    sns.barplot(x='cluster', y='percent', data=df, color=colours[structure], label=structure, zorder=position, ax=ax)
    baselines.update(dict(zip(df['cluster'], df['percent'])))
ax.set_ylim(0, 100)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], title='Structure', bbox_to_anchor=(1.0, 1.0))
plt.xlabel('Cluster')
plt.ylabel('Proportion of mapped structures')
plt.savefig(f'{output_folder}cys_residues_simplestructure.png', dpi=500)
plt.savefig(f'{output_folder}cys_residues_simplestructure.svg')
plt.show()


# asa
calculated_asa = peptide_calculated['compiled_clustered'].copy()

for col in ['pdb_asa', 'dssp_asa']:
    better_boxplot(calculated_asa, col, cluster_col='cluster', cluster_colors=cluster_colors, normed=False, title=f'Peptide _', output_folder=f'{output_folder}')
    layered_distplot(calculated_asa, col, cluster_col='cluster', cluster_colors=cluster_colors, normed=False, title=f'Peptide_', output_folder=f'{output_folder}')


# site info - residue features
residue_features = peptide_calculated['residue_features_chi_observed'].copy().set_index('feature_legend')
residue_features_proportions = residue_features.copy() / residue_features.sum() * 100
residue_features_proportions = pd.melt(residue_features_proportions.reset_index(), id_vars='feature_legend', value_vars=residue_features_proportions.columns, var_name='cluster', value_name='percent').fillna(0)
residue_features_proportions = residue_features_proportions[residue_features_proportions['feature_legend'] == 1]

fig, ax = plt.subplots(figsize=(4, 3))
sns.barplot(x=residue_features_proportions['cluster'], y=residue_features_proportions['percent'], color=residue_features_proportions['cluster'], palette=cluster_colors, order=[1, 2, 3, 4])
plt.ylabel('Proportion of peptides\nwith residue feature')
plt.xlabel('Type')
plt.savefig(f'{output_folder}residue_features.png')
plt.savefig(f'{output_folder}residue_features.svg')
plt.show()

# site info - domain features
domain_features = peptide_calculated['domain_features_chi_observed'].copy().set_index('feature_legend')
domain_features_proportions = domain_features.copy() / domain_features.sum() * 100
domain_features_proportions = pd.melt(domain_features_proportions.reset_index(), id_vars='feature_legend', value_vars=domain_features_proportions.columns, var_name='cluster', value_name='percent').fillna(0)
domain_features_proportions = domain_features_proportions[domain_features_proportions['feature_legend'] == 2]

fig, ax = plt.subplots(figsize=(4, 3))
sns.barplot(x=domain_features_proportions['cluster'], y=domain_features_proportions['percent'], color=domain_features_proportions['cluster'], palette=cluster_colors, order=[1, 2, 3, 4])
plt.ylabel('Proportion of peptides\nwith domain feature')
plt.xlabel('Type')
plt.savefig(f'{output_folder}domain_features.png')
plt.savefig(f'{output_folder}domain_features.svg')
plt.show()

# site info - PFAM domains
domain_features = peptide_calculated['pfam_chi_observed'].copy().set_index('feature_legend')
domain_features_proportions = domain_features.copy() / domain_features.sum() * 100
domain_features_proportions = pd.melt(domain_features_proportions.reset_index(), id_vars='feature_legend', value_vars=domain_features_proportions.columns, var_name='cluster', value_name='percent').fillna(0)
domain_features_proportions = domain_features_proportions[domain_features_proportions['feature_legend'] == 1]

fig, ax = plt.subplots(figsize=(4, 3))
sns.barplot(x=domain_features_proportions['cluster'], y=domain_features_proportions['percent'], color=domain_features_proportions['cluster'], palette=cluster_colors, order=[1, 2, 3, 4])
plt.ylabel('Proportion of residues\nlocated in PFAM domain')
plt.xlabel('Cluster')
plt.savefig(f'{output_folder}pfam_domains.png')
plt.savefig(f'{output_folder}pfam_domains.svg')
plt.show()

# ---------------predicted protein features---------------

if not os.path.exists(f'{output_folder}proteins/'):
    os.makedirs(f'{output_folder}proteins/')

order=[1, 2, 3, 4, 'multiple']

# read in features
protein_predicted = pd.read_excel(f'{protein_folder}predicted_features_anova.xlsx', sheet_name='compiled_data')
protein_predicted.drop([col for col in protein_predicted.copy().columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# Map colours - in the case of peptides, using the 'mixed' col
# as this relates to the specific peptides added to each cluster
cols_to_plot = ['Polar', 'Neutral', 'Hydrophobic', 'Low', 'Medium', 'High', 'charge_positive', 'charge_neutral', 'charge_negative', 'Helix', 'Strand', 'Coil', 'Buried', 'Exposed', 'Intermediate', 'polar_sum', 'hydrophobic_sum']

for title, group_cols in {'Amino acid class': ['Polar', 'Neutral', 'Hydrophobic'], 'Hydrophobicity': ['Low', 'Medium', 'High'], 'Charge': ['charge_positive', 'charge_neutral', 'charge_negative'], 'Secondary structure': ['Helix', 'Strand', 'Coil'], 'Solvent exposure': ['Buried', 'Exposed', 'Intermediate']}.items():
    # generate_heatmap(protein_predicted, group_cols, title=title, cmap='Greys', normed=False, output_folder=output_folder)
    generate_heatmap(protein_predicted, group_cols, cluster_col='unique', title=f'Protein {title}', cmap=sns.diverging_palette(180, 300, sep=80, n=7), center=0, vmin=-0.3, vmax=0.3, normed=True, output_folder=f'{output_folder}proteins/')


protein_disorder = pd.read_excel(f'{protein_folder}disorder_prediction_kruskal.xlsx', sheet_name='compiled_data')
protein_disorder.drop([col for col in protein_disorder.copy().columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
protein_disorder['unique'] = protein_disorder['unique']
better_boxplot(protein_disorder[protein_disorder['disordered'] == 1], 'proportion', cluster_col='unique', cluster_colors=cluster_colors, order=order, normed=False, title=f'protein_disorder')
layered_distplot(protein_disorder[protein_disorder['disordered'] == 1], 'proportion', cluster_col='unique', cluster_colors=cluster_colors, normed=False, title=f'protein_disorder')


# --------------------protein features--------------------
order=[1, 2, 3, 4, 'multiple']

peptides = pd.read_excel(f'{protein_folder}calculated_protein_features_summary.xlsx', sheet_name=None)

molweight = peptides['len_mw_summary'].copy()
molweight['kda'] = molweight['MW'] / 1000

better_boxplot(molweight, 'kda', cluster_col='unique', cluster_colors=cluster_colors, order=order, normed=False, title=f'mol_weight', output_folder=output_folder)

# Molecular weight
fig, ax = plt.subplots(figsize=(4, 3))
# sns.violinplot(x=molweight['unique'], y=molweight['kda'], color=molweight['unique'], palette=cluster_colors, order=[1, 2, 3, 4, 'multiple'])
# plt.setp(ax.collections, alpha=.3)
sns.boxplot(x=molweight['unique'], y=molweight['kda'], color='white', order=order, fliersize=0)
sns.swarmplot(x=molweight['unique'], y=molweight['kda'], color=molweight['unique'], palette=cluster_colors, order=order, alpha=0.5)
plt.ylabel('Molecular weight (kDa)')
plt.xlabel('Cluster')
plt.savefig(f'{output_folder}MolWeight_unique.png')
plt.savefig(f'{output_folder}MolWeight_unique.svg')
plt.show()


better_boxplot(molweight, 'length', cluster_col='unique', cluster_colors=cluster_colors, order=order, normed=False, title=f'length', output_folder=output_folder)

# Length
fig, ax = plt.subplots(figsize=(4, 3))
# sns.violinplot(x=molweight['unique'], y=molweight['kda'], color=molweight['unique'], palette=cluster_colors, order=[1, 2, 3, 4, 'multiple'])
# plt.setp(ax.collections, alpha=.3)
sns.boxplot(x=molweight['unique'], y=molweight['length'], color='white', order=order, fliersize=0)
sns.swarmplot(x=molweight['unique'], y=molweight['length'], color=molweight['unique'], palette=cluster_colors, order=order, alpha=0.5)
plt.ylabel('Protein Length')
plt.xlabel('Cluster')
plt.savefig(f'{output_folder}protein_length_unique.png')
plt.savefig(f'{output_folder}protein_length_unique.svg')
plt.show()


# PFAM domains
fig, ax = plt.subplots(figsize=(4, 3))
# sns.violinplot(x=molweight['unique'], y=molweight['kda'], color=molweight['unique'], palette=cluster_colors, order=[1, 2, 3, 4, 'multiple'])
# plt.setp(ax.collections, alpha=.3)
sns.swarmplot(x=molweight['unique'], y=molweight['pfam_id'], color=molweight['unique'], palette=cluster_colors, order=order, alpha=0.5)
sns.boxplot(x=molweight['unique'], y=molweight['pfam_id'], color='white', order=order, fliersize=0)
plt.ylabel('Number of PFAM domains')
plt.xlabel('Cluster')
plt.savefig(f'{output_folder}protein_pfam_unique.png')
plt.savefig(f'{output_folder}protein_pfam_unique.svg')
plt.show()

# Protein venn diagram - Note cluster colours are slightly off here - need to add manually for now
protein_venn = pd.read_excel(f'{protein_folder}protein_venn.xlsx')
protein_venn.drop([col for col in protein_venn.columns.tolist() if 'Unnamed: ' in str(col)], axis=1, inplace=True)
protein_venn['cluster'] = protein_venn['name'].str.split(' ').str[-1]
labels = protein_venn['name'][0:np.max([len(x) for x in protein_venn['cluster']])]

fig, ax = plt.subplots(figsize=(4, 3))
# venn4(protein_venn['count'], set_labels=labels,  set_colors=[cluster_colors[cluster.split(' ')[-1]] for cluster in labels], ax=ax)
venn4(protein_venn['count'], set_labels=labels,  set_colors=['royalblue', 'rebeccapurple', 'darkorange', 'firebrick'], ax=ax)
plt.savefig(f'{output_folder}protein_clustered_venn.png')
plt.savefig(f'{output_folder}protein_clustered_venn.svg')
plt.show()


# Chaperone venn diagram - Note cluster colours are slightly off here - need to add manually for now
chaperones = pd.read_excel(f'{protein_folder}chaperone_enrichment_summary.xlsx', sheet_name=None)
chaperones.update({key: value.drop([col for col in value.columns.tolist() if 'Unnamed: ' in str(col)], axis=1) for key, value in chaperones.items()})

chaperone_counts = chaperones['overlap']
chaperone_counts['cluster'] = chaperone_counts['name'].str.split(' ').str[-1]
labels = chaperone_counts['name'][0:np.max([len(x) for x in chaperone_counts['cluster']])]

fig, ax = plt.subplots(figsize=(4, 3))
# venn4(chaperone_counts['count'], set_labels=labels,  set_colors=[cluster_colors[cluster.split(' ')[-1]] for cluster in labels], ax=ax)
venn4(chaperone_counts['count'], set_labels=labels,  set_colors=['royalblue', 'rebeccapurple', 'darkorange', 'firebrick'], ax=ax)
plt.savefig(f'{output_folder}chaperone_clustered_venn.png')
plt.savefig(f'{output_folder}chaperone_clustered_venn.svg')
plt.show()

# Generate upset plot
clustered_chaperones = chaperones['clustered_chaperones'].copy()
clustered_chaperones = clustered_chaperones[clustered_chaperones['chaperone'] == 1]
chaperone_membership = {str(cluster): df['Proteins'].unique().tolist() for cluster, df in clustered_chaperones.groupby('cluster')}
counts = from_contents(chaperone_membership)

plot(counts)
plt.savefig(f'{output_folder}chaperone_upset.png')
plt.savefig(f'{output_folder}chaperone_upset.svg')
plt.show()

# Generate heatmap representation
chaperone_details = chaperones['chaperone_details'].copy()[['Proteins', 'mixed', 'Gene names']]
chaperone_details['name'] = chaperone_details['Gene names'].str.split(' ').str[0].str.upper()

chaperone_heatmap = counts.reset_index()
chaperone_heatmap['cluster_num'] = chaperone_heatmap.sum(axis=1)
for col in chaperone_heatmap.set_index(['id', 'cluster_num']).columns.tolist():
    chaperone_heatmap[col] = [int(col) if val else np.nan for val in chaperone_heatmap[col]]
chaperone_heatmap['name'] = chaperone_heatmap['id'].map(dict(chaperone_details[['Proteins', 'name']].values))
chaperone_heatmap.sort_values(['cluster_num', 'name'], inplace=True)
cmap = LinearSegmentedColormap.from_list('clusters', ['darkorange', 'firebrick', 'rebeccapurple', 'royalblue'], 4)

fig, ax = plt.subplots(figsize=(5, 20))
sns.heatmap(chaperone_heatmap.set_index('name')[['1', '2', '3', '4']], cmap=cmap)
plt.xlabel('Cluster')
plt.savefig(f'{output_folder}chaperone_heatmap.png')
plt.savefig(f'{output_folder}chaperone_heatmap.svg')

# Visualise count of chaperones represented in each cluster
chaperone_proportions = chaperones['chaperones_chi_obs'].copy().set_index('chaperone').T.reset_index().rename(columns={'index': 'unique'})
chaperone_proportions['names'] = [sorted(clustered_chaperones[clustered_chaperones['unique'] == cluster]['Proteins'].map(dict(chaperone_details[['Proteins', 'name']].values)).unique().tolist()) for cluster in chaperone_proportions['unique']]
chaperone_proportions['percent'] = round(chaperone_proportions[1] / (chaperone_proportions[0]+chaperone_proportions[1]) * 100, 1)

fig, ax = plt.subplots(figsize=(6, 5))
sns.barplot(x='unique', y=1, data=chaperone_proportions, order=[1, 2, 3, 4, 'multiple'], color='unique', palette=cluster_colors)
# add annotations for proteins in each bar
for x, cluster in enumerate([1, 2, 3, 4, 'multiple']):
    names = chaperone_proportions[chaperone_proportions['unique'] == cluster]['names'].tolist()[0]
    for y, name in enumerate(names):
        plt.annotate(name, (x-0.25, y+0.15), color='white', size=8)
# add percent of total proteins
for x, cluster in enumerate([1, 2, 3, 4, 'multiple']):
    percent = chaperone_proportions[chaperone_proportions['unique'] == cluster]['percent'].tolist()[0]
    max_val = chaperone_proportions[chaperone_proportions['unique'] == cluster][1].tolist()[0]
    plt.annotate(f'{percent}%', (x-0.2, max_val + 1), color='black', size=8)
plt.ylim(-0.9, 20.1)
plt.yticks(np.arange(0, 21, 2))
plt.ylabel('Number of chaperones')
plt.xlabel('Cluster')
plt.savefig(f'{output_folder}chaperone_count.png')
plt.savefig(f'{output_folder}chaperone_count.svg')