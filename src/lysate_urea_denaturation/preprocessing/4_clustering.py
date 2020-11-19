import os, functools
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import skfuzzy as fuzz
from kneed import KneeLocator
from sklearn.decomposition import PCA

from GEN_Utils import FileHandling
from loguru import logger

logger.info("Import ok")

def multiple_PCAs(test_dict):
    """test_dict: dict mapping str(data_type): (df, sample_cols)"""
    pcas = {}
    for data_type, (data, sample_cols) in test_dict.items():
        for_PCA = data[sample_cols].copy().fillna(0)
        pca = PCA(n_components=2)
        principalComponents = pca.fit_transform(for_PCA.values)
        principalDf = pd.DataFrame(data = principalComponents, columns = ['PC1', 'PC2'])
        principalDf.index = data['Sequence']
        pcas[data_type] = principalDf

        logger.info(f'{data_type}: {len(principalDf)}')

        # visualise the PCA
        fig = plt.figure(figsize = (8,8))
        ax = fig.add_subplot(1,1,1) 
        ax.set_xlabel('Principal Component 1', fontsize = 15)
        ax.set_ylabel('Principal Component 2', fontsize = 15)
        ax.set_title(data_type, fontsize = 20)
        ax.scatter(principalDf['PC1'] , principalDf['PC2'], s = 50)
        # plt.savefig(f'{output_folder}{data_type}_PCA.png')
    
    return pcas


def fuzzy_clustering(data_type, data, sample_cols, max_clusters=10):
    alldata = data[sample_cols].fillna(0).T.values
    fpcs = []
    cluster_membership = {}
    cluster_score = {}
    cluster_parameters = {}

    for ncenters in range(2, max_clusters):
        cntr, u, u0, d, jm, p, fpc = fuzz.cluster.cmeans(
            alldata, ncenters, 4, error=0.005, maxiter=1000, init=None)
        cluster_parameters[ncenters] = [cntr, u, u0, d, jm, p, fpc]
        cluster_membership[ncenters] = np.argmax(u, axis=0)
        cluster_score[ncenters] = np.max(u, axis=0)
        # Store fpc values for later
        fpcs.append(fpc)

    # clean cluster parameters
    parameters = pd.DataFrame(cluster_parameters)
    parameters.index = ['cntr', 'u', 'u0', 'd', 'jm', 'p', 'fpc'] 
    
    # clean cluster data
    membership = pd.DataFrame(cluster_membership)
    membership.index = data.index.tolist()
    membership.columns = [f'member_{col}' for col in membership.columns.tolist()]
    score = pd.DataFrame(cluster_score)
    score.index = data.index.tolist()
    score.columns = [f'score_{col}' for col in score.columns.tolist()]
    # Generate merged cluster info
    cluster_membership = pd.merge(score, membership, left_index=True, right_index=True)
    clustered = pd.merge(cluster_membership, data, left_index=True, right_index=True)

    fig2, ax2 = plt.subplots()
    ax2.plot(np.arange(2, max_clusters), fpcs)
    ax2.set_xlabel("Number of centers")
    ax2.set_ylabel("Fuzzy partition coefficient")
    plt.title(data_type)
    plt.show()

    return [cluster_membership, clustered, fpcs, parameters]


def fuzzy_prediction(data_type, data, sample_cols, cntr):
    alldata = data[sample_cols].fillna(0).T.values

    u, u0, d, jm, p, fpc = fuzz.cluster.cmeans_predict(alldata, cntr, 4, error=0.005, maxiter=1000)

    # clean cluster parameters
    parameters = pd.DataFrame([cntr, u, u0, d, jm, p, fpc])
    parameters.index = ['cntr', 'u', 'u0', 'd', 'jm', 'p', 'fpc'] 

    # clean cluster data
    membership = pd.DataFrame(np.argmax(u, axis=0))
    membership.index = data.index.tolist()
    membership.columns = ['member_predicted']

    score = pd.DataFrame(np.max(u, axis=0))
    score.index = data.index.tolist()
    score.columns = ['score_predicted']

    # Generate merged cluster info
    cluster_membership = pd.merge(score, membership, left_index=True, right_index=True)
    return pd.merge(cluster_membership, data, left_index=True, right_index=True)


def multiple_fuzzy_clustering(test_dict, max_clusters=20):
    """test_dict: dict mapping str(data_type): (df, sample_cols)"""

    clustering = {}
    for data_type, (data, sample_cols) in test_dict.items():
        clustering[data_type] = fuzzy_clustering(data_type, data, sample_cols, max_clusters)

    return clustering


def clustered_pca(clustering, pcas, visualise, save=False, output_path=None, palette=False):
    """visualise: dict mapping str(data_type): (cluster_number, sample_cols)"""
    if not palette:
        palette = 'Set1'
    for data_type, (clusters, cols) in visualise.items():
        pca_data = pcas[data_type].reset_index()
        clustered = clustering[data_type][0]
        clustered = clustered[[col for col in clustered.columns.tolist() if f'_{clusters}' in col]]
        for_plotting = pd.merge(clustered, pca_data, left_index=True, right_index=True)

        fig, ax = plt.subplots(figsize=(4, 4))
        sns.scatterplot(data=for_plotting, x='PC1', y='PC2', hue=f'member_{clusters}', palette=palette, size=f'score_{clusters}', legend=None)
        plt.title(data_type)
        if save:
            plt.savefig(f'{output_path}{data_type}_PCA.png')
        plt.tight_layout()
        plt.show()


def clustered_curves(clustering, visualise, range_vals=None, save=False, output_path=None, cluster_colors=None):

    """visualise: dict mapping str(data_type): (cluster_number, sample_cols)"""
    for data_type, (clusters, cols) in visualise.items():
        clustered = clustering[data_type][1]
        cluster_info = [col for col in clustered.columns.tolist() if f'_{clusters}' in str(col)]
        clustered = clustered[cluster_info + cols]
        for_plotting = pd.melt(clustered, id_vars=cluster_info, value_vars=cols, var_name='channel', value_name='ratio').fillna(0)

        for cluster_number, df in for_plotting.groupby(f'member_{clusters}'):
            color = cluster_colors[cluster_number] if cluster_colors else 'Purples'
            fig, axes = plt.subplots()
            sns.lineplot(data=df, x='channel', y='ratio', hue=f'score_{clusters}', palette=color, legend=None)
            plt.title(f'{data_type}_{cluster_number}')
            plt.tight_layout()
            if range_vals:
                plt.ylim(*range_vals)
            if save:
                plt.savefig(f'{output_path}{data_type}_{cluster_number}.png')
            plt.show()


if __name__ == "__main__":


    max_clusters = 10
    input_path = f'results/lysate_denaturation/smoothing/processed_data.xlsx'
    output_folder = f'results/lysate_denaturation/clustering/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # read in raw data
    processed = pd.read_excel(input_path, sheet_name=None)
    info_cols = ['Sequence', 'Proteins', 'replicate', 'coverage']

    # Collect list of datasets to be clustered
    test_dict = {}
    for key, df in processed.items():
        df.drop([col for col in df.columns.tolist() if 'Unnamed: ' in str(col)], axis=1, inplace=True)
        cols = [col for col in df.columns.tolist() if str(col) not in info_cols]
        test_dict[key] = (df, cols)


    # complete PCA on each of the data types
    pcas = multiple_PCAs(test_dict)
    FileHandling.df_to_excel(output_path=f'{output_folder}PCA_data.xlsx', sheetnames=list(pcas.keys()), data_frames=list(pcas.values()))

    ## Test range of clustering numbers, fuzzy factorsusing  unsupervised clustering with the fuzzy c-means algorithm
    clustering = multiple_fuzzy_clustering(test_dict, max_clusters=max_clusters)

    # save cluster data
    dfs = []
    for clustering_df in clustering.values():
        df = clustering_df[1]
        cleaned = df.set_index('Sequence').copy()[[col for col in df.columns.tolist() if "_" in str(col)]]
        dfs.append(cleaned)
    FileHandling.df_to_excel(output_path=f'{output_folder}clusters_data.xlsx', sheetnames=list(clustering.keys()), data_frames=dfs)

    # save clustering fpc dfs
    fpcs = [pd.DataFrame(clustering_df[2]) for clustering_df in clustering.values()]
    FileHandling.df_to_excel(output_path=f'{output_folder}clusters_fpc.xlsx', sheetnames=list(clustering.keys()), data_frames=fpcs)

    # compile cluster, protein, sequence, pca information
    protein_mapper = dict(zip(processed['raw']['Sequence'], processed['raw']['Proteins']))
    summary_dict = {}
    for data_type, (clusters, clustered_data, fpcs, parameters) in clustering.items():
        pca_data = pcas[data_type].reset_index()
        summary_df = pd.merge(clusters, pca_data, left_index=True, right_index=True)
        original_data = test_dict[data_type][0]
        summary_df = pd.merge(summary_df, original_data, on='Sequence')
        summary_df['Proteins'] = summary_df['Sequence'].map(protein_mapper)
        summary_dict[data_type] = summary_df
    FileHandling.df_to_excel(output_path=f'{output_folder}cluster_summary.xlsx', sheetnames=list(summary_dict.keys()), data_frames=list(summary_dict.values()))

    # Determine maximum cluster efficiency, although this method seems to be overly insensitive for our data
    fpcs = pd.read_excel(f'{output_folder}clusters_fpc.xlsx', sheet_name=None)
    for data_type, df in fpcs.items():
        df = df.copy().drop([col for col in df.columns.tolist() if 'Unnamed: ' in str(col)], axis=1)
        y_vals = list(df[0])
        x_vals = list(np.arange(2, len(df)+2))
        kneedle = KneeLocator(x_vals, y_vals, S=0.1, curve='convex', direction='decreasing')

        fig, ax = plt.subplots()
        sns.lineplot(x_vals, y_vals, color='grey')
        try:
            ax.axvline(kneedle.elbow, color='red', linestyle='--')
        except:
            pass
        plt.title(data_type)
        plt.xlabel('Number of clusters')
        plt.ylabel('FPC')
        plt.savefig(f'{output_folder}needleplot_{data_type}.png')
        plt.show()