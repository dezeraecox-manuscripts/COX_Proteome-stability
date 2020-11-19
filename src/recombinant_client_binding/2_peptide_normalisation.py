import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
import itertools

from collections import defaultdict
from loguru import logger
from GEN_Utils import FileHandling, CalcUtils

from fuzzywuzzy import fuzz
from fuzzywuzzy import process

logger.info('Import OK')



def cys_noncys_filter(peptides):
    """ Assigns cys_rank, then separates cys and noncys peptides.
    Collects only peptides from proteins that have at least one cys and one non-cys peptide"""

    # Add cys rank to identify cys peptides
    peptides['cys_rank'] = [
        1 if 'C' in sequence else 0 for sequence in peptides['Sequence']]

    # separate cys and noncys peptides
    cys_peptides = peptides[peptides['cys_rank'] == 1]
    noncys_peptides = peptides[peptides['cys_rank'] == 0]

    # Collect only proteins with at least one cys, one noncys
    common_prots = set(cys_peptides['Proteins']).intersection(
        set(noncys_peptides['Proteins']))
    cys_peptides = cys_peptides[cys_peptides['Proteins'].isin(common_prots)]
    noncys_peptides = noncys_peptides[noncys_peptides['Proteins'].isin(
        common_prots)]

    return cys_peptides, noncys_peptides


def noncys_ci_calculator(noncys_peptides, sample_cols):
    """Calculates relevant info per protein from noncys peptides, including overall mean, SD
    and the per-channel mean (used for cys/noncys calculation)"""

    # for each protein, collect all noncys values and determine mean + SD
    noncys_cis = {}
    for protein, df in noncys_peptides.groupby(['Proteins']):
        values = df[sample_cols].values.flatten()
        # Remove NaN values for calculations
        values = values[~np.isnan(values)]
        noncys_cis[protein] = {'df': df,
                               'noncys_means': df[sample_cols].mean(),
                               'noncys_medians': df[sample_cols].median(),
                               'num_peptides': df.shape[0],
                               'values': values,
                               'pop_mean': np.mean(values),
                               'pop_median': np.median(values),
                               'pop_stdev': np.std(values),
                               'num_vals': len(values)}
    return noncys_cis


def cys_ratio(cys_peptides, noncys_cis, sample_cols):
    """Calculates cys/noncys ratio per cys peptide, using the noncys per protein mean"""
    # Generate cys/noncys ratio
    norm_cys = []
    for protein, df in cys_peptides.groupby('Proteins'):
        data = df.copy()
        noncys_vals = noncys_cis[protein]['noncys_means']
        data[sample_cols] = data[sample_cols] / noncys_vals
        norm_cys.append(data)

    return pd.concat(norm_cys)


def main(sample_name):

    # -----------------------------------Normalisation----------------------------------------------------------
    raw_data = pd.read_excel(
        f'{input_folder}{sample_name}_Compiled.xlsx', sheet_name=None)

    info_cols = ['Sequence', 'Proteins', 'Gene names',
                'Protein names', 'cys_rank', 'replicate']
    peptides = raw_data['Peptides'].copy().drop(
        ['Unique (Groups)', 'Unique (Proteins)'], axis=1)

    # remove any peptides with missed cleavage
    peptides = peptides[peptides['Missed cleavages'].isnull()]  # 165 down to 75?!

    peptides.set_index([col for col in peptides.columns.tolist()
                        if col in info_cols], inplace=True)

    # Adjust column names to be in format label_replicate
    peptides.columns = [('_').join(re.split(r' |_', col)[3:6])
                        for col in peptides.columns.tolist()]

    # collecting and sorting individual replicates
    # initialise as a default dict to allow nested key assignment
    replicate_normalised = defaultdict(dict)
    for replicate in replicates:

        # collect relevant channels
        replicate_cols = [col for col in peptides.columns.tolist()
                        if f'_{replicate}' in col]
        rep_peptides = peptides[replicate_cols]

        # relabel columns with channel number
        sample_cols = [col.split('_')[0] for col in replicate_cols]
        rep_peptides.columns = sample_cols

        rep_peptides['replicate'] = replicate
        rep_peptides['cys_rank'] = [
            1 if 'C' in sequence else 0 for sequence in rep_peptides.reset_index()['Sequence']]

        replicate_normalised[replicate]['raw'] = rep_peptides.reset_index()


    compiled = pd.concat([replicate_normalised[replicate]['raw']
                        for replicate in replicates])

    # scale according to proportion of the sample
    scaling_factor = {
        'P11142': {'1': 2/8, '2': 2/8, '3': 2/8, '4': 2/8, '5': 2/7, '6': 2/7, '7': 2/3, '8': 2/3, '9': 2/2, '10': 2/2, '11': 2/2},
        'P00346': {'1': 5/8, '2': 5/8, '3': 5/8, '4': 5/8, '5': 5/7, '6': 5/7, '7': 0/3, '8': 0/3, '9': 0/2, '10': 0/2, '11': 0/2},
        'P25685': {'1': 1/8, '2': 1/8, '3': 1/8, '4': 1/8, '5': 0/7, '6': 0/7, '7': 1/3, '8': 1/3, '9': 0/2, '10': 0/2, '11': 0/2},
    }
    scaled_vals = compiled.copy()
    scaled = []
    for protein, df in scaled_vals.groupby('Proteins'):
        scaling = scaling_factor[protein]
        for col in sample_cols:
            scaling[col]
            df[col] = df[col] / (scaling[col])
        scaled.append(df)
    scaled = pd.concat(scaled).replace([np.inf, -np.inf], np.nan)

    # separate cys and noncys, generate ratio
    ratios = defaultdict(list)
    for replicate, df in scaled.groupby('replicate'):
        ratios['raw'].append(df)
        df[sample_cols] = np.log10(df[sample_cols])
        cys_peptides, noncys_peptides = cys_noncys_filter(df)
        # Calculate relevant noncys info
        noncys_cis = noncys_ci_calculator(noncys_peptides, sample_cols)

        # Calculate cys/noncys ratio
        cys_noncys_peptides = cys_ratio(cys_peptides, noncys_cis, sample_cols)
        cys_noncys_peptides

        # Calculate noncys/noncys ratio
        noncys_noncys_peptides = cys_ratio(
            noncys_peptides, noncys_cis, sample_cols)

        ratios['cys'].append(cys_peptides)
        ratios['noncys'].append(noncys_peptides)
        ratios['cys_noncys'].append(cys_noncys_peptides)
        ratios['noncys_noncys'].append(noncys_noncys_peptides)

    # regenerate compiled dfs, save to excel
    for key, value in ratios.items():
        ratios[key] = pd.concat(value)

    FileHandling.df_to_excel(output_path=f'{output_folder}normalised_summary_{sample_name}.xlsx', sheetnames=list(
        ratios.keys()), data_frames=list(ratios.values()))

if __name__ == "__main__":
    
    input_folder = f'results/recombinant_client_assay/initial_cleanup/'
    output_folder = f'results/recombinant_client_assay/peptide_normalisation/'
    replicates = ['1', '2', '3', '4']  # (standard replicate must be listed first)

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    for sample_name in ['Heat-A', 'Heat-B']:
        main(sample_name)


