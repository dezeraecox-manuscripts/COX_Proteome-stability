import os, re
import numpy as np
import pandas as pd
from GEN_Utils import FileHandling

from loguru import logger

logger.info('Import OK')


def norm_control_plex(peptides, sample_cols, pooled_col, standard_vals):
    """ Normalises per peptide according to control channel. Returns df with updated sample cols and removes control channel"""
    pooled_factor = pd.merge(pd.DataFrame(peptides[pooled_col]), pd.DataFrame(standard_vals), left_index=True, right_index=True)
    pooled_factor['pooled_factor'] = pooled_factor.iloc[:, 0] / pooled_factor.iloc[:, 1]
    
    # normalise to control channel
    peptides[sample_cols] = peptides[sample_cols].multiply(
        pooled_factor['pooled_factor'], axis=0)

    return peptides


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
                               'num_peptides': df.shape[0],
                               'values': values,
                               'pop_mean': np.mean(values),
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

    cys_noncys_peptides = pd.concat(norm_cys)

    return cys_noncys_peptides


def med_normalise(peptides, control_plex):
    """Calculates correction factor from median of column 1, then normalises all other channels to this factor"""

    medians = peptides.median()
    control_factor = medians[control_plex] / medians

    # normalise to control channel
    peptides = peptides * control_factor

    return peptides


#---------------------------------------------------------------------------------------------

if __name__ == "__main__":

    input_folder = 'results/lysate_denaturation/initial_cleanup/'
    output_folder = 'results/lysate_denaturation/normalised/'
    filter_cols = []
    pooled_col = None
    control_plex = '1'
    replicates = ['1', '2', '3'] #(standard replicate must be listed first)
    quant_threshold = 6

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    raw_data = pd.read_excel(f'{input_folder}Compiled.xlsx', sheet_name=None)

    info_cols = ['Sequence', 'Proteins', 'Gene names', 'Protein names', 'cys_rank', 'replicate']
    peptides = raw_data['Peptides'].copy().set_index([col for col in raw_data['Peptides'].columns.tolist() if col in info_cols])
    urea_conc = {'1': 0.0, '2': 0.5, '3': 1.0, '4': 1.5, '5': 2.0, '6': 2.5, '7': 3.0, '8': 3.5, '9': 4.0, '10': 4.5, '11': 5.0, '12': 5.5, '13': 6.0}

    # remove outlier samples as identified in initial cleanup/batch effects
    peptides = peptides.copy().drop(filter_cols, axis=1)

    # collecting and sorting individual replicates, calculate cy/noncys
    replicate_normalised = {}
    for replicate in replicates:
        replicate

        # collect relevant channels
        replicate_cols = [col for col in peptides.columns.tolist() if f'_{replicate}' in col]
        rep_peptides = peptides[replicate_cols]
        raw = rep_peptides.copy()

        if pooled_col:
            pooled_col = f'{pooled_plex}_{replicate}'
            # drop control channel
            rep_peptides = rep_peptides.drop(pooled_col, axis=1)

        # filter out cols seen to be outliers via PCA
        if filter_cols:
            for col in filter_cols:
                if col in rep_peptides.columns.tolist():
                    # drop control channel
                    rep_peptides = rep_peptides.drop(col, axis=1)

        # relabel columns with channel number
        sample_cols = rep_peptides.columns.tolist()
        new_sample_cols = [col.split('_')[0] for col in sample_cols]
        rep_peptides.columns = new_sample_cols

        # median-normalise across all samples
        rep_peptides = med_normalise(rep_peptides, control_plex)
        med_norm = rep_peptides.copy()

        # separate cys and noncys, generate ratio
        cys_peptides, noncys_peptides = cys_noncys_filter(
            rep_peptides.reset_index())

        # Calculate relevant noncys info
        noncys_cis = noncys_ci_calculator(noncys_peptides, new_sample_cols)

        # Calculate cys/noncys ratio
        cys_noncys_peptides = cys_ratio(cys_peptides, noncys_cis, new_sample_cols)

        # Calculate noncys/noncys ratio
        noncys_noncys_peptides = cys_ratio(
            noncys_peptides, noncys_cis, new_sample_cols)

        # Normalise to control (OM urea) sample such that this is 1 (as there should be no difference with no denaturant)
        norm_cys_noncys_peptides = cys_noncys_peptides.copy()
        norm_cys_noncys_peptides[new_sample_cols] = (
            norm_cys_noncys_peptides[new_sample_cols].T / norm_cys_noncys_peptides[control_plex].T).T

        # remove any peptides quantified in less than quant threshold in this replicate
        dfs = []
        for df in [cys_peptides, noncys_peptides, cys_noncys_peptides, noncys_noncys_peptides, norm_cys_noncys_peptides]:
            dfs.append(df.dropna(subset=new_sample_cols, thresh=quant_threshold))
        df_names = ['cys_peptides', 'noncys_peptides', 'cys_noncys_peptides', 'noncys_noncys_peptides', 'control_norm_cys_ratio']
        replicate_normalised[replicate] = dict(zip(df_names, dfs))

        # add remaining dataframes
        replicate_normalised[replicate].update({'raw': raw, 'med_norm': med_norm})

        # Add column with replicate info
        for key, df in replicate_normalised[replicate].items():
            df['replicate'] = replicate
            replicate_normalised[replicate][key].update(df)

    # save to excel
    for replicate in replicates:
        FileHandling.df_to_excel(output_path=f'{output_folder}normalised_replicate_{replicate}.xlsx', sheetnames=list(
            replicate_normalised[replicate].keys()), data_frames=list(replicate_normalised[replicate].values()))

    # compile replicates
    compiled = {}
    for df_type in replicate_normalised[replicate].keys():
        compiled_df = pd.concat([replicate_normalised[replicate][df_type] for replicate in replicates])
        compiled[df_type] = compiled_df

    FileHandling.df_to_excel(output_path=f'{output_folder}normalised_summary.xlsx', sheetnames=list(
            compiled.keys()), data_frames=list(compiled.values()))
