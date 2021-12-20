import os
import re
from functools import reduce
import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp

from loguru import logger

logger.info('Import OK')

input_folder = 'results/inhibitor_urea_denaturation/initial_cleanup/'
samples = ['VER', 'Control']
pooled_plex = '11'
replicates = ['1', '2', '3']
quant_threshold = 5
urea_conc = {'1': 0.0, '2': 1.0, '3': 2.0, '4': 2.5, '5': 3.0,
             '6': 3.5, '7': 4.0, '8': 4.5, '9': 5.0, '10': 6.0}

output_folder = 'results/inhibitor_urea_denaturation/peptide_normalisation/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


def one_sample_ttest(compiled, sample_cols, group_cols, popmean=1):
    df = compiled.copy()
    ttest_results = []
    for group_key, df in df.groupby(group_cols):
        results = []
        for col in sample_cols:
            test_vals = df[col].values
            if len(test_vals) > 1:
                results.append(tuple(ttest_1samp(test_vals, popmean=popmean, nan_policy='omit')))
            else:
                results.append(tuple([np.nan, np.nan]))
        results = pd.DataFrame(results)
        results.columns = ['t-stat', 'p-val']
        results['sample_col'] = sample_cols
        for x, col in enumerate(group_cols):
            results[col] = group_key[x]
        ttest_results.append(results)
    ttest_results = pd.concat(ttest_results)
    ttest_results[['t-stat', 'p-val']] = ttest_results[['t-stat', 'p-val']].astype(float)

    return ttest_results # 5% of the points detected as significantly different


def z_threshold(data, z_val=1.96):
    """
    Calculate the lower and upper values encompassing a given proportion of the population).

    Common vals:
    ============
    0.5: 38.2%
    1.0: 68.2%
    1.5: 86.6%
    1.96: 95%
    2.0: 95.4%
    2.5: 98.8%
    2.58: 99%
    3.0: 99.8%

    For more explanation: https://upload.wikimedia.org/wikipedia/commons/thumb/2/25/The_Normal_Distribution.svg/1024px-The_Normal_Distribution.svg.png
    """
    range_val = z_val * np.std(data)
    return np.mean(data) - range_val, np.mean(data) + range_val


def pval_smoothing(compiled, sample_cols, group_cols, popmean, penalty_factor=20, complete=False):
    """Scale mean value proportional to pvalue, imposing penalty for variability

    Parameters
    ----------
    compiled : DataFrame
        Longoform pandas df containing descriptive columns (group_cols) and data columns (sample_cols),
        where replicates of each datapoint are stored in columns.
    sample_cols : list[str]
        List of column names where quantitative data can be found. Replicate data points should 
        be contained wholly within single columns
    group_cols : list[str]
        List of column names to group ```compiled``` of, such that grouped df for each group is length of replicates 
    popmean : int
        Hypothesised population mean. Typically, for ratiometric analyses this may be 1 or 0, however 
        can be any value to which samples will be compared
    penalty_factor : int, optional
        Weight to which p-value will be scaled, by default 20. Larger value imposes more severe 
        scaling of the mean value with increased p-value.

    Returns
    -------
    DataFrame
        Smoothed dataframe where replicates have been reduced to the mean value, 
        scaled by p-value smoothing.
    """    
    # Apply t-test to sample
    ttest_results = one_sample_ttest(compiled, sample_cols, group_cols=group_cols, popmean=popmean)
    # Generate scaling factors
    ttest_results['exp_p-val'] = penalty_factor**ttest_results['p-val']
    p_vals = pd.pivot_table(ttest_results, values='exp_p-val', index=group_cols, columns='sample_col')
    p_vals.columns = [f'scalefactor_{col}' for col in p_vals.columns]

    # Calculate mean of the input df
    proportional_pval = compiled.groupby(group_cols).mean()[sample_cols].copy().sort_values(group_cols)
    proportional_pval.columns = [f'mean_{col}' for col in proportional_pval.columns]

    proportional_pval = pd.merge(proportional_pval, p_vals, on=group_cols, how='outer')
    for col in sample_cols:
        proportional_pval[f'scaled_{col}'] = popmean + (proportional_pval[f'mean_{col}'] - popmean) * (1 / proportional_pval[f'scalefactor_{col}'])
    
    # Restoring legacy function to return only scaled values matching input compiled
    smoothed_vals = proportional_pval[[col for col in proportional_pval.columns if 'scaled' in col]].copy()
    smoothed_vals.columns = [col.replace('scaled_', '') for col in smoothed_vals.columns]

    if complete:
        return proportional_pval
    else:
        return smoothed_vals


# ----------------------------------------------------------------------
# Compile all samples to single df
raw_peptide_results = []
for sample_name in samples:

    raw_data = pd.read_excel(
        f'{input_folder}{sample_name}_Compiled.xlsx', sheet_name=None)

    info_cols = ['Sequence', 'Proteins', 'Gene names', 'Protein names', 'cys_rank', 'replicate']
    peptides = raw_data['Peptides'].copy().set_index([col for col in raw_data['Peptides'].columns.tolist() if col in info_cols]).drop(
        ['Unique (Groups)', 'Unique (Proteins)'], axis=1)

    # Adjust column names to be in format label_replicate
    peptides.columns = [('_').join(re.split(r' |_', col)[3:6])
                        for col in peptides.columns.tolist()]

    raw_peptide_results.append(peptides)

merged_peptides = reduce(lambda left, right: pd.merge(left, right, on=['Sequence', 'Proteins', 'Gene names', 'Protein names'], how='outer'), raw_peptide_results)

# Complete sum normalisation for total peptide abundance
scaling_factor = merged_peptides.sum().max() / merged_peptides.sum()
scaled_peptides = (merged_peptides * scaling_factor)

# Melt into longform df
peptides = pd.melt(
    scaled_peptides.reset_index(),
    id_vars=['Sequence', 'Proteins', 'Gene names', 'Protein names'],
    value_vars=scaled_peptides.columns.tolist(),
    value_name='abundance',
    var_name=['sample']
)
peptides[['channel', 'treatment', 'replicate']] = peptides['sample'].str.split('_', expand=True)

# Calculate VER/Control
ver_ratio = pd.pivot_table(
    data=peptides.copy(),
    index=['Sequence', 'Proteins', 'channel', 'replicate'],
    values='abundance',
    columns=['treatment'],
).reset_index()

ver_ratio['VER/Control'] = ver_ratio['VER'] / ver_ratio['Control']
ver_ratio.dropna(subset=['VER/Control'], inplace=True)

# filter for peptides quantified in two or more replicates in that treatment/channel
replicate_filtered = []
for channel, df in ver_ratio.groupby(['channel']):
    replicate_counts = df.groupby('Sequence').count()[
        'VER/Control'].reset_index()
    sequences = replicate_counts[replicate_counts['VER/Control']
                                 > 1]['Sequence'].tolist()
    replicate_filtered.append(df[df['Sequence'].isin(sequences)])
ver_ratio = pd.concat(replicate_filtered)

# Calculate cys/noncys ratio, apply to original table then convert to log
ver_ratio_noncys_av = ver_ratio[~ver_ratio['Sequence'].str.contains(
    'C')].copy().groupby(['Proteins', 'channel', 'replicate']).mean()['VER/Control'].reset_index()
ver_ratio = reduce(lambda left, right: pd.merge(left, right, on=['Proteins', 'channel', 'replicate'], how='outer'), [ver_ratio, ver_ratio_noncys_av.rename(columns={'VER/Control': 'noncys_VER/Control_av'})])
ver_ratio['corrected_VER/Control'] = ver_ratio['VER/Control'] / ver_ratio['noncys_VER/Control_av']
# remove entries for which no noncys was available (only drops cys peptides)
ver_ratio.dropna(subset=['noncys_VER/Control_av'], inplace=True)
# Remove any cys peptides that are now quantified in < 2 replicates
ver_ratio = reduce(lambda left, right: pd.merge(left, right, on=['Sequence', 'channel'], how='outer'), [ver_ratio, ver_ratio.groupby(['Sequence', 'channel']).count()['Proteins'].reset_index().rename(columns={'Proteins': 'rep_count'})])
ver_ratio = ver_ratio[ver_ratio['rep_count'] > 1].copy()
ver_ratio.drop('rep_count', axis=1, inplace=True)

# take log 2
ver_ratio['log2_corrected_VER/Control'] = np.log2(ver_ratio['corrected_VER/Control'])

# perform pVal smoothing
smooth_ver_ratio = pval_smoothing(
    ver_ratio.copy(),
    sample_cols=['log2_corrected_VER/Control'],
    group_cols=['Sequence', 'Proteins', 'channel'], 
    popmean=0,
    penalty_factor=20,
    complete=True
    )

# For single-noncys proteins, corrected val with be all 0 (thus Nan after smoothing)
# Note this could also happen for peptides only identified in 1 replicate
# however these are filtered out above
smooth_ver_ratio['scaled_log2_corrected_VER/Control'] = [val if not np.isnan(val) else (0 if mean_val == 0.0 else np.nan) for val, mean_val in smooth_ver_ratio[['scaled_log2_corrected_VER/Control', 'mean_log2_corrected_VER/Control']].values]

# remove pooled samples
smooth_ver_ratio['urea'] = list(smooth_ver_ratio.reset_index()['channel'].astype(str).map(urea_conc))
smooth_ver_ratio.dropna(subset=['urea'], inplace=True)
smooth_ver_ratio.reset_index(inplace=True)


# save to csv
smooth_ver_ratio.to_csv(f'{output_folder}peptide_ratio_summary.csv')

