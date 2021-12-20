""" Test out pval correction and alternate normalisation methods"""
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from GEN_Utils import FileHandling
from scipy.stats import ttest_1samp

from loguru import logger
logger.info('Import OK')


if not os.path.exists(output_folder):
    os.mkdir(output_folder)


def one_sample_ttest(compiled, sample_cols, popmean=1):
    df = compiled.copy()
    ttest_results = []
    for sequence, df in df.groupby('Sequence'):
        results = []
        for col in sample_cols:
            test_vals = df[col].values
            if len(test_vals) > 1:
                results.append(
                    tuple(ttest_1samp(test_vals, popmean=popmean, nan_policy='omit')))
            else:
                results.append(tuple([np.nan, np.nan]))
        results = pd.DataFrame(results)
        results.columns = ['t-stat', 'p-val']
        results['conc'] = sample_cols
        results['Sequence'] = sequence
        ttest_results.append(results)
    ttest_results = pd.concat(ttest_results)
    ttest_results[['t-stat', 'p-val']
                  ] = ttest_results[['t-stat', 'p-val']].astype(float)

    return ttest_results  # 5% of the points detected as significantly different

def main(sample_name, norm_cols):
    # normalise denatured MDH sample cys_noncys values to corresponding native sample

    # read in raw data for technical replicate means
    raw_data = pd.read_excel(f'{input_folder}/{sample_name}_mean.xlsx', sheet_name='cys_noncys')
    compiled = raw_data.copy().drop([col for col in raw_data.columns.tolist() if 'Unnamed: ' in col], axis=1)
    compiled.drop('cys_rank', inplace=True, axis=1)

    info_cols = ['Sequence', 'Proteins', 'replicate', 'coverage']

    native_norm = compiled.copy()
    for col_1, col_2 in norm_cols.items():
        native_norm[col_1] = native_norm[col_1] / native_norm[col_2]

    # Collect normalised columns
    native_norm = native_norm[['Sequence', 'Proteins', 'replicate']+list(norm_cols.keys())]
    native_norm = native_norm.sort_values(['Proteins', 'Sequence', 'replicate'])

    # test for difference from 1
    native_norm_ttest = one_sample_ttest(native_norm, list(norm_cols.keys()), popmean=1)
    native_norm_ttest.rename(columns={'conc': 'sample'}, inplace=True)
    native_norm_ttest['sign'] = ['*' if pval < 0.05 else np.nan for pval in native_norm_ttest['p-val']]
    native_norm_ttest['sign'] = ['**' if pval < 0.01 else sign for pval, sign in native_norm_ttest[['p-val', 'sign']].values]
    native_norm_ttest['sign'] = ['***' if pval < 0.001 else sign for pval, sign in native_norm_ttest[['p-val', 'sign']].values]

    FileHandling.df_to_excel(
    output_path=f'{output_folder}{sample_name}_native_normed.xlsx',
    sheetnames=['native_norm', 'onesample_ttest'],
    data_frames=[native_norm, native_norm_ttest]
    )

if __name__ == "__main__":
    
    input_folder = f'results/recombinant_client_assay/combine_technical_replicates/'
    output_folder = f'results/recombinant_client_assay/native_normalisation/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    sample_names = ['Heat', 'Urea']
    norm_cols = { #in format denatured/native
        '2': '1',
        '4': '3',
    }

    for sample_name in sample_names:
        main(sample_name, norm_cols)