from GEN_Utils.CalcUtils import sorted_nicely
from scipy.stats import ttest_1samp, ttest_ind
from scipy.optimize import curve_fit
import re
import os
import functools
import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

input_folder = 'results/inhibitor_urea_denaturation/peptide_normalisation/'
output_folder = 'results/inhibitor_urea_denaturation/detect_outliers/'


if not os.path.exists(output_folder):
    os.makedirs(output_folder)


# Read in filtered results
info_cols = ['Sequence', 'Proteins', 'treatment']
smooth_ratios = pd.read_csv(f'{input_folder}peptide_ratio_summary.csv')
smooth_ratios.drop([col for col in smooth_ratios.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)

# Fix data types
smooth_ratios['channel'] = smooth_ratios['channel'].astype(int)
smooth_ratios['urea'] = smooth_ratios['urea'].astype(float)

# determine noncys deviation per protein across all channels, map onto df
noncys_band = smooth_ratios[~smooth_ratios['Sequence'].str.contains('C')].groupby(
    ['Proteins']).agg(['mean', 'std'])[['scaled_log2_corrected_VER/Control']].reset_index()
noncys_band.columns = [f'{i}_{j}' if j != '' else f'{i}' for i,j in noncys_band.columns]

smooth_ratios['noncys_ratio_std'] = smooth_ratios['Proteins'].map(
    dict(noncys_band[['Proteins', 'scaled_log2_corrected_VER/Control_std']].values))
smooth_ratios['noncys_ratio_mean'] = smooth_ratios['Proteins'].map(
    dict(noncys_band[['Proteins', 'scaled_log2_corrected_VER/Control_mean']].values))

# determine outliers as proteins for whom VER/Control for at least one peptide outside non-cys deviation
smooth_ratios['outlier'] = [1 if abs(val) >= thresh else 0 for val, thresh in smooth_ratios[['scaled_log2_corrected_VER/Control', 'noncys_ratio_std']].values]

# Unmelt
smoothed_summary = pd.pivot_table(
    smooth_ratios,
    index=['Proteins', 'Sequence', 'outlier', 'noncys_ratio_std', 'noncys_ratio_mean'],
    values='scaled_log2_corrected_VER/Control',
    columns='urea'
)

smoothed_summary['count'] = smoothed_summary.count(axis=1)
smoothed_summary.reset_index(inplace=True)


# collect only cys peptides quantified in at least 7 out of 10 channels
outlier_summary = smoothed_summary[(smoothed_summary['count'] > 6) & (
    smoothed_summary['Sequence'].str.contains('C'))].copy()

# collect outliers
outliers = outlier_summary[outlier_summary['outlier'] == 1].copy()
outliers['noncys_ratio_std'] = outliers['noncys_ratio_std'].replace(
    0, np.nan)  # remove outiers only detected against a 0 std background
outliers.dropna(subset=['noncys_ratio_std'], inplace=True)


# Save to excel
FileHandling.df_to_excel(
    data_frames=[smooth_ratios.reset_index(), smoothed_summary, outliers],
    sheetnames=['smooth_ratios', 'summary', 'outliers'],
    output_path=f'{output_folder}outlier_summary.xlsx'
    )
