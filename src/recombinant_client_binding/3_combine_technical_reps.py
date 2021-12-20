import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce
from collections import defaultdict

from loguru import logger
from GEN_Utils import FileHandling, CalcUtils

logger.info('Import OK')


def main(sample_name):
    # Read in raw data, collect useful parts
    compiled = defaultdict(list)
    for replicate in technical_reps:
        ratios = pd.read_excel(f'{input_path}normalised_summary_{sample_name}-{replicate}.xlsx', sheet_name=None)
        for info_type, df in ratios.items():
            df.drop([col for col in df.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
            df['technical_replicate'] = replicate
            compiled[info_type].append(df)
    compiled_ratios = {key: pd.concat(value) for key, value in compiled.items()}

    FileHandling.df_to_excel(
        output_path=f'{output_folder}{sample_name}_compiled.xlsx',
        data_frames=list(compiled_ratios.values()),
        sheetnames=list(compiled_ratios.keys())
    )

    # generate dataset for easy transfer to graphpad
    for info_type, df in compiled_ratios.items():
        info_cols = ['Sequence', 'Proteins', 'replicate', 'cys_rank']
        compiled_ratios.update({info_type: df.groupby([col for col in df.columns.tolist() if col in info_cols]).mean().reset_index()})
    
    
    FileHandling.df_to_excel(
        output_path=f'{output_folder}{sample_name}_mean.xlsx',
        data_frames=list(compiled_ratios.values()),
        sheetnames=list(compiled_ratios.keys())
    )

if __name__ == "__main__":
    input_path = f'results/recombinant_client_assay/peptide_normalisation/'
    output_folder = f'results/recombinant_client_assay/combine_technical_replicates/'

    technical_reps = ['A', 'B']
    sample_names = ['Heat', 'Urea']

    # Define info variables
    info_cols = ['Sequence', 'Proteins', 'Gene names',
                'Protein names', 'replicate', 'cys_rank']
    quant_cols = ['1', '2', '3', '4']

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for sample_name in sample_names:
        main(sample_name)
