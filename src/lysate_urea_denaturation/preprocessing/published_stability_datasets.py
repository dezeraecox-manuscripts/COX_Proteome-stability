import shutil
import os, re

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

input_folder = 'results/published_stability_correlations/'
output_folder = f'results/lysate_denaturation/published_datasets/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)


# -----------------------------------Copy published published summaries-----------------------------------
shutil.copy2(f'{input_folder}published_stability_datasets.xlsx', f'{output_folder}compiled_dfs.xlsx')
shutil.copy2(f'{input_folder}correlation/dataset_details.xlsx', f'{output_folder}dataset_details.xlsx')

