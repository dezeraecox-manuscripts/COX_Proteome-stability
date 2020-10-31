import os, re
import zipfile
from shutil import copyfile
from loguru import logger

from utilities import database_collection
from utilities.database_map_and_filter import gz_unzipper, tar_file_to_folder

logger.info('Import OK')


if __name__ == "__main__":

    url = 'https://' # update to repository address
    folder_name = 'lysate_urea_denaturation/'
    output_folder = 'raw_data/lysate_urea_denaturation/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Download file from repository
    # database_collection.download_resources(filename=f'{folder_name}.zip', url=url, resource_folder=output_folder) 
    # with zipfile.ZipFile(f'{output_folder}{folder_name}.zip', 'r') as zip_ref:
    #     zip_ref.extractall(f'{output_folder}{folder_name}')

    # Current solution 
    input_folder = 'C:/Users/Dezerae/Documents/Current Writing/180701_Manuscipt_Urea denaturation in lysate/data_and_analysis/experiments/200831_Analysis_Experiment 20 Refactoring/raw_data/'
    output_folder = 'raw_data/lysate_urea_denaturation/'
    file_list = [filename for filename in os.listdir(input_folder) if '.xlsx' in filename]
    for filename in file_list:
        copyfile(f'{input_folder}{filename}', f'{output_folder}{filename}')

    # Download useful databases to resources folder
    database_collection.main(tax_ids=['10090', '9606'], resource_folder='resources/bioinformatics_databases/', caller_id = "www.github.com/dezeraecox")
    


