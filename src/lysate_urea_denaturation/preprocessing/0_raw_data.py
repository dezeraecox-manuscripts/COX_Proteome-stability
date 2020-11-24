import os, re
import zipfile
from shutil import copyfile
from loguru import logger

from utilities import database_collection
from utilities.database_map_and_filter import gz_unzipper, tar_file_to_folder

logger.info('Import OK')


if __name__ == "__main__":

    url = 'https://zenodo.org/record/4280621/files/lysate_denaturation.zip?download=1' # update to repository address
    folder_name = 'lysate_denaturation'
    output_folder = 'raw_data/'
    resource_folder='resources/bioinformatics_databases/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Download file from repository
    database_collection.download_resources(filename=f'{folder_name}.zip', url=url, resource_folder=output_folder) 
    with zipfile.ZipFile(f'{output_folder}{folder_name}.zip', 'r') as zip_ref:
        zip_ref.extractall(f'{output_folder}')
    