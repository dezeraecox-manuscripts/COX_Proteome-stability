import os, re
import zipfile
from shutil import copyfile
from loguru import logger

from utilities import database_collection
from utilities.database_map_and_filter import gz_unzipper, tar_file_to_folder

logger.info('Import OK')


if __name__ == "__main__":

    url = 'https://' # update to repository address
    folder_name = 'published_stability_correlations/'
    output_folder = 'raw_data/published_stability_correlations/'
    resource_folder='resources/bioinformatics_databases/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Download supplementary info files from various repositories
    # Note some datasets may require login credentials to access beyond journal paywalls. Manually downloaded datasets should be placed in the raw_data/published_stability_correlations/ folder

    data_paths = {
        'Jarzab_2020' : f'https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-020-0801-4/MediaObjects/41592_2020_801_MOESM4_ESM.xlsx',
        'Jarzab_details': f'https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-020-0801-4/MediaObjects/41592_2020_801_MOESM3_ESM.xlsx',
        'Ball_2020': f'https://europepmc.org/articles/PMC7021718/bin/42003_2020_795_MOESM2_ESM.xlsx',
        'Sridharan_2019': 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6411743/bin/41467_2019_9107_MOESM4_ESM.xlsx',
        'Walker_2019': 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6442572/bin/pnas.1819851116.sd04.xlsx', 
        'Miettinen_2018': f'https://www-embopress-org.eu1.proxy.openathens.net/action/downloadSupplement?doi=10.15252%2Fembj.201798359&file=embj201798359-sup-0003-Tableev2.xlsx', 
        'Savitski_2018': f'https://data.mendeley.com/public-files/datasets/8pzhg2tdyb/files/9dabf032-c1a1-4910-a9b7-2a32e0b5aeed/file_downloaded', #unzip --> Fig5_SD6_reference_melting_curves
        'Ogburn_2017A': 'https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.7b00442/suppl_file/pr7b00442_si_003.xlsx', 
        'Ogburn_2017B': 'https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.7b00442/suppl_file/pr7b00442_si_004.xlsx',
        'Leuenberger_2017': 'https://science.sciencemag.org/highwire/filestream/690833/field_highwire_adjunct_files/2/aai7825_Leuenberger_Table-S3.xlsx', 
        'Roberts_2016': 'https://pubs.acs.org/doi/suppl/10.1021/acs.jproteome.6b00927/suppl_file/pr6b00927_si_003.xlsx', 
        'Becher_2016': f'https://static-content.springer.com/esm/art%3A10.1038%2Fnchembio.2185/MediaObjects/41589_2016_BFnchembio2185_MOESM256_ESM.xlsx', 
        'Franken_2015': f'https://static-content.springer.com/esm/art%3A10.1038%2Fnprot.2015.101/MediaObjects/41596_2015_BFnprot2015101_MOESM411_ESM.xlsx', 
        'Savitski_2014A': 'https://science.sciencemag.org/highwire/filestream/595533/field_highwire_adjunct_files/8/Table_S4_Thermal_Profiling_Staurosporine_cell_extract.xlsx',
        'Savitski_2014B': 'https://science.sciencemag.org/highwire/filestream/595533/field_highwire_adjunct_files/8/Table_S3_Thermal_Profiling_ATP_cell_extract.xlsx',
        'dataset_annotations' : ''
    }

    for dataset, url in data_paths.items():
        database_collection.download_resources(filename=f'{url.split("/")[-1]}', url=url,resource_folder=output_folder)
