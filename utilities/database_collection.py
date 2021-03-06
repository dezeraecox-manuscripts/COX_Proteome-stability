from pyorthomap import FindOrthologs
from pybiomart import Server
import gzip, shutil, os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pyorthomap import findOrthologsHsMm, findOrthologsMmHs
import time
import requests
import xmltodict
from collections import OrderedDict
import json
from contextlib import closing
import urllib.request as request
import time

from utilities.decorators import ProgressBar
from utilities.database_map_and_filter import convert_geneids, gz_unzipper

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')


def download_resources(filename, url, resource_folder):
    """
    Worker function to download and save file from URL.
    
    inputs
    ======
    filename: (str) name of output file (including extension)
    url: (str) complete location of file to be downloaded
    output_path: (str) relative or complete path to directory where folder will be saved.

    returns:
    ======
    None

    """
    if not os.path.exists(resource_folder):
        os.makedirs(resource_folder)

    try:
        with closing(request.urlopen(url)) as r:
            with open(f'{resource_folder}{filename}', 'wb') as f:
                shutil.copyfileobj(r, f)
        logger.info(f'Downloaded {filename}')
    except:
        logger.info(f'Downloaded failed for {filename}.')


# ----------------------------------Collect/generate static databases--------------------------

# ASA from PyMol
def retrieve_asa(protein, chains=False, resource_folder='resources/bioinformatics_databases/'):
    
    """
    Using pymol command to calculate relative per-residue solvent accessible surface area https://pymol.org/pymol-command-ref.html#get_sasa_relative
    
    Calculates the relative per-residue solvent accessible surface area and optionally labels and colors residues. The value is relative to full exposure of the residue, calculated by removing all other residues except its two next neighbors, if present.
    
    Loads a value beteween 0.0 (fully buried) and 1.0 (fully exposed) into the b-factor property, available in "iterate", "alter" and "label" as "b".
    """
    output_folder=f'{resource_folder}pymol_asa/'
    
    if not os.path.exists(f'{output_folder}'):
        os.mkdir(f'{output_folder}')

    cmd.delete('all')
    structure = cmd.fetch(protein)
    if not chains:
        chains = cmd.get_chains()
    chain_info = {}
    for chain in chains:
        fasta = cmd.get_fastastr(f"chain {chain}").replace(f'>{protein}_{chain}', '').replace('\n', '').replace('?', '')
        test_asa =  cmd.get_sasa_relative(f"chain {chain}")
        keys = pd.DataFrame([key for key in list(test_asa.keys())])
        values = pd.DataFrame([value for value in list(test_asa.values())])
        asa = pd.merge(keys, values, left_index=True, right_index=True)
        asa.columns = ['pdb', 'aa', 'chain', 'position', 'asa']
        chain_info[chain] = asa
    summary = pd.concat(list(chain_info.values()))
    summary.to_csv(f'{output_folder}{protein}.csv')

    return summary

def asa_calculator(pdb_info, resource_folder='resources/bioinformatics_databases/'):

    pymol_asa = []
    nonmapped = []
    for structure, seqs in pdb_info.groupby('pdb_id'):
        chains = seqs['chain'].unique().tolist()
        try:
            pymol_asa.append(retrieve_asa(structure, chains, resource_folder))
            yield
        except:
            nonmapped.append(structure)
    asa_summary = pd.concat(pymol_asa)
    asa_summary.to_csv(f'{resource_folder}asa_summary.csv')
    return asa_summary


# ----------------------------------API/live search functions----------------------------------

def uniprot_download(genes, output_folder):

    if not os._exists(output_folder):
        os.makedirs(output_folder)

    for gene in genes:
        url = f'https://www.uniprot.org/uniprot/{gene}.xml'
        response = requests.get(url)
        with open(f'{output_folder}{gene}.xml', 'w') as f:
            f.write(response.text)


def uniprot_features(genes, resource_folder='resources/bioinformatics_databases/'):

    def collect_feature(feature_dict, feature_name):
        try:
            return feature_dict[feature_name]
        except:
            return np.nan

    if not os.path.exists(f'{resource_folder}uniprot_download/'):
        os.makedirs(f'{resource_folder}uniprot_download/')

    protein_features = {}
    genes_unprocessed = []
    for gene in genes:
        try:
            with open(f'{resource_folder}uniprot_download/{gene}.xml') as fd:
                gene_dict = xmltodict.parse(fd.read())
        except FileNotFoundError:
            uniprot_download(genes, f'{resource_folder}uniprot_download/')
            with open(f'{resource_folder}uniprot_download/{gene}.xml') as fd:
                gene_dict = xmltodict.parse(fd.read())
        try:
            feature_names = ['location', '@type', '@description', '@evidence']
            new_features = []
            for feature in gene_dict['uniprot']['entry']['feature']:
                new_dict = {}
                for feature_name in feature_names:
                    new_dict[feature_name] = collect_feature(feature, feature_name)
                new_features.append(new_dict)
            features = pd.DataFrame(new_features)
            features = pd.merge(features[['@type', '@description', '@evidence']], features['location'].apply(pd.Series), left_index=True, right_index=True)
            for location_col in ['begin', 'end', 'position']:
                if location_col in features.columns.tolist():
                    features[f'{location_col}_key'] = [list(pos.keys())[0] if type(pos) == OrderedDict else np.nan for pos in features[location_col] ]
                    features = features[~features[f'{location_col}_key'].astype(str).str.contains('@status')]
                    features[location_col] = [int(pos['@position']) if type(pos) == OrderedDict else np.nan for pos in features[location_col] ]
            if len(features) < 1:
                raise Exception(f"No features found for {gene}.")
            features['entry'] = gene
            cols = ['@type', '@description', '@evidence', 'begin', 'end', 'position', 'entry']
            protein_features[gene] = features[[col for col in features.columns.tolist() if col in cols]].rename(columns={'@type': 'feature_type', '@description': 'feature_description', '@evidence': 'evidence'})
        except:
            logger.info(f'{gene} not processed.')
            genes_unprocessed.append(gene)
    
    return protein_features, genes_unprocessed


def uniprot_function(genes, resource_folder='resources/bioinformatics_databases/'):

    protein_function = {}
    genes_unprocessed = []
    for gene in genes:
        try:
            with open(f'{resource_folder}uniprot_download/{gene}.xml') as fd:
                gene_dict = xmltodict.parse(fd.read())
        except FileNotFoundError:
            uniprot_download(genes, f'{resource_folder}uniprot_download/')
            with open(f'{resource_folder}uniprot_download/{gene}.xml') as fd:
                gene_dict = xmltodict.parse(fd.read())
        try:
            protein_function[gene] = {words['@id'] : words['#text'] for words in gene_dict['uniprot']['entry']['keyword']}
        except:
            logger.info(f'{gene} not processed.')
            genes_unprocessed.append(gene)

    return protein_function


def uniprot_sequences(genes, resource_folder='resources/bioinformatics_databases/'):

    protein_sequence = {}
    genes_unprocessed = []
    for gene in genes:
        try:
            with open(f'{resource_folder}uniprot_download/{gene}.xml') as fd:
                gene_dict = xmltodict.parse(fd.read())
        except:
            uniprot_download(genes, f'{resource_folder}uniprot_download/')
            with open(f'{resource_folder}uniprot_download/{gene}.xml') as fd:
                gene_dict = xmltodict.parse(fd.read())
        try:
            protein_sequence[gene] = gene_dict['uniprot']['entry']['sequence']['#text']
        except:
            logger.info(f'{gene} not processed.')
            genes_unprocessed.append(gene)

    return protein_sequence



def string_geneid_mapper(genes, tax_id, id_limit=1, echo_query=1, string_api_url="https://string-db.org/api", output_format='tsv', caller_id="www.github.com/dezeraecox"):

    # set parameters
    method = "get_string_ids"
    params = {
        "identifiers" : "\r".join(genes), # your protein list
        "species" : tax_id, # species NCBI identifier 
        "limit" : id_limit, # only one (best) identifier per input protein
        "echo_query" : echo_query, # see your input identifiers in the output
        "caller_identity" : caller_id # your app name
    }

    # Construct URL
    request_url = "/".join([string_api_url, output_format, method])
    # Call STRING
    response = requests.post(request_url, data=params)

    # format response into dataframe
    results = pd.DataFrame(response.text.strip().split("\n"))
    results = results[0].str.strip().str.split("\t", expand=True).T.set_index(0).T
    results['string'] = results['stringId'].str.split('.').str[1]

    return results


def network_interactions(genes, tax_id, id_type='string', string_api_url="https://string-db.org/api", output_format='tsv', caller_id="www.github.com/dezeraecox"):

    """
    id_type: STRING=premapped ids to Esembl, UNIPROT ids will be premapped to STRING
    and converted back after interaction mapping"""

    if not id_type == 'string':
        try:
            id_map = string_geneid_mapper(genes, tax_id, id_limit=1, echo_query=1)
            genes = id_map['string'].tolist()
        except:
            logger.info('Unable to map genes to STRING format.')
    
    # set parameters
    method = "network"
    params = {
        "identifiers" : "%0d".join(genes), # your proteins
        "species" : tax_id, # species NCBI identifier 
        "caller_identity" : caller_id # your app name
    }

    # Construct URL
    request_url = "/".join([string_api_url, output_format, method])
    # Call STRING
    response = requests.post(request_url, data=params)

    results = pd.DataFrame(response.text.strip().split("\n"))
    results = results[0].str.strip().str.split("\t", expand=True).T.set_index(0).T
    if not id_type == 'string':
        results['originalId_A'] = results['stringId_A'].map(dict(zip(id_map['stringId'], id_map['queryItem'])))
        results['originalId_B'] = results['stringId_B'].map(dict(zip(id_map['stringId'], id_map['queryItem'])))

    return results


def all_interactions( genes, tax_id, max_partners=1000, id_type='string', string_api_url="https://string-db.org/api", output_format='tsv', caller_id="www.github.com/dezeraecox", confidence_threshold=0.7):

    """
    id_type: STRING=premapped ids to Esembl, UNIPROT ids will be premapped to STRING
    and converted back after interaction mapping"""

    if not id_type == 'string':
        try:
            id_map = string_geneid_mapper(genes, tax_id, id_limit=1, echo_query=1)
            genes = id_map['string'].tolist()
        except:
            logger.info('Unable to map genes to STRING format.')
    
    def interactors_search(genes):
        # set parameters
        method = "interaction_partners"
        params = [
            f'identifiers={"%0D".join(genes)}', # your proteins
            f'species={tax_id}', # species NCBI identifier 
            f'limit={max_partners}',
            f'caller_identity={caller_id}' # your app name
        ]

        # Construct URL
        request_url = "/".join([string_api_url, output_format, method])
        request_url = f'{request_url}?{"&".join(params)}'

        # Call STRING
        response = requests.get(request_url)

        results = pd.DataFrame(response.text.strip().split("\n"))
        results = results[0].str.strip().str.split("\t", expand=True).T.set_index(0).T
    
        return results
    
    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    interactions = []
    for number, gene_list in enumerate(chunks(genes, 10)):
        try:
            try:
                interactions.append(interactors_search(gene_list))
            except:
                time.sleep(5)
                interactions.append(interactors_search(gene_list))
            logger.info(f'Processed chunk number {number}')
        except:
            logger.info(f'{gene_list} not processed.')    

    interactions = pd.concat(test)
    interactions = interactions[interactions['score'].astype(float) >= confidence_threshold]
    interactions = interactions[interactions['stringId_A'].str.strip(f'{tax_id}.').isin(genes) | interactions['stringId_B'].str.strip(f'{tax_id}.').isin(genes)]

    genes_to_map = interactions['stringId_A'].unique().tolist() + interactions['stringId_B'].unique().tolist()
    gene_map = convert_geneids(genes_to_map, tax_id=10090, id_from='STRING', id_to='uniprot')
    gene_map['string_id'] = gene_map['ID'].str.split('.').str[1]
    gene_map = dict(zip(gene_map['string_id'], gene_map['UniProtKB-AC']))
    gene_map.update(dict(zip(id_map['string'], id_map['queryItem'])))

    interactions['uniprot_A'] = interactions['stringId_A'].str.strip(f'{tax_id}.').map(gene_map)
    interactions['uniprot_B'] = interactions['stringId_B'].str.strip(f'{tax_id}.').map(gene_map)
    
    return interactions


def interaction_enrichment( genes, tax_id, id_type='string', string_api_url="https://string-db.org/api", output_format='tsv', caller_id="www.github.com/dezeraecox"):

    if not id_type == 'string':
        try:
            id_map = string_geneid_mapper(genes, tax_id, id_limit=1, echo_query=1)
            genes = id_map['string'].tolist()
        except:
            logger.info('Unable to map genes to STRING format.')
    # set parameters
    method = "ppi_enrichment"
    params = {
    "identifiers" : "%0d".join(genes), # your proteins
    "species" : tax_id, # species NCBI identifier 
    "caller_identity" : caller_id # your app name
    }

    # Construct URL
    request_url = "/".join([string_api_url, output_format, method])
    # Call STRING
    response = requests.post(request_url, data=params)

    results = pd.DataFrame(response.text.strip().split("\n"))
    results = results[0].str.strip().str.split("\t", expand=True).T.set_index(0).T

    return results


def disorder_prediction_iupred(accession_ids, output_folder):

    
    disorder = []
    for x, accession in enumerate(accession_ids):
        if x % 10 == 0:
            logger.info(f'Processing protein number {x}')
        # Call IUPreD Database
        request_url = f'http://iupred2a.elte.hu/iupred2a/{accession}.json'

        response = requests.get(request_url)
        results = json.loads(response.text)

        # format disorder into dataframe
        results = pd.DataFrame([results['iupred2'], results['sequence']], index=['disorder_probability', 'sequence']).T.reset_index().rename(columns={'index': 'residue'})
        results['disorder_probability'] = results['disorder_probability'].astype(float)
        results['disordered'] = [1 if probability > 0.5 else 0 for probability in results['disorder_probability']]
        results['Proteins'] = accession
        disorder.append(results)

    return pd.concat(disorder)


def pfam_domains(accession_ids, resource_folder='resources/bioinformatics_databases/'):

    if not os.path.exists(f'{resource_folder}pfam_domains/'):
        os.makedirs(f'{resource_folder}pfam_domains/')

    pfam = []
    for x, accession in enumerate(accession_ids):
        if x % 10 == 0:
            logger.info(f'Processing protein number {x}')
        if not os.path.exists(f'{resource_folder}pfam_domains/{accession}.xml'):
            # Call Database
            request_url = f'https://pfam.xfam.org/protein/{accession}?output=xml'
            response = requests.get(request_url)
            with open(f'{resource_folder}pfam_domains/{accession}.xml', 'w') as f:
                f.write(response.text)

        with open(f'{resource_folder}pfam_domains/{accession}.xml') as fd:
            gene_dict = xmltodict.parse(fd.read())
        try:
            # format disorder into dataframe
            matches = gene_dict['pfam']['entry']['matches']['match']
            if type(matches) == list:
                for match in matches:
                    pfam_acc = match['@accession']
                    pfam_id = match['@id']
                    pfam_type = match['@type']
                    start, end = (match['location']['@start'], match['location']['@end'])
                    pfam.append(pd.DataFrame([accession, pfam_acc, pfam_id, pfam_type, start, end], index=['Proteins', 'pfam_acc', 'pfam_id', 'pfam_type', 'pfam_start', 'pfam_end']).T)
            else:
                pfam_acc = matches['@accession']
                pfam_id = matches['@id'] 
                pfam_type = matches['@type'] 
                start, end = (matches['location']['@start'], matches['location']['@end'])
                pfam.append(pd.DataFrame([accession, pfam_acc, pfam_id, pfam_type, start, end], index=['Proteins', 'pfam_acc', 'pfam_id', 'pfam_type', 'pfam_start', 'pfam_end']).T)
        except:
            logger.info(f'{accession} not processed.')

    return pd.concat(pfam).reset_index(drop=True)


def main(tax_ids={'MOUSE': '10090', 'HUMAN': '9606'}, resource_folder='resources/bioinformatics_databases/', caller_id = "www.github.com/dezeraecox", fasta_codes={'10090': 'UP000000589', '9606': 'UP000005640'}):

    # download useful databases if they don't already exist
    string_api_url = "https://string-db.org/api"
    output_format = "tsv"

    if not os.path.exists(resource_folder):
        os.makedirs(resource_folder)

    for species, tax_id in tax_ids.items():
        if not os.path.exists(f'{resource_folder}{tax_id}_idmapping.tab.gz'):
            download_resources(
                filename=f'{tax_id}_idmapping.tab.gz',
                url=f'https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/{species}_{tax_id}_idmapping_selected.tab.gz',
                resource_folder=resource_folder)
        gz_unzipper(f'{tax_id}_idmapping.tab', input_path=resource_folder, output_path=resource_folder)


        if not os.path.exists(f'{resource_folder}{tax_id}.tab.gz'):
            download_resources(
                filename=f'{tax_id}.tab.gz',
                url=f'https://www.uniprot.org/uniprot/?query={tax_id}&%28{species}%29+%5B%22&fil=&offset=0&compress=yes&format=tab',
                resource_folder=resource_folder)
        gz_unzipper(f'{tax_id}.tab', input_path=resource_folder, output_path=resource_folder)
    
        if not os.path.exists(f'{resource_folder}{tax_id}.fasta.gz'):
            code = fasta_codes[tax_id]
            download_resources(
                filename=f'{tax_id}.fasta.gz',
                url=f'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/{code}_{tax_id}.fasta.gz',
                resource_folder=resource_folder)
        gz_unzipper(f'{tax_id}.fasta', input_path=resource_folder, output_path=resource_folder)
    
        if not os.path.exists(f'{resource_folder}{tax_id}_idmapping.dat.gz'):
            download_resources(
                filename=f'{tax_id}_idmapping.dat.gz',
                url=f'https://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/{species}_{tax_id}_idmapping.dat.gz',
                resource_folder=resource_folder)
                

    if not os.path.exists(f'{resource_folder}PANTHERGOslim.obo'):
        download_resources(
            filename=f'PANTHERGOslim.obo',
            url=f'http://data.pantherdb.org/PANTHER15.0/ontology/PANTHERGOslim.obo',
            resource_folder=resource_folder)

    if not os.path.exists(f'{resource_folder}ss_dis.txt'):
        download_resources(
            filename=f'ss_dis.txt.gz',
            url=f'https://cdn.rcsb.org/etl/kabschSander/ss_dis.txt.gz',
            resource_folder=resource_folder)
        gz_unzipper('ss_dis.txt', input_path=resource_folder, output_path=resource_folder)
        
    if not os.path.exists(f'{resource_folder}pdb_seqres.txt'):
        download_resources(
            filename=f'pdb_seqres.txt.gz',
            url=f'https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz',
            resource_folder=resource_folder)
        gz_unzipper('pdb_seqres.txt', input_path=resource_folder, output_path=resource_folder)

    if not os.path.exists(f'{resource_folder}HOM_MouseHumanSequence.txt'):
        download_resources(
            filename=f'HOM_MouseHumanSequence.txt',
            url=f'http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt',
            resource_folder=resource_folder)

    if not os.path.exists(f'{resource_folder}hmmer_pdb_all.txt'):
        download_resources(
            filename=f'hmmer_pdb_all.txt',
            url=f'https://www.rcsb.org/pdb/rest/hmmer?file=hmmer_pdb_all.txt',
            resource_folder=resource_folder)    
        

if __name__ == "__main__":
    main()




