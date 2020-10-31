
import os
import os, re
import pandas as pd
import numpy as np
import json
import networkx
import obonet
import pandas as pd
import tarfile
import time
from pyorthomap import findOrthologsHsMm, findOrthologsMmHs

from GEN_Utils import FileHandling
from loguru import logger
from utilities.decorators import ProgressBar

logger.info(f'Import OK')

resource_folder = 'resources/bioinformatics_databases/'


def gz_unzipper(filename, input_path=resource_folder, output_path=resource_folder):
    with gzip.open(f'{input_path}{filename}.gz', 'rb') as f_in:
        with open(f'{output_path}{filename}', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def tar_file_to_folder(input_path, output_path):
    tar = tarfile.open(f'{input_path}', 'r')
    tar.extractall(f'{output_path}')
    tar.close()



def go_lineage_tracer(go_term, obo_path, alt_ids=False, direct=False):
    """Return all nodes underneath (i.e. all children) of go_term of interest. 
    Default is to collect entire tree; if stopping at 'direct',
    then only direct decendants are collected."""

    # Read the ontology
    graph = obonet.read_obo(obo_path)
    # Collect all nodes into child:parents dataframe
    children = {}
    for node in graph.nodes:
        try:
            children[node] = graph.nodes[node]['is_a']
        except:
            children[node] = []

    child_terms = pd.DataFrame()
    child_terms['child_term'] = list(children.keys())
    child_terms['parents'] = [';'.join(terms) for terms in list(children.values())]

    # collect any alternate ids for go term
    search_list = []
    if alt_ids:
        try:
            search_list.append(graph.nodes[go_term]['alt_id'])
        except:
            pass
    search_list = [go_term] + [item for sublist in search_list for item in sublist]

    # Collect all children where term of interest is parent
    def search_terms(search_list, term_family=[], direct=False):
        term_family.append(search_list)
        search_list = '|'.join(search_list)
        terms = child_terms[child_terms['parents'].str.contains(search_list)]['child_term'].tolist()
        if direct:
            return terms
        if len(terms) > 0:
            search_terms(terms, term_family)
        return [item for sublist in term_family for item in sublist]

    # collect list of terms of interest
    family = search_terms(search_list, direct=False)

    return family


def uniprot_go_genes(tax_id, go_term, resource_folder=resource_folder, child_terms=True, direct=False, output='list'):
    """Collect all genes from annotated (reviewed) uniprot database containing the GO term of interest.
    tax_id: uniprot id corresponding to saved databases in resources folder e.g. '10090', '9606'
    go_term: term id e.g. 'GO:0032991'
    resource_folder: directory to where stored databases are
    child_terms: default(True) collects all terms for which the term if interest is a parent
    direct: default(false) limits child terms to direct descendents i.e. child term 'is_a' go_term
    output: default(list) choose type of output from 'list' ('Entry' ids), 'df' (complete genes df) or directory (save)"""

    # read in uniprot database for the species, with xref details
    uniprot_database = uniprot_summary(tax_id=tax_id, resource_folder=resource_folder, reviewed=True)
    genes = uniprot_database.dropna(subset=['GO']) # 16525/17474 = 95% have annotated GO terms

    # collect search terms according to optional child terms and direct lineage
    if child_terms:
        search_terms = go_lineage_tracer(go_term, obo_path='PANTHERGOslim.obo', alt_ids=False, direct=direct)
        search_terms = '|'.join(search_terms)
    else:
        search_terms = go_term

    # Collect all genes with ontology_id
    gene_list = genes[genes['GO'].str.contains(search_terms)]

    # generate output
    if output == 'list':
        return gene_list['Entry'].tolist()
    elif output == 'df':
        return gene_list
    else:
        logger.info('Output format not detected. Attempting output to path.')
        gene_list.to_csv(output)


def ontology_wordfinder(words, obo_path='PANTHERGOslim.obo', resource_folder=resource_folder):
    """Retrieves all ontology terms containing 'words' in the name.
    words: list of words to search
    return: df of term, name matches"""

    # Read the ontology
    graph = obonet.read_obo(obo_path)

    # Collect term_ids, term names
    terms = [node for node in graph.nodes]
    names = [graph.nodes[node]['name'] for node in terms]
    terms = pd.DataFrame([terms, names], index=['go_term', 'go_name']).T

    # Collect only terms containing name of interest
    search_words = '|'.join(words)

    return terms[terms['go_name'].str.contains(search_words)]


def go_term_details(go_terms, obo_path='PANTHERGOslim.obo', resource_folder=resource_folder):
    """Retrieves details for all go terms as df.
    go_terms: list of term_ids to search
    return: df of term, name matches"""

    # Read the ontology
    graph = obonet.read_obo(obo_path)

    # Generate df for go_term
    cleaned_terms = []
    for go_term in go_terms:
        try:
            test_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in graph.nodes[go_term].items() ]))
            df = pd.DataFrame([';'.join(test_df[col].dropna()) for col in test_df.columns.tolist()], index=test_df.columns.tolist()).T
            df['go_id'] = go_term
            cleaned_terms.append(df)
        except:
            pass
    cleaned_terms = pd.concat(cleaned_terms)

    return cleaned_terms


def go_uniprot_proteins(protein_names, tax_id, resource_folder=resource_folder, name_type= 'Entry', output='list'):
    """For any given gene, pull out the associated GO terms annotated in uniprot as a list"""

    # read in uniprot database for the species, with xref details
    uniprot_database = uniprot_summary(tax_id=tax_id, resource_folder=resource_folder, reviewed=True)
    gene_details = uniprot_database[uniprot_database[name_type].isin(protein_names)]
    
   # generate output
    if output == 'list':
        return gene_details['GO'].tolist()
    elif output == 'df':
        return gene_details
    else:
        logger.info('Output format not detected')


def taxonomy_id(uniprot_tax_ids, resource_folder):
    species = pd.read_table(f'{resource_folder}orthodb_v10.1/odb10v1_species.tab.gz', compression='gzip', header=None)
    species.columns = ['ncbi_tax_id', 'ortho_tax_id', 'organism_name', 'genome_assembly_id', 'ortho_gene_count', 'ortho_group_count', 'mapping_type']
    species_dict = dict(zip(species['ncbi_tax_id'], species['ortho_tax_id']))

    return [species_dict[ncbi_id] for ncbi_id in uniprot_tax_ids]


@ProgressBar(step=1/41)
def create_genes(resource_folder, ortho_tax_ids):
    try:
        genes = pd.read_excel(f'{resource_folder}{"_".join(ortho_tax_ids)}_genes.xlsx')
        genes.drop([col for col in genes.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
        yield
    except:
        logger.info('File not found. Processing database.')

        gene_chunks = pd.read_table(f'{resource_folder}orthodb_v10.1/odb10v1_genes.tab.gz', compression='gzip', chunksize=1000000, header=None)
        genes = []

        for df in gene_chunks:
            genes.append(df[df[1].isin(ortho_tax_ids)]) # note this relies on the order of the columns - see ortho README
            yield
        genes = pd.concat(genes)
        genes.columns = ['ortho_gene_id', 'ortho_organism_id', 'original_id', 'synonyms', 'mapped_uniprot_id', 'mapped_ensembl_ids', 'ncbi_gene_name', 'mapped_description']

        for tax_id in ortho_tax_ids:
            logger.info(f'{len(genes[genes["ortho_organism_id"] == tax_id])} genes found for {tax_id}')
        
        FileHandling.df_to_excel(f'{resource_folder}{"_".join(ortho_tax_ids)}_genes.xlsx', sheetnames=['all_genes'], data_frames=[genes])

    return genes


@ProgressBar(step=1/196)
def create_og2genes(resource_folder, ortho_tax_ids):
    try:
        og2genes = pd.read_excel(f'{resource_folder}{"_".join(ortho_tax_ids)}_go2genes.xlsx')
        og2genes.drop([col for col in og2genes.columns.tolist() if 'Unnamed: ' in col], axis=1, inplace=True)
        yield
    except:
        logger.info('File not found. Processing database.')

        og2genes_chunks = pd.read_table(f'{resource_folder}orthodb_v10.1/odb10v1_OG2genes.tab.gz', compression='gzip', header=None, chunksize=1000000)

        search_ids = '|'.join(ortho_tax_ids)
        og2genes = []
        for df in og2genes_chunks:
            og2genes.append(df[df[1].str.contains(search_ids)])  # Takes a while!
            yield
        og2genes = pd.concat(og2genes)
        # note this relies on the order of the columns - see ortho README
        og2genes.columns = ['og_id', 'ortho_gene_id']

        FileHandling.df_to_excel(f'{resource_folder}{"_".join(ortho_tax_ids)}_go2genes.xlsx', sheetnames=[search_ids], data_frames=[og2genes])

    return og2genes


@ProgressBar(step=1/16)
def create_uniprot_db(input_path, gene_ids=[]):
    # Testing collect uniprot entries
    from_uniprot_db_chunks = pd.read_table(f'{input_path}', compression='gzip', chunksize=10000)

    from_uniprot_db = []
    for df in from_uniprot_db_chunks:
        if len(gene_ids):
            search_ids = '|'.join(gene_ids)
            from_uniprot_db.append(df[df['Entry'].str.contains(search_ids)])  # Takes a while!
        else:
            from_uniprot_db.append(df)  # Takes a while!
        yield
            
    from_uniprot_db = pd.concat(from_uniprot_db)

    return from_uniprot_db


@ProgressBar(step=0.05)
def create_uniprot_map(input_path, gene_ids=[]):
    # collect ensembl ids from uniprot mapped file
    uniprot_chunks = pd.read_table(f'{input_path}', compression='gzip', chunksize=10000, header=None)

    uniprot_map = []
    for df in uniprot_chunks:
        if len(gene_ids) > 0:
            search_ids = '|'.join(gene_ids)
            uniprot_map.append(df[df[0].str.contains(search_ids)])  # Takes a while!
        else:
            uniprot_map.append(df)  # Takes a while!
        yield

    uniprot_map = pd.concat(uniprot_map)
    uniprot_map.columns = ['UniProtKB-AC', 'UniProtKB-ID', 'GeneID(EntrezGene)', 'RefSeq', 'GI', 'PDB', 'GO', 'UniRef100', 'UniRef90', 'UniRef50', 'UniParc', 'PIR', 'NCBI-taxon', 'MIM', 'UniGene', 'PubMed', 'EMBL', 'EMBL-CDS', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO', 'Additional PubMed']

    return uniprot_map


@ProgressBar(step=0.05)
def create_uniprot_xref(input_path, tax_id, gene_ids=[], id_type=None):
    # collect uniprot mapped file for specific genes
    uniprot_chunks = pd.read_table(f'{input_path}{tax_id}_idmapping.dat.gz', compression='gzip', chunksize=10000, header=None)

    uniprot_xref = []
    for df in uniprot_chunks:
        if len(gene_ids) > 0:
            search_ids = '|'.join(gene_ids)
            uniprot_xref.append(df[df[0].str.contains(search_ids)]) # Takes a while!
        else:
            uniprot_xref.append(df)  # Takes a while!
        yield

    uniprot_xref = pd.concat(uniprot_xref)
    uniprot_xref.columns = ['UniProtKB-AC', 'ID_type', 'ID']

    if id_type:
        uniprot_xref = uniprot_xref[uniprot_xref['ID_type'] == id_type]

    return uniprot_xref


@ProgressBar(step=0.05)
def convert_geneids(gene_ids, tax_id, id_from, id_to, resource_folder=resource_folder):
    # collect uniprot mapped file for specific genes
    uniprot_chunks = pd.read_table(f'{resource_folder}{tax_id}_idmapping.dat.gz', compression='gzip', chunksize=10000, header=None)

    uniprot_xref = []
    for df in uniprot_chunks:
        df.columns = ['UniProtKB-AC', 'ID_type', 'ID']
        search_ids = '|'.join(gene_ids)
        if id_from == 'uniprot':
            df = df[df['UniProtKB-AC'].str.contains(search_ids)] # Takes a while!
        else:
            df = df[df['ID'].str.contains(search_ids)]
        uniprot_xref.append(df)  # Takes a while!
        yield

    uniprot_xref = pd.concat(uniprot_xref)
    if not id_to == 'uniprot':
        uniprot_xref = uniprot_xref[uniprot_xref['ID_type'] == id_to]
    if not id_from == 'uniprot':
        uniprot_xref = uniprot_xref[uniprot_xref['ID_type'] == id_from]

    return uniprot_xref


def uniprot_genename_mapper(uniprot_tax_ids, gene_ids, reviewed=False, orthodb=False, cols=[]):

    # Collect Uniprot info
    compiled = {}
    for direction, tax in uniprot_tax_ids.items():
        if direction == 'from':
            genes = gene_ids
        else:
            genes = []
        uniprot_db = create_uniprot_db(f'{resource_folder}{tax}.tab.gz', genes)
        uniprot_map = create_uniprot_map(f'{resource_folder}{tax}_idmapping.tab.gz', genes)
        merged = pd.merge(uniprot_db, uniprot_map, left_on='Entry', right_on='UniProtKB-AC', how='outer')
        merged['name'] = merged['Entry name'].str.split('_').str[0]
        compiled[direction] = merged

    # Join based on gene name
    to_from_compiled = pd.merge(compiled['from'], compiled['to'], on='name', how='inner', suffixes=('_from', '_to'))

    # annotate matching orthoDB references
    to_from_compiled['orthologous'] = [1 if val_1 == val_2 else 0 for val_1, val_2 in to_from_compiled[[
        'Cross-reference (OrthoDB)_from', 'Cross-reference (OrthoDB)_to']].values]

    # Apply filters
    if reviewed:
            merged = merged[merged['Status'] == 'reviewed']
    if orthodb:
        to_from_compiled = to_from_compiled[to_from_compiled['orthologous'] == 1]
    if cols:
        cols = [[f'{col}_from', f'{col}_to'] for col in cols]
        cols = ['name', 'orthologous'] + [item for sublist in cols for item in sublist]

    # collect unmapped genes
    unmapped_genes = set(gene_ids) - set(to_from_compiled['Entry_from'].tolist())

    return to_from_compiled, unmapped_genes


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def biomart_h2m(gene_list, identifier_type='link_ensembl_gene_id'):
    mapped_ids = []
    for chunk in chunks(gene_list, 100):
        mapped_ids.append(findOrthologsHsMm(from_filters=identifier_type,
                                            from_values=chunk).map())
        time.sleep(10)
    return pd.concat(mapped_ids)


def biomart_m2h(gene_list, identifier_type='link_ensembl_gene_id'):
    mapped_ids = []
    for chunk in chunks(gene_list, 100):
        mapped_ids.append(findOrthologsMmHs(from_filters=identifier_type,
                                            from_values=chunk).map())
        time.sleep(10)
    return pd.concat(mapped_ids)


def ortholog_map(gene_ids, direction, output_folder, resource_folder=resource_folder):

    # Generate mapped ids using orthologsBioMART
    # Collect appropriate id's - note this is not necessarily going to generate unique id's
    gene_identifiers_from = create_uniprot_xref(input_path=f'{resource_folder}9606_idmapping.dat.gz', gene_ids=gene_ids, id_type='Ensembl')
    gene_identifiers_from = dict(zip(gene_identifiers_from['UniProtKB-AC'], gene_identifiers_from['ID']))
    # Collect orthologous genes
    if direction == 'h2m':
        mapped_ids = biomart_h2m(list(gene_identifiers_from.values()), identifier_type='link_ensembl_gene_id')
        uniprot_tax_ids = {'from': '9606', 'to': '10090'}
    elif direction == 'm2h':
        mapped_ids = biomart_m2h(list(gene_identifiers_from.values()), identifier_type='link_ensembl_gene_id')
        uniprot_tax_ids = {'from': '10090', 'to': '9606'}
    # map ensembl back to uniprot
    gene_identifiers_to = create_uniprot_xref(input_path=f'{resource_folder}10090_idmapping.dat.gz', gene_ids=[], id_type='Ensembl')
    gene_identifiers_to = dict(zip(gene_identifiers_to['ID'], gene_identifiers_to['UniProtKB-AC']))
    # Add Uniprot info back to the mapped ids
    mapped_ids['Entry_to'] = mapped_ids['mouse_ensembl_gene_id'].map(gene_identifiers_to)

    # generate mapped ids using UniProt mapper
    cols = ['Entry', 'Entry name', 'Status', 'GeneID(EntrezGene)', 'Ensembl', 'Cross-reference (OrthoDB)']
    compiled, unmapped_genes = uniprot_genename_mapper(
        uniprot_tax_ids, gene_ids, reviewed=True, orthodb=False)


    FileHandling.df_to_excel(
        output_path=f'{output_folder}mapped_ids.xlsx',
        sheetnames=['BioMART_map', 'UniProt_map', 'UniProt_unmapped'],
        data_frames=[mapped_ids, compiled, pd.DataFrame(unmapped_genes)]
    )

    return compiled, mapped_ids


def uniprot_summary(tax_id, resource_folder, genes=[], reviewed=False):
    uniprot_db = create_uniprot_db(f'{resource_folder}{tax_id}.tab.gz', genes)
    uniprot_map = create_uniprot_map(f'{resource_folder}{tax_id}_idmapping.tab.gz', genes)
    merged = pd.merge(uniprot_db, uniprot_map, left_on='Entry', right_on='UniProtKB-AC', how='outer')
    merged['name'] = merged['Entry name'].str.split('_').str[0]
    if reviewed:
        merged = merged[merged['Status'] == 'reviewed']
    return merged

