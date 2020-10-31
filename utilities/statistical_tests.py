import functools
import gzip
# Investigate the sequences belonging to each cluster
# Attempt enrichment analysis
import os
import re
import shutil
import time
from collections import OrderedDict, defaultdict

import matplotlib.pyplot as plt
import networkx
import numpy as np
import obonet
import pandas as pd
import requests
import scikit_posthocs as sp
import seaborn as sns
import statsmodels.api as sa
import statsmodels.formula.api as sfa
import xmltodict
from bs4 import BeautifulSoup
from GEN_Utils import FileHandling
from loguru import logger
from pybiomart import Server
from pyorthomap import FindOrthologs, findOrthologsHsMm, findOrthologsMmHs
from scipy.stats import fisher_exact, pearsonr, spearmanr, kruskal, chi2_contingency
from scipy.linalg import helmert
from statsmodels.stats.multitest import multipletests

from scikit_posthocs import posthoc_dunn
logger.info('Import OK')


def enrichment(organism, refOrganism, annotDataSet, enrichmentTestType, correction, genes, obo_path, reference=[]):

    """
    geneInputList * - string (query)
    ================================
    Each identifier to be delimited by comma i.e. ','. Maximum of 100000 Identifiers can be any of the following: Ensemble gene identifier, Ensemble protein identifier, Ensemble transcript identifier, Entrez gene id, gene symbol, NCBI GI, HGNC Id, International protein index id, NCBI UniGene id, UniProt accession andUniProt id

    organism * - string (query)	
    ===================================
    One taxon id required. Refer to 'supportedgenomes' service to retrieve supported ids

    refInputList - string  (query)	
    ===================================
    If not specified, the system will use all the genes for the specified organism. Each identifier to be delimited by comma i.e. ','. Maximum of 100000 Identifiers can be any of the following: Ensemble gene identifier, Ensemble protein identifier, Ensemble transcript identifier, Entrez gene id, gene symbol, NCBI GI, HGNC Id, International protein index id, NCBI UniGene id, UniProt accession andUniProt id

    refOrganism - string (query)
    =====================================	
    This parameter is only required, if parameter 'refInputList' has been specified. Only one taxon id can be specified. Retrieve the list of supported genomes and corresponding taxon ids from the 'supportedgenomes' service.

    annotDataSet * - string (query)
    =====================================	
    One of the supported PANTHER annotation data types. Use the 'supportedannotdatasets' service to retrieve list of supported annotation data types

    enrichmentTestType - string (query)
    =====================================	
    Fisher's Exact test will be used by default

    Available values : FISHER, BINOMIAL

    Default value : FISHER

    correction - string (query)
    =====================================	
    Available values : FDR, BONFERRONI, NONE
    Default value : FDR

    Note: this works, however there are some discrepancies compared the PantherDB online search in all except the top-most level. This is because "If a term is a parent of more than one term in the results table, it is shown only under its first descendant." To deal with this, then we would need to sort the terms by p-value and assign parent terms that occur in more than one family only to the top-most level 0 term
    """

    base_url = 'http://pantherdb.org/services/oai/pantherdb/enrich/overrep?'

    request_parameters = {"organism": organism,
                        "refOrganism": refOrganism,
                        "annotDataSet": annotDataSet,
                        "enrichmentTestType": enrichmentTestType,
                        "correction": correction}

    # Format gene list to comma-delimited string
    gene_list = []
    for l in genes:
        gene_list.append(l.rstrip())
    gene_list = ",".join(gene_list)
    request_parameters["geneInputList"] = gene_list

    if len(reference)>0:
        reference_list = []
        for l in reference:
            reference_list.append(l.rstrip())
        reference_list = ",".join(reference_list)
        request_parameters["refInputList"] = reference_list

    # Format parameters for url
    parameter_strings = []
    for parameter, query in request_parameters.items():
        if (" ") in parameter:
            # Format special characters for URL. Ex: "GO:0008150" -> "GO%3A0008150"
            query = quote(query)
        parameter_string = f"{parameter}={query}"
        parameter_strings.append(parameter_string)
    url = base_url + "&".join(parameter_strings)
    logger.info(f"Request URL - {url}")

    # Collect response using query url
    response = requests.get(url)

    try:
        results = response.json()['results']['result']
    except:
        print("ERROR: Parsing response failed. Full response:\n{}".format(response))
        return response

    all_results = pd.DataFrame(results)
    sig_results = all_results[all_results['pValue'] < 0.05]
    # handle "unclassified" term
    sig_results['term_id'] = [terms['id'] if 'id' in terms.keys() else 'None' for terms in sig_results['term']]
    sig_results['term_label'] = [terms['label'] if 'label' in terms.keys() else 'None' for terms in sig_results['term']]

    # read in geneontology database
    goslim = obonet.read_obo(obo_path)

    # Collect parent terms from ontolgy
    # assign top-level terms (i.e. with no child-terms)
    parent_terms = []
    for term_id in sig_results['term_id']:
        if 'UNCLASSIFIED' in term_id:
            parent_terms.append(['None'])
        else:
            terms = []
            try:
                terms.append([term for term in goslim.nodes[term_id]['is_a']])
            except:
                pass
            try:
                # to capture related terms that signal a term is part of something else
                terms.append([related.split(' ')[-1] for related in goslim.nodes[term_id]['relationship'] if 'part_of' in related])
            except:
                pass
            parent_terms.append([item for sublist in terms for item in sublist])

    sig_results['related_terms'] = parent_terms
    all_parents = [item for sublist in parent_terms for item in sublist]
    top_terms = [term for term in sig_results['term_id'] if term not in all_parents]

    # Generate level and family maps


    def level_mapper(term_id, level, level_map):
        level_map[term_id] = level
        parents = []
        try:
            parents.append([term for term in goslim.nodes[term_id]['is_a']])
        except:
            pass
        try:
            # to capture related terms that signal a term is part of something else
            parents.append([related.split(' ')[-1] for related in goslim.nodes[term_id]['relationship'] if 'part_of' in related])
        except:
            pass
        parents = [item for sublist in parents for item in sublist]
        for term in parents:
            if term in term_list:
                new_level = level + 1
                level_mapper(term, new_level, level_map)
        return level_map


    level_maps = []
    term_list = list(sig_results['term_id'])
    for term_id in top_terms:
        if 'UNCLASSIFIED' not in term_id:
            level_maps.append(level_mapper(term_id, level=0, level_map={}))
    levels = {}
    families = {}
    for number, family in enumerate(level_maps):
        levels.update(family)
        families[number] = list(family.keys())

    # invert families dictionary for mapping
    family_mapper = {}
    for key, value in families.items():
        for val in value:
            family_mapper[val] = key

    # Add family and level info to original df
    sig_results['family'] = sig_results['term_id'].map(family_mapper)
    sig_results['level'] = sig_results['term_id'].map(levels)

    return sig_results


def apply_enrichment(df, searches=None, obo_path='PANTHERGOslim.obo', organism='10090', refOrganism='10090', enrichmentTestType='FISHER', correction='BONFERONNI', min_proteins=5, reference=None):
    """
    Worker function to apply GO enrichment for df of proteins against the reference list (background) on a per-column basis. Column labels are returned in the 'column' key
    """

    if not searches:
        searches = {'Bio Process': 'ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP', 'Mol Function': 'ANNOT_TYPE_ID_PANTHER_GO_SLIM_MF', 'Cell Component': 'ANNOT_TYPE_ID_PANTHER_GO_SLIM_CC'}

    enrichment_dict = {}
    for col in df.columns.tolist():
        if len(df[col].dropna()) >= min_proteins:
            for abbr, search_type in searches.items():
                enrichment_dict[(col, abbr)] = enrichment(organism=organism, refOrganism=refOrganism, annotDataSet=search_type, enrichmentTestType=enrichmentTestType, correction=correction, genes=df[col].dropna(), obo_path=obo_path, reference=reference)

    summary = []
    for (col, search_type), result in enrichment_dict.items():
        result['column'] = col
        result['search_type'] = search_type
        summary.append(result)
    summary = pd.concat(summary)

    return summary


def fischers_test(go_genes, hit_genes, background, output_folder=None):

    genes = {'GO': go_genes,'hits': hit_genes,'background': background}
    
    # make sure all lists are unique
    for list_name, gene_list in genes.items():
        logger.info(f'Number of genes passed for {list_name}: {len(gene_list)}')
        genes[list_name] = set(gene_list)
        logger.info(
            f'Number of unique genes for {list_name}: {len(set(gene_list))}')

    assert genes['hits'] & genes['background'] == genes['hits'] # check that all hits are in background

    # Initialize variables for the contingency table
    not_hits = genes['background'] - genes['hits']
    genes['not_hits'] = not_hits
    logger.info(f'Number of unique non-hit genes: {len(not_hits)}')

    # Determine elements of the contingency table:
    hit_in_go = genes['hits'] & genes['GO']
    hit_not_go = genes['hits'] - genes['GO']
    not_hit_in_go = genes['not_hits'] & genes['GO']
    not_hit_not_go = genes['not_hits'] - genes['GO']

    # Summarise genes in each category
    gene_summary = pd.DataFrame(
        [hit_in_go, hit_not_go, not_hit_in_go, not_hit_not_go],
        index=['hit_in_go', 'hit_not_go', 'not_hit_in_go', 'not_hit_not_go']).T

    # build contingency table
    contingency_df = pd.DataFrame(
        {'Hits': [len(hit_in_go), len(hit_not_go)],
         'Not Hits': [len(not_hit_in_go), len(not_hit_not_go)]}, index=['In GO', 'Not in GO'])
    contingency_df['Total'] = contingency_df.sum(axis=1) # column total
    contingency_df.loc['Total'] = contingency_df.sum() # row total

    # perform test!
    odd_ratio, pval = fisher_exact([[len(hit_in_go), len(not_hit_in_go)], [len(hit_not_go), len(not_hit_not_go)]])

    if output_folder:
        test_df = pd.DataFrame({odd_ratio, pval}, index=['odds_ratio', 'pval'])
        FileHandling.df_to_excel(output_path=f'{output_folder}fischer_result.xlsx', data_frames=[gene_summary, contingency_df, test_df], sheetnames=['gene_summary', 'contingency_df', 'test_df'])

    return odd_ratio, pval, contingency_df


def apply_fisher_test(hit_summary, background_genes, ontology_genes, test_name, output_folder=False):
    """Worker function for completing tests on predefined datasets
    hit_summary: df, each column is tested for enrichment against the background for ontology_genes
    background_genes: list, background geneset. Must contain all hit genes + non-hit genes
    ontology_genes: list, benchmark genes e.g. specific ontology term
    test_name: str appended to filename of output
    """
    contingencies = {}
    test_results = {}
    for column_name in hit_summary.columns.tolist():
        hit_genes = hit_summary[column_name].dropna()
        oddsratio, pval, test_result = fischers_test(ontology_genes, hit_genes, background_genes, output_folder=None)
        test_results[column_name] = [oddsratio, pval]
        contingencies[column_name] = test_result
    contingencies['test_summary'] = pd.DataFrame(test_results, index=['oddsratio', 'pval']).T
    if output_folder:
        FileHandling.df_to_excel(output_path=f'{output_folder}{test_name}_fisher_summary.xlsx', data_frames=list(contingencies.values()), sheetnames=list(contingencies.keys()))

    return contingencies


def stability_summary_calculator(compilation, id_col):
    """Generate a summary stability measure for published or measured stability datasets,
    where dataset contains mixture of protein/peotide level measurements with or without replicates - preference order is (1) mean/median of replicate protein values, (2) protein value, (2) mean/median stability from peptides, (3) peptide stability or (5) NaN"""
    stability_summary = []
    for (resource, id_col), df in compilation.groupby(['resource', id_col]):
        protein = df['protein_stability'].dropna()
        peptide = df['peptide_stability'].dropna()
        if len(protein) > 0:
            df['mean_stability'] = protein.mean()
            df['median_stability'] = protein.median()
        elif len(peptide) > 0:
            df['mean_stability'] = peptide.mean()
            df['median_stability'] = peptide.median()
        else:
            df['mean_stability'] = np.nan
            df['median_stability'] = np.nan
        stability_summary.append(df)
    return pd.concat(stability_summary)


def r_squared_calculator(x_vals, y_vals, function, popt):
    residuals = y_vals - function(x_vals, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y_vals-np.mean(y_vals))**2)
    return 1 - (ss_res / ss_tot)


def correlation_df(x, y, corr_type=None):
    data = pd.DataFrame()
    data['y'] = y
    data['x'] = x

    data.dropna(inplace=True)
    if len(data) < 5:
        return pd.DataFrame([np.nan, np.nan, len(data)], index=['pearsons_r', 'pearsons_pval', 'count'])
    if corr_type == 'pearsons':
        (r, p) = pearsonr(data['x'], data['y'])
        count = len(data)
        return pd.DataFrame([r, p, count], index=['pearsons_r', 'pearsons_pval', 'count'])
    elif corr_type == 'spearmans':
        (r, p) = spearmanr(data['x'], data['y'])
        count = len(data)
        return pd.DataFrame([r, p, count], index=['spearmans_r', 'spearmans_pval', 'count'])
    else:
        (pearsons_r, pearsons_p) = pearsonr(data['x'], data['y'])
        (spearmans_r, spearmans_p) = spearmanr(data['x'], data['y'])
        count = len(data)
        return pd.DataFrame([pearsons_r, pearsons_p, spearmans_r, spearmans_p, count], index=['pearsons_r', 'pearsons_pval', 'spearmans_r', 'spearmans_pval', 'count'])


def apply_oneway_anova(df, xcol, group_col):
    ols_model = f'{xcol} ~ C({group_col})'
    lm = sfa.ols(ols_model, data=df).fit()
    anova = sa.stats.anova_lm(lm)

    # complete post-hoc test with bonferroni correction
    pairs = lm.t_test_pairwise(f"C({group_col})", method='bonferroni')

    return anova, pairs.result_frame

def anova():
    pass

def apply_kruskal(df, val_col, group_col, p_adjust='bonferroni', sig_cutoff=0.05):
    test_val, pval = kruskal(*[data[val_col].tolist() for group, data in df.groupby(group_col)], nan_policy='omit')

    if pval < sig_cutoff:
        post_hocs = posthoc_dunn(df, val_col=val_col, group_col=group_col, p_adjust=p_adjust, sort=True)
    else:
        post_hocs = pd.DataFrame()

    return test_val, pval, post_hocs

def ilr_transformation(data):

    """Isometric logratio transformation (not vectorized).
    Adapted from https://github.com/ofgulban/compoda/blob/master/compoda/core.py
    Parameters
    ----------
    data : 2d numpy array, shape [n_samples, n_coordinates]
        Barycentric coordinates (closed) in simplex space.
    Returns
    -------
    out : 2d numpy array, shape [n_samples, n_coordinates-1]
        Coordinates in real space.
    Reference
    ---------
    [1] Pawlowsky-Glahn, V., Egozcue, J. J., & Tolosana-Delgado, R.
        (2015). Modelling and Analysis of Compositional Data, pg. 37.
        Chichester, UK: John Wiley & Sons, Ltd.
        DOI: 10.1002/9781119003144
    """
    dims = data.shape
    out = np.zeros((dims[0], dims[1]-1))
    helmertian = helmert(dims[1]).T
    for i in range(data.shape[0]):
        out[i, :] = np.dot(np.log(data[i, :]), helmertian)
    return out




def apply_chisquare(df, p_adjust='bonferroni', sig_cutoff=0.05):
    test_val, pval, dof, contingency = chi2_contingency(df.fillna(0).values)

    if pval < sig_cutoff:
        posthoc_test = defaultdict(defaultdict)
        for col_1 in df.columns.tolist():
            for col_2 in df.columns.tolist():
                if col_1 == col_2:
                    continue
                ph_test_val, ph_pval, ph_dof, ph_contingency = chi2_contingency(df.fillna(0)[[col_1, col_2]].values)
                posthoc_test[col_1][col_2] = ph_pval
        posthoc_test = pd.melt(pd.DataFrame(posthoc_test).reset_index().rename(columns={'index': 'col_1'}), id_vars='col_1', value_vars=None, var_name='col_2', value_name='pval')
        posthoc_test['corrected_pval'] = multipletests(posthoc_test['pval'], alpha=0.05, method='bonferroni', is_sorted=False, returnsorted=False)[1]
    else:
        posthoc_test = pd.DataFrame()
    # add initial pvalue for summary purposes
    posthoc_test = pd.concat([posthoc_test, pd.DataFrame(['all', 'all', pval, np.nan], index=['col_1', 'col_2', 'pval', 'corrected_pval']).T]).reset_index(drop=True)

    return test_val, pval, dof, contingency, posthoc_test
