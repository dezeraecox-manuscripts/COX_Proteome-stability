import os
import pandas as pd
from utilities.statistical_tests import apply_enrichment, apply_fisher_test

from loguru import logger

logger.info('Import OK')

input_path = 'results/inhibitor_urea_denaturation/detect_outliers/outlier_summary.xlsx'
output_folder = 'results/inhibitor_urea_denaturation/enrichment_analysis/'
ontology_path = 'results/inhibitor_urea_denaturation/protein_interactions/HSPAs_interactors.csv'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


# read in raw results
compiled = pd.read_excel(input_path, sheet_name=None)
compiled.update({key: df.drop([col for col in df.columns.tolist(
) if 'Unnamed: ' in str(col)], axis=1) for key, df in compiled.items()})

# Collect background as all smoothed proteins
background = compiled['smooth_ratios'].copy()
background = background['Proteins'].unique().tolist()

# Collect outliers
outliers = compiled['outliers'].copy()

# optional: annotate direction??
# note - only one 'up' protein was unique to the 'up' group (P32921), the rest were also found in the down list
info_cols = ['Proteins', 'Sequence', 'outlier', 'noncys_ratio_std', 'noncys_ratio_mean', 'count']
proteins_up = outliers.set_index(info_cols)[outliers.set_index(
    info_cols) > outliers['noncys_ratio_std']].dropna(how='all').reset_index()['Proteins'].unique().tolist()

proteins_down = outliers.set_index(info_cols)[outliers.set_index(
    info_cols) < outliers['noncys_ratio_std']].dropna(how='all').reset_index()['Proteins'].unique().tolist()

# Collect list of combinations to test
combinations = {
    'outliers_all': outliers['Proteins'].unique().tolist(),
    'outliers_up': proteins_up,
    'outliers_down': proteins_down,
    'outliers_mixed': list(set(proteins_down).intersection(set(proteins_up))),
}
combinations = pd.DataFrame(
    dict([(k, pd.Series(v)) for k, v in combinations.items()]))

# perform standard Panther enrichment test
enrichment = apply_enrichment(combinations, searches=None, obo_path='resources/bioinformatics_databases/PANTHERGOslim.obo',organism='10090', refOrganism='10090', enrichmentTestType='FISHER', correction='FDR', min_proteins=2, reference=background)

# Save all to excel
enrichment.to_csv(f'{output_folder}go_enrichment.csv')

# Look for overrepresentation of HSP70 interactors
ontology_map = pd.read_csv(f'{ontology_path}')
ontology_genes = ontology_map['Proteins'].unique().tolist()
len(ontology_genes)

apply_fisher_test(
    hit_summary=combinations, 
    background_genes=background, 
    ontology_genes=ontology_genes, 
    test_name='outliers',
    output_folder=output_folder)
