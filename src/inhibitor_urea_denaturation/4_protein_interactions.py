""" NOTE: current versions of the py4cytoscape require cytoscape to be open and running."""

import time
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import py4cytoscape as cyto

from loguru import logger

logger.info('Import OK')

input_path = 'results/inhibitor_urea_denaturation/detect_outliers/outlier_summary.xlsx'
output_folder = 'results/inhibitor_urea_denaturation/protein_interactions/'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Check cytoscape has been started
# start_cyto = subprocess.Popen(
#     r"C:/Program Files/Cytoscape_v3.9.0/cytoscape.exe", shell=True)
cyto.cytoscape_ping()

font = {'family': 'normal',
        'weight': 'normal',
        'size': 14}
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'


def cytoscape_STRING_import(proteins, species='Mus Musculus', confidence=0.4, add_interactors=0):
    proteins = ','.join(proteins)
    # Collect network from STRING
    string_cmd = f'string protein query species="{species}" query="{proteins}" cutoff="{confidence}" limit="{add_interactors}"'
    cyto.commands_run(string_cmd)

    # Get table of mapped nodes
    nodes = cyto.get_table_columns()

    return nodes


def cytoscape_map_table(df, df_key='Proteins', node_key='query term'):
    # Map onto node table
    cyto.load_table_data(df, data_key_column=df_key,
                         table_key_column=node_key)
    # nodes var is not automatically updated with new info
    nodes = cyto.get_table_columns()
    return nodes


def cytoscape_create_style(style_name, color_col, map_vals=[-1, 0, 1], map_colors=['#05008a', '#adadad', '#a80000']):

    # Create new style with some mappings
    defaults = {'NODE_SHAPE': "circle", 'NODE_SIZE': 10,
                'EDGE_TRANSPARENCY': 70, 'NODE_LABEL_POSITION': "W,E,c,0.00,0.00"}
    node_labels = cyto.map_visual_property(
        visual_prop='node label', table_column='display name', mapping_type='p')  # 'p' means 'passthrough' mapping
    node_color = cyto.map_visual_property(
        visual_prop='node fill color', table_column=color_col, mapping_type='c',
        table_column_values=map_vals, visual_prop_values=map_colors)  # 'p' means 'passthrough' mapping, col_vals are the extremes to map, vis_prop are the colors
    cyto.create_visual_style(style_name, defaults, [node_labels, node_color])


def find_interactors(edges, proteins_of_interest):
    """Finds a list of unique, direct interactions with the proteins of interest

    Parameters
    ----------
    edges : DataFrame
        edges dataframe collected from cytoscape with 'shared name' column mapping ppis
    proteins_of_interest : list(str)
        List of display names for proteins of interest to which returned proteins share a direct interaction

    Returns
    -------
    list
        list of proteins with direct interaction to proteins of interest
    """
    test_edges = edges['shared name'].copy().reset_index()
    test_edges[['Protein A', 'Protein B']] = test_edges['shared name'].str.split(
        r' \(pp\) ', expand=True)
    test_edges['A_interactor'] = [
        1 if protein in proteins_of_interest else 0 for protein in test_edges['Protein A']]
    test_edges['B_interactor'] = [
        1 if protein in proteins_of_interest else 0 for protein in test_edges['Protein B']]
    test_edges['interactor'] = [1 if sum_val == 1 else 0 for sum_val in test_edges[['A_interactor', 'B_interactor']].sum(
        axis=1)]  # note - removes interactions BETWEEN proteins of interest (where sum would be 2)
    interactors = [item for sublist in test_edges[test_edges['interactor'] == 1][[
        'Protein A', 'Protein B']].values for item in sublist]
    return list({protein for protein in interactors if protein not in proteins_of_interest})


def edge_distance(edges, seed_proteins):
    interactors = find_interactors(edges, seed_proteins)
    interactions = [seed_proteins, interactors]
    while len(interactors) > 0:
        interactors = find_interactors(edges, interactors)
        interactors = [protein for protein in interactors if protein not in [
            item for sublist in interactions for item in sublist]]
        interactions.append(interactors)
    interactions = dict(zip(range(len(interactions)), interactions))
    interactions = pd.DataFrame(
        dict([(k, pd.Series(v)) for k, v in interactions.items()]))
    return pd.melt(
        interactions,
        var_name='edge_distance',
        value_vars=interactions.columns.tolist(),
        value_name='display_name'
    ).dropna()[['display_name', 'edge_distance']]


def interactor_count(edges):
    test_edges = edges.copy()
    test_edges[['Protein A', 'Protein B']] = test_edges['shared name'].str.split(
        r' \(pp\) ', expand=True)
    unique_interactors = list(
        set(test_edges['Protein A'].tolist() + test_edges['Protein B'].tolist()))
    degree = {
        protein: len(find_interactors(edges, [protein]))
        for protein in unique_interactors
    }
    return pd.DataFrame([degree.keys(), degree.values()], index=['display_name', 'interactor_count']).T

# -------------------------------Prepare quantitative data-------------------------------
# Read in summary data
compiled = pd.read_excel(input_path, sheet_name=None)
compiled.update({key: df.drop([col for col in df.columns.tolist(
) if 'Unnamed: ' in str(col)], axis=1) for key, df in compiled.items()})

# Collect background as all smoothed proteins
background = compiled['smooth_ratios']['Proteins'].unique().tolist()

# Collect outliers
outliers = compiled['outliers'].copy()
outliers['max_change'] = outliers[[col for col in outliers.columns if type(col) != str]].abs().max(axis=1)

# Add HSPAs to list of ID's for mapping
hspas = ['P63017', 'Q61696', 'P17879'] # 3 major isoforms

# ---------------------------Prepare STRING PPI map---------------------------
# Get STRING map for proteins of interest
protein_ids = outliers['Proteins'].unique().tolist()

cytoscape_STRING_import(protein_ids+hspas, species='Mus Musculus',
                        confidence=0.4, add_interactors=0)

# Map maximum deviation per protein 
cytoscape_map_table(
    df=outliers[['Proteins', 'max_change']].groupby('Proteins').max().reset_index(),
    )


# Transfer to cytoscape, manually organise layout according to distance from mapped HSPAs

# Get data from table and remap onto node table (this should be automatic but appears broken)
nodes = cyto.get_table_columns('node')
nodes.to_csv(f'{output_folder}cytoscape_outlier_node_summary.csv')

# Determine degree of interaction with HSPAs
edges = cyto.get_table_columns('edge')
edges.to_csv(f'{output_folder}cytoscape_edge_summary.csv')

display_map = dict(nodes[['query term', 'display name']].values)
hspas = [display_map[entry] for entry in hspas]
interaction_distance = edge_distance(
    edges, hspas)
interaction_distance.to_csv(f'{output_folder}HSPA_edge_distance.csv')
cytoscape_map_table(
    df=interaction_distance,
    df_key='display_name',
    node_key='display name'
)

# determine connectivity of each node within the network
interactor_degree = interactor_count(edges)
interactor_degree['interactor_count'] = interactor_degree['interactor_count'].astype(int)
interactor_degree.to_csv(f'{output_folder}interactor_counts.csv')

cytoscape_map_table(
    df=interactor_degree,
    df_key='display_name',
    node_key='display name'
)
 
# Create a new style 
defaults = {'NODE_SHAPE': "circle", 'EDGE_TRANSPARENCY': 70,
            'NODE_LABEL_POSITION': "W,E,c,0.00,0.00", "NODE_BORDER_WIDTH": 2, "NODE_FILL_COLOR": '#000000', "NODE_SIZE": 20}
node_labels = cyto.map_visual_property(
    visual_prop='node label', table_column='display name', mapping_type='p')  # 'p' means 'passthrough' mapping

node_color = cyto.map_visual_property(
    visual_prop='node fill color',
    table_column='max_change',
    mapping_type='c',
    table_column_values=[0, 1.5],
    visual_prop_values=['#ffffff', '#000000']
)

cyto.create_visual_style('max_change', defaults, [node_labels, node_color])
cyto.style_mappings.set_node_size_mapping('interactor_count', table_column_values=[0, 10], sizes=[5, 50], style_name='max_change', mapping_type='c')
cyto.set_visual_style('max_change')

cyto.notebook_export_show_image()
cyto.networks.rename_network('outlier_degree')

cyto.export_image(f'outlier_degree_change', type='SVG',
                  network=f'outlier_degree')
# copy image file to Notebook directory
cyto.sandbox_get_from(
    f'outlier_degree_change.svg', f'{output_folder}outlier_degree_change.svg')


cyto.clone_network('outlier_degree')
time.sleep(10)
cyto.networks.rename_network(f'outlier_plain')

# Create a new style 
defaults = {'NODE_SHAPE': "circle", 'EDGE_TRANSPARENCY': 70,
            'NODE_LABEL_POSITION': "W,E,c,0.00,0.00", "NODE_BORDER_WIDTH": 2, "NODE_FILL_COLOR": '#000000', "NODE_SIZE": 20}
node_labels = cyto.map_visual_property(
    visual_prop='node label', table_column='display name', mapping_type='p')  # 'p' means 'passthrough' mapping
cyto.create_visual_style('plain_labelled', defaults, [node_labels, node_color])
cyto.set_visual_style('plain_labelled')

cyto.export_image(f'outlier_plain', type='SVG', network=f'outlier_plain')
# copy image file to Notebook directory
cyto.sandbox_get_from(
    f'outlier_plain.svg', f'{output_folder}outlier_plain.svg')

# Save session file - copy session file from Notebook directory to workstation
cyto.save_session('outliers')
cyto.sandbox_get_from('outliers.cys', f'{output_folder}outliers.cys')


# ----------------Generate first-shell interactors for HSPAs----------------
# for each HSPA, import string network with up to 1000 extra first-shell interactors (all are below that anyway)

hspas = ['P63017', 'Q61696', 'P17879']  # 3 major isoforms
for protein in hspas:
    protein
    cytoscape_STRING_import([protein], species='Mus Musculus',
                            confidence=0.4, add_interactors=1000)
    cyto.networks.rename_network(f'{protein}')

# Collect node and edge dataframes, extract interactors
edge_dfs = []
node_dfs = []
for protein in hspas:
    cyto.set_current_network(protein)
    time.sleep(10)
    edge_dfs.append(cyto.get_table_columns('edge'))
    node_dfs.append(cyto.get_table_columns('node'))
edges = pd.concat(edge_dfs)
nodes = pd.concat(node_dfs)

# compile list of unique interactors with any of the HSPAs
display_map = dict(nodes[['query term', 'display name']].values)
hspa_names = [display_map[entry] for entry in hspas]
interactors = []
for protein in hspa_names:
    interactor_distance = edge_distance(edges, seed_proteins=[protein])
    interactor_distance = interactor_distance[interactor_distance['edge_distance'] == 1]
    interactor_distance['HSPA'] = protein
    interactors.append(interactor_distance)
interactors = pd.concat(interactors)
interactors['Proteins'] = interactors['display_name'].map(
    dict(nodes[['display name', 'stringdb::canonical name']].values))

# Save to excel for enrichment testing
interactors.to_csv(f'{output_folder}HSPAs_interactors.csv')

# Save session file - copy session file from Notebook directory to workstation
cyto.save_session('HSPA_interactors')
cyto.sandbox_get_from('HSPA_interactors.cys', f'{output_folder}HSPA_interactors.cys')
