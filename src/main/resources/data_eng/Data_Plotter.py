"""
This python script is used to plot the data for the Reactome immport web app manuscript. The script uses scanpy
for umap and other related packages for violin plot. See the required packages in the immport list.
"""
import math
import random

import numpy as np
import scanpy as sc
import pandas as pd
import scanpy.external as sce
from os.path import exists
import seaborn as sns
import matplotlib.pyplot as plt
import textwrap

data_file_name = '../data/gene_expression_matrix.csv'
meta_file_name = '../data/biosample_metadata.csv'
h5ad_file_name = '../data/gene_expression.h5ad'

random_state = 1256082


def load_data_to_sc(data_file_name: str = data_file_name,
                    meta_file_name: str = meta_file_name,
                    h5ad_file_name: str = h5ad_file_name) -> sc.AnnData:
    if exists(h5ad_file_name):
        return sc.read_h5ad(h5ad_file_name)
    # Use DataFrame to read to avoid NA issue
    df = pd.read_csv(data_file_name)
    df = df.set_index('GeneGene')
    df = df.transpose()
    # To make UMAP work, replace NA as 0
    df.fillna(0, inplace=True)
    adata = sc.AnnData(X=df)
    # Attach the meta file
    meta = pd.read_csv(meta_file_name)
    meta.set_index(keys='gsm', inplace=True)
    adata.obs = adata.obs.join(meta)
    # Keep it for further loading
    adata.write_h5ad(h5ad_file_name)
    return adata


def run_umap(adata: sc.AnnData,
             batch_key: str = 'gse'):
    # Do a little bit pre-process
    # These parameters are quite arbitry
    adata_umap = adata.copy()
    sc.pp.filter_cells(adata_umap, min_genes=200)
    sc.pp.filter_genes(adata_umap, min_cells=100)
    sc.pp.normalize_total(adata_umap, 1E+4)
    sc.pp.log1p(adata_umap)
    sc.pp.highly_variable_genes(adata_umap)
    sc.pp.pca(adata_umap, random_state=random_state)
    # sc.pp.neighbors(adata_umap, random_state=random_state)
    if batch_key is not None:
        sce.pp.bbknn(adata_umap, batch_key=batch_key)
    else:
        # bbknn do neighbor analysis. If no bbknn, we need run neighbors here.
        sc.pp.neighbors(adata_umap, random_state=random_state)
    # sc.pp.neighbors(adata_umap, use_rep='X_pca_harmony', random_state=random_state)
    sc.tl.umap(adata_umap, random_state=random_state)
    return adata_umap


def plot_umap(adata: sc.AnnData,
              obs_color: list = ['gse', 'gpl', 'vaccine']):
    # Do a little bit update
    if 'vaccine' in obs_color:
        _add_vaccine_short(adata)
    sc.pl.umap(adata, color=obs_color, wspace=0.5)


def _add_vaccine_short(adata: sc.AnnData):
    if 'vaccine_short' not in adata.obs.keys():
        adata.obs['vaccine_full'] = adata.obs['vaccine']
        adata.obs['vaccine'] = adata.obs['vaccine'].map(lambda x: x[0:12])


def plot_expression(adata: sc.AnnData,
                    sample_genes: int = 1000,
                    group_by: str = 'vaccine',
                    vaccines_for_plot: list = None) -> pd.DataFrame:
    adata_copy = adata.copy()
    # Need to pick the top variable genes for this plot
    sc.pp.normalize_total(adata_copy, 1E+4)
    sc.pp.log1p(adata_copy)
    sc.pp.highly_variable_genes(adata_copy)
    highly_variable_genes = adata.var_names[adata_copy.var_vector('highly_variable')]
    # Randomly sample 1,000 genes for plot
    random.seed(random_state)
    highly_variable_genes = random.sample(list(highly_variable_genes), sample_genes)
    adata_copy = adata.copy()
    adata_copy = adata_copy[:, highly_variable_genes]
    if group_by == 'vaccine':
        _add_vaccine_short(adata_copy)
        # group_by = 'vaccine_short'
    # Create a dataframe for plot
    violin_df = adata_copy.to_df()
    violin_df[group_by] = adata_copy.obs[group_by]
    violin_df = violin_df.melt(id_vars=[group_by], var_name='gene', value_name='expression_value')
    # Reset 0 as nan
    violin_df['expression_value'] = violin_df['expression_value'].replace(0, np.nan)
    violin_df.dropna(inplace=True)
    # Do some slicing if needed
    if vaccines_for_plot is not None:
        which_rows = violin_df['vaccine'].isin(vaccines_for_plot)
        violin_df = violin_df[which_rows]
        violin_df.vaccine.cat = vaccines_for_plot
    plt.subplots_adjust(bottom=0.30)
    # Apparently boxplot gives us a better view in this case
    g = sns.boxplot(x='vaccine', y='expression_value', data=violin_df)
    # g = sns.violinplot(x='vaccine_short',
    #                    y='expression_value', data=violin_df)
    g.set_xticklabels(g.get_xticklabels(), rotation=90)
    g.set_xlabel('vaccine')
    # For selected two vaccines we, may need to adjust the font sizes
    return violin_df, g


def plot_pathways_in_use_case_1(need_end_results: bool = True,
                                need_temp_result: bool = True):
    """
    This function is used to plot a pre-selected list of pathways in use case 1.
    :return:
    """
    pathways = [
        "Immune System",
        "Adaptive Immune System",
        "Class I MHC mediated antigen processing & presentation",
        "Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell",
        "Innate Immune System",
        "Neutrophil degranulation",
        "Cytokine Signaling in Immune system",
        "Interferon Signaling",
        "Interferon alpha/beta signaling",
        "Interferon gamma signaling",
        "Antiviral mechanism by IFN-stimulated genes",
        "Signaling by Interleukins",
        "Interleukin-4 and Interleukin-13 signaling",
        "Interleukin-10 signaling",
        "Cell Cycle"
    ]
    dir_name = '/Volumes/ssd/results/immport-ws/Manuscript_UseCase_1_Results'
    # File names for the end results
    weeks = [1, 2, 4]
    week_names = ['OneWeek', 'TwoWeek', 'FourWeek']
    print('Check numbers of genes selected for pathway enrichment analysis:')
    # Get the total genes used in the pathway analysis
    for i in range(len(weeks)):
        # End results
        file_name = dir_name + '/' + week_names[i] + '_GeneExpressionAnalysis.csv'
        print(week_names[i], ' ', file_name)
        _print_selected_genes(file_name)
        # Get the results between
        if i > 0:
            file_name = dir_name + '/' + week_names[i - 1] + 'To' + week_names[i] + '_GeneExpressionAnalysis.csv'
            print(week_names[i - 1] + 'To' + week_names[i], ' ', file_name)
            _print_selected_genes(file_name)
    # Generate data frames for plot
    types = ['AllGenes', 'PosGenes', 'NegGenes']
    if need_end_results:
        end_results_df = _plot_end_results(dir_name, pathways, types, week_names, weeks)
    if need_temp_result:
        temp_results_df = _plot_temporal_results(dir_name, pathways, types, week_names, weeks)
    return temp_results_df


def _plot_temporal_results(dir_name, pathways, types, week_names, weeks):
    temp_results_df = None;
    print('Plot temporal results:')
    for i in range(len(weeks)):
        for type in types:
            file_name = None
            if i == 0:
                file_name = '{}/{}_{}_PathwayEnrichmentAnalysis.csv'.format(dir_name,
                                                                            week_names[i],
                                                                            type)
            else:
                file_name = '{}/{}To{}_{}_PathwayEnrichmentAnalysis.csv'.format(dir_name,
                                                                                week_names[i - 1],
                                                                                week_names[i],
                                                                                type)
            temp_results_df = _load_pathway_df(temp_results_df, file_name, pathways, type, weeks[i])
            temp_results_df['Pathway'] = temp_results_df['Pathway Name']
    # Print out the end results
    fig, ax = plt.subplots(len(types), 2)
    fig.tight_layout(pad=.05)
    plt.subplots_adjust(bottom=0.05, top=0.95)
    plt.rcParams.update({'font.size': 16})  # This apply to titles only
    # Default pathways into two groups
    pathways_1 = pathways.copy()
    pathways_2 = pathways_1[6:14].copy() # Cytokine signaling
    for pathway in pathways_2:
        pathways_1.remove(pathway)
    pathways_list = [pathways_1, pathways_2]
    for j in range(len(pathways_list)):
        pathways = pathways_list[j]
        which_rows = temp_results_df['Pathway'].isin(pathways)
        temp_results_df_pathways = temp_results_df[which_rows]
        for i in range(len(types)):
            temp_results_df_filter = temp_results_df_pathways[temp_results_df_pathways['Type'] == types[i]]
            for pathway in pathways:
                temp_results_df_filter = temp_results_df_filter.append({'Pathway': pathway,
                                                                        '-Log10(FDR)': 0,
                                                                        'Week': 0,
                                                                        'Type': types[i]},
                                                                       ignore_index=True)
            # Get the following markers from https://github.com/mwaskom/seaborn/issues/1513
            filled_markers = ('o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X')
            # Make sure we have the same order of pathways across all
            temp_results_df_filter = temp_results_df_filter.sort_values(by=['Pathway', 'Week'])
            g = sns.lineplot(data=temp_results_df_filter,
                             x='Week',
                             y='-Log10(FDR)',
                             hue='Pathway',
                             markers=filled_markers,
                             style='Pathway',
                             dashes=False,
                             markersize=7,
                             alpha=0.8,
                             ax=ax[i][j])
            g.set_title(types[i])
            g.set(ylim=(0, 15))
            g.set_ylabel(g.get_ylabel(), size=15)
            g.set_xlabel(g.get_xlabel(), size=15)
            g.set_xticks([0, 1, 2, 4])
            if i < len(types) - 1:
                g.legend().remove()
            else:
                handles, labels = g.get_legend_handles_labels()
                order = [labels.index(pathway) for pathway in pathways]
                # The first is the title: pathway
                order.insert(0, 0)
                labels = wrap_pathway_names(labels, width=45)
                g.legend([handles[i] for i in order], [labels[i] for i in order],
                         loc='upper left',
                         fontsize=12)
            if i < len(types) - 1:
                g.set_xticklabels([])
                g.set_xlabel(None)
            else:
                g.set_xticklabels(g.get_xticks(), size=15)
    return temp_results_df_filter


def _plot_end_results(dir_name, pathways, types, week_names, weeks):
    end_results_df = None
    print('Plot end results:')
    for i in range(len(weeks)):
        for type in types:
            file_name = '{}/{}_{}_PathwayEnrichmentAnalysis.csv'.format(dir_name,
                                                                        week_names[i],
                                                                        type)
            end_results_df = _load_pathway_df(end_results_df, file_name, pathways, type, weeks[i])
    # Print out the end results
    fig, ax = plt.subplots(len(types), 1)
    fig.tight_layout(pad=.1)
    plt.subplots_adjust(bottom=0.2, top=0.95)
    plt.rcParams.update({'font.size': 16})  # This apply to titles only
    for i in range(len(types)):
        end_results_df_type = end_results_df[end_results_df['Type'] == types[i]]
        # Need to sort based on pathway index
        end_results_df_type['X Order'] = end_results_df_type['Pathway Name'].map(lambda n: pathways.index(n))
        end_results_df_type.sort_values(by='Pathway Name', inplace=True)
        print(end_results_df_type.shape)
        g = sns.barplot(data=end_results_df_type,
                        x='X Order',
                        y='-Log10(FDR)',
                        hue='Week',
                        ax=ax[i],
                        ci=None)
        g.set_title(types[i])
        g.set(ylim=(0, 14))
        g.set_ylabel(g.get_ylabel(), size=15)
        if i > 0:
            g.legend().remove()
        # Force to use the sorted pathways
        if i < len(types) - 1:
            g.set_xticklabels([], rotation=90)
            g.set_xlabel(None)
        else:
            g.set_xticklabels(wrap_pathway_names(pathways, width=18), rotation=90, size=15)
            g.set_xlabel('Pathway')
    return end_results_df


def _load_pathway_df(end_results_df, file_name, pathways, type, week):
    df = pd.read_csv(file_name)
    filter = df['Pathway Name'].isin(pathways)
    df_filter = df[filter]
    # Check to make sure all pathways there
    for pathway in pathways:
        if pathway not in df_filter['Pathway Name'].to_list():
            # print("Need to add " + pathway)
            df_filter = df_filter.append({'Pathway Name': pathway,
                                          'Entities FDR': 1.0},
                                         ignore_index=True)
    # Add a new column
    df_filter['Type'] = type
    df_filter['Week'] = week
    if end_results_df is None:
        end_results_df = df_filter
    else:
        end_results_df = pd.concat([end_results_df, df_filter])
    end_results_df['-Log10(FDR)'] = end_results_df['Entities FDR'].map(lambda x: -math.log10(x))
    return end_results_df


def _print_selected_genes(file_name):
    df = pd.read_csv(file_name)
    # All genes: log2FC > 0.2 or <-0.2 and p-value < 0.05
    filters = ((df['Log2FC'] > 0.2) | (df['Log2FC'] < -0.2)) & (df['pValue'] < 0.05)
    df_filtered = df[filters]
    print('Pos and Neg genes: {}'.format(df_filtered.shape[0]))
    filters = (df['Log2FC'] > 0.2) & (df['pValue'] < 0.05)
    df_filtered = df[filters]
    print('Pos genes: {}'.format(df_filtered.shape[0]))
    filters = (df['Log2FC'] < -0.2) & (df['pValue'] < 0.05)
    df_filtered = df[filters]
    print('Neg genes: {}'.format(df_filtered.shape[0]))
    print()


def wrap_pathway_names(pathways, width=10):
    """
    This function is modified from https://medium.com/dunder-data/automatically-wrap-graph-labels-in-matplotlib-and-seaborn-a48740bc9ce
    :param pathways:
    :return:
    """
    labels = []
    for name in pathways:
        name = textwrap.fill(name, width=width, break_long_words=False)
        labels.append(name)
    return labels


# plot_pathways_in_use_case_1()