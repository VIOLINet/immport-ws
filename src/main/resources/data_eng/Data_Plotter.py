"""
This python script is used to plot the data for the Reactome immport web app manuscript. The script uses scanpy
for umap and other related packages for violin plot. See the required packages in the immport list.
"""
import random

import numpy as np
import scanpy as sc
import pandas as pd
import scanpy.external as sce
from os.path import exists
import seaborn as sns
import matplotlib.pyplot as plt

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
    sce.pp.bbknn(adata_umap, batch_key=batch_key)
    # sc.pp.neighbors(adata_umap, use_rep='X_pca_harmony', random_state=random_state)
    sc.tl.umap(adata_umap, random_state=random_state)
    return adata_umap


def plot_umap(adata: sc.AnnData,
              obs_color: list=['gse', 'gpl', 'vaccine']):
    # Do a little bit update
    if 'vaccine' in obs_color:
        _add_vaccine_short(adata)
        v_index = obs_color.index('vaccine')
        obs_color[v_index] = 'vaccine_short'
    sc.pl.umap(adata, color=obs_color, wspace=0.5)


def _add_vaccine_short(adata:sc.AnnData):
    if not 'vaccine_short' in adata.obs.keys():
        adata.obs['vaccine_short'] = adata.obs['vaccine'].map(lambda x: x[0:12])


def plot_expression(adata: sc.AnnData,
                    sample_genes: int=1000,
                    group_by: str='vaccine') -> pd.DataFrame:
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
        group_by = 'vaccine_short'
    # Create a dataframe for plot
    violin_df = adata_copy.to_df()
    violin_df[group_by] = adata_copy.obs[group_by]
    violin_df = violin_df.melt(id_vars=[group_by], var_name='gene', value_name='expression_value')
    # Reset 0 as nan
    violin_df['expression_value'] = violin_df['expression_value'].replace(0, np.nan)
    violin_df.dropna(inplace=True)
    plt.subplots_adjust(bottom=0.30)
    # Apparently boxplot gives us a better view in this case
    g = sns.boxplot(x='vaccine_short', y='expression_value', data=violin_df)
    # g = sns.violinplot(x='vaccine_short',
    #                    y='expression_value', data=violin_df)
    g.set_xticklabels(g.get_xticklabels(), rotation=90)
    g.set_xlabel('vaccine')
    return violin_df, g

