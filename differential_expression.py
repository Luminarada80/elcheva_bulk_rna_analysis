import os
import re
from typing import Tuple, Union

import gffutils
import gseapy
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sanbomics.tools import id_map

results_dir = "/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA"

metadata_file = "/gpfs/Labs/Uzun/METADATA/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/Sample_Annotation.BATCH_2025_01.txt"
comparison_file = "/gpfs/Labs/Uzun/METADATA/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/Comparisons.BATCH_2025_01.txt"
annotation_file = "/gpfs/Labs/Uzun/RESULTS/GRANT_APPS/2025.NYNRIN.DOD.ELCHEVA/gene_annotation_lookup.parquet"

gene_count_dir = os.path.join(results_dir, "gene_counts")

comparison_df = pd.read_csv(comparison_file, sep="\t", header=0, index_col=0)
metadata = pd.read_csv(metadata_file, sep="\t", header=0, index_col=0)

def build_sample_expression_matrix(gene_count_dir: str):
    """
    Creates a sample x Ensembl gene ID count matrix for each sample in `gene_count_dir`.
    
    Args:
        gene_count_dir (str): 
            Directory containing the gene_counts.txt files for each sample. Assumes the file
            name follows the <sample_name>_gene_counts.txt format.

    Returns:
        pd.DataFrame: 
            DataFrame object of the gene count matrix, with sample names as the index and the 
            Ensembl gene ids as the columns.
    """
    gene_count_df_list = []
    for file in os.listdir(gene_count_dir):
        sample_name = file.replace("_gene_counts.txt", "")     
        filepath = os.path.join(gene_count_dir, file)

        test_gene_count_df = pd.read_csv(
            filepath, 
            sep="\t", 
            header=0, 
            skiprows=1,
            )
        test_gene_count_df = test_gene_count_df.rename(columns={
            test_gene_count_df.columns[-1] : sample_name
            })
        test_gene_count_df["Geneid"] = test_gene_count_df["Geneid"].str.replace(r'\.\d+$', '', regex=True)
        test_gene_count_df = test_gene_count_df[["Geneid", sample_name]].set_index("Geneid")
        gene_count_df_list.append(test_gene_count_df)
    
    counts = pd.concat(gene_count_df_list, axis=1, join="inner")
    
    # Filter to only keep rows with at least one gene expressed
    counts = counts[counts.sum(axis=1) > 0]
    
    # Transpose to have genes as columns and samples as rows
    counts = counts.T
    
    return counts

def check_expression_matrix_samples_match_metadata(
    expr_matrix_df: pd.DataFrame, 
    metadata: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Checks to ensure that the index for the gene expression DataFrame and the metadata contain the
    same samples. Reindexes the expr_matrix_df matrix to use the samples in the metadata index.

    Args:
        expr_matrix_df (pd.DataFrame): 
            DataFrame of gene expression counts where the index contains the sample 
            names and the columns contain the Ensembl gene IDs.
        metadata (pd.DataFrame): 
            DataFrame of the samples and their corresponding groups, with the sample name as the 
            index.

    Returns:
        expr_matrix_df,metadata: 
            Returns the expr_matrix_df and metadata DataFrames back, with the expr_matrix_df DataFrame reindexed
            according the the sample names in the metadata DataFrame.
    """
    
    expr_matrix_df.index = expr_matrix_df.index.str.strip()
    metadata.index = metadata.index.str.strip()

    missing_in_meta   = expr_matrix_df.index.difference(metadata.index)
    missing_in_counts = metadata.index.difference(expr_matrix_df.index)
    if len(missing_in_meta) or len(missing_in_counts):
        print("In counts but not metadata:", list(missing_in_meta))
        print("In metadata but not counts:", list(missing_in_counts))
        common = expr_matrix_df.index.intersection(metadata.index)
        expr_matrix_df   = expr_matrix_df.loc[common]
        metadata = metadata.loc[common]
    expr_matrix_df = expr_matrix_df.reindex(metadata.index)
    assert expr_matrix_df.index.equals(metadata.index)
    
    return expr_matrix_df, metadata

def convert_ensembl_id_to_gene_symbol(deseq_stat_df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds a "Symbol" column to the DataFrame containing the gene symbol corresponding to
    the Ensembl ID in the index. Assumes the species is human.

    Args:
        deseq_stat_df (pd.DataFrame): DataFrame where the index contains Ensembl gene IDs

    Returns:
        pd.DataFrame: 
            DataFrame containing a new "Symbol" column with the gene symbol corresponding to the
            Ensembl gene ID for the row.
    """
    mapper = id_map(species = 'human')
    mapper.mapper
    deseq_stat_df['Symbol'] = deseq_stat_df.index.map(mapper.mapper)
    
    deseq_stat_df = deseq_stat_df.sort_values("pvalue").drop_duplicates("Symbol", keep="first")
    
    return deseq_stat_df

def filter_out_non_protein_coding_genes(expr_matrix_df: pd.DataFrame) -> pd.DataFrame:
    """
    Uses a preprocessed gene annotation DataFrame to remove non-protein coding genes from the 
    gene expression DataFrame. 

    Args:
        expr_matrix_df (pd.DataFrame): Ensembl gene ID x sample gene expression counts matrix.

    Returns:
        pd.DataFrame: Filtered DataFrame with non-protein coding genes removed.
    """
    annotation_df = pd.read_parquet(annotation_file)

    # Get the set of protein-coding Ensembl IDs (no version)
    protein_ids = (
        annotation_df.loc[annotation_df["gene_type"].str.lower() == "protein_coding", :]
        .assign(ensg_nover=annotation_df["ensg_nover"].astype(str))
        .drop_duplicates("ensg_nover")["ensg_nover"]
    )

    # Subset the expr_matrix_df matrix to protein-coding genes only
    # (expr_matrix_df has samples as rows, genes as columns)
    expr_matrix_df = expr_matrix_df.loc[:, expr_matrix_df.columns.intersection(protein_ids)]

    return expr_matrix_df

def plot_significant_de_gene_heatmap(graph_df):
    mat = graph_df.copy()

    if isinstance(mat.columns, pd.MultiIndex):
        mat.columns = mat.columns.get_level_values(-1)

    # Split up the column names for the samples to get the cell type, gene target, group, and replicate num
    cols = mat.columns.astype(str)
    meta = cols.to_series().str.extract(
        r'^(?P<cell>[^_]+)_(?P<group>[^_]+)_(?P<rep>Rep\d+)$'
    )
    meta.index = cols
    meta['gene_target'] = meta['group']
    meta['group'] = meta['group'].replace({'shCNTL': 'Control'})  # nicer label
    meta['rep_num'] = meta['rep'].str.extract(r'Rep(\d+)').astype(int)

    # MultiIndex columns: level0 = Group, level1 = Rep
    mi = pd.MultiIndex.from_arrays([meta['group'], meta['rep']],
                                names=['Group','Rep'])
    mat.columns = mi

    # Set the order to display the samples based on their groups
    group_order = ['Control', meta['group'][1]]
    rep_order   = [f'Rep{i}' for i in sorted(meta['rep_num'].unique())]

    # Create the multiindex for the columns and reindex mat
    wanted = pd.MultiIndex.from_product([group_order, rep_order])
    mat = mat.reindex(columns=[c for c in wanted if c in mat.columns])

    # Safeguard to make sure groups exist before doing mat.xs
    present_groups = set(mat.columns.get_level_values('Group'))
    missing = set(group_order) - present_groups
    if missing:
        raise ValueError(f"Missing groups in matrix: {missing}")


    mu_ctrl  = mat.xs(group_order[0],   axis=1, level='Group').mean(axis=1)
    mu_treat = mat.xs(group_order[1], axis=1, level='Group').mean(axis=1)

    delta = (mu_treat - mu_ctrl)
    delta = delta.fillna(0)

    N = 15 # Number of top up- and down-regulated genes to plot
    top_up = delta.nlargest(N)
    top_dn = delta.nsmallest(N)

    keep_rows = pd.Index(top_up.index.tolist() + top_dn.index.tolist())
    mat_top = mat.loc[keep_rows]

    order = delta.loc[keep_rows].sort_values(ascending=False).index
    mat_sorted  = mat.loc[order]

    g = sns.clustermap(mat_sorted, z_score=0, cmap='RdYlBu_r',
                    row_cluster=False, col_cluster=False,
                    figsize=(6,8))

    # --- get group/rep labels from columns (MultiIndex OR strings) ---
    cols = g.data2d.columns

    if isinstance(cols, pd.MultiIndex):
        # use named levels if present, else 0/1
        group_level = 'Group' if 'Group' in cols.names else 0
        rep_level   = 'Rep'   if 'Rep'   in cols.names else 1
        groups = pd.Series(cols.get_level_values(group_level)).reset_index(drop=True)
        reps   = pd.Series(cols.get_level_values(rep_level)).reset_index(drop=True)
    else:
        cols_s = pd.Index(cols.astype(str))
        groups = cols_s.to_series().str.extract(r'^[^_]+_([^_]+)_')[0].reset_index(drop=True)
        reps   = cols_s.to_series().str.extract(r'(Rep\d+)$')[0].reset_index(drop=True)

    # optional rename for readability
    groups = groups.replace({'shCNTL': 'Control'})

    # if columns were clustered, honor new order
    if getattr(g, "dendrogram_col", None) is not None:
        order = g.dendrogram_col.reordered_ind
        groups = groups.iloc[order].reset_index(drop=True)
        reps   = reps.iloc[order].reset_index(drop=True)

    # set per-replicate tick labels

    # --- add group labels centered under their blocks ---
    def add_group_labels_under_xticks(grid, group_series, pad=0.07, fontsize=13):
        ax = grid.ax_heatmap
        n = len(group_series)
        trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        start = 0
        for i in range(1, n + 1):
            if i == n or group_series[i] != group_series[start]:
                center = (start + i) / 2.0   # centers at 0..n-1
                ax.text(center, -pad, str(group_series[start]),
                        ha='center', va='top', transform=trans,
                        fontsize=fontsize, clip_on=False)
                start = i

    add_group_labels_under_xticks(g, groups)

    g.cax.set_position([0.05, 0.20, 0.05, 0.55])
    g.cax.tick_params(labelsize=10, length=0)

    # x-label and x-tick settings
    g.ax_heatmap.set_xticklabels(reps, rotation=0, fontsize=11)
    g.ax_heatmap.set_xlabel("")

    # y-label and y-tick settings
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=9)
    g.ax_heatmap.set_ylabel("")

    # Colorbar settings
    x0, _y0, _w, _h = g.cbar_pos
    g.ax_cbar.set_title('Differential\nExpression\n(log1p)')

    # Turn off dendrograms
    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)

    g.ax_heatmap.set_title(
        f"{meta['cell'][0]}\n{meta['gene_target'][0]} vs {meta['gene_target'][1]}\nSignificant DE Genes", 
        fontsize=14)

    out = os.path.join(
        results_dir,
        f"differential_expr_plots/{meta['cell'][0]}_{meta['gene_target'][0]}_vs_{meta['gene_target'][1]}.png"
    )
    g.figure.savefig(out, dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(g.figure)
    
if __name__ == "__main__":
    expr_matrix_df = build_sample_expression_matrix(gene_count_dir)
    
    expr_matrix_df = filter_out_non_protein_coding_genes(expr_matrix_df)
    
    expr_matrix_df, metadata = check_expression_matrix_samples_match_metadata(expr_matrix_df, metadata)
    
    # Only keep genes that are expressed in more than 2 samples
    group = metadata.loc[expr_matrix_df.index, "Group_Name"]
    keep = (expr_matrix_df >= 10).groupby(group, axis=0).any().sum(axis=0) >= 2
    expr_matrix_df = expr_matrix_df.loc[:, keep.index[keep]]
    
    # Build a DESeq dataset using the counts matrix and the metadata groups
    dds = DeseqDataSet(
        counts=expr_matrix_df, 
        metadata=metadata,
        design="Group_Name"
        )
    
    dds.deseq2()
    
    # Use the metadata file to find the different comparisons we want to run
    comparisons = [comparison_df.iloc[i, :].to_list() for i in range(len(comparison_df))]
    
    for i, comparison in enumerate(comparisons):
    
        deseq_stat_results = DeseqStats(dds, n_cpus=8, contrast = ('Group_Name', comparisons[i][0], comparisons[i][1]))
        deseq_stat_results.summary()
        deseq_stat_df = deseq_stat_results.results_df
        
        deseq_stat_df = convert_ensembl_id_to_gene_symbol(deseq_stat_df)
        
        # Isolate genes with a significant differential expression
        lfc_cut = np.log2(1.5)
        sig_de_df = deseq_stat_df[(deseq_stat_df.padj < 0.05) & (deseq_stat_df.log2FoldChange.abs() > lfc_cut)]
        
        # Rank the genes based on their stat value
        ranking = deseq_stat_df[['Symbol', 'stat']].dropna().sort_values('stat', ascending = False).drop_duplicates()

        # Log1p normalize the normalized counts, save as a separate layer
        dds.layers['log1p'] = np.log1p(dds.layers['normed_counts'])
        
        # Extract only the significant genes from the DeseqDataSet
        dds_sigs = dds[:, sig_de_df.index]
        
        # Create a DataFrame of the log1p expression values for significant DE genes for a comparison
        grapher = pd.DataFrame(dds_sigs.layers['log1p'].T,
                            index=dds_sigs.var_names, columns=dds_sigs.obs_names)
        
        # Prep the graph DataFrame
        mapper = id_map(species = 'human')
        mapper.mapper
        grapher.index = grapher.index.map(mapper.mapper)
        base = grapher.columns.astype(str).str.split('_Rep', n=1).str[0]
        mask = base.isin(comparisons[i])
        grapher = grapher.loc[:, mask]
        grapher = grapher.sort_index(axis=0)
        
        # plot a heatmap of the top 15 significant DE genes for the coomparison
        plot_significant_de_gene_heatmap(grapher)
    
    