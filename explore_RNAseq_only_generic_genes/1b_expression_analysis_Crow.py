# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.9.1+dev
#   kernelspec:
#     display_name: Python [conda env:generic_expression] *
#     language: python
#     name: conda-env-generic_expression-py
# ---

# # Expression of Crow data
#
# This notebook tests the hypothesis that the RNA-seq generic genes are those that are not well captured on microarray technology. [Zhao et. al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3894192/) found that one possible reason for differences in differentially expressed genes detected on different platforms is due to how genes are captured. RNA-Seq is more likely to detect the changes at two different conditions for genes with very low expression or very high expression compared to arrays.
#
# This data was generated by running `download_Crow_data.R` script that downloads expression data from https://github.com/PavlidisLab/gemmaAPI.R

# %load_ext autoreload
import os
import random
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from ponyo import utils
from generic_expression_patterns_modules import ranking

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general.tsv")
)

params = utils.read_config(config_filename)
local_dir = params["local_dir"]
dataset_name = params["dataset_name"]
project_id = params["project_id"]
col_to_rank_genes = params["rank_genes_by"]
mapped_compendium_filename = params["mapped_compendium_filename"]

if_single_experiment = False
# -

# Read in recount2 expression compendium
recount2_expression = pd.read_csv(
    mapped_compendium_filename, sep="\t", index_col=0, header=0
).T

# ## Format Crow expression data
#
# * Include only genes that were used in our original analysis

# Read in Crow expression data
crow_expression_filename = os.path.join(local_dir, "Crow_expression_data_union.tsv")
crow_expression_data = pd.read_csv(
    crow_expression_filename, sep="\t", index_col=0, header=0
).T

crow_expression_data.shape

crow_expression_data.head()

# +
# Load gene_summary_filename
gene_summary_filename = os.path.join(
    base_dir, dataset_name, f"generic_gene_summary_{project_id}.tsv"
)

summary_gene_ranks = pd.read_csv(gene_summary_filename, sep="\t", index_col=0, header=0)
# -

summary_gene_ranks.head()

# +
# Subset genes
our_gene_ids = list(summary_gene_ranks.index)
crow_gene_ids = list(crow_expression_data.index)

shared_gene_ids = set(crow_gene_ids).intersection(our_gene_ids)

expression_data = crow_expression_data.loc[shared_gene_ids]
# -

print(expression_data.shape)
expression_data.head()

# ## (optional) Select gene subset of samples
#
# Select samples from a specific experiment to examine local gene expression behavior within a single experiment (local) in addition to across all samples (global)
#
# We would actually like to do this for crow data but, we do not have metadata mapping experiment ids to sample ids. So this analysis option is only available for recount2 for now

recount2_metadata_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", "recount2_metadata.tsv"
)


# Function scraped from ponyo since we're already using a different version of ponyo in this repo
def get_sample_ids_random_experiment(
    metadata_filename, delimiter, experiment_colname, sample_id_colname, rn_seed
):
    """
    Returns sample ids (found in gene expression df) associated with
    a given list of experiment ids (found in the metadata)

    Arguments
    ----------
    metadata_filename: str
        Metadata file path. An example metadata file can be found
        here: https://github.com/greenelab/ponyo/blob/master/human_tests/data/metadata/recount2_metadata.tsv

    delimiter: str
        Delimiter for metadata file

    experiment_colname: str
        Column header that contains the experiment ids

    sample_id_colname: str
        Column header that contains sample id that maps expression data
        and metadata

    """
    random.seed(rn_seed)

    # Read in metadata
    metadata = pd.read_csv(metadata_filename, header=0, sep=delimiter, index_col=None)

    # Set index column to experiment id column
    metadata.set_index(experiment_colname, inplace=True)

    # Select random experiment
    rn_experiment_id = random.choice(list(np.unique(metadata.index)))

    # Select samples associated with experiment id
    selected_metadata = metadata.loc[rn_experiment_id]
    sample_ids = list(selected_metadata[sample_id_colname])

    return sample_ids


if if_single_experiment:
    # Get sample ids for random experiment
    recount2_sample_ids = get_sample_ids_random_experiment(
        recount2_metadata_filename, "\t", "project", "run", 1
    )

    # Subset expression data
    recount2_expression = recount2_expression.loc[recount2_sample_ids]

# ## Get uncorrelated genes

# +
# Get generic genes identified by Crow et. al.
DE_prior_filename = params["reference_gene_filename"]
ref_gene_col = params["reference_gene_name_col"]
ref_rank_col = params["reference_rank_col"]

figure_filename = f"gene_ranking_{col_to_rank_genes}_tmp.svg"

corr, shared_ranking = ranking.compare_gene_ranking(
    summary_gene_ranks,
    DE_prior_filename,
    ref_gene_col,
    ref_rank_col,
    figure_filename,
)
# -

shared_ranking.head()

# +
# Get uncorrelated gene ids
uncorrelated_ranking = shared_ranking[
    (shared_ranking["Percentile (simulated)"] > 80)
    & (shared_ranking["DE_Prior_Rank"] < 20)
]

uncorrelated_genes = uncorrelated_ranking["Gene_Name"]
print(len(uncorrelated_genes))

# +
# Get correlated gene ids
correlated_ranking = shared_ranking[
    (shared_ranking["Percentile (simulated)"] > 80)
    & (shared_ranking["DE_Prior_Rank"] > 80)
]

correlated_genes = correlated_ranking["Gene_Name"]
print(len(correlated_genes))
# -

# Save uncorrelated genes
uncorrelated_genes.to_csv("uncorrelated_genes.tsv", sep="\t")

# ## Plot average expression

# Get average expression of SOPHIE trained recount2 dataset
recount2_expression_mean = recount2_expression.mean(axis=1)

recount2_expression_mean.head()

# Get average expression of Crow dataset
crow_expression_mean = crow_expression_data.mean(axis=1)

crow_expression_mean.head()

# Check that we selecting the correct genes
uncorrelated_genes = list(uncorrelated_genes.values)
uncorrelated_genes[0:5]

recount2_expression_mean[uncorrelated_genes].head()

crow_expression_mean.reindex(uncorrelated_genes).head()

recount2_expression_mean.head()

# +
# Format df for plotting
recount2_expression_mean_toplot = pd.DataFrame(
    data={
        "all genes": np.log10(recount2_expression_mean),
        "RNAseq and array generic": np.log10(
            recount2_expression_mean[correlated_genes]
        ),
        "only RNAseq generic": np.log10(recount2_expression_mean[uncorrelated_genes]),
    }
)


recount2_expression_mean_toplot.head()
# -

# Violin plot of average recount2 expression highlighing uncorrelated genes
print(
    f"Number of uncorrelated gene data available: {len(recount2_expression_mean[uncorrelated_genes])}"
)
f = sns.violinplot(
    data=recount2_expression_mean_toplot,
    palette=["lightgrey", "#2c7fb8", "#add8e6"],
    orient="h",
)
f.set_title("Average recount2 expression")
f.set_xlabel("log10(average expression)")

# +
# Format df for plotting
crow_expression_mean_toplot = pd.DataFrame(
    data={
        "all genes": np.log10(crow_expression_mean),
        "RNAseq and array generic": np.log10(crow_expression_mean[correlated_genes]),
        "only RNAseq generic": np.log10(crow_expression_mean[uncorrelated_genes]),
    }
)


crow_expression_mean_toplot.head()
# -

# Violin plot of average array expression highlighing uncorrelated genes
print(
    f"Number of uncorrelated gene data available: {len(crow_expression_mean.reindex(uncorrelated_genes))}"
)
g = sns.violinplot(
    data=crow_expression_mean_toplot,
    palette=["lightgrey", "#2c7fb8", "#add8e6"],
    orient="h",
)
g.set_title("Average array (Crow et. al.) expression")
g.set_xlabel("log10(average expression)")

g.get_figure().savefig(
    "array_expression_dist_gene_groups_highlight.svg",
    format="svg",
    bbox_inches="tight",
    transparent=True,
    pad_inches=0,
    dpi=300,
)

# **Takeaway:**
# * Our hypothesis is that these RNA-seq generic genes are those that are not well captured on microarray technology.
# * Based on the distribution of the array data, it looks like these genes are fairly lowly expressed, but based on the density of the violin plot there appear to be many genes that have a similar range of expression. So these RNA-seq generic genes are as well captured as RNA-seq/array generic genes.
# * This hypothesis does not seem to hold

# **Other thoughts:**
#
# Looking to characterize _who_ these RNA-seq generic genes are, we used https://academic.oup.com/nar/article/48/D1/D174/5588346 to lookup the RNA-seq generic genes to determine if they have 3' end processing (i.e. polyadenylation sites)
#
# Manually looking up individual genes (since there doesn't seem to be a way to do this in batch), we found that most genes have at least one polyadenylated site.
