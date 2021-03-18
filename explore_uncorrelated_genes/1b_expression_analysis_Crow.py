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
# This notebook tests the hypothesis that the uncorrelated genes from
#
# This data was generated by running `download_Crow_data.R` script that downloads expression data from https://github.com/PavlidisLab/gemmaAPI.R

# %load_ext autoreload
import os
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
# -

# ## Format data
#
# * Drop extra rows so that the data frame only includes gene samples x gene symbol
# * Include only genes that were used in our original analysis

# Read in combined expression data
crow_expression_filename = os.path.join(local_dir, "Crow_expression_data.tsv")
crow_expression_data = pd.read_csv(
    crow_expression_filename, sep="\t", index_col=0, header=0
)

crow_expression_data.shape

crow_expression_data.head(20)

# +
# TO DO
# Drop extra rows

# +
# Load gene_summary_filename
gene_summary_filename = os.path.join(
    base_dir, dataset_name, f"generic_gene_summary_{project_id}.tsv"
)

summary_gene_ranks = pd.read_csv(gene_summary_filename, sep="\t", index_col=0, header=0)
# -

# TO DO
# Subset genes
"""our_gene_ids = list(summary_gene_ranks.index())
crow_gene_ids = list(crow_expression_data.columns)

shared_gene_ids = set(crow_gene_ids).intersection(our_gene_ids)

expression_data = crow_expression_data.loc[shared_gene_ids]"""

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
# -

# ## Get average expression data

# +
# Get average expression of SOPHIE trained recount2 dataset
recount2_expression = pd.read_csv(
    mapped_compendium_filename, sep="\t", index_col=0, header=0
)

recount2_expression_mean = recount2_expression.mean()

# +
# Get average expression of Crow dataset
# crow_expression_mean = crow_expression_data.mean()

# +
# TO DO
# Normalize or scale data?
# Coloring
# -

# Violin plot of average recount2 expression highlighing uncorrelated genes
f = sns.violinplot(recount2_expression_mean)
f = sns.swarmplot(recount2_expression_mean.loc[uncorrelated_genes], color="red")
plt.xscale("log")
f.set_title("Average recount2 expression with uncorrelated genes highlighted")

# Violin plot of average array expression highlighing uncorrelated genes
"""g = sns.violinplot(crow_expression_mean)
g = sns.stripplot(crow_expression_mean.loc[uncorrelated_genes], color="red")
plt.xscale("log")
g.set_title("Average array expression with uncorrelated genes highlighted")"""

# Distribution plot with distribution of Crow vs uncorrelated genes
sns.distplot(recount2_expression_mean)
sns.distplot(recount2_expression_mean.loc[uncorrelated_genes])