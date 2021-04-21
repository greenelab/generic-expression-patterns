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

# # Comparisons against template
#
# Get some baselines by comparing gene ranking of:
# * Template experiment vs Crow et al
# * Template experiment vs SOPHIE

# +
# %load_ext autoreload
# %autoreload 2

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from ponyo import utils
from generic_expression_patterns_modules import process, stats, ranking

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general.tsv")
)

params = utils.read_config(config_filename)
# -

# Load params
local_dir = params["local_dir"]
project_id = params["project_id"]
dataset_name = params["dataset_name"]

# +
# Read in summary results
gene_summary_filename = os.path.join(
    base_dir, dataset_name, f"generic_gene_summary_{project_id}_noise_model.tsv"
)

gene_summary = pd.read_csv(gene_summary_filename, sep="\t", index_col=0, header=0)
# -

gene_summary.head()

# ## Template vs Crow et al

# +
# Scale ranking of genes in template experiment
gene_summary_scaled = ranking.scale_reference_ranking(gene_summary, "Rank (Real)")

gene_summary_scaled.head()

# +
# Get scaled Crow et al ranking from summary df
DE_prior_filename = params["reference_gene_filename"]
ref_gene_col = params["reference_gene_name_col"]
ref_rank_col = params["reference_rank_col"]

figure_filename = "tmp.svg"

corr, shared_ranking = ranking.compare_gene_ranking(
    gene_summary, DE_prior_filename, ref_gene_col, ref_rank_col, figure_filename
)
shared_ranking.head()

# +
# Merge template and Crow et al ranking using gene names
gene_summary_merged = gene_summary_scaled.merge(
    shared_ranking, left_on="Gene ID", right_on="Gene_Name"
)

gene_summary_merged.head()

# +
# Make joint plot
fig = sns.jointplot(
    data=gene_summary_merged,
    x="Rank (Real)",
    y=ref_rank_col,
    kind="hex",
    marginal_kws={"color": "white", "edgecolor": "white"},
)

cbar_ax = fig.fig.add_axes([0.9, 0.25, 0.05, 0.4])  # x, y, width, height
cb = plt.colorbar(cax=cbar_ax)
cb.set_label("Number of genes")

fig.set_axis_labels(
    "Template experiment",
    "DE prior (Crow et. al. 2019)",
    fontsize=14,
    fontname="Verdana",
)
# -

# Looks like template experiment already expresses generic genes, so adding a small amount of noise (Normal(0,2)) will still find these generic results.
#
# This generic pattern does seem to be disrupted when we add enough noise (Normal(0.1000))

# ## Template vs SOPHIE

# +
# Make joint plot
fig = sns.jointplot(
    data=gene_summary_merged,
    x="Rank (Real)",
    y="Percentile (simulated)_x",
    kind="hex",
    marginal_kws={"color": "white", "edgecolor": "white"},
)

cbar_ax = fig.fig.add_axes([0.9, 0.25, 0.05, 0.4])  # x, y, width, height
cb = plt.colorbar(cax=cbar_ax)
cb.set_label("Number of genes")

fig.set_axis_labels("Template experiment", "SOPHIE", fontsize=14, fontname="Verdana")
# -

# There is a distinct sigmoid shape -- highly changed genes in the template experiment are pushed to having a higher rank in SOPHIE while low/unchanged genes in template experiment are pushed to having a low rank in SOPHIE
#
# SOPHIE pushes genes to extreme ranks -- a characteristic of the VAE
