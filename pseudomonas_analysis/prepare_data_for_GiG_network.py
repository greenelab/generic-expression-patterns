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

# # Prepare data for G-G network
#
# This notebook prepares input data to create gene-gene co-expression network using correlation amongst eADAGE latent variables. We are using the [visualize_gene_network function from ADAGEpath](https://rdrr.io/github/greenelab/ADAGEpath/man/visualize_gene_network.html) to create and plot this gene-gene network.
#
# The goal is to create a gene-gene network, highlighting the generic genes identified by SOPHIE. In order to highlight these genes we need to provide a dataframe mapping each gene to a label (generic vs not generic). This notebook is creating `gene_color_value` argument that will be passed into the `visualize_gene_network` function found in [make_GiG_network.R](../gene-expression-modules/make_GiG_network.R)

# +
# %load_ext autoreload
# %load_ext rpy2.ipython
# %autoreload 2

import os
import pandas as pd

from ponyo import utils, simulate_expression_data

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_pseudomonas_33245.tsv")
)
params = utils.read_config(config_filename)

# +
# Load params
normalized_compendium_filename = params["normalized_compendium_filename"]

generic_genes_filename = os.path.join("data", "SOPHIE_GAPE_generic.tsv")
# -

# Create dataframe with two columns: gene id, label=1 if generic, 0 otherwise

# Read expression data
expression_data = pd.read_csv(
    normalized_compendium_filename, sep="\t", index_col=0, header=0
).T
expression_data.head()

# Read generic gene data
generic_genes_data = pd.read_csv(
    generic_genes_filename, sep="\t", index_col=0, header=0
)
generic_gene_ids = list(generic_genes_data["gene id"])

# Map generic genes
expression_data["label"] = 0
expression_data.loc[generic_gene_ids, "label"] = 1

# Truncate df
annot_df = expression_data["label"].to_frame()
annot_df.head()

# Save
annot_df.to_csv("annot_df.tsv", sep="\t")
