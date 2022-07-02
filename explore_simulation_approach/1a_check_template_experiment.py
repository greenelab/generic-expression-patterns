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

# # Check template experiment
#
# We want to verify that we have a good differential expression signal - i.e. control and perturb samples separate based on specific and common genes

# +
# %load_ext autoreload
# %autoreload 2
import os
import pickle
import numpy as np
import pandas as pd
import random
import seaborn as sns
from ponyo import utils

np.random.seed(1)

# +
i = 1

generic_gene_ids_filename = "generic_gene_ids.pickle"
# -

template_filename = f"/home/alexandra/Documents/Data/Generic_expression_patterns/reviewer_experiment/raw_template_{i}.tsv"

template_specific_gene_ids_filename = f"/home/alexandra/Documents/Data/Generic_expression_patterns/reviewer_experiment/template_specific_gene_ids_{i}.pickle"

template_experiment = pd.read_csv(template_filename, sep="\t", index_col=0, header=0)

# +
# Load saved specific (template-specific) and generic gene ids (same across all template experiments)
with open(template_specific_gene_ids_filename, "rb") as specific_fh:
    specific_gene_ids = pickle.load(specific_fh)

with open(generic_gene_ids_filename, "rb") as generic_fh:
    generic_gene_ids = pickle.load(generic_fh)
# -

# Get NA genes
all_gene_ids = template_experiment.columns
all_gene_ids_tmp = all_gene_ids.difference(specific_gene_ids)
na_gene_ids = all_gene_ids_tmp.difference(generic_gene_ids)

# Template data subsets
template_specific_df = template_experiment[specific_gene_ids]
template_common_df = template_experiment[generic_gene_ids]
template_na_df = template_experiment[na_gene_ids]

print(template_specific_df.shape)
template_specific_df

print(template_common_df.shape)
template_common_df

print(template_na_df.shape)
template_na_df

f = sns.clustermap(template_specific_df.T, cmap="viridis")
f.fig.suptitle("Template experiment specific genes")

f = sns.clustermap(template_common_df.T, cmap="viridis")
f.fig.suptitle("Template experiment common genes")

f = sns.clustermap(template_na_df.T, cmap="viridis")
f.fig.suptitle("Template experiment NA genes")
