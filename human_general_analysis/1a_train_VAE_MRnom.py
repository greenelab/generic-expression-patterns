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

# ## Train VAE on MR normalized data

# %load_ext autoreload
# %load_ext rpy2.ipython
# %autoreload 2
import os
import pandas as pd
import seaborn as sns
from ponyo import utils, train_vae_modules
from generic_expression_patterns_modules import process

# Set seeds to get reproducible VAE trained models
process.set_all_seeds()

# +
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_general_MRnorm.tsv")
)

params = utils.read_config(config_filename)

local_dir = params["local_dir"]
dataset_name = params["dataset_name"]

# File that contains gene ranks identified by Crow et. al.
DE_prior_filename = params["reference_gene_filename"]

# Template experiment ID
project_id = params["project_id"]

# Output files of recount2 template experiment data
# processed_template_filename = params["processed_template_filename"]

# Output files of recount2 compendium data
MRnormalized_compendium_filename = params["MRnormalized_compendium_filename"]
normalized_compendium_filename = params["normalized_compendium_filename"]


# Output file: pickled scaler (generated during compendium normalization)
scaler_filename = params["scaler_filename"]

# Output: size factor for MR normalization
sf_filename = "data/metadata/MR_norm_compendium_size_factor.tsv"

# +
# Create VAE directories if needed
output_dirs = [
    os.path.join(base_dir, dataset_name, "models"),
    os.path.join(base_dir, dataset_name, "logs"),
]

NN_architecture = params["NN_architecture"]

# Check if NN architecture directory exist otherwise create
for each_dir in output_dirs:
    sub_dir = os.path.join(each_dir, NN_architecture)
    os.makedirs(sub_dir, exist_ok=True)
# -

# Train VAE on new compendium data
train_vae_modules.train_vae(config_filename, normalized_compendium_filename)
