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
#     display_name: Python [conda env:generic_expression_new] *
#     language: python
#     name: conda-env-generic_expression_new-py
# ---

# # Process and train vae
# This notebook inputs a compendium of gene expression data and trains a VAE model that will be used to simulate new experiments.
#
# The output of this notebook, which includes the vae models and template experiment, will be used in the next notebook

# +
# %load_ext autoreload
# %autoreload 2

import os
import sys
import pandas as pd
import numpy as np
from sklearn import preprocessing
import pickle

from ponyo import utils, train_vae_modules, simulate_expression_data
from generic_expression_patterns_modules import process
# -

# ## User inputs needed
#
# User needs to define the following in the [config file](../configs/config_new_model_experiment.tsv):
#
# 1. Directory on your local machine to store intermediate and output data files generated (`local_dir`). Make sure to end with `\`.
# 2. Template experiment (`raw_template_filename`). This is the experiment you are interested in studying. This experiment is expected to be a matrix with samples as row and genes as columns (tab-delimited).
# 3. Training compendium used to train VAE (`processed_compendium_filename`). This dataset is expected to be a matrix with samples as row and genes as columns (tab-delimited). Note: if using human gene ids from ensembl and you want to convert these to HGNC symbols, functions are available to do this in `generic_expression_patterns_modules/process_names.R` and `generic_expression_patterns_modules/process.py`. See [example](../human_general_analysis/1_process_recount2_data.ipynb)
# 4. Scaler transform (`scaler_filename`) used to normalize the training compendium. This can be found in the `data/` directory within the analysis folder.
# 5. Directory (`vae_model_dir`) containing trained VAE model (.h5 files) from the previous notebook.
# 6. Size of the latent dimension (`latent_dim`).
# 7. File that maps experiment ids to the associated sample ids (`experiment_to_sample_filename`)
# 8. The delimiter used in the 'experiment_to_sample_filename' file (`metadata_delimiter`)
# 9. The column header/name that contains the experiment ids (`experiment_id_colname`)
# 10. Experiment id (`project_id`) to label newly create simulated experiments.
# 11. The column header/name in the metadatathat contains the sample ids (`sample_id_colname`)
#
# The remaining parameters within the `config` file specify values needed to run the next notebook or filenames that are intermediate data files that will be generated when SOPHIE runs.

# Set seeds to get reproducible VAE trained models
process.set_all_seeds()

# ## Set parameters for data processing

# +
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_new_model_experiment.tsv")
)

params = utils.read_config(config_filename)

local_dir = params["local_dir"]
dataset_name = params["dataset_name"]

# Column header containing sample ids
metadata_colname = params["metadata_colname"]

# Template experiment ID
project_id = params["project_id"]

# Output file: pickled list of shared genes(generated during gene ID mapping)
shared_genes_filename = params["shared_genes_filename"]

# Output files of pseudomonas template experiment data
raw_template_filename = params["raw_template_filename"]
processed_template_filename = params["processed_template_filename"]

# Output files of compendium data
processed_compendium_filename = params["processed_compendium_filename"]
normalized_compendium_filename = params["normalized_compendium_filename"]

# Output file: pickled scaler (generated during compendium normalization)
scaler_filename = params["scaler_filename"]

# Load metadata file with mapping between experiments and associated samples
metadata_filename = params["experiment_to_sample_filename"]
metadata_delimiter = params["metadata_delimiter"]
experiment_id_colname = params["experiment_id_colname"]
# -

# ## Normalize compendium
# Here we will 0-1 normalize expression data

process.normalize_compendium(
    processed_compendium_filename,
    normalized_compendium_filename,
    scaler_filename,
)

# ## Get raw template experiment

# +
# Get sample ids associated with selected project id
sample_ids = simulate_expression_data.get_sample_ids(
    metadata_filename,
    metadata_delimiter,
    experiment_id_colname,
    project_id,
    sample_id_colname,
)

# Get samples from experiment id
processed_compendium = pd.read_csv(
    processed_compendium_filename, header=0, index_col=0, sep="\t"
)
template_data = processed_compendium.loc[sample_ids]

template_data.to_csv(raw_template_filename, sep="\t")
# -

# ### Train VAE

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
