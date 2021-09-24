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

# # Process pseudomonas data
# This notebook trains a VAE on the PAO1 and PA14 RNA-seq compendia to support the analysis in https://github.com/greenelab/core-accessory-interactome
#
# 1. Selects template experiment from the compendium
# 2. Normalizes the gene expression data from the RNA-seq Pseudomonas compendium
# 3. Train VAE on the normalized data

# +
# %load_ext autoreload
# %autoreload 2

import os
import sys
import pandas as pd
import numpy as np
from sklearn import preprocessing
import pickle

from ponyo import utils, train_vae_modules
from generic_expression_patterns_modules import process
# -

# Set seeds to get reproducible VAE trained models
process.set_all_seeds()

# ### Set parameters for data processing
#
# Most parameters are read from `config_filename`. We manually selected bioproject [GEOD-33245](https://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-33245/?s_sortby=col_8&s_sortorder=ascending), as the template experiment, which contains multiple different comparisons including WT vs *crc* mutants, WT vs *cbr* mutants in different conditions.

# +
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_pseudomonas_pa14_rnaseq.tsv")
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

# Output files of pseudomonas compendium data
# raw_compendium_filename = params['raw_compendium_filename']
processed_compendium_filename = params["processed_compendium_filename"]
normalized_compendium_filename = params["normalized_compendium_filename"]

# Output file: pickled scaler (generated during compendium normalization)
scaler_filename = params["scaler_filename"]

# Load metadata file with grouping assignments for samples
sample_id_metadata_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", f"{project_id}_process_samples.tsv"
)
# -

# ### Normalize compendium
# The compendium is MR normalized, but here we will 0-1 normalize it

process.normalize_compendium(
    processed_compendium_filename,
    normalized_compendium_filename,
    scaler_filename,
)
"""process.process_raw_compendium_pseudomonas(
    raw_compendium_filename,
    processed_compendium_filename,
    normalized_compendium_filename,
    scaler_filename,
)"""

test = pd.read_csv(normalized_compendium_filename, sep="\t", index_col=0, header=0)

test.head()


# ## Get raw pseudomonas template experiment

# +
# The function to pull out the template experiment from the compendium in this environment's version of ponyo
# doesn't allow us to pass in a metadata file and instead assumes a fixed set of metadata files.
# To manually work around this, we will locally define the functions here
def get_sample_ids_tmp(
    metadata_filename, delimiter, experiment_colname, experiment_id, sample_id_colname
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

    experiment_id: str
        Experiment id selected to retrieve sample ids for

    sample_id_colname: str
        Column header that contains sample id that maps expression data
        and metadata

    """

    # Read in metadata
    metadata = pd.read_csv(metadata_filename, header=0, sep=delimiter, index_col=None)

    # Set index column to experiment id column
    metadata.set_index(experiment_colname, inplace=True)

    # Select samples associated with experiment id
    selected_metadata = metadata.loc[experiment_id]
    sample_ids = list(selected_metadata[sample_id_colname])

    return sample_ids


def process_raw_template_pseudomonas_tmp(
    processed_compendium_filename,
    project_id,
    metadata_filename,
    metadata_delimiter,
    experiment_colname,
    metadata_colname,
    raw_template_filename,
):
    """
    Create processed pseudomonas template data file based on
    processed compendium file (`compendium_filename`),
    drop sample rows if needed, and save updated
    template data on disk.
    """

    # Get sample ids associated with selected project id
    sample_ids = get_sample_ids_tmp(
        metadata_filename,
        metadata_delimiter,
        experiment_colname,
        project_id,
        metadata_colname,
    )

    # Get samples from experiment id
    processed_compendium = pd.read_csv(
        processed_compendium_filename, header=0, index_col=0, sep="\t"
    )
    template_data = processed_compendium.loc[sample_ids]

    template_data.to_csv(raw_template_filename, sep="\t")


# +
metadata_filename = "data/metadata/SraRunTable.csv"

process_raw_template_pseudomonas_tmp(
    processed_compendium_filename,
    project_id,
    metadata_filename,
    ",",
    "SRA_study",
    metadata_colname,
    raw_template_filename,
)
# -

test2 = pd.read_csv(raw_template_filename, sep="\t", index_col=0, header=0)

test2.head()

test2.shape

# **Note:**
# * We are training our VAE model using ALL the data in the compendium.
# * The template experiment is using a subset of the samples in the real experiment and using those in the DE analysis in order to ensure the comparison of samples with consistent backgrounds (i.e. some experiments have samples with 3 different biological conditions and for now our statistical test is doing a binary comparison).
# * Simulated experiments are generated by shifting the template experiment (using ALL samples in the real experiment) in the latent space. Then dropping the samples to match the template experiment and perform DE analysis.
#
#
# So there is an inconsistency in the samples used to learn a low-dimensional representation and those used to calculate DE statistics. This inconsistency should not not change the simulated experiments since all samples in the template experiment are moved the same amount in the latent space. The only way for this inconsistency to effect the simulated experiments is if the low dimensional space is significantly different including all the experiment samples vs only including a subset. However, we believe that such few samples will likely not effect the space. Furthermore, the dataset used to train the VAE should be a general representation of gene expression patterns and shouldn't have to be include the template experiment.

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
