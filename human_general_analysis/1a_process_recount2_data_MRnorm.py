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

# # Process recount2 data
# This notebook does the following:
#
# 1. Select template experiment. This template experiment will be used in the next [notebook](2_identify_generic_genes_pathways.ipynb) to simulate experiments with the same experimental design but testing a different biological process.
#
#
# 2. Uses pre-downloaded data from [notebook]()
#
# 3. Normalizes data using MRnorm
# Check what the distribution looks like
#
# 4. Train VAE on recount2 data

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

# ### Set parameters for data processing
#
# Most parameters are read from `config_filename`. We manually selected bioproject [SRP012656](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37764) as the template experiment, which contains primary non-small cell lung adenocarcinoma tumors and adjacent normal tissues of 6 never-smoker Korean female patients with 2 replicates each.

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
normalized_compendium_filename = params["normalized_compendium_filename"]

# Output file: pickled scaler (generated during compendium normalization)
scaler_filename = params["scaler_filename"]
# -

# ### Load template data file

mapped_template_filename = params["mapped_template_filename"]

mapped_template = pd.read_csv(mapped_template_filename, sep="\t", index_col=0, header=0)
mapped_template.head()

# ### Load recount2

mapped_compendium_filename = params["mapped_compendium_filename"]

# +
# mapped_compendium = pd.read_csv(mapped_compendium_filename, sep="\t", index_col=0, header=0)
# -

# ### MR normalize

metadata = pd.DataFrame(
    data=[i for i in range(mapped_template.shape[0])],
    index=mapped_template.index,
    columns=["group"],
)
metadata.head()

metadata_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", "MRnorm_metadata.tsv"
)
metadata.to_csv(metadata_filename, sep="\t")

# + magic_args="-i base_dir -i mapped_template_filename -i metadata_filename -i normalized_compendium_filename -i scaler_filename" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/normalize_counts.R'))
#
# MRnorm_expression(mapped_template_filename, metadata_filename, normalized_compendium_filename, scaler_filename)
# -

# Read in MR normalized data
MRnorm_expression_data = pd.read_csv(
    normalized_compendium_filename, sep="\t", index_col=0, header=0
)

scale_factor = pd.read_csv(scaler_filename, sep="\t", index_col=0, header=0)

scale_factor

MRnorm_expression_data.head()

sns.distplot(MRnorm_expression_data["SRR493937"])
sns.distplot(MRnorm_expression_data["SRR493938"])

# +
# Check that we're normalizing by gene
# Check that dist looks similar across genes
# Determine how to re-scale using size factors
# May need to do a hack job for ponyo
# -

# ### Train VAE
# Performed exploratory analysis of compendium data [here](../explore_data/viz_recount2_compendium.ipynb) to help interpret loss curve.

"""# Create VAE directories if needed
output_dirs = [
    os.path.join(base_dir, dataset_name, "models"),
    os.path.join(base_dir, dataset_name, "logs"),
]

NN_architecture = params["NN_architecture"]

# Check if NN architecture directory exist otherwise create
for each_dir in output_dirs:
    sub_dir = os.path.join(each_dir, NN_architecture)
    os.makedirs(sub_dir, exist_ok=True)"""

# +
# Train VAE on new compendium data
# train_vae_modules.train_vae(config_filename, normalized_compendium_filename)
