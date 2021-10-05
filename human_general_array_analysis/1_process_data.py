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

# # Process Crow et al. data
# This notebook does the following:
#
# 1. Select template experiment. This template experiment will be used in the next [notebook](2_identify_generic_genes_pathways.ipynb) to simulate experiments with the same experimental design but testing a different biological process.
#
# 2. Crow et al. data was downloaded using `explore_RNAseq_only_generic_genes/download_Crow_data.R` script that downloads expression data from https://github.com/PavlidisLab/gemmaAPI.R
#
# Note: For the analysis exploring the RNA-seq only common DEGs we used the union of genes per experiment, which resulted in some samples having NaNs for some samples. For this analysis we are taking the intersection so that we can remove all NaNs to train.
#
# 3. Train VAE on processed data.

# %load_ext autoreload
# %load_ext rpy2.ipython
# %autoreload 2
import os
import pandas as pd
import pickle
import numpy as np
import seaborn as sns
from ponyo import utils, train_vae_modules, simulate_expression_data
from generic_expression_patterns_modules import process

# Set seeds to get reproducible VAE trained models
process.set_all_seeds()

# ### Set parameters for data processing
#
# Most parameters are read from `config_filename`. We manually selected bioproject [GSE11352](https://www.ncbi.nlm.nih.gov/gds/?term=GSE11352[Accession]) as the template experiment, which contains breast cell lines treated with estradiol at 12H, 24H and 48H.

# +
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

# Read in config variables
config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_human_Crow.tsv")
)

params = utils.read_config(config_filename)

local_dir = params["local_dir"]
dataset_name = params["dataset_name"]

# File that contains gene ranks identified by Crow et. al.
DE_prior_filename = params["reference_gene_filename"]

# Template experiment ID
project_id = params["project_id"]

# Column name containing sample ids
metadata_colname = params["metadata_colname"]

# Output file: pickled list of shared genes(generated during gene ID mapping)
shared_genes_filename = params["shared_genes_filename"]

# Output files of template experiment data
mapped_template_filename = params["mapped_template_filename"]
processed_template_filename = params["processed_template_filename"]

# Output files of Rani's compendium data
raw_compendium_filename = params["raw_compendium_filename"]
mapped_compendium_filename = params["mapped_compendium_filename"]
normalized_compendium_filename = params["normalized_compendium_filename"]

# Output file: pickled scaler (generated during compendium normalization)
scaler_filename = params["scaler_filename"]
# -

# ### Load compendium data
#
# [Gemma](https://pavlidislab.github.io/Gemma/) contains nearly 4K expression profiling studies/datasets. The data comes from different platforms, including array and RNA-seq platforms. To enable comparisons across platforms, we perform sequence analysis and gene assignment based on the current genome annotations: https://pubmed.ncbi.nlm.nih.gov/16237126/
#
# Data was normalized using RMA, which outputs log2 transformed data: https://pubmed.ncbi.nlm.nih.gov/12925520/. Positive logFC indicate the logarithmic foldness of UPregulation. Negative logFC indicate the logarithmic foldness of DOWNregulation

raw_compendium = pd.read_csv(raw_compendium_filename, sep="\t", header=0, index_col=0)
print(raw_compendium.shape)
raw_compendium.head()

# Are there any negative values in the raw downloaded data
# Expression data was downloaded from GEMMA (https://gemma.msl.ubc.ca/home.html)
# Description for how the data was processed can be found: https://pavlidislab.github.io/Gemma/curation.html
# Looks like there are some negative values in the raw downloaded data.
# Possibly from the log transform
tmp = raw_compendium[(raw_compendium < 0).any(axis=1)]
tmp2 = tmp.loc["GSE12108_Biomat_8___BioAssayId=46616Name=SchuS4infectedD8_S4"]
tmp2[tmp2 < 0]

# ### Process compendium data
#
# 1. Drop probe column
# 2. Transpose
# 3. Get only shared genes from Crow et. al.
# 4. Normalize

# I manually looked up GSE experiment ids with ~3K - ~9K NA genes using gemma.msl.ubc.ca. These experiments were all using the  GPL570 (Affymetrix Human Genome U133 Plus 2.0 Array) platform described in [Crow et al.](https://www.pnas.org/content/116/13/6491).
#
# Some examples:
# * https://gemma.msl.ubc.ca/expressionExperiment/showExpressionExperiment.html?id=3195
# * https://gemma.msl.ubc.ca/expressionExperiment/showExpressionExperiment.html?id=7470
# * https://gemma.msl.ubc.ca/expressionExperiment/showExpressionExperiment.html?id=7934
#
# If i look up the genes that have NA, it says there is "no data" though I'm not sure what is causing this. So for now I will try to remove as many of the NAs as I can with filtering.

# +
# Matrix is sample x gene
# Note, there are NaNs in the matrix. I'm not exactly sure what is causing this since the data
# downloaded was based on what Crow et al used in their analysis and therefore should have been
# filtered by platform.

# Filter out samples that are all NaNs ***
samples_to_drop = raw_compendium[
    raw_compendium.isna().sum(axis=1) == raw_compendium.shape[1]
].index
processed_compendium = raw_compendium.drop(samples_to_drop)

# All genes have at least 1 NaN so dropping all genes with NaN removes all the data
# Instead we will move genes if they are NaN in _most_ samples (>90%)
x = processed_compendium.isna().sum()
processed_compendium = processed_compendium[x.index[x < 100]]

# Remove remaining samples with NaNs
processed_compendium = processed_compendium[processed_compendium.isna().sum(axis=1) < 1]

assert processed_compendium.isna().sum().sum() == 0

# +
# Log transformed the data since the spread is very large
# If we do not log10 transform the data, the few very samples with very large expression values
# will be the dominant signal that will drown out the small signals in the data
# Depending on the biological interest of these large expression signals we can choose to process the
# data differently. In this case, since we are looking at the data wholistically we will use the
# log10 transform of the data here.
processed_compendium_transform = np.log10(processed_compendium)

# Replace -inf (from 0 input) and nans (from negative values) with 0
# Negatives probably created from taking the log of very small values
# Normalization performed using RMA
# Note: the input data is already log2 transformed from the RMA bioconductor library
# so using pseudocounts will not work in this case
processed_compendium_transform = processed_compendium_transform.replace(
    -np.inf, 0.0
).replace(np.nan, 0.0)

assert processed_compendium_transform.isna().sum().sum() == 0

# +
# Get only gene expression data for genes in Crow et. al.
our_gene_ids_hgnc = list(processed_compendium_transform.columns)

published_generic_genes = process.get_published_generic_genes(DE_prior_filename)
shared_genes_hgnc = list(set(our_gene_ids_hgnc).intersection(published_generic_genes))

# In Python, the order of elements in a list that is converted from a set
# is non-deterministic, so it is sorted here to have reproducible result.
shared_genes_hgnc.sort()

# Pickle `shared_genes_hgnc` and save as `shared_genes_filename`
if not os.path.exists(shared_genes_filename):
    with open(shared_genes_filename, "wb") as pkl_fh:
        pickle.dump(shared_genes_hgnc, pkl_fh)

mapped_compendium = processed_compendium_transform[shared_genes_hgnc]

assert mapped_compendium.isna().sum().sum() == 0

print(mapped_compendium.shape)
mapped_compendium.head(10)
# -

# Save
mapped_compendium.to_csv(mapped_compendium_filename, sep="\t")

# Normalize data
process.normalize_compendium(
    mapped_compendium_filename, normalized_compendium_filename, scaler_filename
)

# Check normalized data values
normalized_compendium = pd.read_csv(
    normalized_compendium_filename, sep="\t", index_col=0
)

normalized_compendium.head(10)

# +
# Check that there are no negative values
# -

sns.displot(normalized_compendium["AA06"])

# ### Select and process template data
#
# 1. Get gene expression associated with `project_id`, which was manually selected by the user and specified in the config file.
#
# Note: The data is not normalized so that we can perform DE analysis in next notebook

# +
# metadata file with mapping from experiment to sample
experiment_to_sample_metadata_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", "experiment_sample_annotations.csv"
)

sample_ids = simulate_expression_data.get_sample_ids(
    experiment_to_sample_metadata_filename,
    ",",
    "Experiment id",
    project_id,
    metadata_colname,
)

# Get expression data
template_mapped = mapped_compendium.loc[sample_ids]
print(template_mapped.shape)

# Save
template_mapped.to_csv(mapped_template_filename, sep="\t")
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
