# -*- coding: utf-8 -*-
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

# # Identify generic genes and pathways
#
# This notebook is meant to identify common differentially expressed genes in the PAO1 and PA14 compendia

# +
# %load_ext autoreload
# %load_ext rpy2.ipython
# %autoreload 2
# %matplotlib inline

import os
import sys
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import pickle
import scipy.stats as ss
from keras.models import load_model
from rpy2.robjects import pandas2ri
from ponyo import utils, simulate_expression_data
from generic_expression_patterns_modules import process, stats, ranking

pandas2ri.activate()

np.random.seed(123)

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = os.path.abspath(
    os.path.join(base_dir, "configs", "config_pseudomonas_pa14_rnaseq.tsv")
)
params = utils.read_config(config_filename)

# +
# Load params
local_dir = params["local_dir"]
dataset_name = params["dataset_name"]
NN_architecture = params["NN_architecture"]
num_runs = params["num_simulated"]
project_id = params["project_id"]
metadata_col_id = params["metadata_colname"]
raw_template_filename = params["raw_template_filename"]
processed_template_filename = params["processed_template_filename"]
normalized_compendium_filename = params["normalized_compendium_filename"]
scaler_filename = params["scaler_filename"]
col_to_rank_genes = params["rank_genes_by"]
logFC_name = params["DE_logFC_name"]
pvalue_name = params["DE_pvalue_name"]
latent_dim = params["latent_dim"]

# Load metadata file with grouping assignments for samples
sample_id_metadata_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", f"{project_id}_process_samples.tsv"
)

# Load metadata file with grouping assignments for samples
metadata_filename = os.path.join(
    base_dir, dataset_name, "data", "metadata", f"{project_id}_groups.tsv"
)

# Load pickled file
scaler = pickle.load(open(scaler_filename, "rb"))

# Percentile threshold to identify generic genes
percentile_threshold = 80.0

metadata_simulate_filename = "data/metadata/SraRunTable.csv"
metadata_delimiter = ","
experiment_id_colname = "SRA_study"
# -

# Output files
gene_summary_filename = os.path.join(
    base_dir, dataset_name, f"generic_gene_summary_{project_id}.tsv"
)

# ## Need to customize code from ponyo
#
# The current simulation-related function in ponyo, `get_sample_ids` assumes that the user is using one of two different metadata files (one associated with the pseudomonas compendium and another associated with recount2). The compendium dataset we are using here has a slightly different format for their metadata file.
#
# Here we are temporarily writing our own function customized for this Pa RNA-seq compendia

# ### Simulate experiments using selected template experiment
# Workflow:
#
# 1. Get the gene expression data for the selected template experiment
# 2. Encode this experiment into a latent space using the trained VAE model
# 3. Linearly shift the encoded template experiment in the latent space
# 4. Decode the samples. This results in a new experiment
# 5. Repeat steps 1-4 to get multiple simulated experiments

# Simulate multiple experiments
# This step creates the following files in "<local_dir>/pseudo_experiment/" directory:
#   - selected_simulated_data_SRP012656_<n>.txt
#   - selected_simulated_encoded_data_SRP012656_<n>.txt
#   - template_normalized_data_SRP012656_test.txt
# in which "<n>" is an integer in the range of [0, num_runs-1]
os.makedirs(os.path.join(local_dir, "pseudo_experiment"), exist_ok=True)
simulate_expression_data.shift_template_experiment(
    normalized_compendium_filename,
    NN_architecture,
    latent_dim,
    dataset_name,
    scaler,
    metadata_simulate_filename,
    metadata_delimiter,
    experiment_id_colname,
    metadata_col_id,
    project_id,
    local_dir,
    base_dir,
    num_runs,
)

# +
simulated_filename = os.path.join(
    local_dir, "pseudo_experiment", f"selected_simulated_data_{project_id}_1.txt"
)

test = pd.read_csv(simulated_filename, sep="\t", index_col=0, header=0)
# -

test.head()

# ### Process template and simulated data
#
# * Remove samples not required for comparison.
# * Make sure ordering of samples matches metadata for proper comparison

# +
if not os.path.exists(sample_id_metadata_filename):
    sample_id_metadata_filename = None

stats.process_samples_for_DESeq(
    raw_template_filename,
    metadata_filename,
    processed_template_filename,
    sample_id_metadata_filename,
)

for i in range(num_runs):
    simulated_filename = os.path.join(
        local_dir, "pseudo_experiment", f"selected_simulated_data_{project_id}_{i}.txt"
    )
    stats.process_samples_for_DESeq(
        simulated_filename,
        metadata_filename,
        None,
        sample_id_metadata_filename,
    )

# +
# Quick check
template_data = pd.read_csv(
    processed_template_filename, header=0, index_col=0, sep="\t"
)

assert template_data.shape[0] == 6
# -

template_data.head()

# ### Differential expression analysis

# Create subdirectory: "<local_dir>/DE_stats/"
os.makedirs(os.path.join(local_dir, "DE_stats"), exist_ok=True)

# + magic_args="-i metadata_filename -i project_id -i processed_template_filename -i local_dir -i base_dir" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/DE_analysis.R'))
#
# get_DE_stats_DESeq(metadata_filename,
#                    project_id,
#                    processed_template_filename,
#                    "template",
#                    local_dir,
#                    "real")

# + magic_args="-i metadata_filename -i project_id -i base_dir -i local_dir -i num_runs -o num_sign_DEGs_simulated" language="R"
#
# source(paste0(base_dir,'/generic_expression_patterns_modules/DE_analysis.R'))
#
# num_sign_DEGs_simulated <- c()
#
# for (i in 0:(num_runs-1)){
#     simulated_data_filename <- paste(
#         local_dir,
#         "pseudo_experiment/selected_simulated_data_",
#         project_id,
#         "_",
#         i,
#         ".txt",
#         sep=""
#     )
#
#     run_output <- get_DE_stats_DESeq(
#         metadata_filename,
#         project_id,
#         simulated_data_filename,
#         "simulated",
#         local_dir,
#         i
#     )
#     num_sign_DEGs_simulated <- c(num_sign_DEGs_simulated, run_output)
# }
# -

# ### Rank genes

analysis_type = "DE"
template_DE_stats_filename = os.path.join(
    local_dir, "DE_stats", f"DE_stats_template_data_{project_id}_real.txt"
)
template_DE_stats, simulated_DE_summary_stats = ranking.process_and_rank_genes_pathways(
    template_DE_stats_filename,
    local_dir,
    num_runs,
    project_id,
    analysis_type,
    col_to_rank_genes,
    logFC_name,
    pvalue_name,
)

# ### Gene summary table

# +
summary_gene_ranks = ranking.generate_summary_table(
    template_DE_stats_filename,
    template_DE_stats,
    simulated_DE_summary_stats,
    col_to_rank_genes,
    local_dir,
    "gene",
    params,
)

summary_gene_ranks.sort_values(by="Z score", ascending=False).head()
# -

# Add gene name as column to summary dataframe
summary_gene_ranks = ranking.add_pseudomonas_gene_name_col(summary_gene_ranks, base_dir)
summary_gene_ranks.sort_values(by="Z score", ascending=False).head()

summary_gene_ranks.sort_values(by="Percentile (simulated)", ascending=False).head()

# Check if there is an NaN values, there should not be
summary_gene_ranks.isna().any()

# Create `gene_summary_filename`
summary_gene_ranks.to_csv(gene_summary_filename, sep="\t")

# ## Compare gene ranking
#
# We can only compare the ranking between the PAO1 RNA-seq compendium vs GAPE, where we still see good concordance as expected.
#
# When we look for common genes, we do this based on percentiles generated by SOPHIE for both the PAO1 and PA14 compendia to be consistent.

if "pao1" in config_filename:
    # Get generic genes identified by Crow et. al.
    GAPE_filename = params["reference_gene_filename"]
    ref_gene_col = params["reference_gene_name_col"]
    ref_rank_col = params["reference_rank_col"]

    figure_filename = f"gene_ranking_{col_to_rank_genes}.svg"

    corr, shared_ranking = ranking.compare_gene_ranking(
        summary_gene_ranks, GAPE_filename, ref_gene_col, ref_rank_col, figure_filename
    )
