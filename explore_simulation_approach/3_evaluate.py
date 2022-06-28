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

# # Evaluation
#
# How often do regular differential expression analysis vs sophie prioritize the specific vs generic genes?
#
# 1. Simulate 1 template perturbation experiment using the technique above
# 2. Apply SOPHIE to get ranking of specific and generic genes based on their z-score.
# 3. Apply traditional DE analysis and get ranking of specific and generic genes based on their log fold change value
# 4. Compare the difference in ranking between specific and generic genes using SOPHIE vs traditional metrics.

# +
# %load_ext autoreload
# %autoreload 2
# %load_ext rpy2.ipython
import os
import pickle
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from rpy2.robjects import pandas2ri
from ponyo import utils, train_vae_modules, simulate_expression_data
from generic_expression_patterns_modules import (
    process,
    new_experiment_process,  # REMOVE
    stats,
    ranking,
)

np.random.seed(1)

# +
# Read in config variables
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))

config_filename = "config_sophie_vs_trad.tsv"

params = utils.read_config(config_filename)

# +
# Load config params

# Local directory to store intermediate files
local_dir = params["local_dir"]

#
dataset_name = params["dataset_name"]

# File containing un-normalized template experiment
raw_template_filename = params["raw_template_filename"]

# Un-normalized compendium filename
raw_compendium_filename = params["raw_compendium_filename"]

# Normalized compendium filename
normalized_compendium_filename = params["normalized_compendium_filename"]

# ID for template experiment to be selected
project_id = params["project_id"]

# Number of simulated experiments to generate
num_runs = params["num_simulated"]

# Directory containing trained VAE model
vae_model_dir = params["vae_model_dir"]

# Size of the latent dimension
latent_dim = params["latent_dim"]

# Scaler transform used to scale compendium data into 0-1 range for training
scaler_filename = params["scaler_filename"]

# Which DE method to use
# We recommend that if data is RNA-seq then use DESeq2
# If data is microarray then use Limma
de_method = params["DE_method"]

# If using DE-seq, setting this parameter will
# remove genes below a certain threshold
count_threshold = params["count_threshold"]

# Metadata file that specifies which samples to keep for DE analysis (Optional)
template_process_samples_filename = params["template_process_samples_filename"]

# Metadata file that specifies sample grouping for DE analysis
template_DE_grouping_filename = params["template_DE_grouping_filename"]

# Statistic to use to rank genes or pathways by
# Choices are {} FILL IN
col_to_rank_genes = params["rank_genes_by"]

# Pickle files saving specific and generic gene ids
template_specific_gene_ids_filename = params["template_specific_gene_ids_filename"]
generic_gene_ids_filename = "generic_gene_ids.pickle"

# +
# Files generated by this notebook

# File storing template experiment with gene ids mapped to compendium gene ids
mapped_template_filename = params["mapped_template_filename"]

# File storing normalized template experiment
normalized_template_filename = params["normalized_template_filename"]

# File storing processed template experiment,
# after samples have been selected for comparison in DE analysis
processed_template_filename = params["processed_template_filename"]

# Output summary file
output_filename = params["output_filename"]
# -

# ## SOPHIE

# Process template
new_experiment_process.process_template_experiment(
    raw_template_filename,
    raw_compendium_filename,
    scaler_filename,
    mapped_template_filename,
    normalized_template_filename,
)

# +
# Simulate multiple experiments UPDATE COMMENT
# This step creates the following files in "<local_dir>/pseudo_experiment/" directory:
#   - selected_simulated_data_SRP012656_<n>.txt
#   - selected_simulated_encoded_data_SRP012656_<n>.txt
#   - template_normalized_data_SRP012656_test.txt
# in which "<n>" is an integer in the range of [0, num_runs-1]

# REMOVE LATER
# dataset_name = "pre_model_unseen_template"
# Load pickled file
scaler = pickle.load(open(scaler_filename, "rb"))

# Update simulated dir
os.makedirs(os.path.join(local_dir, "pseudo_experiment"), exist_ok=True)

# Update to take in file to be consisten
normalized_compendium = pd.read_csv(
    normalized_compendium_filename, header=0, sep="\t", index_col=0
)
normalized_template = pd.read_csv(
    normalized_template_filename, header=0, sep="\t", index_col=0
)
# ------------
# Update call when new version of ponyo
for run_id in range(num_runs):
    new_experiment_process.embed_shift_template_experiment(
        normalized_compendium,
        normalized_template,
        vae_model_dir,
        project_id,
        scaler_filename,
        local_dir,
        latent_dim,
        run_id,
    )

# +
## Update simulated dir
if not os.path.exists(template_process_samples_filename):
    template_process_samples_filename = None

if de_method == "deseq":
    # Process template data
    stats.process_samples_for_DESeq(
        raw_template_filename,
        template_DE_grouping_filename,
        processed_template_filename,
        count_threshold,
        template_process_samples_filename,
    )

    # Process simulated data
    for i in range(num_runs):
        simulated_filename = os.path.join(
            local_dir,
            "pseudo_experiment",
            f"selected_simulated_data_{project_id}_{i}.txt",
        )
        out_simulated_filename = os.path.join(
            local_dir,
            "pseudo_experiment",
            f"selected_simulated_data_{project_id}_{i}_processed.txt",
        )
        stats.process_samples_for_DESeq(
            simulated_filename,
            template_DE_grouping_filename,
            out_simulated_filename,
            count_threshold,
            template_process_samples_filename,
        )
else:
    stats.process_samples_for_limma(
        raw_template_filename,
        template_DE_grouping_filename,
        processed_template_filename,
        template_process_samples_filename,
    )

    for i in range(num_runs):
        simulated_filename = os.path.join(
            local_dir,
            "pseudo_experiment",
            f"selected_simulated_data_{project_id}_{i}.txt",
        )
        stats.process_samples_for_limma(
            simulated_filename,
            template_DE_grouping_filename,
            None,
            template_process_samples_filename,
        )
# -

# Create subdirectory: "<local_dir>/DE_stats/"
os.makedirs(os.path.join(local_dir, "DE_stats"), exist_ok=True)

# + magic_args="-i template_DE_grouping_filename -i project_id -i processed_template_filename -i local_dir -i base_dir -i de_method" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/DE_analysis.R'))
#
# # File created: "<local_dir>/DE_stats/DE_stats_template_data_<project_id>_real.txt"
# if (de_method == "deseq"){
#     get_DE_stats_DESeq(
#         template_DE_grouping_filename,
#         project_id,
#         processed_template_filename,
#         "template",
#         local_dir,
#         "real"
#     )
# }
# else{
#     get_DE_stats_limma(
#         template_DE_grouping_filename,
#         project_id,
#         processed_template_filename,
#         "template",
#         local_dir,
#         "real"
#     )
# }

# + magic_args="-i template_DE_grouping_filename -i project_id -i base_dir -i local_dir -i num_runs -i de_method" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/DE_analysis.R'))
#
# # Files created: "<local_dir>/DE_stats/DE_stats_simulated_data_<project_id>_<n>.txt"
# for (i in 0:(num_runs-1)){
#     simulated_data_filename <- paste(
#         local_dir,
#         "pseudo_experiment/selected_simulated_data_",
#         project_id,
#         "_",
#         i,
#         "_processed.txt",
#         sep = ""
#     )
#     if (de_method == "deseq"){
#         get_DE_stats_DESeq(
#             template_DE_grouping_filename,
#             project_id,
#             simulated_data_filename,
#             "simulated",
#             local_dir,
#             i
#             )
#     }
#     else {
#         get_DE_stats_limma(
#             template_DE_grouping_filename,
#             project_id,
#             simulated_data_filename,
#             "simulated",
#             local_dir,
#             i
#             )
#         }
#     }

# +
analysis_type = "DE"
template_DE_stats_filename = os.path.join(
    local_dir, "DE_stats", f"DE_stats_template_data_{project_id}_real.txt"
)

# Added
if de_method == "deseq":
    logFC_name = "log2FoldChange"
    pvalue_name = "padj"
else:
    logFC_name = "logFC"
    pvalue_name = "adj.P.Val"

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

# +
# Get summary table
summary_gene_ranks = ranking.generate_summary_table(
    template_DE_stats_filename,
    template_DE_stats,
    simulated_DE_summary_stats,
    col_to_rank_genes,
    local_dir,
    "gene",
    params,
)

summary_gene_ranks.sort_values(by="Z score", ascending=False).head(10)
# -

summary_gene_ranks_sorted = summary_gene_ranks.sort_values(
    by="Z score", ascending=False
)

# Add ranking based on Z-score
summary_gene_ranks_sorted["rank"] = summary_gene_ranks_sorted["Z score"].rank(
    ascending=True
)

summary_gene_ranks_sorted.head(10)

# ## Traditional DE

# + magic_args="-i template_DE_grouping_filename -i project_id -i processed_template_filename -i local_dir -i base_dir -i de_method" language="R"
#
# source(paste0(base_dir, '/generic_expression_patterns_modules/DE_analysis.R'))
#
# # File created: "<local_dir>/DE_stats/DE_stats_template_data_<project_id>_real.txt"
# if (de_method == "deseq"){
#     get_DE_stats_DESeq(
#         template_DE_grouping_filename,
#         project_id,
#         processed_template_filename,
#         "template",
#         local_dir,
#         "real"
#     )
# }
# else{
#     get_DE_stats_limma(
#         template_DE_grouping_filename,
#         project_id,
#         processed_template_filename,
#         "template",
#         local_dir,
#         "real"
#     )
# }

# +
# Load DE statistics file
trad_de_stats_filename = os.path.join(
    local_dir, "DE_stats", f"DE_stats_template_data_{project_id}_real.txt"
)

trad_de_stats = pd.read_csv(trad_de_stats_filename, sep="\t", index_col=0, header=0)
# -

# Sort by log fold change
trad_de_stats_sorted = trad_de_stats.sort_values(by="log2FoldChange", ascending=False)

# Add ranking based on log2FoldChange
trad_de_stats_sorted["rank"] = trad_de_stats_sorted["log2FoldChange"].rank(
    ascending=True
)

trad_de_stats_sorted.head(10)

# ## Compare
#
# 1. mean rank of specific genes - mean rank of generic genes for template experiment
# 2. We will need to re-run this notebook for each template experiment and then plot the distribution of difference scores
#
# We want to compare the mean ranking of specific genes vs the mean ranking of generic genes. If the mean difference is large then yes, that would indicate that there is a difference between the specific and generic genes that we can detect. In addition to the difference we want the specific genes to be higher ranked compared to the generic ones, so we want to see a large positive value if the method is performing better.

# +
# Load pickled file
with open(template_specific_gene_ids_filename, "rb") as specific_fh:
    specific_gene_ids = pickle.load(specific_fh)

with open(generic_gene_ids_filename, "rb") as generic_fh:
    generic_gene_ids = pickle.load(generic_fh)
# -

# Get mean of specific gene ranks
sophie_specific_mean = summary_gene_ranks_sorted.loc[specific_gene_ids, "rank"].mean()
trad_specific_mean = trad_de_stats_sorted.loc[specific_gene_ids, "rank"].mean()

print(sophie_specific_mean)
print(trad_specific_mean)

summary_gene_ranks_sorted.loc[specific_gene_ids, "rank"]

trad_de_stats_sorted.loc[specific_gene_ids, "rank"]

# Get mean of generic gene ranks
sophie_generic_mean = summary_gene_ranks_sorted.loc[generic_gene_ids, "rank"].mean()
trad_generic_mean = trad_de_stats_sorted.loc[generic_gene_ids, "rank"].mean()

print(sophie_generic_mean)
print(trad_generic_mean)

summary_gene_ranks_sorted.loc[generic_gene_ids, "rank"]

trad_de_stats_sorted.loc[generic_gene_ids, "rank"]

# +
# Difference
diff_sophie = sophie_specific_mean - sophie_generic_mean
diff_trad = trad_specific_mean - trad_generic_mean

print("sophie difference: ", diff_sophie)
print("traditional difference: ", diff_trad)
